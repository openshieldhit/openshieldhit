#include "gemca/parse/osh_gemca2_parse_zone.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/parse/osh_gemca2_parse_keys.h"
#include "gemca/parse/osh_gemca2_parse_stack.h"

// static char* _skip_zone_op(char *s);
static int _key_is_zone_continuation(const char *key);
static struct body *_body_from_name(char *bname, struct gemca_workspace *g);

static struct cgnode *_new_node_comp(struct stack *st, char operator);
static struct cgnode *_new_node_body(struct body *b);
static struct cgnode *_build_ast(struct zone *z, struct gemca_workspace *g);

static size_t _reformat(char const *input, char **output);
static int _tokenizer(char const *input, char ***t);
static int _reverse_tokens(char **tokens, int ntokens);
static int _is_operator(char o);
static void _concat(char **a, char const *b);

/**
 * @brief Initialize a zone.
 *
 * @param[in,out] **body - body to be initialized, memory will be allocated and zeroed,
 *    and memory will also be allocated to the sub-bodies in this zone.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca_zone_init(struct zone **zone) {

    struct zone *z;

    z = calloc(1, sizeof(struct zone));

    if (z == NULL) {
        osh_error(EX_UNAVAILABLE, "osh_gemca_zone_init: could not allocate memory\n");
    }

    z->medium = 0;
    z->id = 0;
    // z->node = NULL;
    z->name = NULL;

    *zone = z;
    return 1;
}

/**
 * @brief Count the number of zones in the given geo.dat file.
 *
 * @param[in] shf - a oshfile struct which also is aware of the line number
 *
 * @returns The number of zones found in the geo.dat file.
 *
 * @author Niels Bassler
 */
size_t osh_gemca_parse_count_zones(struct oshfile *shf) {
    int lineno;
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;

    int nzone = 0;
    int in_block = 0;

    rewind(shf->fp);

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(OSH_GEMCA_KEY_END, key) != 0) {
            if (in_block == 2) {
                if (!_key_is_zone_continuation(key)) {
                    /* new zone found */
                    ++nzone;
                }
            }
        } else {
            /* proceed to next block */
            in_block++;
        }

        if (in_block == 0)
            /* were done reading the header */
            in_block = 1;

        free(line);
    }
    free(line);
    printf("Found %i zones in geo.dat file\n", nzone);
    return nzone;
}

/**
 * @brief Parse zone information
 *
 * @details This function parses the second part of the geo.dat file.
 * (1st part is body description, 2nd is zone description, 3rd is material description)
 *
 * @param[in] *fp - file pointer to file open for reading.
 * @param[in,out] *g - gemca workspace
 *
 * @returns 1 if format is OK, 0 otherwise
 *
 * @author Niels Bassler
 */
int osh_gemca_parse_zones(struct oshfile *shf, struct gemca_workspace *g) {

    char *key = NULL;
    char *args = NULL;
    char *line = NULL;

    char *bstr = NULL;    /* entire body string as given by user */
    char *tstr = NULL;    /* entire body string, formatted parser-friendly */
    char **tokens = NULL; /* list of tokens */

    size_t izone;
    int zone_active = 0;

    int lineno;  /* current line number */
    int ntokens; /* number of tokens */
    int i;
    int len;

    /* move to the second end statement */
    rewind(shf->fp);

    /* forward to the first line after the first END statement */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(OSH_GEMCA_KEY_END, key) == 0) {
            break;
        }
        free(line);
    }
    free(line);

    /* prepare the temprary bstr buffer, which will hold all the user given logic for a single zone */
    /* important: it must hold a NULL byte, so the _concat() function will work. */
    bstr = calloc(1, sizeof(char));
    tstr = calloc(1, sizeof(char));

    /* now proceed parsing */
    zone_active = 0;
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        // printf("LINE %i: KEY='%s' ARGS='%s'\n", lineno, key, args);

        /**
         * Checks if the KEY is an existing body name and if the args start with a valid operator.
         * If KEY matches a body name AND args begins with '+', '-', or '|' operator,
         * it is recognized as a body reference. The 'OR' operator is not supported at the
         * beginning of a line; use '|' instead.
         * If these conditions are not met, the input is treated as a new zone description.
         */
        /* check is KEY an existing body name. */
        /* If not, then assume that this is a new zone description. */
        /* A key is only a body if it matches a body name AND args starts with a valid operator prefix. */
        if (_key_is_zone_continuation(key)) {
            /* Continuation before any zone header is illegal */
            if (!zone_active) {
                osh_error(EX_CONFIG, "zone continuation before first zone at line %d", lineno);
            }
            _concat(&bstr, key);
        } else {
            // printf("Found ZONE name or END: '%s'\n", key);
            /* key is not a body, then it must be a new zone name or the END card.
               1) Process the old zone and attach it to the zone */
            if (zone_active && strlen(bstr) > 0) {

                printf("\n");
                printf("------------------------------------------------------------------------------\n");
                printf("ZONE: #%3lli - '%s'\n", (long long int) g->zones[izone]->id, g->zones[izone]->name);
                printf("USERGIVEN STRING: '%s'\n", bstr);
                _reformat(bstr, &tstr);
                printf("PRE-TOKEN STRING: '%s'\n", tstr);
                ntokens = _tokenizer(tstr, &tokens);
                _reverse_tokens(tokens, ntokens);

                printf("number of tokens: %i\n", ntokens);

                for (i = 0; i < ntokens; i++) {
                    printf("token #%i '%s'\n", i, tokens[i]);
                }
                /* attach the tokens of this zone to the geometry */
                g->zones[izone]->tokens = tokens;
                g->zones[izone]->ntokens = ntokens;

                /* build the abstract syntax tree from the tokens */
                _build_ast(g->zones[izone], g);

                bstr[0] = '\0'; /* reset bstr */
                tstr[0] = '\0'; /* reset tstr */
            }

            /* 2) check if we have reached the second END statement */
            if (strcasecmp(key, OSH_GEMCA_KEY_END) == 0) {
                /* there is nothing more to read, bail out */
                break; /* break out of while loop */
            }

            /* 3) new zone: increment zone index and allocate memory for the name, and copy it into the placeholder.*/
            if (!zone_active) {
                zone_active = 1;
                izone = 0;
            } else {
                izone++;
            }
            /* we may have crossed the END card, which will exceed the number of zones */
            len = strlen(key);
            g->zones[izone]->name = calloc(len + 1, sizeof(char));
            g->zones[izone]->id = izone + 1;  /* zone IDs start at 1, not at 0 */
            g->zones[izone]->lineno = lineno; /* save the line number where this zone was defined */
            strncpy(g->zones[izone]->name, key, len + 1);
        }

        /* now parse whatever is left in the args string, these are always logic and bodies */
        _concat(&bstr, args);
        free(line);
    } /* end while loop */

    free(line);
    free(bstr);
    free(tstr);

    return 0;
}

/* @brief Check if the given key is a zone continuation operator */
static int _key_is_zone_continuation(const char *key) {
    return key && (key[0] == '+' || key[0] == '-' || key[0] == '|' || key[0] == '(' || key[0] == ')');
}

/**
 * @brief lookup in g->bodies for a body with the name bname, and return a pointer to this body, if found.
 *
 * @param[in] *bname - character string holding the body name to be looked up
 * @param[in] *g - gemca workspace, where the bodies with their names have been loaded already.
 * @param[out] **return_body - pointer to the body found.
 *
 * @returns body or NULL if not found or some problem occurred.
 *
 * @author Niels Bassler
 */
static struct body *_body_from_name(char *bname, struct gemca_workspace *g) {

    size_t i = 0;

    // printf("_body_from_name()  '%s'\n", bname);
    for (i = 0; i < g->nbodies; i++) {
        if (strcmp(bname, g->bodies[i]->name) == 0) {
            // printf("_body_from_name()  found '%s'\n", g->bodies[i]->name);
            return g->bodies[i];
        }
    }
    return NULL;
}

/**
 * @brief Allocate and return pointer to a new cgnode of the type _OSH_GEMCA_CGNODE_COMPOSITE.
 *
 * @details This node holds two cgnodes and an operator.
 *          Memory will be allocated and stack will be popped by two items.
 *
 * @param[in] *st - pointer to a stack holding at least 2 cgnodes which will be used to build this node.
 * @param[in] op - the operator which will apply for the two cgnodes in this node: Left node - operator - right node.
 *
 * @returns pointer to a cgnode of type _OSH_GEMCA_CGNODE_COMPOSITE
 *
 * @author Niels Bassler
 */
static struct cgnode *_new_node_comp(struct stack *st, char operator) {

    struct cgnode *node;
    struct stackitem *si = NULL;

    node = calloc(1, sizeof(struct cgnode));
    node->type = _OSH_GEMCA_CGNODE_COMPOSITE;
    node->op = operator;

    // osh_gemca_stack_print(st);

    /* first element popped must go on the left leg */
    si = osh_gemca_stack_pop(st);
    node->left = si->v.cgnode;
    free(si);

    /* second element poppes must go on the right leg */
    si = osh_gemca_stack_pop(st);
    node->right = si->v.cgnode;
    free(si);

    return node;
}

/**
 * @brief Allocate and return pointer to a new cgnode of the type _OSH_GEMCA_CGNODE_BODY
 *
 * @details Memory will be allocated to this node, and the body pointer will be copied into the node.
 *          The node type will be set accordingly.
 *
 * @param[in] *b - pointer to a body which will be initialized
 *
 * @returns pointer to a cgnode of type _OSH_GEMCA_CGNODE_BODY
 *
 * @author Niels Bassler
 */
static struct cgnode *_new_node_body(struct body *b) {

    struct cgnode *node;

    node = calloc(1, sizeof(struct cgnode));
    node->type = _OSH_GEMCA_CGNODE_BODY;
    node->body = b;
    return node;
}

/**
 * @brief Build the abstract syntax tree.
 *
 * @details This uses the shunting yard algorithm, and undertands paranthesis, +, -, | operators.
 *
 * @param[in] *z - input struct zone from which the AST will be generated.
 * @param[in] *g - input struct gemca (needed to lookup the bodies from the token names)
 *
 * @returns pointer to the top level cgnode.
 *
 * @see https://en.wikipedia.org/wiki/Shunting-yard_algorithm (retrieved medio June 2020).
 *      It is based on the "detailed algorithm description" presented there.
 *
 * @author Niels Bassler
 */
static struct cgnode *_build_ast(struct zone *z, struct gemca_workspace *g) {

    size_t i;
    char *token;
    struct stackitem *si = NULL;
    struct cgnode *node = NULL;
    struct stack *opst = NULL;
    struct stack *st = NULL;
    struct body *b;

    for (i = 0; i < z->ntokens; i++) {

        token = z->tokens[i];

        /* check for regular operators which will be pushed onto the operator stack */
        if ((token[0] == '+') || (token[0] == '-') || (token[0] == '|') || (token[0] == '(')) {

            /* allocate memory for a new stack item */
            si = calloc(1, sizeof(struct stackitem));
            si->type = _OSH_GEMCA_STACKITEM_OPERATOR;
            si->v.op = token[0];
            /* push to operator stack. This also allocates any memory */
            osh_gemca_stack_push(&opst, si);
        } else if (token[0] == ')') {
            /* pop the entire operator stack until we see a matching opening '(' */
            while (opst->ni > 0) {
                si = osh_gemca_stack_pop(opst);
                if (si->v.op == '(') {
                    free(si);
                    break;
                } else {
                    /* make a csgnode from the popped operator and push the new node to the csgstack */
                    node = _new_node_comp(st, si->v.op);
                    free(si);
                    si = calloc(1, sizeof(struct stackitem));
                    si->type = _OSH_GEMCA_STACKITEM_CGNODE;
                    si->v.cgnode = node;
                    osh_gemca_stack_push(&st, si);
                }
            }
        } else {
            /* this is a simple body / leaf node. Push it to the stack as such, */
            /* make csgnode from popped and push to csgstack.  */
            /* But first, lookup body from token name. */
            b = _body_from_name(token, g);
            if (b == NULL) {
                osh_error(EX_CONFIG, "_build_ast() coudn't find body names '%s'", token);
            }
            node = _new_node_body(b);

            si = calloc(1, sizeof(struct stackitem));
            si->type = _OSH_GEMCA_STACKITEM_CGNODE;
            si->v.cgnode = node;
            osh_gemca_stack_push(&st, si);
        }
    } /* end of for loop over tokens */

    /* check if there are any remaining operators left on the operator stack. If so, pop them. */
    while (opst->ni > 0) {
        si = osh_gemca_stack_pop(opst);

        if ((si->v.op == '(') || (si->v.op == ')')) {
            osh_error(EX_CONFIG, "unbalanced paranthesis in zone description");
        } else {
            /* make csgnode from popped and push to csgstack */
            node = _new_node_comp(st, si->v.op);
            free(si);
            si = calloc(1, sizeof(struct stackitem));
            si->type = _OSH_GEMCA_STACKITEM_CGNODE;
            si->v.cgnode = node;
            osh_gemca_stack_push(&st, si);
        }
    }
    /* what is left on the stack is the AST object. */
    si = osh_gemca_stack_pop(st);
    z->node = *si->v.cgnode;
    free(si);
    return &z->node;
}

/**
 * @brief Reformat a user-given zone string into a new string which is easier to tokenize.
 *
 * @details - This reformatter makes SHIELD-HIT12A compatible with the free-format from FLUKA,
 *            without breaking the fixed format notation, and the peculiaraity that SH12A allows single 'OR' use
 *            in user strings.
 *
 * @example
 *            '+1'                                 ->       '(1)'
 *            'OR +foobar   OR +6'                 ->       '(foobar)|(6)'
 *            '+1     -4     -foobar     -6'       ->       '(1-4-foobar-6)'
 *
 * @param[in] *input - user-given zone description string, only containing bodies, operators and any whitespaces.
 * @param[in] **output - pointer to a new string. Memory will be reallocated as needed.
 *
 * @returns length of new string
 *
 * @author Niels Bassler
 */
static size_t _reformat(char const *input, char **output) {

    size_t i; /* index for input string */
    size_t j; /* index for output string */
    size_t il;
    size_t ol;

    il = strlen(input);

    /* a conservative estimate how much space we need for the output string.
       normally it is less than the input string, but here we will just allocate at least 0xFF bytes or
       twice the size of the input string */
    ol = 0x100; /* start with just 256 bytes for the output string */
    if (il > ol)
        ol = 2 * il;
    *output = realloc(*output, ol * sizeof(char));
    if (*output == NULL) {
        osh_error(EX_UNAVAILABLE, "_reformat(): cannot malloc");
    }
    memset(*output, 0, ol);

    i = 0;
    j = 0;

    if (input[0] == '\0') {
        return 0; /* string is empty, nothing is to be done */
    }

    /*   add a leading paranthesis */
    (*output)[0] = '(';
    j++;

    /* check if first two bytes are OR or | */
    if ((input[0] == 'O') && (input[1] == 'R')) {
        i += 2; /* skip first two bytes */
    } else if (input[0] == '|') {
        i++; /* skip first byte */
    }

    /* scan input */
    while ((input[i] != '\0') && ((j + 4) < ol)) {
        /* do not add '+' or '-' if there was a leading '(' */
        if ((*output)[j - 1] == '(') {
            if (input[i] == '-') {
                osh_error(EX_CONFIG, "leading body cannot be a '-' body, only '+'");
            } else if (input[i] == '+') {
                /* do not add unary '+' to output stream, i.e. (+2-1) -> (2-1) */
                i++;
                continue;
            }
        }

        /* check if we have a real OR or | */
        if ((input[i] == 'O') && (input[i + 1] == 'R')) {
            i++; /* move an extra byte of the input stream  */
            /* substitute with ')|()' */
            (*output)[j++] = ')';
            (*output)[j++] = '|';
            (*output)[j++] = '(';
        }

        else if (input[i] == '|') {
            /* substitute with ')|()' */
            (*output)[j++] = ')';
            (*output)[j++] = '|';
            (*output)[j++] = '(';
        }

        /* remove any whitespaces and otherwise just copy the next character to the output stream */
        else if (!(isspace(input[i]))) {
            (*output)[j++] = input[i];
        }
        i++;

        /* check if output string will need more memory */
        if ((j + 4) > ol) {
            ol += 0x100; /* add another 256 bytes */
            *output = realloc(*output, ol * sizeof(char));
            if (*output == NULL) {
                osh_error(EX_UNAVAILABLE, "_reformat(): cannot realloc #2");
            }
        }

    } /* end while */

    (*output)[j++] = ')';
    (*output)[j] = '\0';

    /* trim size of output string to the memory actually needed */
    *output = realloc(*output, (j + 1) * sizeof(char));
    if (*output == NULL) {
        osh_error(EX_UNAVAILABLE, "_reformat(): cannot realloc #3");
    }

    return j;
}

/**
 * @brief Tokenize the output from _reformat()
 *
 * @details - ptokens must be a tripple pointer, since memory is allocated here and the address needs to be returned.

 * @param[in] *input - pointer to a char string prepared with _reformat().
 * @param[in] ***ptokens - pointer to an arraypointer of pointers to strings
 *
 * @returns Number of tokens found and parsed from the given string.
 *
 * @author Niels Bassler
 */
static int _tokenizer(char const *input, char ***ptokens) {
    size_t ilen;
    size_t n;
    size_t i;
    size_t j;

    char **tokens = *ptokens;

    /* Instead of realloc() every tooken, we scan first how many tokens there are, and do a single calloc instead. */
    n = 0;
    ilen = strlen(input);

    /* Count the number of tokens needed. */
    for (i = 0; i < ilen; i++) {
        if (_is_operator(input[i])) {
            n++; /* increment token counter for the operator token we just found */
        } else {
            /* if not an operator, it is a body name. Keep scanning characters until we find the next operator */
            while (!(_is_operator(input[i]))) {
                i++;
            }
            i--; /* we scanned already into the operator, so step one back */
            n++; /* increment token counter for the body token we just found */
        }
    }

    /* allocate memory for n tokens for the token list */
    tokens = calloc(n, sizeof(char *));

    /* Now repeat the loop, while also adding the tokens to the token list. */
    n = 0; /* notice, n is used as an index first, since the increment happens after the assignments */
    for (i = 0; i < ilen; i++) {
        if (_is_operator(input[i])) {
            tokens[n] = calloc(2, sizeof(char));
            tokens[n][0] = input[i];
            tokens[n][1] = '\0';
        } else {
            /* This must be a body. Scan ahead until we find the next operator. */
            j = 0;
            while (!(_is_operator(input[i + j]))) {
                j++;
            }

            /* bodies need variable amounts of memoery, since they represented by character strings with var. size. */
            tokens[n] = calloc(j + 1, sizeof(char)); /* carefully including the terminal null byte */
            strncpy(tokens[n], input + i, j);

            /* skip to next operator position */
            i += j - 1;
        }
        n++; /* increment for either operator or for a body */
    }
    /* We are now done making the token list. Attach it to the ptokens pointer so it can be returned */
    *ptokens = tokens;
    /* n is now no longer an index, but holds the actual number of tokens  due to the last icreement */
    return n;
}

/**
 * @brief Reverses the list of tokens, and flips the paranthesises.
 *
 * @param[in,out] **ptokens - array of pointers to strings
 * @param[in] *input - pointer to a char string prepared with _reformat().
 *
 * @returns Number of tokens found and parsed from the given string.
 *
 * @author Niels Bassler
 */
static int _reverse_tokens(char **tokens, int ntokens) {

    char *c; /* temporary placeholder for a pointer to string */
    int i;
    size_t len, j;

    for (i = 0; i < ntokens; i++) {
        // printf("Tokens: '%s'\n", tokens[i]);
    }

    for (i = 0; i < ntokens / 2; i++) {
        c = tokens[i]; /* temporarily save the first pointer */
        tokens[i] = tokens[ntokens - i - 1];
        tokens[ntokens - i - 1] = c;
    }

    printf("\n");
    for (i = 0; i < ntokens; i++) {
        // printf("Tokens: '%s'\n", tokens[i]);
    }

    for (i = 0; i < ntokens; i++) {
        /* reverse any paranthesises */
        len = strlen(tokens[i]);
        // printf("len : %li\n", len);
        for (j = 0; j < len; j++) {
            // printf("%s '%c'\n", tokens[i], tokens[i][j]);
            if (tokens[i][j] == '(') {
                // printf("flipped a (\n");
                tokens[i][j] = ')';
            } else if (tokens[i][j] == ')') {
                // printf("flipped a )\n");
                tokens[i][j] = '(';
            }
        }
    }
    return 1;
}

/**
 * @brief Checks if the input character is an operator.
 *
 * @param[in] o - character to be checked against "+-()|"
 *
 * @returns 1 if o is an operator. 0 If o is not an operator.
 *
 * @author Niels Bassler
 */
static int _is_operator(char o) {

    char const ops[] = "+-()|"; /* available operators in a formatted string for the tokenizer */
    int const nops = 5;         /* number of operators */
    int i;

    for (i = 0; i < nops; i++) {
        if (o == ops[i])
            return 1;
    }
    return 0;
}

/**
 * @brief Copy string b to a, and increase memoery of a accordingly.
 *
 * @details Both a and b must hold a '\0' byte so the initial string lengths can be determined.
 *
 * @param[out] *a - pointer to destination string. Memory will be reallocated if needed.
 * @param[in] b - input string.
 *
 * @returns 1 if o is an operator. 0 If o is not an operator.
 *
 * @author Niels Bassler
 */
static void _concat(char **a, char const *b) {
    int la, lb;

    la = strlen(*a);
    lb = strlen(b);

    *a = realloc(*a, sizeof(char) * (la + lb + 1));
    if (*a == NULL) {
        osh_error(EX_SOFTWARE, "_concat(): cannot reallocate memory");
    }

    strncat(*a, b, lb);
    return;
}
