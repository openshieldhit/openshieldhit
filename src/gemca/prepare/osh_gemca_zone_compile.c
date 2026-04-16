#include "gemca/prepare/osh_gemca_zone_compile.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "gemca/osh_gemca2.h"
#include "gemca/prepare/osh_gemca_stack.h"

static struct body *_body_from_name(char *bname, struct gemca_workspace *g);

static struct cgnode *_new_node_comp(struct stack *st, char operator);
static struct cgnode *_new_node_body(struct body *b);
static struct cgnode *_build_ast(struct zone *z, struct gemca_workspace *g);

static size_t _reformat(char const *input, char **output);
static int _tokenizer(char const *input, char ***t);
static void _reverse_tokens(char **tokens, int ntokens);
static int _is_operator(char o);

/**
 * @brief Initialize a zone.
 *
 * @param[out] zone Receives the newly allocated internal zone object.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_gemca_zone_init(struct zone **zone) {

    struct zone *z;

    z = calloc(1, sizeof(struct zone));

    if (z == NULL) {
        osh_error("osh_gemca_zone_init: could not allocate memory");
        return OSH_ENOMEM;
    }

    z->material_idx = (size_t) -1;
    z->name = NULL;
    z->material_name = NULL;

    *zone = z;
    return OSH_OK;
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

    for (i = 0; i < g->nbodies; i++) {
        if (strcmp(bname, g->bodies[i]->name) == 0) {
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
 * @details This uses the shunting yard algorithm, and undertands parentheses, +, -, | operators.
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
            while (opst != NULL && opst->ni > 0) {
                si = osh_gemca_stack_pop(opst);
                if (si->v.op == '(') {
                    free(si);
                    si = NULL;
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
            if (opst == NULL) {
                osh_error("%s:%zu: unbalanced parenthesis in zone description", g->filename, z->lineno);
            }
        } else {
            /* this is a simple body / leaf node. Push it to the stack as such, */
            /* make csgnode from popped and push to csgstack.  */
            /* But first, lookup body from token name. */
            b = _body_from_name(token, g);
            if (b == NULL) {
                osh_error("_build_ast() coudn't find body names '%s'", token);
            }
            node = _new_node_body(b);

            si = calloc(1, sizeof(struct stackitem));
            si->type = _OSH_GEMCA_STACKITEM_CGNODE;
            si->v.cgnode = node;
            osh_gemca_stack_push(&st, si);
        }
    } /* end of for loop over tokens */

    /* check if there are any remaining operators left on the operator stack. If so, pop them. */
    while (opst != NULL && opst->ni > 0) {
        si = osh_gemca_stack_pop(opst);

        if ((si->v.op == '(') || (si->v.op == ')')) {
            osh_error("%s:%zu: unbalanced parenthesis in zone description", g->filename, z->lineno);
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
    if (si == NULL) {
        osh_error("empty zone description");
    }
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
    if (il > ol) {
        ol = 2 * il;
    }
    *output = realloc(*output, ol * sizeof(char));
    if (*output == NULL) {
        osh_error("_reformat(): cannot malloc");
    }
    memset(*output, 0, ol);

    i = 0;
    j = 0;

    if (input[0] == '\0') {
        return 0; /* string is empty, nothing is to be done */
    }

    /*   add a leading parenthesis */
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
                osh_error("leading body cannot be a '-' body, only '+'");
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
                osh_error("_reformat(): cannot realloc #2");
            }
        }

    } /* end while */

    (*output)[j++] = ')';
    (*output)[j] = '\0';

    /* trim size of output string to the memory actually needed */
    *output = realloc(*output, (j + 1) * sizeof(char));
    if (*output == NULL) {
        osh_error("_reformat(): cannot realloc #3");
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
    if (n == 0) {
        *ptokens = NULL;
        return 0;
    }

    *ptokens = (char **) calloc(n, sizeof(char *));
    if (*ptokens == NULL) {
        osh_error("_tokenizer(): cannot allocate token list");
    }

    /* Now repeat the loop, while also adding the tokens to the token list. */
    n = 0; /* notice, n is used as an index first, since the increment happens after the assignments */
    for (i = 0; i < ilen; i++) {
        if (_is_operator(input[i])) {
            (*ptokens)[n] = (char *) calloc(2, sizeof(char));
            (*ptokens)[n][0] = input[i];
            (*ptokens)[n][1] = '\0';
        } else {
            /* This must be a body. Scan ahead until we find the next operator. */
            j = 0;
            while (!(_is_operator(input[i + j]))) {
                j++;
            }

            /* bodies need variable amounts of memoery, since they represented by character strings with var. size. */
            (*ptokens)[n] = (char *) calloc(j + 1, sizeof(char)); /* carefully including the terminal null byte */
            strncpy((*ptokens)[n], input + i, j);

            /* skip to next operator position */
            i += j - 1;
        }
        n++; /* increment for either operator or for a body */
    }
    /* We are now done making the token list. Attach it to the ptokens pointer so it can be returned */
    /* n is now no longer an index, but holds the actual number of tokens  due to the last icreement */
    return n;
}

/**
 * @brief Reverses the list of tokens, and flips the parentheses.
 *
 * @param[in,out] **ptokens - array of pointers to strings
 * @param[in] *input - pointer to a char string prepared with _reformat().
 *
 * @returns Number of tokens found and parsed from the given string.
 *
 * @author Niels Bassler
 */
static void _reverse_tokens(char **tokens, int ntokens) {

    char *c; /* temporary placeholder for a pointer to string */
    int i;
    size_t len;
    size_t j;

    for (i = 0; i < ntokens / 2; i++) {
        c = tokens[i]; /* temporarily save the first pointer */
        tokens[i] = tokens[ntokens - i - 1];
        tokens[ntokens - i - 1] = c;
    }

    for (i = 0; i < ntokens; i++) {
        /* reverse any parentheses */
        len = strlen(tokens[i]);
        for (j = 0; j < len; j++) {
            if (tokens[i][j] == '(') {
                tokens[i][j] = ')';
            } else if (tokens[i][j] == ')') {
                tokens[i][j] = '(';
            }
        }
    }
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
        if (o == ops[i]) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Compile a zone boolean expression into a cgnode AST.
 *
 * @details Exposes the _reformat → _tokenizer → _reverse_tokens → _build_ast
 *          pipeline as a public entry point so that callers who build cold
 *          struct zone objects without going through the file parser can still
 *          produce a ready-to-use cgnode tree.
 *
 *          Bodies must already be set up in @p g before this is called because
 *          _build_ast() looks up body pointers by name.
 *
 * @param[in,out] z     Zone to populate.
 * @param[in]     expr  Raw zone expression (e.g. "+BODY1 -BODY2 | +BODY3").
 * @param[in]     g     Gemca workspace with bodies[] already initialised.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
/* Returns 1 if any leaf node in the cgnode tree has a NULL body pointer. */
static int _has_null_body(struct cgnode const *node) {
    if (!node) {
        return 0;
    }
    if (node->type == _OSH_GEMCA_CGNODE_BODY) {
        return node->body == NULL;
    }
    return _has_null_body(node->left) || _has_null_body(node->right);
}

enum osh_status osh_gemca_zone_compile_expr(struct zone *z, char const *expr, struct gemca_workspace *g) {
    char *tstr = NULL;
    char **tokens = NULL;
    int ntokens;

    if (!z || !expr || !g) {
        return OSH_EINVAL;
    }

    /* _reformat() uses realloc(), so the base pointer must be non-NULL and
     * point to at least a NUL byte so that strlen() is safe on it. */
    tstr = calloc(1, sizeof(char));
    if (!tstr) {
        return OSH_ENOMEM;
    }

    _reformat(expr, &tstr);
    ntokens = _tokenizer(tstr, &tokens);
    _reverse_tokens(tokens, ntokens);
    free(tstr);

    z->tokens = tokens;
    z->ntokens = (size_t) ntokens;

    _build_ast(z, g);

    /* _build_ast() logs errors for unresolved body names but does not abort.
     * Check the compiled tree: any NULL body pointer means the zone expr
     * referenced a body name that does not exist in the workspace. */
    if (_has_null_body(&z->node)) {
        return OSH_EPARSE;
    }

    return OSH_OK;
}
