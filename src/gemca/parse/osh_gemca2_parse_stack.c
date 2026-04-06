#include "osh_gemca2_parse_stack.h"

#include <stdlib.h>

#include "common/osh_logger.h"

struct stackitem *osh_gemca_stack_pop(struct stack *s) {
    struct stackitem *si;

    (s->ni)--;
    si = s->si[s->ni];
    return si;
}

/* push a given stackitem pointer onto stack. */
/* Stack will be initilized, if it does not exist */
/* returns the new stack size in number of elements */
size_t osh_gemca_stack_push(struct stack **ps, struct stackitem *i) {

    struct stack *s;
    struct stackitem **tmp_si;
    size_t new_n;

    /* check if stack exists, if not, allocate memory */
    if (*ps == NULL) {
        /* allocate memory to the stack */
        s = (struct stack *) calloc(1, sizeof(struct stack));
        if (s == NULL) {
            osh_alloc_failed("osh_gemca_stack_push(): stack");
        }

        /* allocate memory to 32 stack pointers */
        s->n = 32;
        s->si = (struct stackitem **) calloc(s->n, sizeof(struct stackitem *));
        if (s->si == NULL) {
            free(s);
            osh_alloc_failed("osh_gemca_stack_push(): stack items");
        }

        /* number of elements in stack is so far just 0 */
        s->ni = 0;
    } else {
        s = *ps;
    }

    /* increase number of element counter by one */
    (s->ni)++;

    /* check if we need more memory */
    if ((s->ni) > (s->n)) {
        new_n = s->n + 32; /* allocate memory for another 32 elements */
        tmp_si = (struct stackitem **) realloc((void *) s->si, new_n * sizeof(struct stackitem *));
        if (tmp_si == NULL) {
            /* Keep existing stack untouched to avoid dangling pointers at call sites. */
            osh_alloc_failed("osh_gemca_stack_push(): stack grow");
        }
        s->si = tmp_si;
        s->n = new_n;
    }

    /* push pointer onto stack */
    s->si[(s->ni) - 1] = i;
    *ps = s;

    return s->ni;
}

void osh_gemca_stack_free(struct stack **ps) {
    size_t i;
    if ((ps == NULL) || (*ps == NULL))
        return;
    for (i = 0; i < (*ps)->n; i++) {
        free((*ps)->si[i]);
    }
    free((*ps)->si);
    free(*ps);
    *ps = NULL;
}

void osh_gemca_stack_print(struct stack *s) {
    size_t i;
    osh_debug(OSH_LOG_HLINE);
    osh_debug("STACK : %p\n", (void *) s);
    osh_debug("NELEM : %llu\n", (unsigned long long) s->ni);
    osh_debug("NMEM  : %llu\n", (unsigned long long) s->n);

    for (i = 0; i < s->ni; i++) {
        if (s->si[i]->type == _OSH_GEMCA_STACKITEM_OPERATOR) {
            osh_debug("    StackITEM: %llu: %p  type: %i  '%c'\n",
                      (unsigned long long) i,
                      (void *) s->si[i],
                      s->si[i]->type,
                      s->si[i]->v.op);
        } else {
            osh_debug("    StackITEM: %llu: %p  type: %i\n", (unsigned long long) i, (void *) s->si[i], s->si[i]->type);
        }
    }
}
