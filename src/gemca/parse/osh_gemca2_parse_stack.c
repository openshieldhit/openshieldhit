#include "osh_gemca2_parse_stack.h"

#include <stdlib.h>

struct stackitem *osh_gemca_stack_pop(struct stack *s) {
    struct stackitem *si;

    (s->ni)--;
    // printf("READ FROM INDEX: %li\n", s->ni);
    si = s->si[s->ni];
    // printf("POPPED: %p from %p\n", si, s);
    // osh_gemca_stack_print(s);
    return si;
}

/* push a given stackitem pointer onto stack. */
/* Stack will be initilized, if it does not exist */
/* returns the new stack size in number of elements */
size_t osh_gemca_stack_push(struct stack **ps, struct stackitem *i) {

    struct stack *s;

    /* check if stack exists, if not, allocate memory */
    if (*ps == NULL) {
        /* allocate memory to the stack */
        // printf("Allocate memory to stack\n");

        s = calloc(1, sizeof(struct stack));

        /* allocate memory to 32 stack pointers */
        s->n = 32;
        s->si = calloc(s->n, sizeof(struct stackitem *));

        /* number of elements in stack is so far just 0 */
        s->ni = 0;
    } else {
        s = *ps;
    }

    /* increase number of element counter by one */
    (s->ni)++;

    /* check if we need more memory */
    if ((s->ni) > (s->n)) {
        s->n += 32; /* allocate memory for another 32 elements */
        s->si = realloc(s->si, s->n * sizeof(struct stackitem *));
    }

    // printf("WRITE TO INDEX: %li\n", s->ni - 1);
    // printf("PUSHED: %p on %p\n", i, s);
    /* push pointer onto stack */
    s->si[(s->ni) - 1] = i;
    *ps = s;

    // osh_gemca_stack_print(s);
    return s->ni;
}

void osh_gemca_stack_free(struct stack **ps) {
    size_t i;
    for (i = 0; i < (*ps)->n; i++) {
        free((*ps)->si[i]);
    }
    free(*ps);
}

void osh_gemca_stack_print(struct stack *s) {
    size_t i;
    printf("\n");
    printf("-------------\n");
    printf("STACK : %p\n", (void *) s);
    printf("NELEM : %llu\n", (unsigned long long) s->ni);
    printf("NMEM  : %llu\n", (unsigned long long) s->n);

    for (i = 0; i < s->ni; i++) {
        printf("    StackITEM: %llu: %p  type: %i   ", (unsigned long long) i, (void *) s->si[i], s->si[i]->type);
        if (s->si[i]->type == _OSH_GEMCA_STACKITEM_OPERATOR) {
            printf("'%c'\n", s->si[i]->v.op);
        } else
            printf("\n");
    }

    printf("\n");
}
