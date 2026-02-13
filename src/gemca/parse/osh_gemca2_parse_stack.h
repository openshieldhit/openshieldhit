#ifndef _OSH_GEMCA2_PARSE_STACK
#define _OSH_GEMCA2_PARSE_STACK

#include "gemca/osh_gemca2.h"

#define _OSH_GEMCA_STACKITEM_OPERATOR 0
#define _OSH_GEMCA_STACKITEM_CGNODE 1

struct stack {
    struct stackitem **si; /* list of pointers to stack items */
    size_t ni;             /* number of last item on stack. It has index i = ni-1 */
    size_t n;              /* current memory allocated to stack in terms of elements */
};

struct stackitem {
    union {
        struct cgnode *cgnode; /* stack item holds a pointer to a cgnodes */
        char op;               /* stack item holds a single operator */
    } v;
    int type; /* stack type  _OSH_GEMCA_STACKITEM_* */
};

size_t osh_gemca_stack_push(struct stack **ps, struct stackitem *i);
struct stackitem *osh_gemca_stack_pop(struct stack *s);
void osh_gemca_stack_free(struct stack **ps);
void osh_gemca_stack_print(struct stack *s);

#endif /* _OSH_GEMCA2_PARSE_STACK */
