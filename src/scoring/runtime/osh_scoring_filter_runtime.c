#include "scoring/runtime/osh_scoring_filter_runtime.h"

#include <math.h>

#include "common/osh_const.h"
#include "scoring/runtime/osh_scoring_runtime.h"

static int filter_compare_int(unsigned int lhs, enum osh_scoring_filter_op op, unsigned int rhs);
static int filter_compare_double(double lhs, enum osh_scoring_filter_op op, double rhs);
static int
filter_rule_passes(struct osh_scoring_filter_runtime_rule const *rule, struct particle const *part, double energy_mev);

int osh_scoring_page_passes_filters(struct osh_scoring_runtime const *rt,
                                    struct osh_scoring_page_runtime const *page,
                                    struct particle const *part,
                                    double energy_mev) {
    size_t i;
    size_t j;
    struct osh_scoring_filter_runtime const *filter;

    for (i = 0; i < page->nfilters; ++i) {
        if (page->filters[i].filter_idx >= rt->nfilters) {
            return 0;
        }
        filter = &rt->filters[page->filters[i].filter_idx];
        for (j = 0; j < filter->nrules; ++j) {
            if (!filter_rule_passes(&filter->rules[j], part, energy_mev)) {
                return 0;
            }
        }
    }
    return 1;
}

static int filter_compare_int(unsigned int lhs, enum osh_scoring_filter_op op, unsigned int rhs) {
    switch (op) {
    case OSH_SCORING_FILTER_OP_LT:
        return lhs < rhs;
    case OSH_SCORING_FILTER_OP_LE:
        return lhs <= rhs;
    case OSH_SCORING_FILTER_OP_GT:
        return lhs > rhs;
    case OSH_SCORING_FILTER_OP_GE:
        return lhs >= rhs;
    case OSH_SCORING_FILTER_OP_EQ:
        return lhs == rhs;
    case OSH_SCORING_FILTER_OP_NE:
        return lhs != rhs;
    default:
        return 0;
    }
}

static int filter_compare_double(double lhs, enum osh_scoring_filter_op op, double rhs) {
    double eps = 1.0e-10;

    switch (op) {
    case OSH_SCORING_FILTER_OP_LT:
        return lhs < rhs;
    case OSH_SCORING_FILTER_OP_LE:
        return lhs <= rhs;
    case OSH_SCORING_FILTER_OP_GT:
        return lhs > rhs;
    case OSH_SCORING_FILTER_OP_GE:
        return lhs >= rhs;
    case OSH_SCORING_FILTER_OP_EQ:
        return fabs(lhs - rhs) < eps;
    case OSH_SCORING_FILTER_OP_NE:
        return fabs(lhs - rhs) >= eps;
    default:
        return 0;
    }
}

static int
filter_rule_passes(struct osh_scoring_filter_runtime_rule const *rule, struct particle const *part, double energy_mev) {
    double mass_amu;

    switch (rule->field) {
    case OSH_SCORING_FILTER_FIELD_ID:
        return filter_compare_int((unsigned int) part->pdg, rule->op, (unsigned int) rule->value);
    case OSH_SCORING_FILTER_FIELD_Z:
        return filter_compare_int((unsigned int) part->z, rule->op, (unsigned int) rule->value);
    case OSH_SCORING_FILTER_FIELD_A:
        return filter_compare_int((unsigned int) part->a, rule->op, (unsigned int) rule->value);
    case OSH_SCORING_FILTER_FIELD_AMASS:
        return filter_compare_double(part->mass, rule->op, rule->value);
    case OSH_SCORING_FILTER_FIELD_AMU:
        mass_amu = part->mass / OSH_AMU;
        return filter_compare_double(mass_amu, rule->op, rule->value);
    case OSH_SCORING_FILTER_FIELD_E:
        return filter_compare_double(energy_mev, rule->op, rule->value);
    case OSH_SCORING_FILTER_FIELD_ENUC:
        if (part->a != 0u) {
            return filter_compare_double(energy_mev / (double) part->a, rule->op, rule->value);
        }
        return filter_compare_double(energy_mev, rule->op, rule->value);
    case OSH_SCORING_FILTER_FIELD_EAMU:
        mass_amu = part->mass / OSH_AMU;
        if (mass_amu > 0.0) {
            return filter_compare_double(energy_mev / mass_amu, rule->op, rule->value);
        }
        return filter_compare_double(energy_mev, rule->op, rule->value);
    case OSH_SCORING_FILTER_FIELD_GEN:
        return filter_compare_int(part->gen, rule->op, (unsigned int) rule->value);
    case OSH_SCORING_FILTER_FIELD_NPRIM:
        return filter_compare_int(part->nprim, rule->op, (unsigned int) rule->value);
    default:
        return 0;
    }
}
