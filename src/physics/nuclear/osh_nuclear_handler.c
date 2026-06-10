#include "physics/nuclear/osh_nuclear_handler.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/material.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_pp.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"

/* ---- Handler lifecycle --------------------------------------------------- */

enum osh_status osh_nuclear_handler_compile(struct osh_material_workspace const *ws,
                                            struct osh_nuclear_handler *out) {
    size_t total_elems;
    size_t i;
    size_t j;
    struct osh_nuclear_elem *ep;

    if (!ws || !out) {
        return OSH_EINVAL;
    }

    memset(out, 0, sizeof(*out));
    out->nmaterials = ws->nmaterials;

    out->elem_offset = (size_t *) malloc(ws->nmaterials * sizeof(size_t));
    out->elem_count  = (size_t *) malloc(ws->nmaterials * sizeof(size_t));
    if (!out->elem_offset || !out->elem_count) {
        osh_nuclear_handler_free(out);
        return OSH_ENOMEM;
    }

    total_elems = 0u;
    for (i = 0u; i < ws->nmaterials; ++i) {
        out->elem_offset[i] = total_elems;
        out->elem_count[i]  = ws->materials[i].nelements;
        total_elems += ws->materials[i].nelements;
    }

    if (total_elems > 0u) {
        out->elem_pool = (struct osh_nuclear_elem *) malloc(total_elems * sizeof(struct osh_nuclear_elem));
        if (!out->elem_pool) {
            osh_nuclear_handler_free(out);
            return OSH_ENOMEM;
        }
    }

    ep = out->elem_pool;
    for (i = 0u; i < ws->nmaterials; ++i) {
        for (j = 0u; j < ws->materials[i].nelements; ++j) {
            struct osh_material_element const *src = &ws->materials[i].elements[j];
            ep->z            = src->z;
            ep->a            = src->a;
            ep->mass_fraction = (float) src->mass_fraction;
            ++ep;
        }
    }

    return OSH_OK;
}

void osh_nuclear_handler_free(struct osh_nuclear_handler *h) {
    if (!h) {
        return;
    }
    free(h->elem_pool);
    free(h->elem_offset);
    free(h->elem_count);
    memset(h, 0, sizeof(*h));
}

/* ---- Handler step -------------------------------------------------------- */

void osh_nuclear_handler_step(struct osh_nuclear_handler const *handler,
                              double rate_energy_mev,
                              double final_energy_mev,
                              double const incident_dir[3],
                              size_t material_idx,
                              double ds_gcm2,
                              struct particle const *projectile,
                              struct osh_transport_params const *params,
                              struct osh_rng *rng,
                              struct osh_nuclear_event *event_out) {
    struct osh_nuclear_elem const *elems;
    size_t nelem;
    size_t i;
    double at;
    double e_per_nucleon;
    double sigma_inel;
    double lambda_inel;
    double rate_inel;
    double hydrogen_mf;
    double sigma_el;
    double lambda_pp;
    double rate_pp;
    double rate_tot;
    double p_event;
    double cos_cm;
    double cos_phi, sin_phi;
    double cos1, e1, cos2, e2;
    double sin1, sin2;
    double a_proj;

    event_out->kind        = OSH_NUCLEAR_EVENT_NONE;
    event_out->n_secondaries = 0u;

    if (ds_gcm2 <= 0.0) {
        return;
    }
    if (!params->nuclear_inelastic && !params->nuclear_elastic) {
        return;
    }
    if (material_idx >= handler->nmaterials) {
        return;
    }

    elems = handler->elem_pool + handler->elem_offset[material_idx];
    nelem = handler->elem_count[material_idx];

    /*
     * Tripathi inelastic rate at rate_energy_mev.
     * at = effective target mass number (= z_mean / z_over_a).
     * Compute effective Z and A from the element list.
     */
    {
        double z_sum = 0.0, a_sum = 0.0, wsum = 0.0;
        for (i = 0u; i < nelem; ++i) {
            z_sum += (double) elems[i].z * (double) elems[i].mass_fraction;
            a_sum += (double) (elems[i].a > 0u ? elems[i].a : elems[i].z * 2u) * (double) elems[i].mass_fraction;
            wsum  += (double) elems[i].mass_fraction;
        }
        at = (wsum > 0.0) ? (a_sum / wsum) : 1.0;
        {
            double zt = (wsum > 0.0) ? (z_sum / wsum) : 1.0;
            a_proj = (projectile->a > 0u) ? (double) projectile->a : 1.0;
            e_per_nucleon = rate_energy_mev / a_proj;
            sigma_inel = params->nuclear_inelastic
                ? osh_nuclear_tripathi_sigma(projectile->z, projectile->a, zt, at, e_per_nucleon)
                : 0.0;
        }
    }
    lambda_inel = (sigma_inel > 0.0) ? osh_nuclear_lambda_gcm2(at, sigma_inel) : 1.0e30;
    rate_inel   = 1.0 / lambda_inel;

    /* pp elastic rate — only for proton projectile */
    rate_pp = 0.0;
    lambda_pp = 1.0e30;
    if (params->nuclear_elastic && projectile->pdg == OSH_PART_PDG_PROTON) {
        hydrogen_mf = 0.0;
        for (i = 0u; i < nelem; ++i) {
            if (elems[i].z == 1u) {
                hydrogen_mf += (double) elems[i].mass_fraction;
            }
        }
        if (hydrogen_mf > 0.0) {
            sigma_el  = osh_nuclear_pp_sigma_el(rate_energy_mev);
            lambda_pp = osh_nuclear_pp_lambda_gcm2(hydrogen_mf, sigma_el);
            rate_pp   = 1.0 / lambda_pp;
        }
    }

    rate_tot = rate_inel + rate_pp;
    if (rate_tot <= 0.0) {
        return;
    }

    p_event = 1.0 - exp(-ds_gcm2 * rate_tot);
    if (osh_rng_double(rng) >= p_event) {
        return;
    }

    /* Event fired: select channel by exact rate ratio */
    if (rate_pp > 0.0 && osh_rng_double(rng) < rate_pp / rate_tot) {
        /* pp elastic — use final_energy_mev for kinematics */
        cos_cm = osh_nuclear_pp_sample_cos_theta_cm(rate_energy_mev, rng);
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);

        osh_kinematics_elastic_equal_mass_lab(final_energy_mev, OSH_PART_MASS_PROTON, cos_cm,
                                             &cos1, &e1, &cos2, &e2);

        sin1 = sqrt(1.0 - cos1 * cos1);
        sin2 = sqrt(1.0 - cos2 * cos2);

        osh_kinematics_rotate_dir_cos(incident_dir, event_out->primary_dir,
                                      cos1, sin1, cos_phi, sin_phi);
        osh_kinematics_rotate_dir_cos(incident_dir, event_out->secondaries[0].dir,
                                      cos2, sin2, -cos_phi, -sin_phi);

        event_out->kind           = OSH_NUCLEAR_EVENT_ELASTIC_PP;
        event_out->primary_energy = e1;
        event_out->n_secondaries  = 1u;
        event_out->secondaries[0].energy  = e2;
        event_out->secondaries[0].species = projectile; /* proton = same species */
    } else {
        /* Inelastic absorption */
        event_out->kind           = OSH_NUCLEAR_EVENT_ABSORB;
        event_out->primary_energy = 0.0;
    }
}
