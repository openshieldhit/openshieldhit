#include "physics/nuclear/osh_nuclear_handler.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/material.h"
#include "particle/osh_isotope_db_generated.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_abrasion.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"
#include "physics/nuclear/osh_nuclear_pp.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"

/* p+p has no inelastic channel below the single-pion production threshold
 * (~280 MeV lab kinetic energy).  Without this gate, the Tripathi σ_R
 * evaluated on hydrogen would fire the ABSORB fallback and unphysically
 * swallow the full proton energy.  Above threshold the fallback remains in
 * effect until a pion-production model exists. */
#define H_INELASTIC_THRESHOLD_MEV 280.0

/** 1 if the inelastic channel must be skipped for this target element:
 *  hydrogen target, proton projectile, below the pion-production threshold. */
static inline int
_skip_hydrogen_inelastic(struct particle const *projectile, unsigned int target_z, double rate_energy_mev) {
    return target_z == 1u && projectile->pdg == OSH_PART_PDG_PROTON && rate_energy_mev < H_INELASTIC_THRESHOLD_MEV;
}

/* ---- Handler lifecycle --------------------------------------------------- */

enum osh_status osh_nuclear_handler_compile(struct osh_material_workspace const *ws, struct osh_nuclear_handler *out) {
    size_t total_elems;
    size_t i;
    size_t j;
    struct osh_nuclear_elem *ep;
    enum osh_status rc;

    if (!ws || !out) {
        return OSH_EINVAL;
    }

    /* Start from a clean object so failure paths can call free() safely. */
    memset(out, 0, sizeof(*out));
    out->nmaterials = ws->nmaterials;

    /* Per-material index arrays point into one flat element pool. */
    out->elem_offset = (size_t *) malloc(ws->nmaterials * sizeof(size_t));
    out->elem_count = (size_t *) malloc(ws->nmaterials * sizeof(size_t));
    if (!out->elem_offset || !out->elem_count) {
        osh_nuclear_handler_free(out);
        return OSH_ENOMEM;
    }

    /* Build prefix offsets, expanding natural (A=0) elements into per-isotope
     * entries so each slot in elem_pool carries a concrete (Z, A) nuclide and a
     * correctly weighted mass fraction. */
    total_elems = 0u;
    for (i = 0u; i < ws->nmaterials; ++i) {
        size_t expanded = 0u;
        for (j = 0u; j < ws->materials[i].nelements; ++j) {
            struct osh_material_element const *src = &ws->materials[i].elements[j];
            if (src->a == 0u) {
                unsigned int idx = osh_isotopes_idx[src->z];
                unsigned int len = osh_isotopes_len[src->z];
                unsigned int k;
                size_t n_nat = 0u;
                for (k = 0u; k < len; ++k) {
                    if (osh_isotope_db[idx + k].abund > 0.0) {
                        ++n_nat;
                    }
                }
                expanded += (n_nat > 0u) ? n_nat : 1u; /* 1u fallback: element has no natural isotope data */
            } else {
                ++expanded;
            }
        }
        out->elem_offset[i] = total_elems;
        out->elem_count[i] = expanded;
        total_elems += expanded;
    }

    if (total_elems > 0u) {
        out->elem_pool = (struct osh_nuclear_elem *) malloc(total_elems * sizeof(struct osh_nuclear_elem));
        if (!out->elem_pool) {
            osh_nuclear_handler_free(out);
            return OSH_ENOMEM;
        }
    }

    /* Fill the pool: A=0 → expand into isotopes; A>0 → copy directly. */
    ep = out->elem_pool;
    for (i = 0u; i < ws->nmaterials; ++i) {
        for (j = 0u; j < ws->materials[i].nelements; ++j) {
            struct osh_material_element const *src = &ws->materials[i].elements[j];
            if (src->a == 0u) {
                unsigned int z = src->z;
                unsigned int idx = osh_isotopes_idx[z];
                unsigned int len = osh_isotopes_len[z];
                unsigned int k;
                double m_nat = 0.0;
                for (k = 0u; k < len; ++k) {
                    m_nat += osh_isotope_db[idx + k].abund * (double) osh_isotope_db[idx + k].a;
                }
                if (m_nat <= 0.0) {
                    /* No natural abundance data: keep one entry as-is. */
                    ep->z = z;
                    ep->a = 0u;
                    ep->mass_fraction = (float) src->mass_fraction;
                    ++ep;
                } else {
                    /* One entry per naturally occurring isotope; mf_i = f × abund_i × a_i / M_nat */
                    for (k = 0u; k < len; ++k) {
                        const struct isotope *iso = &osh_isotope_db[idx + k];
                        if (iso->abund <= 0.0) {
                            continue;
                        }
                        ep->z = z;
                        ep->a = iso->a;
                        ep->mass_fraction = (float) (src->mass_fraction * iso->abund * (double) iso->a / m_nat);
                        ++ep;
                    }
                }
            } else {
                ep->z = src->z;
                ep->a = src->a;
                ep->mass_fraction = (float) src->mass_fraction;
                ++ep;
            }
        }
    }

    /* Fermi break-up channel table for prefragment de-excitation. */
    rc = osh_nuclear_fermi_breakup_compile(&out->fbu);
    if (rc != OSH_OK) {
        osh_nuclear_handler_free(out);
        return rc;
    }

    /* Persistent (Z,A) ion species table for recoil / fragment transport and
     * scoring.  Built once here (const during stepping): a heavy recoil or FBU
     * fragment needs a stable species pointer to be injected into the ion pool
     * or attributed in a point deposit. */
    {
        size_t const ndense = (size_t) (OSH_FERMI_BREAKUP_ZMAX + 1) * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1);
        unsigned int z;
        unsigned int a;

        out->recoil_species = (struct particle *) calloc(ndense, sizeof(struct particle));
        if (!out->recoil_species) {
            osh_nuclear_handler_free(out);
            return OSH_ENOMEM;
        }
        for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
            for (a = z; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
                size_t widx = (size_t) z * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1) + (size_t) a;
                int pdg = OSH_PART_PDG_HIBASE + (int) z * 10000 + (int) a * 10;
                struct particle p;
                /* osh_particle_from_pdg fills mass from the isotope DB; it fails
                 * for (Z,A) that are not real isotopes, which we leave zeroed
                 * (a==0) to mark absent. */
                if (osh_particle_from_pdg(&p, pdg)) {
                    out->recoil_species[widx] = p;
                }
            }
        }
    }

    return OSH_OK;
}

void osh_nuclear_handler_free(struct osh_nuclear_handler *h) {
    if (!h) {
        return;
    }
    osh_nuclear_fermi_breakup_free(&h->fbu);
    free(h->elem_pool);
    free(h->elem_offset);
    free(h->elem_count);
    free(h->recoil_species);
    memset(h, 0, sizeof(*h));
}

struct particle const *
osh_nuclear_handler_recoil_species(struct osh_nuclear_handler const *h, unsigned int z, unsigned int a) {
    size_t widx;

    if (!h || !h->recoil_species) {
        return NULL;
    }
    if (z == 0u || a == 0u || z > (unsigned int) OSH_FERMI_BREAKUP_ZMAX || a > (unsigned int) OSH_FERMI_BREAKUP_AMAX) {
        return NULL;
    }
    widx = (size_t) z * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1) + (size_t) a;
    if (h->recoil_species[widx].a == 0u) {
        return NULL; /* absent isotope */
    }
    return &h->recoil_species[widx];
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

    /* Projectile energy used by Tripathi tables: total KE converted to MeV/u. */
    double a_proj;
    double e_per_nucleon;

    /* Inelastic process hazard, summed over material elements in g^-1 cm^2. */
    double lambda_inel;
    double rate_inel;

    /* Concrete target nucleus selected after the inelastic channel fires. */
    double selected_a;
    double selected_z;
    double selected_sigma_inel;

    /* pp elastic hazard for hydrogen-containing materials. */
    double hydrogen_mf;
    double sigma_el;
    double lambda_pp;
    double rate_pp;

    /* Competing-process event sampling. */
    double rate_tot;
    double p_event;

    /* pp elastic final-state kinematics. */
    double cos_cm;
    double cos_phi;
    double sin_phi;
    double cos1;
    double e1;
    double cos2;
    double e2;
    double sin1;
    double sin2;

    event_out->kind = OSH_NUCLEAR_EVENT_NONE;
    event_out->n_secondaries = 0u;
    event_out->n_fragments = 0u;

    if (ds_gcm2 <= 0.0) {
        return;
    }
    if (!params->nuclear_inelastic && !params->nuclear_elastic) {
        return;
    }
    if (material_idx >= handler->nmaterials) {
        return;
    }

    nelem = handler->elem_count[material_idx];
    if (nelem == 0u) {
        return;
    }
    elems = handler->elem_pool + handler->elem_offset[material_idx];

    a_proj = (projectile->a > 0u) ? (double) projectile->a : 1.0;
    e_per_nucleon = rate_energy_mev / a_proj;

    /*
     * Inelastic rate: sum per-element hazards, then sample the struck element
     * proportional to its hazard when the inelastic channel wins. This keeps
     * compound materials physically meaningful for fragmentation: water is
     * struck on oxygen (abrasion + break-up), while hydrogen is excluded
     * below the pion-production threshold (no p+p inelastic channel there;
     * pp elastic is handled separately).
     */
    rate_inel = 0.0;
    if (params->nuclear_inelastic) {
        for (i = 0u; i < nelem; ++i) {
            double ai = (double) (elems[i].a > 0u ? elems[i].a : elems[i].z * 2u);
            double zi = (double) elems[i].z;
            double sigma_i;
            if (_skip_hydrogen_inelastic(projectile, elems[i].z, rate_energy_mev)) {
                continue;
            }
            sigma_i = osh_nuclear_tripathi_sigma(projectile->z, projectile->a, zi, ai, e_per_nucleon);
            if (sigma_i > 0.0) {
                lambda_inel = osh_nuclear_lambda_gcm2(ai, sigma_i);
                rate_inel += (double) elems[i].mass_fraction / lambda_inel;
            }
        }
    }

    /* pp elastic rate — only for proton projectile */
    rate_pp = 0.0;
    if (params->nuclear_elastic && projectile->pdg == OSH_PART_PDG_PROTON) {
        hydrogen_mf = 0.0;
        for (i = 0u; i < nelem; ++i) {
            if (elems[i].z == 1u) {
                hydrogen_mf += (double) elems[i].mass_fraction;
            }
        }
        if (hydrogen_mf > 0.0) {
            sigma_el = osh_nuclear_pp_sigma_el(rate_energy_mev);
            if (sigma_el > 0.0) {
                lambda_pp = osh_nuclear_pp_lambda_gcm2(hydrogen_mf, sigma_el);
                rate_pp = 1.0 / lambda_pp;
            }
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

        osh_kinematics_elastic_equal_mass_lab(final_energy_mev, OSH_PART_MASS_PROTON, cos_cm, &cos1, &e1, &cos2, &e2);

        /* Møller convention: primary keeps the higher-energy particle.
         * For backward CM scatter (cos_CM < 0) e2 > e1; swap so gen=0 stays
         * with the forward-going proton, not the low-energy sideways one. */
        if (e2 > e1) {
            double tmp;
            tmp = e1;
            e1 = e2;
            e2 = tmp;
            tmp = cos1;
            cos1 = cos2;
            cos2 = tmp;
            cos_phi = -cos_phi;
            sin_phi = -sin_phi;
        }

        sin1 = sqrt(fmax(0.0, 1.0 - (cos1 * cos1)));
        sin2 = sqrt(fmax(0.0, 1.0 - (cos2 * cos2)));

        osh_kinematics_rotate_dir_cos(incident_dir, event_out->primary_dir, cos1, sin1, cos_phi, sin_phi);
        osh_kinematics_rotate_dir_cos(incident_dir, event_out->secondaries[0].dir, cos2, sin2, -cos_phi, -sin_phi);

        event_out->kind = OSH_NUCLEAR_EVENT_ELASTIC_PP;
        event_out->primary_energy = e1;
        event_out->n_secondaries = 1u;
        event_out->secondaries[0].energy = e2;
        event_out->secondaries[0].species = projectile; /* proton = same species */
    } else {
        /* Select the actual struck element, not a compound-average nucleus. */
        double threshold;
        double cumulative;

        selected_a = 0.0;
        selected_z = 0.0;
        selected_sigma_inel = 0.0;
        threshold = osh_rng_double(rng) * rate_inel;
        cumulative = 0.0;

        for (i = 0u; i < nelem; ++i) {
            double ai = (double) (elems[i].a > 0u ? elems[i].a : elems[i].z * 2u);
            double zi = (double) elems[i].z;
            double sigma_i;
            /* Must mirror the rate_inel sum above so hazards stay consistent. */
            if (_skip_hydrogen_inelastic(projectile, elems[i].z, rate_energy_mev)) {
                continue;
            }
            sigma_i = osh_nuclear_tripathi_sigma(projectile->z, projectile->a, zi, ai, e_per_nucleon);
            if (sigma_i > 0.0) {
                double elem_rate = (double) elems[i].mass_fraction / osh_nuclear_lambda_gcm2(ai, sigma_i);
                selected_a = ai;
                selected_z = zi;
                selected_sigma_inel = sigma_i;
                cumulative += elem_rate;
                if (threshold <= cumulative) {
                    break;
                }
            }
        }

        if (selected_sigma_inel <= 0.0 || projectile->pdg != OSH_PART_PDG_PROTON || selected_a <= 1.5) {
            event_out->kind = OSH_NUCLEAR_EVENT_ABSORB;
            event_out->primary_energy = 0.0;
            event_out->n_secondaries = 0u;
            event_out->n_fragments = 0u;
            return;
        }

        /* Abrasion (fast stage) followed by Fermi break-up de-excitation of
         * the surviving prefragment (ablation stage).  The break-up consumes
         * the fragment slot when it fires and appends its products to the
         * same event; un-broken residues stay for the fragment pool. */
        osh_nuclear_abrasion_step(
            final_energy_mev, incident_dir, selected_a, selected_z, selected_sigma_inel, rng, event_out);
        if (event_out->n_fragments > 0u) {
            osh_nuclear_fermi_breakup_step(&handler->fbu, &event_out->fragments[0], rng, event_out);
        }
    }
}
