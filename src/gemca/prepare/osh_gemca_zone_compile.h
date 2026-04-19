#ifndef OSH_GEMCA_ZONE_COMPILE_H
#define OSH_GEMCA_ZONE_COMPILE_H

#include "gemca/osh_gemca2.h"
#include "openshieldhit/logger.h"

/**
 * @file osh_gemca_zone_compile.h
 * @brief Internal GEMCA zone-compilation helpers used by the prepare layer.
 *
 * @details
 * These helpers compile the raw boolean expression stored in the public cold
 * @ref osh_geometry_zone into the internal pointer-linked CSG tree used by
 * the compatibility workspace @ref osh_gemca_prepared.
 */

enum osh_status osh_gemca_zone_init(struct zone **zone);

/**
 * @brief Compile a zone boolean expression into a cgnode AST.
 *
 * @details Runs the _reformat → _tokenizer → _reverse_tokens → _build_ast
 *          pipeline on @p expr and stores the result in @p z.  Bodies must
 *          already be set up in @p g (osh_gemca_body_setup() called) so that
 *          name → body pointer resolution works.
 *
 * @param[in,out] z     Zone to populate (name and material_name already set).
 * @param[in]     expr  Raw zone expression string (e.g. "+BODY1 -BODY2").
 * @param[in]     g     Gemca workspace with all bodies already initialised.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_gemca_zone_compile_expr(struct zone *z,
                                            char const *expr,
                                            struct osh_gemca_prepared *g,
                                            struct osh_diag_sink const *diag);

#endif /* OSH_GEMCA_ZONE_COMPILE_H */
