#ifndef OSH_GEOMETRY_PREPARED_H
#define OSH_GEOMETRY_PREPARED_H

/**
 * @file osh_geometry_prepared.h
 * @brief Private prepared state for a cold osh_geometry_workspace.
 *
 * @details
 * This type is forward-declared in the public geometry.h and stored via
 * pointer in osh_geometry_workspace.prepared so that callers never see the
 * internal GEMCA types.
 *
 * During migration the prepared state is implemented as a thin wrapper
 * around the existing struct gemca_workspace.  Once the analytic path has
 * been fully migrated away from gemca_workspace, this struct can be
 * replaced with a direct reference to the compiled runtime representation.
 *
 * Ownership: allocated by osh_geometry_workspace_prepare() and freed by
 * osh_geometry_workspace_free().  No caller outside src/gemca/ should
 * allocate or free this struct directly.
 */

struct gemca_workspace; /* defined in src/gemca/osh_gemca2.h */

struct osh_geometry_prepared {
    struct gemca_workspace *gemca; /**< Internal compiled GEMCA workspace (compatibility layer). */
};

#endif /* OSH_GEOMETRY_PREPARED_H */
