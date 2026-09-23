#ifndef OSH_VECT_H
#define OSH_VECT_H

#define OSH_VECT_DIM 3
#define OSH_VECT_EPS 1e-10

#include <stddef.h>

/**
 * @brief  Add p to q.
 *
 * @param[in] p - first point or vector
 * @param[in] q - second point or vector
 * @param[out] u - the sum <u> = <p> + <q>
 *
 * @author Niels Bassler
 */
void osh_vect_add(double const *p, double const *q, double *u);

/**
 * @brief  Add p to q * scalar d.
 *
 * @param[in] p - first point or vector
 * @param[in] q - second point or vector
 * @param[in] d - scalar which will be multiplied to q
 * @param[out] u - the sum <u> = <p> + <q> * d
 *
 * @author Niels Bassler
 */
void osh_vect_addmul(double const *p, double const *q, double d, double *u);

/**
 * @brief  Subtracts p from q.
 *
 * @param[in] p - first point or vector
 * @param[in] q - second point or vector
 * @param[out] u - the difference <u> = <p> - <q>
 *
 * @author Niels Bassler
 */
void osh_vect_sub(double const *p, double const *q, double *u);

/**
 * @brief  Deep copy of u into v.
 *
 * @param[in] u - first point or vector
 * @param[out] v - copy point or vector
 *
 * @author Niels Bassler
 */
void osh_vect_copy(double const *u, double *v);

/**
 * @brief  Returns a vector in the opposite direction.
 *
 * @param[in] u - first vector
 * @param[out] v = -<u>
 *
 * @author Niels Bassler
 */
void osh_vect_reverse(double const *u, double *v);

/**
 * @brief  Return the length^2 of a given vector
 *
 * @param[in] u - a vector
 *
 * @returns the length of the vector u squared, |u|^2.
 *
 * @author Niels Bassler
 */
double osh_vect_len2(double const *u);

/**
 * @brief  Return the dot product of two input vectors.
 *
 * @param[in] u - first vector
 * @param[out] v - second vector
 *
 * @returns <u> dot <v>
 *
 * @author Niels Bassler
 */
double osh_vect_dot(double const *u, double const *v);

/**
 * @brief Return the cross product of two input vectors.
 *
 * @param[in] u First input vector (3D)
 * @param[in] v Second input vector (3D)
 * @param[out] w Resulting cross product vector (u × v), must be allocated
 *
 * @author Niels Bassler
 */
void osh_vect_cross(double const *u, double const *v, double *w);

/**
 * @brief Scalar projection of vector u onto vector v.
 *
 * @param[in] u - first vector
 * @param[in] v - second vector where vector u will be projected on.
 *
 * @returns scalar projection of u onto v: <u>,<v>/|v|
 *
 * @author Niels Bassler
 */
double osh_vect_sproj(double const *u, double const *v);

/**
 * @brief Normalize a vector, updates u to hold <u>/|u|
 *
 * @param[in,out] u - vector which will be normalized
 *
 * @author Niels Bassler
 */
void osh_vect_norm(double *u);

/**
 * @brief Normalize a vector, return as a new vector.
 *
 * @param[in] u - vector to be normalized
 * @param[out] v - normalized vector
 *
 * @author Niels Bassler
 */
void osh_vect_norm2(double const *u, double *v);

/**
 * @brief For a given vector v, calculate two vectors which are perpendicular to
 * it and each other.
 *
 * @details Simple replacement for GVEC90 mode 2
 *  Strategy: take x,y,z axis unit vectors and dot each of them to the input
 * vector. The one with the lowest abs(dot product), must be "most orthogonal".
 *  This leads to more stable numerical situations for calculating the cross
 * product later on.
 *
 * (u,v,w) is right-handed, but u and v are not unit vectors -- the caller must
 * osh_vect_norm() them if it needs an orthonormal triad, or call
 * osh_vect_orthonormal_basis(), which returns one directly.
 *
 * @param[in]  w - input unit vector
 * @param[out] u - vector orthogonal to <w> and <v>
 * @param[out] v - vector orthogonal to <w> and <u>
 *
 * @author Niels Bassler, Leszek Grzanka
 */
void osh_vect_orthogonal_basis(double const *w, double *u, double *v);

/**
 * @brief Like osh_vect_orthogonal_basis(), but the outputs are unit vectors,
 * so the caller needs no follow-up osh_vect_norm().
 *
 * @details Projects out the least-aligned Cartesian axis (Gram-Schmidt) rather
 * than crossing with e3, which is numerically better conditioned and has no
 * degenerate case to special-case.  (u,v,w) is right-handed for every @p w.
 * Preferred over osh_vect_orthogonal_basis() for new code.
 *
 * TODO: migrate gemca callers of osh_vect_orthogonal_basis to this function.
 *
 * @param[in]  w  Input unit direction.
 * @param[out] u  First transverse unit vector, perpendicular to @p w.
 * @param[out] v  Second transverse unit vector, perpendicular to both.
 */
void osh_vect_orthonormal_basis(double const *w, double *u, double *v);

/**
 * @brief
 * Computes the coefficients (A, B, C, D) of the plane equation Ax + By + Cz + D
 * = 0, given a point `p` in the plane and a normal vector `u` orthogonal to it.
 *
 * The output is not normalized: the vector (A, B, C) has the same direction and
 * magnitude as `u`, and D is calculated as -dot(u, p).
 *
 * @param[in]  p     A point lying in the plane (3D vector).
 * @param[in]  u     A vector orthogonal to the plane (not required to be
 * normalized).
 * @param[out] pp[4] Output array for the plane coefficients (A, B, C, D).
 *
 * @author Niels Bassler
 */
void osh_vect_eqpln(double const *p, double const *u, double *pp);

/**
 * @brief For a given vector u, rotate it clockwise by a given angle around OY
 * axis. u is updated.
 *
 * @param[in] alpha - rotation angle in [radians]
 * @param[in,out] u - input vector.
 *
 * @author Niels Bassler, Leszek Grzanka
 */
void osh_vect_rot_y(double alpha, double *u);

/**
 * @brief For a given vector u, rotate it clockwise by a given angle around OZ
 * axis. u is updated.
 *
 * @param[in] phi - rotation angle in [radians]
 * @param[in,out] u - input vector.
 *
 * @author Niels Bassler, Leszek Grzanka
 */
void osh_vect_rot_z(double alpha, double *u);

/**
 * @brief Prints the contents of a 3D vector to stdout.
 *
 * @author Niels Bassler
 */
void osh_vect_print(double const *v);

/**
 * @brief Prints the contents of the 4x4 transformation matrix
 *
 * @author Niels Bassler
 */
void osh_vect_matrix4_print(double const *tm);

/**
 * @brief Builds the UNIVERSE --> frame transformation matrix, where "frame" is
 * the rotated local system of a body or beam -- OSH_COORD_BZALIGN and
 * OSH_COORD_PZALIGN, whose consumers apply the full 3x4 affine.
 *
 * Not for OSH_COORD_BCALIGN: that tag means a pure translation, and its
 * consumers add t[3], t[7], t[11] while ignoring the rotation rows entirely
 * (see the BCALIGN branches in osh_gemca_runtime.c).  A matrix built here
 * would have its rotation silently dropped.
 *
 * @details @p origin_universe maps to (0,0,0) and @p axis maps onto local +z.
 *
 * With M = [S|T|R] the local basis expressed in UNIVERSE, this stores M^T:
 * row i is the basis vector itself, and tm[i*4+3] = -<basis_i,P>, so that the
 * standard affine rule p_out = M_rows * p_in + t evaluates
 *   p_frame = M^T (p_universe - P)
 *
 * The basis (S,T,R) is always right-handed.  The inverse direction is built by
 * osh_vect_tmatrix_frame_to_universe(); note that the two take their origin in
 * different systems -- each takes it in the system it transforms *from*, so
 * that they compose as exact inverses when P_universe = M * P_frame.
 *
 * @param[in] origin_universe[3] - frame origin, in UNIVERSE coordinates
 * @param[in] axis[3] - vector for calculating the rotation matrix (not the
 * rotation axis!); need not be normalised
 * @param[out] tm[16] - transformation matrix, which maps axis onto e3 (0,0,1)
 *
 * @author Niels Bassler
 */
void osh_vect_tmatrix_universe_to_frame(double const *origin_universe, double const *axis, double *tm);

/**
 * @brief Builds the frame --> UNIVERSE transformation matrix, i.e. the opposite
 * direction of osh_vect_tmatrix_universe_to_frame() above.
 *
 * @details The columns of M are the local basis vectors (S,T,R) expressed in
 * UNIVERSE, so the standard affine rule p_out = M * p_in + t evaluates
 *   p_universe = M * (p_frame + origin_frame)
 *
 * Following the convention above, the origin is given in the system this
 * matrix transforms *from* -- here the local frame, not UNIVERSE.  The stored
 * translation is therefore derived as t = M * origin_frame.  Callers holding a
 * UNIVERSE-space origin must convert it first, or use
 * osh_vect_tmatrix_universe_to_frame() for the other direction.
 *
 * @param[in] origin_frame[3] - frame origin, in local frame coordinates
 * @param[in] axis_universe[3] - local +z direction, in UNIVERSE coordinates
 * @param[out] tm[16] - transformation matrix, which maps e3 (0,0,1) onto axis
 */
void osh_vect_tmatrix_frame_to_universe(double const *origin_frame, double const *axis_universe, double *tm);

/**
 * @brief Apply a standard affine matrix to a point.
 */
void osh_vect_trans_point(double const *p, double *pt, double const *tm);

/**
 * @brief Apply only the rotation part of a standard affine matrix to a vector.
 */
void osh_vect_trans_vector(double const *v, double *vt, double const *tm);

#endif /* OSH_VECT_H */
