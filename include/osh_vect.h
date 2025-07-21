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
 * @brief For a given vector v, calculate two vectors which are perpendicular to it and each other.
 *
 * @details Simple replacement for GVEC90 mode 2
 *  Strategy: take x,y,z axis unit vectors and dot each of them to the input vector.
 *  The one with the lowest abs(dot product), must be "most orthogonal".
 *  This leads to more stable numerical situations for calculating the cross product later on.
 *
 * @param[in] u - input unit vector
 * @param[out] v - unit vector orthogonal to <u> and <w>
 * @param[out] w - unit vector orthogonal to <u> and <v>
 *
 * @author Niels Bassler, Leszek Grzanka
 */
void osh_vect_orthogonal_basis(double const *w, double *u, double *v);

/**
 * @brief
 * Computes the coefficients (A, B, C, D) of the plane equation Ax + By + Cz + D = 0,
 * given a point `p` in the plane and a normal vector `u` orthogonal to it.
 *
 * The output is not normalized: the vector (A, B, C) has the same direction and
 * magnitude as `u`, and D is calculated as -dot(u, p).
 *
 * @param[in]  p     A point lying in the plane (3D vector).
 * @param[in]  u     A vector orthogonal to the plane (not required to be normalized).
 * @param[out] pp[4] Output array for the plane coefficients (A, B, C, D).
 *
 * @author Niels Bassler
 */
void osh_vect_eqpln(double const *p, double const *u, double *pp);


/**
 * @brief For a given vector u, rotate it clockwise by a given angle around OY axis. u is updated.
 *
 * @param[in] alpha - rotation angle in [radians]
 * @param[in,out] u - input vector.
 *
 * @author Niels Bassler, Leszek Grzanka
 */
void osh_vect_rot_y(double alpha, double *u);

/**
 * @brief For a given vector u, rotate it clockwise by a given angle around OZ axis. u is updated.
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
 * @brief Builds the transformation matrix osh_COORD_BZALIGN --> osh_COORD_UNIVERSE system.
 *
 * @details p will be at (0,0,0) and r will be along Z in the osh_COORD_BZALIGN system.
 *
 * @param[in] p[3] - translation vector
 * @param[in] r[3] - vector for calculating the rotation matrix (not the rotation axis!)
 * @param[out] tm[16] - transformation matrix, which maps e3 (0,0,1) into R
 *
 * @returns
 *
 * @author Niels Bassler
 */
void osh_vect_setup_tmatrix_bzalign(double *p, double *r, double *tm);

#endif /* OSH_VECT_H */

