//
// Created by zhaox on 2021/10/4.
//

#ifndef OPT_WOLFE_SIMPLEX_H
#define OPT_WOLFE_SIMPLEX_H
/**
 * \file include/shg/mathprog.h
 * Mathematical programming.
 * \date Created on 16 July 2014.
 */
#include <except.h>
#include <matrix.h>

namespace SHG {

/**
 * Wolfe method for quadratic programming problem.
 *
 * Minimizes \f$px + x^TCx\f$ subject to \f$Ax = b\f$, \f$x \geq 0\f$.
 * \f$p\f$ and \f$x\f$ are \f$n \times 1\f$ vectors, \f$b\f$ is an
 * \f$m \times 1\f$ vector, \f$A\f$ is an \f$m \times n\f$ matrix,
 * \f$C\f$ is an \f$n \times n\f$ matrix. The matrix \f$C\f$ must be
 * symmetric (only upper right triangle is used) and should be
 * positive definite. The number of linear equality constraints
 * \f$m\f$ and the number of variables \f$n\f$ must be greater than 0.
 *
 * \retval 0 ok, solution found
 * \retval 1 invalid arguments
 * \retval 2 not enough memory
 * \retval 3 the maximum number of simplex iteration steps exceeded
 * \retval 4 no feasible solution in initiation phase
 * \retval 5 no feasible solution in recursion phase
 * \retval 6 unbounded solution in initiation phase
 * \retval 7 unbounded solution in recursion phase
 *
 * If the returned value is 0, \a x contains the solution, \a f
 * contains the minimum value, otherwise they are udefined.
 *
 * \note The parameter \a C should be a vector of size n x (n + 1) / 2
 * which contains the upper right triangle of the matrix \f$C\f$
 * written by rows.
 *
 * \date Written on 10 May 2009.
 *
 * \internal
 * Implementation.
 * ===============
 *
 * Problem.
 * --------
 *
 * Let A be an m x n matrix, b >= 0 an m x 1 vector, p an 1 x n
 * vector, C an n x n symmetric, positive semidefinite matrix, and x
 * an n x 1 vector. The quadratic programming problem is:
 *
 *     minimize px + x^(T)Cx
 *     subject to Ax = b,
 *     x >= 0.
 *
 * The Wolfe method is based on the theorem that x >= 0 solves the
 * quadratic programming problem if Ax = b and there exist n x 1
 * vector v >= 0 and m x 1 vector u such that
 *
 *     v^(T)x = 0,
 *     2Cx - v + A^(T)u + p^(T) = 0.
 *
 * Convergence of the "short form" of the method, used here, requires
 * that either p = 0 or C is positive definite.
 *
 * Algorithm.
 * ----------
 *
 * Let z1, z2 be n x 1 vectors and w be an m x 1 vector.
 *
 * 1. Initiation. Use the simplex method to minimize w(1) + ... + w(m)
 * to zero subject to
 *
 *     Ax + w = b,
 *     2Cx - v + A^(T)u + z1 - z2 = -p^(T),
 *     x, v, z1, z2, w >= 0,
 * keeping u and v zero. Discard w and the unused components of z1 and
 * z2. Let the remaining n components be denoted by z, and their
 * coefficients by E. We have a solution of the system
 *
 *     Ax = b,
 *     2Cx - v + A^(T)u + Ez = -p^(T),
 *     x, v, z >= 0.
 *
 * 2. Recursion. Take u = u1 - u2, u1, u2 >= 0. Use the simplex method
 * to minimize z(1) + z(2) + ... + z(n) to zero subject to
 *
 *     Ax = b,
 *     2Cx - v + A^(T)u1 - A^(T)u2 + Ez = -p^(T),
 *     x, v, u1, u2, z >= 0,
 * under the side condition: for k = 1, ..., n: if x(k) is in the
 * basis, do not admit v(k); if v(k) is in the basis, do not admit
 * x(k).
 *
 * 3. Termination. If p = 0 or C is positive definite, the objective
 * function z(1) + ... + z(n) will vanish in at most 3n over n (Newton
 * symbol) iterations, yielding z = 0. The x-part of the terminal
 * basic solution is a solution of the quadratic programming problem.
 *
 * Implementation.
 * ---------------
 *
 * The initial simplex table for the initiation phase is
 *
 * |Basis| x | v | u | z | w |RHS |
 * |:---:|:-:|:-:|:-:|:-:|:-:|:--:|
 * |  w  | A | 0 | 0 | 0 | I | b  |
 * |  z  | 2C| -I|A^T| I | 0 |-p^T|
 * |     |   | 0 | 0 | 0 | 0 |    |
 *
 * The n variables z are choosen to have -p^T >= 0. The same table may
 * be used in the recursion phase. The numberes of basis variables
 * will be kept in a special vector, so there is no need to keep their
 * coefficients in the table. The table in both initiation and
 * recursion will be
 *
 * |Basis| x | v | u |RHS |
 * |:---:|:-:|:-:|:-:|:-: |
 * |  w  | A | 0 | 0 | b  |
 * |  z  | 2C| -I|A^T|-p^T|
 * |     |   | 0 | 0 |    |
 *
 * In the initiation phase u and v will not enter the basis and in the
 * recursion phase w will not enter the basis and the side condition
 * will be used. Numbers of variables will be as follows:
 *
 *     x(1) ... x(n) - 0 ... n - 1
 *     v(1) ... v(n) - n ... 2n - 1
 *     u(1) ... u(m) - 2n ... 2n + m - 1
 *     z(1) ... z(n) - 2n + m ... 3n + m - 1
 *     w(1) ... w(m) - 3n + m ... 3n + 2m - 1
 *
 * The simplex table will be kept in the matrix
 *
 *     t(m + n + 1, 2n + m + 1),
 *
 * the m + n-th row and 2n + m-th column being added for coefficients.
 * The numbers of basis and non-basis variables will be kept in the
 * integer vectors B(m + n) and N(2n + m).
 *
 * Both phases will be performed within the same loop. The flow is
 * arranged to affect whole rows of the table t in order to make
 * eventual version using scratch file easier to write.
 *
 * References.
 * -----------
 *
 * \cite wolfe-1959, \cite gass-1980, p. 276-282,
 * \cite grabowski-1980, p. 245-251,
 * \cite bronsztejn-siemiendiajew-musiol-muhlig-2004, p. 943-945
 */
    int wolfe(const Vecdouble& p, const Vecdouble& C, const Matdouble& A,
              const Vecdouble& b, Vecdouble& x, double& f);

/** \} */ /* end of group mathematical_programming */

}  // namespace SHG

#endif //OPT_WOLFE_SIMPLEX_H
