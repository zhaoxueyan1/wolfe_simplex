//
// Created by zhaox on 2021/10/4.
//
#include "wolfe_simplex.h"
namespace SHG {

int wolfe(const Vecdouble& p, const Vecdouble& C, const Matdouble& A,
          const Vecdouble& b, Vecdouble& x, double& f) {
    static const double eps = 1e-11;
    static const int maxiter = 10000;
    const int m = A.nrows();
    const int n = A.ncols();
    if (m <= 0 || n <= 0)
        return 1;  // invalid parameter m or n
    {
        const size_t ms = m;
        const size_t ns = n;
        if (p.size() != ns || C.size() != ns * (ns + 1) / 2 ||
            A.nrows() != ms || A.ncols() != ns || b.size() != ms ||
            x.size() != ns)
            return 1;  // invalid dimensions
    }
    const int g = m + n;      // added row in simplex table
    const int h = 2 * n + m;  // added column in simplex table
    const int                 // intervals of variable numbers
    // x:  0 .. nv - 1
    nv = n,             // v: nv .. nu - 1
    nu = 2 * n,         // u: nu .. nz - 1
    nz = 2 * n + m,     // z: nz .. nw - 1
    nw = 3 * n + m;     // w: nw .. 3 * n + 2 * m - 1
    Matdouble t;             // simplex table
    Vecint B, N;             // basis and non-basis variables
    int q = m;               // last row to leave the basis
    bool initiation = true;  // true = initiation, false = recursion
    int iter = 0;
    int i, j, k, l = 0;  // Dummy initialization for compiler.
    double d, s;

    try {
        t.resize(g + 1, h + 1);
        B.resize(g);
        N.resize(h);
    } catch (const std::bad_alloc&) {
        return 2;  // not enough memory
    }

    // Initialize basis and non-basis variable numbers.
    for (i = 0; i < g; i++)
        B[i] = i < m ? nw + i : nz + (i - m);
    for (i = 0; i < h; i++)
        N[i] = i;

    // Initialize first m rows of the simplex table.
    for (i = 0; i < m; i++) {
        if (b[i] < 0.0) {
            for (j = 0; j < n; j++)
                t[i][j] = -A[i][j];
            t[i][h] = -b[i];
        } else {
            for (j = 0; j < n; j++)
                t[i][j] = A[i][j];
            t[i][h] = b[i];
        }
        for (j = n; j < h; j++)
            t[i][j] = 0.0;
    }

    // Initialize next n rows of the simplex table.
    k = 0;
    for (i = m; i < g; i++) {
        t[i][i - m] = 2.0 * C[k++];
        for (j = i - m + 1; j < n; j++)
            t[i][j] = t[m + j][i - m] = 2.0 * C[k++];
        for (j = n; j < 2 * n; j++)
            t[i][j] = 0.0;
        t[i][n + (i - m)] = -1.0;
        for (j = 2 * n; j < h; j++)
            t[i][j] = A[j - 2 * n][i - m];
        t[i][h] = -p[i - m];
        if (t[i][h] < 0.0) {
            for (j = 0; j <= h; j++)
                t[i][j] = -t[i][j];
        }
    }

    // Initialize the g-th row of the simplex table.
    for (j = 0; j < n; j++) {
        s = 0.0;
        for (i = 0; i < m; i++)
            s += t[i][j];
        t[g][j] = s;
    }
    for (j = n; j < h; j++)
        t[g][j] = 0.0;
    s = 0.0;
    for (i = 0; i < m; i++)
        s += t[i][h];
    t[g][h] = s;

    for (;;) {  // for initiation and recursion
        if (++iter > maxiter)
            return 3;  // too many iteration steps

        // Find variable to enter the basis.
        s = 0.0;
        if (initiation) {
            for (j = 0; j < h; j++)
                // exclude u, v, z
                if (N[j] < nv || N[j] >= nw)
                    if (t[g][j] > s)
                        s = t[g][l = j];
        } else {
            for (j = 0; j < h; j++)
                if (N[j] < nw) {
                    // exclude w
                    if (N[j] >= nv && N[j] < nu) {
                        // side condition for v
                        k = N[j] - n;
                        for (i = 0; i < g && B[i] != k; i++)
                            ;
                        if (i < g)
                            continue;
                    } else if (N[j] < nv) {
                        // side condition for x
                        k = N[j] + n;
                        for (i = 0; i < g && B[i] != k; i++)
                            ;
                        if (i < g)
                            continue;
                    } else if (N[j] >= nu && N[j] < nz &&
                               t[g][j] < 0.0) {
                        // adjust u, as we have u1 and u2
                        for (i = 0; i <= g; i++)
                            t[i][j] = -t[i][j];
                    }
                    if (t[g][j] > s)
                        s = t[g][l = j];
                }
        }

        // If s < eps, there are three possibilities:
        // - no feasible solution,
        // - the method finishes successfully in recursion phase,
        // - the phase is changed from initiation to recursion.

        if (s < eps) {
            if (t[g][h] > eps)
                // no feasible solution
                return initiation ? 4 : 5;
            if (!initiation)
                break;  // finished successfully

            // Change the phase.
            initiation = false;
            q = g;

            // Discard unused components of z.
            for (i = m; i < g; i++)
                if (t[i][h] < 0.0) {
                    for (j = 0; j <= h; j++)
                        t[i][j] = -t[i][j];
                }

            // Calculate last row for the new objective function.
            for (j = 0; j <= h; j++) {
                s = 0.0;
                for (i = 0; i < g; i++)
                    if (B[i] >= nz && B[i] < nw)
                        s += t[i][j];
                if (j < h && N[j] >= nz && N[j] < nw)
                    s -= 1.0;
                t[g][j] = s;
            }
            continue;
        }

        // Find variable to leave the basis.
        k = -1;
        for (i = 0; i < q; i++)
            // discard w in recursion
            if ((initiation || B[i] < nw) && t[i][l] >= eps) {
                d = t[i][h] / t[i][l];
                if (k < 0 || d < s) {
                    k = i;
                    s = d;
                }
            }

        if (k < 0)
            // unbounded solution
            return initiation ? 6 : 7;

        // Change variables.
        i = B[k];
        B[k] = N[l];
        N[l] = i;

        // Transform the simplex table around t[k][l].
        s = t[k][l];
        t[k][l] = 1.0;
        for (j = 0; j <= h; j++)
            t[k][j] /= s;
        for (i = 0; i <= g; i++)
            if (i != k) {
                d = t[i][l];
                t[i][l] = 0.0;
                for (j = 0; j <= h; j++)
                    t[i][j] -= t[k][j] * d;
            }
    }  // for initiation and recursion

    // Put the solution.
    for (i = 0; i < n; i++)
        x[i] = 0.0;

    for (i = 0; i < g; i++)
        if (B[i] < n)
            x[B[i]] = t[i][h];

    // Calculate minimum value.
    f = 0.0;
    k = 0;
    for (i = 0; i < n; i++) {
        s = 0.0;
        d = C[k++];
        for (j = i + 1; j < n; j++)
            s += C[k++] * x[j];
        f += x[i] * (2.0 * s + d * x[i] + p[i]);
    }
    return 0;
}
}