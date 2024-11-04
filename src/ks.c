/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2012   The R Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* ks.c
   Compute the asymptotic distribution of the one- and two-sample
   two-sided Kolmogorov-Smirnov statistics, and the exact distributions
   in the two-sided one-sample and two-sample cases.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

static double K(int n, double d);
static void m_multiply(double *A, double *B, double *C, int m);
static void m_power(double *A, int eA, double *V, int *eV, int m, int n);

static void
pkstwo(int n, double *x, double tol)
{
    double new, old, s, w, z;
    int i, k, k_max;

    k_max = (int) sqrt(2 - log(tol));

    for(i = 0; i < n; i++) {
	if(x[i] < 1) {
	    z = - (M_PI_2 * M_PI_4) / (x[i] * x[i]);
	    w = log(x[i]);
	    s = 0;
	    for(k = 1; k < k_max; k += 2) {
		s += exp(k * k * z - w);
	    }
	    x[i] = s / M_1_SQRT_2PI;
	}
	else {
	    z = -2 * x[i] * x[i];
	    s = -1;
	    k = 1;
	    old = 0;
	    new = 1;
	    while(fabs(old - new) > tol) {
		old = new;
		new += 2 * s * exp(z * k * k);
		s *= -1;
		k++;
	    }
	    x[i] = new;
	}
    }
}

static double psmirnov2x(double *x, int m, int n)
{
    double md, nd, q, *u, w;
    int i, j;

    if(m > n) {
	i = n; n = m; m = i;
    }
    md = (double) m;
    nd = (double) n;
    q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
    u = (double *) R_alloc(n + 1, sizeof(double));

    for(j = 0; j <= n; j++) {
	u[j] = ((j / nd) > q) ? 0 : 1;
    }
    for(i = 1; i <= m; i++) {
	w = (double)(i) / ((double)(i + n));
	if((i / md) > q)
	    u[0] = 0;
	else
	    u[0] = w * u[0];
	for(j = 1; j <= n; j++) {
	    if(fabs(i / md - j / nd) > q) 
		u[j] = 0;
	    else
		u[j] = w * u[j] + u[j - 1];
	}
    }
    return u[n];
}

static double
K(int n, double d)
{
    int k, m, i, j, g, eH, eQ;
    double h, s, *H, *Q;

    k = (int) (n * d) + 1;
    m = 2 * k - 1;
    h = k - n * d;
    H = R_Calloc(m * m, double);
    Q = R_Calloc(m * m, double);
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++)
	    if(i - j + 1 < 0)
		H[i * m + j] = 0;
	    else
		H[i * m + j] = 1;
    for(i = 0; i < m; i++) {
	H[i * m] -= pow(h, i + 1);
	H[(m - 1) * m + i] -= pow(h, (m - i));
    }
    H[(m - 1) * m] += ((2 * h - 1 > 0) ? pow(2 * h - 1, m) : 0);
    for(i = 0; i < m; i++)
	for(j=0; j < m; j++)
	    if(i - j + 1 > 0)
		for(g = 1; g <= i - j + 1; g++)
		    H[i * m + j] /= g;
    eH = 0;
    m_power(H, eH, Q, &eQ, m, n);
    s = Q[(k - 1) * m + k - 1];
    for(i = 1; i <= n; i++) {
	s = s * i / n;
	if(s < 1e-140) {
	    s *= 1e140;
	    eQ -= 140;
	}
    }
    s *= pow(10., eQ);
    R_Free(H);
    R_Free(Q);
    return(s);
}

static void
m_multiply(double *A, double *B, double *C, int m)
{
    int i, j, k;
    double s;
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++) {
	    s = 0.;
	    for(k = 0; k < m; k++)
		s+= A[i * m + k] * B[k * m + j];
	    C[i * m + j] = s;
	}
}

static void
m_power(double *A, int eA, double *V, int *eV, int m, int n)
{
    double *B;
    int eB, i;

    if(n == 1) {
	for(i = 0; i < m * m; i++)
	    V[i] = A[i];
	*eV = eA;
	return;
    }
    m_power(A, eA, V, eV, m, n / 2);
    B = R_Calloc(m * m, double);
    m_multiply(V, V, B, m);
    eB = 2 * (*eV);
    if((n % 2) == 0) {
	for(i = 0; i < m * m; i++)
	    V[i] = B[i];
	*eV = eB;
    }
    else {
	m_multiply(A, B, V, m);
	*eV = eA + eB;
    }
    if(V[(m / 2) * m + (m / 2)] > 1e140) {
	for(i = 0; i < m * m; i++)
	    V[i] = V[i] * 1e-140;
	*eV += 140;
    }
    R_Free(B);
}

SEXP pSmirnov2x(SEXP statistic, SEXP snx, SEXP sny)
{
    int nx = asInteger(snx), ny = asInteger(sny);
    double st = asReal(statistic);
    return ScalarReal(psmirnov2x(&st, nx, ny));
}

SEXP pKS2(SEXP statistic, SEXP stol)
{
    int n = LENGTH(statistic);
    double tol = asReal(stol);
    SEXP ans = duplicate(statistic);
    pkstwo(n, REAL(ans), tol);
    return ans;
}

SEXP pKolmogorov2x(SEXP statistic, SEXP sn)
{
    int n = asInteger(sn);
    double st = asReal(statistic), p;
    p = K(n, st);
    return ScalarReal(p);
}

