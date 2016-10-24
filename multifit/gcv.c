/* multifit/gcv.c
 * 
 * Copyright (C) 2016 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * References:
 *
 * [1] P. C. Hansen, "Discrete Inverse Problems: Insight and Algorithms,"
 * SIAM Press, 2010.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_min.h>

typedef struct
{
  const gsl_vector * S;
  const gsl_vector * UTy;
  double delta0;
  size_t np;
  gsl_vector * workp;
} gcv_params;

static int gcvmin(const gsl_vector * reg_param,
                  const gsl_vector * G,
                  gcv_params * params,
                  double * lambda,
                  gsl_multifit_linear_workspace * work);
static double gcv_func(double lambda, void * params);

/*
gsl_multifit_linear_gcv()
  Calculate Generalized Cross Validation curve for a set
of regularization parameters

Inputs: y         - right hand side vector
        reg_param - (output) regularization parameters
        G         - (output) GCV curve values
        work      - workspace
*/

int
gsl_multifit_linear_gcv(const gsl_vector * y,
                        gsl_vector * reg_param,
                        gsl_vector * G,
                        gsl_multifit_linear_workspace * work)
{
  const size_t n = y->size;
  const size_t N = G->size; /* number of points on GCV curve */

  if (n != work->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else if (reg_param->size != N)
    {
      GSL_ERROR ("size of reg_param and G vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t p = work->p;
      size_t i;

      gsl_matrix_view U = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view xt = gsl_vector_subvector(work->xt, 0, p);
      gsl_vector_view workp = gsl_matrix_subcolumn(work->QSI, 0, 0, p);

      const double smax = gsl_vector_get(&S.vector, 0);
      const double smin = gsl_vector_get(&S.vector, p - 1);
      gcv_params params;

      double dr; /* residual error from projection */
      double delta0 = 0.0;

      double normy = gsl_blas_dnrm2(y);
      double normUTy;
      double regmin;

      /* compute projection xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, &U.matrix, y, 0.0, &xt.vector);
      normUTy = gsl_blas_dnrm2(&xt.vector);

      /* dr = ||y||^2 - ||U^T y||^2 */
      dr = (normy + normUTy) * (normy - normUTy);

      /* calculate regularization parameters */
      gsl_multifit_linear_lreg(smin, smax, reg_param);

      if (n > p && dr > 0.0)
        delta0 = dr;

      params.S = &S.vector;
      params.UTy = &xt.vector;
      params.delta0 = delta0;
      params.np = n - p;
      params.workp = &workp.vector;

      for (i = 0; i < N; ++i)
        {
          double lambdai = gsl_vector_get(reg_param, i);
          double Gi = gcv_func(lambdai, &params);

          gsl_vector_set(G, i, Gi);
        }

      status = gcvmin(reg_param, G, &params, &regmin, work);

      return status;
    }
}

/*
gcvmin()
  Calculate Generalized Cross Validation curve for a set
of regularization parameters

Inputs: reg_param - regularization parameters
        G         - GCV curve values
        params    - gcv parameters
        lambda    - (output) regularization parameter which minimizes
                    GCV curve
        work      - workspace
*/

static int
gcvmin(const gsl_vector * reg_param,
       const gsl_vector * G,
       gcv_params * params,
       double * lambda,
       gsl_multifit_linear_workspace * work)
{
  const size_t npts = reg_param->size;
  const size_t max_iter = 500;
  const double tol = 1.0e-4;
  int status;
  int idxG = (int) gsl_vector_min_index(G);
  double a = gsl_vector_get(reg_param, GSL_MIN(idxG + 1, (int) npts - 1));
  double b = gsl_vector_get(reg_param, GSL_MAX(idxG - 1, 0));
  double m = gsl_vector_get(reg_param, idxG);
  size_t iter = 0;

  gsl_min_fminimizer *min_s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  gsl_function F;

  F.function = gcv_func;
  F.params = params;

  gsl_min_fminimizer_set(min_s, &F, m, a, b);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate(min_s);

      a = gsl_min_fminimizer_x_lower(min_s);
      b = gsl_min_fminimizer_x_upper(min_s);

      status = gsl_min_test_interval(a, b, 0.0, tol);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status == GSL_SUCCESS)
    *lambda = gsl_min_fminimizer_minimum(min_s);
  else
    status = GSL_EMAXITER;

  gsl_min_fminimizer_free(min_s);

  return status;
}

static double
gcv_func(double lambda, void * params)
{
  gcv_params * par = (gcv_params *) params;
  const gsl_vector *S = par->S;
  const gsl_vector *UTy = par->UTy;
  double delta0 = par->delta0;
  size_t np = par->np;
  gsl_vector *workp = par->workp;
  const size_t p = S->size;
  size_t i;
  double lambda_sq = lambda * lambda;
  double G, d, norm;
  double sumf = 0.0;

  /* compute workp = 1 - filter_factors */
  for (i = 0; i < p; ++i)
    {
      double si = gsl_vector_get(S, i);
      double fi = lambda_sq / (si * si + lambda_sq);
      gsl_vector_set(workp, i, fi);
      sumf += fi;
    }

  d = (double)np + sumf;

  gsl_vector_mul(workp, UTy);
  norm = gsl_blas_dnrm2(workp);

  G = (norm*norm + delta0) / (d * d);

  return G;
}
