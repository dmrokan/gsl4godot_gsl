/* movstat/test_mean.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_movstat.h>

/* compute moving mean by explicitely constructing window and computing mean */
int
slow_movmean(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
             const int H, const int J)
{
  const int n = (int) x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  int i, j;

  for (i = 0; i < n; ++i)
    {
      int idx1, idx2, wsize;
      double mean;

      if (etype == GSL_MOVSTAT_END_TRUNCATE)
        {
          idx1 = GSL_MAX(i - H, 0);
          idx2 = GSL_MIN(i + J, n - 1);
        }
      else
        {
          idx1 = i - H;
          idx2 = i + J;
        }

      wsize = idx2 - idx1 + 1;

      /* fill sliding window */
      for (j = idx1; j <= idx2; ++j)
        {
          int widx = j - idx1;

          if (j < 0)
            {
              /* initial condition */
              if (etype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (etype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, 0);
            }
          else if (j >= (int) n)
            {
              if (etype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (etype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, n - 1);
            }
          else
            {
              window[widx] = gsl_vector_get(x, j);
            }
        }

      mean = gsl_stats_mean(window, 1, wsize);
      gsl_vector_set(y, i, mean);
    }

  free(window);

  return GSL_SUCCESS;
}

static void
test_mean_proc(const double tol, const size_t n, const size_t H, const size_t J,
               const gsl_movstat_end_t etype, gsl_rng * rng_p)
{
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * z = gsl_vector_alloc(n);
  char buf[2048];

  random_vector(x, rng_p);

  /* y = mean(x) with slow brute force method */
  slow_movmean(etype, x, y, H, J);

  /* y = mean(x) with fast method */
  gsl_movstat_mean(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u mean random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = mean(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_mean(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u mean random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_movstat_free(w);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
}

static void
test_mean(gsl_rng * rng_p)
{
  const double eps = 1.0e-10;

  test_mean_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 1000, 0, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 1000, 5, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mean_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_mean_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 1000, 0, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 1000, 5, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mean_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_mean_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 1000, 0, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 1000, 5, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mean_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
