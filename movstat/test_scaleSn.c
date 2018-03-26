/* movstat/test_scaleSn.c
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
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_test.h>

static void
test_scaleSn(void)
{
  const size_t n = 19;
  const double tol = GSL_DBL_EPSILON;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *work = gsl_vector_alloc(n);
  double Sn;
  size_t i;

  for (i = 0; i < n; ++i)
    gsl_vector_set(x, i, i + 1.0);

  Sn = gsl_movstat_full_scaleSn(x);
  gsl_test_rel(Sn, 5.0, tol, "S_n, x = [1:19]");

  /* add some outliers */
  for (i = 10; i < n; ++i)
    gsl_vector_set(x, i, i + 91.0);

  /*Sn = gsl_movstat_full_scaleSn(x);*/
  gsl_sort_vector(x);
  Sn = gsl_stats_sn_from_sorted_data(x->data, x->stride, x->size, work->data);
  gsl_test_rel(Sn, 9.0, tol, "S_n, x = [1:10, 100 + 1:9]");

  gsl_vector_free(x);
  gsl_vector_free(work);
}
