/* movstat/test_medacc.c
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
#include <gsl/gsl_movstat.h>

void
test_medacc(void)
{
  const double tol = 1.0e-12;

  /* first test medacc on a full dataset without a sliding window */
  {
    const double data[] = { -1.0, 5.7, 3.4, 1.1, 9.5, -23.7, -5.6, 0.2 };

    /* expected_median[i] = median[ data(1:i+1) ] */
    const double expected_median[] = { -1.0, 2.35, 3.4, 2.25, 3.4, 2.25, 1.1, 0.65 };
    size_t k;

    /* compute the median of each data(1:k) up to the full array */
    for (k = 1; k <= 8; ++k)
      {
        gsl_movstat_medacc_workspace *w = gsl_movstat_medacc_alloc(k);
        double median;
        size_t i;

        for (i = 0; i < k; ++i)
          gsl_movstat_medacc_insert(data[i], w);

        median = gsl_movstat_medacc_median(w);
        gsl_test_rel(median, expected_median[k - 1], tol, "mediator fixed k=%zu");

        gsl_movstat_medacc_free(w);
      }
  }

  /* test reset function */
  {
    gsl_movstat_medacc_workspace *w = gsl_movstat_medacc_alloc(111);
    double median;

    gsl_movstat_medacc_insert(1.0, w);
    gsl_movstat_medacc_insert(2.0, w);
    gsl_movstat_medacc_insert(3.0, w);
    median = gsl_movstat_medacc_median(w);
    gsl_test_rel(median, 2.0, tol, "mediator [1 2 3]");

    gsl_movstat_medacc_reset(w);

    gsl_movstat_medacc_insert(4.0, w);
    gsl_movstat_medacc_insert(5.0, w);
    gsl_movstat_medacc_insert(6.0, w);
    median = gsl_movstat_medacc_median(w);
    gsl_test_rel(median, 5.0, tol, "mediator [4 5 6]");

    gsl_movstat_medacc_free(w);
  }
}
