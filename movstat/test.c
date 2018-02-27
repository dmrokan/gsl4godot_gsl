/* movstat/test.c
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

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_ieee_utils.h>

#define SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

/* compare two vectors */
static void
compare_vectors(const double tol, const gsl_vector * v, const gsl_vector * expected,
                const char * desc)
{
  const size_t n = v->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double vi = gsl_vector_get(v, i);
      double ui = gsl_vector_get(expected, i);

      gsl_test_rel(vi, ui, tol, "%s i=%zu", desc, i);
    }
}

/* generate random vector with elements in [-1,1] */
static void
random_vector(gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double vi = 2.0 * gsl_rng_uniform(r) - 1.0; /* in [-1,1] */
      gsl_vector_set(v, i, vi);
    }
}

/* generate random vector of integers with elements in [a,b] */
static void
random_vector_int(const int a, const int b, gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      int vi = (int) ((b - a) * gsl_rng_uniform(r) + (double) a); /* in [a,b] */
      gsl_vector_set(v, i, (double) vi);
    }
}

/* find median of array z of length n by sorting */
static double
median_select(const long n, double *z)
{
  long low, high;
  long median;
  long middle, ll, hh;

  low = 0;
  high = n-1;
  median = (low + high) / 2;
  for (;;)
    {
      /* One element only */
      if (high <= low)
        return z[median];

      /* Two elements only */
      if (high == low + 1)
        {
          if (z[low] > z[high])
              SWAP(z[low], z[high]);
          return z[median];
        }

      /* Find median of low, middle and high items; swap to low position */
      middle = (low + high) / 2;
      if (z[middle] > z[high])
        SWAP(z[middle], z[high]);
      if (z[low] > z[high])
        SWAP(z[low], z[high]);
      if (z[middle] > z[low])
        SWAP(z[middle], z[low]);

      /* Swap low item (now in position middle) into position (low+1) */
      SWAP(z[middle], z[low+1]);

      /* Work from each end towards middle, swapping items when stuck */
      ll = low + 1;
      hh = high;
      for (;;)
        {
          do
            ll++;
          while (z[low] > z[ll]);
          do
            hh--;
          while (z[hh] > z[low]);

          if (hh < ll)
            break;

          SWAP(z[ll], z[hh]);
        }

      /* Swap middle item (in position low) back into correct position */
      SWAP(z[low], z[hh]);

      /* Reset active partition */
      if (hh <= median)
        low = ll;
      if (hh >= median)
        high = hh - 1;
    }
}

#include "test_mad.c"
#include "test_medacc.c"
#include "test_median.c"
#include "test_minmaxacc.c"

int
main()
{
  test_medacc();
  test_median();
  test_minmaxacc();
  test_mad();

  exit (gsl_test_summary());
}
