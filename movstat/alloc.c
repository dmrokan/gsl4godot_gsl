/* movstat/alloc.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>

#include "medacc.c"
#include "mmacc.c"
#include "mvacc.c"
#include "qnacc.c"
#include "snacc.c"
#include "sumacc.c"

/*
gsl_movstat_alloc()
  Allocate a workspace for moving window statistics.
The window around sample x_i is defined as:

W_i^{H,J} = {x_{i-H},...,x_i,...x_{i+J}}

The total window size is:

K = H + J + 1

Inputs: K - total samples in window (H = J = K / 2)

Return: pointer to workspace

Notes:
1) If K is even, it is rounded up to the next odd
*/

gsl_movstat_workspace *
gsl_movstat_alloc(const size_t K)
{
  const size_t H = K / 2;
  return gsl_movstat_alloc2(H, H);
}
 
/*
gsl_movstat_alloc2()
  Allocate a workspace for moving window statistics.
The window around sample x_i is defined as:

W_i^{H,J} = {x_{i-H},...,x_i,...x_{i+J}}

The total window size is:

K = H + J + 1

Inputs: H - number of samples before current sample
        J - number of samples after current sample

Return: pointer to workspace
*/

gsl_movstat_workspace *
gsl_movstat_alloc2(const size_t H, const size_t J)
{
  gsl_movstat_workspace *w;
  size_t state_size = 0;

  w = calloc(1, sizeof(gsl_movstat_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = H;
  w->J = J;
  w->K = H + J + 1;

  /*
   * determine maximum number of bytes needed for the various accumulators;
   * the accumulators will all share the same workspace
   */
  state_size = GSL_MAX(state_size, mvacc_size(w->K));  /* mean/variance accumulator */
  state_size = GSL_MAX(state_size, mmacc_size(w->K));  /* min/max accumulator */
  state_size = GSL_MAX(state_size, sumacc_size(w->K)); /* sum accumulator */
  state_size = GSL_MAX(state_size, medacc_size(w->K)); /* median accumulator */
  state_size = GSL_MAX(state_size, qnacc_size(w->K));  /* Q_n accumulator */
  state_size = GSL_MAX(state_size, snacc_size(w->K));  /* S_n accumulator */

  w->state = malloc(state_size);
  if (w->state == 0)
    {
      gsl_movstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for accumulator state", GSL_ENOMEM);
    }

  w->work = malloc(w->K * sizeof(double));
  if (w->work == 0)
    {
      gsl_movstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for work", GSL_ENOMEM);
    }

  return w;
}

void
gsl_movstat_free(gsl_movstat_workspace * w)
{
  if (w->work)
    free(w->work);

  if (w->state)
    free(w->state);

  free(w);
}
