/* movstat/movstat_common.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

/*
movstat_fill_window()
  Fill window for sample 'idx' from x using given end conditions

Inputs: endtype - how to handle end points
        idx     - index of center sample in window
        H       - number of samples left of center to include
        J       - number of samples right of center to include
        x       - input vector
        window  - (output) window of samples centered on x_{idx}

Return: number of samples in window
*/

static int
movstat_fill_window(const gsl_movstat_end_t endtype, const int idx, const int H, const int J,
                    const gsl_vector * x, double * window)
{
  const int n = x->size;
  int idx1, idx2, j;
  int wsize;

  if (endtype == GSL_MOVSTAT_END_TRUNCATE)
    {
      idx1 = GSL_MAX(idx - H, 0);
      idx2 = GSL_MIN(idx + J, n - 1);
    }
  else
    {
      idx1 = idx - H;
      idx2 = idx + J;
    }

  wsize = idx2 - idx1 + 1;

  /* fill sliding window */
  for (j = idx1; j <= idx2; ++j)
    {
      int widx = j - idx1;

      if (j < 0)
        {
          /* initial condition */
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            window[widx] = 0.0;
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            window[widx] = gsl_vector_get(x, 0);
        }
      else if (j >= n)
        {
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            window[widx] = 0.0;
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            window[widx] = gsl_vector_get(x, n - 1);
        }
      else
        {
          window[widx] = gsl_vector_get(x, j);
        }
    }

  wsize = idx2 - idx1 + 1;

  return wsize;
}
