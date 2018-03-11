/* movstat/minmaxacc.c
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

/*
 * This module contains routines for tracking minimum/maximum values of a
 * moving fixed-sized window. It is based on the algorithm of:
 *
 * [1] Daniel Lemire, Streaming Maximum-Minimum Filter Using No More than Three Comparisons per Element,
 *     Nordic Journal of Computing, Volume 13, Number 4, pages 328-339, 2006
 *
 * Also available as a preprint here: https://arxiv.org/abs/cs/0610046
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

#include "deque.c"

/*
gsl_movstat_minmaxacc_alloc()
  Allocate workspace for tracking min/max values of a sliding window

Inputs: n - size of window

Return: pointer to workspace
*/

gsl_movstat_minmaxacc_workspace *
gsl_movstat_minmaxacc_alloc(const size_t n)
{
  gsl_movstat_minmaxacc_workspace *w;

  if (n < 1)
    {
      GSL_ERROR_NULL ("window size must be at least 1", GSL_EBADLEN);
    }

  w = calloc(1, sizeof(gsl_movstat_minmaxacc_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->window = malloc(n * sizeof(double));
  if (w->window == NULL)
    {
      gsl_movstat_minmaxacc_free(w);
      GSL_ERROR_NULL ("failed to allocate space for window", GSL_ENOMEM);
    }

  w->maxque = deque_alloc(n + 1);
  if (w->maxque == NULL)
    {
      gsl_movstat_minmaxacc_free(w);
      GSL_ERROR_NULL ("failed to allocate space for maxque", GSL_ENOMEM);
    }

  w->minque = deque_alloc(n + 1);
  if (w->minque == NULL)
    {
      gsl_movstat_minmaxacc_free(w);
      GSL_ERROR_NULL ("failed to allocate space for minque", GSL_ENOMEM);
    }

  w->n = n;

  gsl_movstat_minmaxacc_reset(w);

  return w;
}

/*
gsl_movstat_minmaxacc_free()
  Free workspace
*/

void
gsl_movstat_minmaxacc_free(gsl_movstat_minmaxacc_workspace * w)
{
  if (w->window)
    free(w->window);

  if (w->maxque)
    deque_free(w->maxque);

  if (w->minque)
    deque_free(w->minque);

  free(w);
}

int
gsl_movstat_minmaxacc_reset(gsl_movstat_minmaxacc_workspace * w)
{
  /* empty the queues */
  deque_empty(w->minque);
  deque_empty(w->maxque);

  w->idx = 0;
  w->init = 0;
  w->xprev = 0.0;

  return GSL_SUCCESS;
}

/*
gsl_movstat_minmaxacc_insert()
  Insert a new sample into sliding window. If window is saturated, the oldest
sample is discarded and the new sample takes its place.

Inputs: x - value of sample
        w - workspace

Return: success/error

Notes:
1) On entrance, w->idx is in [0,K-1] and points to the index in window[] where to
store the next sample. On exit, w->idx is incremented, wrapping back to 0 if necessary.
*/

int
gsl_movstat_minmaxacc_insert(const double x, gsl_movstat_minmaxacc_workspace * w)
{
  if (w->init == 0)
    {
      /* first sample */
      w->window[0] = x;
      deque_push_front(0, w->maxque);
      deque_push_front(0, w->minque);
      w->idx = 1;
      w->init = 1;
    }
  else
    {
      const size_t window_size = w->n;

      if (x > w->xprev)
        {
          deque_pop_back(w->maxque);

          while (!deque_is_empty(w->maxque))
            {
              if (x <= w->window[deque_peek_back(w->maxque) % window_size])
                break;

              deque_pop_back(w->maxque);
            }
        }
      else
        {
          deque_pop_back(w->minque);

          while (!deque_is_empty(w->minque))
            {
              if (x >= w->window[deque_peek_back(w->minque) % window_size])
                break;

              deque_pop_back(w->minque);
            }
        }

      /* store new sample and its index in window[] */
      w->window[w->idx % window_size] = x;
      deque_push_back(w->idx, w->maxque);
      deque_push_back(w->idx, w->minque);

      if (w->idx == window_size + deque_peek_front(w->maxque))
        deque_pop_front(w->maxque);
      else if (w->idx == window_size + deque_peek_front(w->minque))
        deque_pop_front(w->minque);

      ++(w->idx);
    }

  w->xprev = x;

  return GSL_SUCCESS;
}

int
gsl_movstat_minmaxacc_minmax(double * min, double * max, const gsl_movstat_minmaxacc_workspace * w)
{
  if (w->init == 0)
    {
      GSL_ERROR ("no samples yet added to workspace", GSL_EINVAL);
    }
  else
    {
      *min = w->window[deque_peek_front(w->minque) % w->n];
      *max = w->window[deque_peek_front(w->maxque) % w->n];
      return GSL_SUCCESS;
    }
}

double
gsl_movstat_minmaxacc_min(const gsl_movstat_minmaxacc_workspace * w)
{
  if (w->init == 0)
    {
      GSL_ERROR_VAL ("no samples yet added to workspace", GSL_EINVAL, 0.0);
    }
  else
    {
      return (w->window[deque_peek_front(w->minque) % w->n]);
    }
}

double
gsl_movstat_minmaxacc_max(const gsl_movstat_minmaxacc_workspace * w)
{
  if (w->init == 0)
    {
      GSL_ERROR_VAL ("no samples yet added to workspace", GSL_EINVAL, 0.0);
    }
  else
    {
      return (w->window[deque_peek_front(w->maxque) % w->n]);
    }
}
