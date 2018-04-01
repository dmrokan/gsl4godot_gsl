/* movstat/mvacc.c
 *
 * Moving window mean/variance accumulator
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

typedef double ringbuf_type;
#include "ringbuf.c"

typedef struct
{
  size_t n;      /* window size */
  size_t k;      /* number of samples currently in window */
  double mean;   /* current window mean */
  ringbuf *rbuf; /* ring buffer storing current window */
} mvacc_state_t;

static size_t
mvacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(mvacc_state_t);
  size += ringbuf_size(n);

  return size;
}

static int
mvacc_init(const size_t n, void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  state->n = n;
  state->k = 0;
  state->mean = 0.0;

  state->rbuf = vstate + sizeof(mvacc_state_t);
  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
mvacc_add(const double x, void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  if (ringbuf_is_full(state->rbuf))
    {
      /* remove oldest window element and add new one */
      state->mean += (x - ringbuf_peek_back(state->rbuf)) / (double) state->n;
    }
  else
    {
      state->mean += (x - state->mean) / (state->k + 1.0);
      ++(state->k);
    }

  /* add new element to ring buffer */
  ringbuf_add(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
mvacc_delete(void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf) && state->k > 0)
    {
      double old = ringbuf_peek_back(state->rbuf);
      state->mean += (state->mean - old) / (state->k - 1.0);
      ringbuf_pop_back(state->rbuf);
      --(state->k);
    }

  return GSL_SUCCESS;
}

static double
mvacc_mean(const void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;
  return state->mean;
}
