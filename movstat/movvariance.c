/* movstat/movvar.c
 *
 * Routines related to a moving window variance
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

#include "movstat_common.c"
#include "mvacc.c"

/*
gsl_movstat_variance()
  Apply moving variance to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector, size n
        w       - workspace
*/

int
gsl_movstat_variance(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = movstat_apply(endtype, x, y, mvacc_init, mvacc_add, mvacc_delete, mvacc_variance, w);
  return status;
}

/*
gsl_movstat_sd()
  Apply moving standard deviation to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector, size n
        w       - workspace
*/

int
gsl_movstat_sd(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = movstat_apply(endtype, x, y, mvacc_init, mvacc_add, mvacc_delete, mvacc_sd, w);
  return status;
}
