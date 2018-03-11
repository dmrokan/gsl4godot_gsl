/* filter/gsl_filter.h
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

#ifndef __GSL_FILTER_H__
#define __GSL_FILTER_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef enum
{
  GSL_FILTER_END_PADZERO = GSL_MOVSTAT_END_PADZERO,
  GSL_FILTER_END_PADVALUE = GSL_MOVSTAT_END_PADVALUE,
  GSL_FILTER_END_TRUNCATE = GSL_MOVSTAT_END_TRUNCATE
} gsl_filter_end_t;

/* workspace for recursive median filter */
typedef struct
{
  size_t H;        /* window half-length (K / 2) */
  size_t K;        /* window size */
  gsl_movstat_minmaxacc_workspace *minmaxacc_workspace_p;
} gsl_filter_rmedian_workspace;

gsl_filter_rmedian_workspace *gsl_filter_rmedian_alloc(const size_t K);
void gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w);
int gsl_filter_rmedian(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w);
int gsl_filter_rmedian2(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w);

typedef struct
{
  size_t K;   /* window size */
  gsl_movstat_workspace *movstat_workspace_p;
} gsl_filter_hampel_workspace;

gsl_filter_hampel_workspace *gsl_filter_hampel_alloc(const size_t K);
void gsl_filter_hampel_free(gsl_filter_hampel_workspace * w);
int gsl_filter_hampel(const gsl_filter_end_t endtype, const double nsigma, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian,
                      gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_hampel_workspace * w);
int gsl_filter_hampel2(const gsl_filter_end_t endtype, const double epsilon, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian,
                       gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_hampel_workspace * w);

__END_DECLS

#endif /* __GSL_FILTER_H__ */
