/* movstat/gsl_movstat.h
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

#ifndef __GSL_MOVSTAT_H__
#define __GSL_MOVSTAT_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

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
  GSL_MOVSTAT_END_PADZERO,
  GSL_MOVSTAT_END_PADVALUE,
  GSL_MOVSTAT_END_TRUNCATE
} gsl_movstat_end_t;

/* workspace for moving window statistics */

typedef struct
{
  size_t H;     /* number of previous samples in window */
  size_t J;     /* number of after samples in window */
  size_t K;     /* window size K = H + J + 1 */

  union
    {
      void *state;  /* state workspace for various accumulators */
      double *work; /* workspace, size K */
    };
} gsl_movstat_workspace;

/* alloc.c */

gsl_movstat_workspace *gsl_movstat_alloc(const size_t K);
gsl_movstat_workspace *gsl_movstat_alloc2(const size_t H, const size_t J);
void gsl_movstat_free(gsl_movstat_workspace * w);

int gsl_movstat_mean(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_variance(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_sd(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_median(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_min(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_max(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_minmax(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max, gsl_movstat_workspace * w);
int gsl_movstat_mad(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian,
                    gsl_vector * xmad, gsl_movstat_workspace * w);
int gsl_movstat_qqr(const gsl_movstat_end_t endtype, const gsl_vector * x, const double q,
                    gsl_vector * xqqr, gsl_movstat_workspace * w);
int gsl_movstat_Sn(const gsl_movstat_end_t endtype, const gsl_vector * x,
                   gsl_vector * xscale, gsl_movstat_workspace * w);
int gsl_movstat_Qn(const gsl_movstat_end_t endtype, const gsl_vector * x,
                   gsl_vector * xscale, gsl_movstat_workspace * w);
int gsl_movstat_sum(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);

__END_DECLS

#endif /* __GSL_MOVSTAT_H__ */
