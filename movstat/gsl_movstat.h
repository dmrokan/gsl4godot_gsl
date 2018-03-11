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

/****** Accumulators *******/

typedef double gsl_movstat_medacc_t;

/* workspace for tracking median of a fixed-size moving window */
typedef struct
{
  gsl_movstat_medacc_t *data; /* circular queue of values, size k */
  int *pos;                    /* index into `heap` for each value, size 2*k */
  int *heap;                   /* max/median/min heap holding indices into `data`. */
  int k;                       /* window size */
  int idx;                     /* position in circular queue */
  int ct;                      /* count of items in queue */
} gsl_movstat_medacc_workspace;

gsl_movstat_medacc_workspace *gsl_movstat_medacc_alloc(const size_t k);
void gsl_movstat_medacc_free(gsl_movstat_medacc_workspace * w);
int gsl_movstat_medacc_reset(gsl_movstat_medacc_workspace * w);
int gsl_movstat_medacc_insert(const gsl_movstat_medacc_t v, gsl_movstat_medacc_workspace * w);
gsl_movstat_medacc_t gsl_movstat_medacc_median(const gsl_movstat_medacc_workspace * w);

/* workspace for tracking minimum/maximum of a fixed-size moving window */
typedef struct
{
  size_t n;           /* window size */
  double *window;     /* samples in current window, size n */
  size_t idx;         /* index where to store next sample */
  int init;           /* 0 if no samples have been added yet */
  double xprev;       /* previous sample added to window */
  void *minque;       /* double-ended queue of min values (L) */
  void *maxque;       /* double-ended queue of max values (U) */
} gsl_movstat_minmaxacc_workspace;

gsl_movstat_minmaxacc_workspace *gsl_movstat_minmaxacc_alloc(const size_t n);
void gsl_movstat_minmaxacc_free(gsl_movstat_minmaxacc_workspace * w);
int gsl_movstat_minmaxacc_reset(gsl_movstat_minmaxacc_workspace * w);
int gsl_movstat_minmaxacc_insert(const double x, gsl_movstat_minmaxacc_workspace * w);
int gsl_movstat_minmaxacc_minmax(double * min, double * max, const gsl_movstat_minmaxacc_workspace * w);
double gsl_movstat_minmaxacc_min(const gsl_movstat_minmaxacc_workspace * w);
double gsl_movstat_minmaxacc_max(const gsl_movstat_minmaxacc_workspace * w);

/********* Routines accepting full vector inputs/outputs ***********/

/* workspace for moving window statistics */

typedef struct
{
  size_t H;     /* number of previous samples in window */
  size_t J;     /* number of after samples in window */
  size_t K;     /* window size K = H + J + 1 */
  double *work; /* workspace, size K */
  gsl_movstat_medacc_workspace *medacc_workspace_p;
  gsl_movstat_minmaxacc_workspace *minmaxacc_workspace_p;
} gsl_movstat_workspace;

/* alloc.c */

gsl_movstat_workspace *gsl_movstat_alloc(const size_t K);
gsl_movstat_workspace *gsl_movstat_alloc2(const size_t H, const size_t J);
void gsl_movstat_free(gsl_movstat_workspace * w);

int gsl_movstat_median(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w);
int gsl_movstat_minmax(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max, gsl_movstat_workspace * w);
int gsl_movstat_mad(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian,
                    gsl_vector * xmad, gsl_movstat_workspace * w);

__END_DECLS

#endif /* __GSL_MOVSTAT_H__ */
