/* movstat/medacc.c
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
 *
 * Original copyright notice:
 * Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>
 */
 
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>
 
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)
 
/*--- Helper Functions ---*/
 
#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap 

static int mmless(gsl_movstat_medacc_workspace* m, int i, int j);
static int mmexchange(gsl_movstat_medacc_workspace* w, int i, int j);
static int mmCmpExch(gsl_movstat_medacc_workspace* w, int i, int j);
static void minSortDown(gsl_movstat_medacc_workspace* w, int i);
static void maxSortDown(gsl_movstat_medacc_workspace* w, int i);
static int minSortUp(gsl_movstat_medacc_workspace* w, int i);
static int maxSortUp(gsl_movstat_medacc_workspace* w, int i);

/*
gsl_movstat_medacc_alloc()
  Allocate a workspace for median filtering. The workspace
is set up to calculate a running median with a given window size.

Inputs: k - number of samples in window

Return: pointer to workspace
*/

gsl_movstat_medacc_workspace *
gsl_movstat_medacc_alloc(const size_t k)
{
  gsl_movstat_medacc_workspace *w;

  w = calloc(1, sizeof(gsl_movstat_medacc_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->data = malloc(k * sizeof(gsl_movstat_medacc_t));
  if (w->data == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for data", GSL_ENOMEM);
    }

  w->pos = malloc(2 * k * sizeof(int));
  if (w->pos == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for pos", GSL_ENOMEM);
    }

  w->heap = w->pos + k + (k/2); /* points to middle of storage */
  w->k = k;

  /* initialize heap */
  gsl_movstat_medacc_reset(w);

  return w;
}

void
gsl_movstat_medacc_free(gsl_movstat_medacc_workspace * w)
{
  if (w->data)
    free(w->data);

  if (w->pos)
    free(w->pos);

  free(w);
}

/*
gsl_movstat_medacc_reset()
  Reset mediator to initial state
*/

int
gsl_movstat_medacc_reset(gsl_movstat_medacc_workspace * w)
{
  int kk = (int) w->k;

  w->ct = 0;
  w->idx = 0;

  /* set up initial heap fill pattern: median,max,min,max,... */
  while (kk--)
    {
      w->pos[kk] = ((kk + 1)/2) * ((kk & 1) ? -1 : 1);
      w->heap[w->pos[kk]] = kk;
    }

  return GSL_SUCCESS;
}

//Inserts item, maintains median in O(lg k)
int
gsl_movstat_medacc_insert(const gsl_movstat_medacc_t v, gsl_movstat_medacc_workspace * w)
{
  int isNew = (w->ct < (int) w->k);
  int p = w->pos[w->idx];
  gsl_movstat_medacc_t old = w->data[w->idx];

  w->data[w->idx] = v;
  w->idx = (w->idx+1) % w->k;
  w->ct += isNew;

  if (p>0)         //new item is in minHeap
    {
      if (!isNew && ItemLess(old, v))
        minSortDown(w, p * 2);
     else if (minSortUp(w, p))
       maxSortDown(w, -1);
    }
  else if (p<0)   //new item is in maxheap
    {
      if (!isNew && ItemLess(v, old))
        maxSortDown(w, p * 2);
     else if (maxSortUp(w, p))
       minSortDown(w, 1);
    }
  else            //new item is at median
    {
      if (maxCt(w))
        maxSortDown(w, -1);
      if (minCt(w))
        minSortDown(w, 1);
    }

  return GSL_SUCCESS;
}
 
//returns median item (or average of 2 when item count is even)
gsl_movstat_medacc_t
gsl_movstat_medacc_median(const gsl_movstat_medacc_workspace * w)
{
  gsl_movstat_medacc_t v = w->data[w->heap[0]];

  if ((w->ct & 1) == 0)
    v = ItemMean(v, w->data[w->heap[-1]]);

  return v;
}
 
//returns 1 if heap[i] < heap[j]
static int
mmless(gsl_movstat_medacc_workspace * w, int i, int j)
{
  return ItemLess(w->data[w->heap[i]], w->data[w->heap[j]]);
}
 
//swaps items i&j in heap, maintains indexes
static int
mmexchange(gsl_movstat_medacc_workspace* w, int i, int j)
{
  int t = w->heap[i];
  w->heap[i]=w->heap[j];
  w->heap[j]=t;
  w->pos[w->heap[i]]=i;
  w->pos[w->heap[j]]=j;
  return 1;
}
 
//swaps items i&j if i<j;  returns true if swapped
static int
mmCmpExch(gsl_movstat_medacc_workspace* w, int i, int j)
{
  return (mmless(w,i,j) && mmexchange(w,i,j));
}
 
//maintains minheap property for all items below i/2.
static void
minSortDown(gsl_movstat_medacc_workspace * w, int i)
{
  for (; i <= minCt(w); i*=2)
    {
      if (i>1 && i < minCt(w) && mmless(w, i+1, i))
        ++i;

      if (!mmCmpExch(w,i,i/2))
        break;
   }
}
 
//maintains maxheap property for all items below i/2. (negative indexes)
static void
maxSortDown(gsl_movstat_medacc_workspace* w, int i)
{
  for (; i >= -maxCt(w); i*=2)
    {
      if (i<-1 && i > -maxCt(w) && mmless(w, i, i-1))
        --i;

      if (!mmCmpExch(w,i/2,i))
        break;
   }
}
 
//maintains minheap property for all items above i, including median
//returns true if median changed
static int
minSortUp(gsl_movstat_medacc_workspace* w, int i)
{
   while (i>0 && mmCmpExch(w,i,i/2)) i/=2;
   return (i==0);
}
 
//maintains maxheap property for all items above i, including median
//returns true if median changed
static int
maxSortUp(gsl_movstat_medacc_workspace* w, int i)
{
   while (i<0 && mmCmpExch(w,i/2,i))  i/=2;
   return (i==0);
}
