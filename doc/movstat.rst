.. index::
   single: moving window statistics
   single: statistics, moving window

************************
Moving Window Statistics
************************

This chapter describes routines for computing *moving window
statistics* (also called rolling statistics and running statistics),
using a window around a sample which is used to calculate various
local statistical properties of an input data stream. The window is
then slid forward by one sample to process the next data point and so on.

The functions described in this chapter are declared in the header file
:file:`gsl_movstat.h`.

Introduction
============

This chapter is concerned with calculating various statistics from
subsets of a given dataset. The main idea is to compute statistics
in the vicinity of a given data sample by defining a *window* which
includes the sample itself as well as some specified number of samples
before and after the sample in question. For a sample :math:`x_i`, we
define a window :math:`W_i^{H,J}` as

.. math:: W_i^{H,J} = \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+J} \right\}

The parameters :math:`H` and :math:`J` are non-negative integers specifying
the number of samples to include before and after the sample :math:`x_i`.
Statistics such as the mean and standard deviation of the window :math:`W_i^{H,J}`
may be computed, and then the window is shifted forward by one sample to
focus on :math:`x_{i+1}`. The total number of samples in the window is
:math:`K = H + J + 1`. To define a symmetric window centered on :math:`x_i`,
one would set :math:`H = J = K / 2`.

Handling Endpoints
==================

When processing samples near the ends of the input signal, there will not
be enough samples to fill the window :math:`W_i^{H,J}` defined above.
Therefore the user must specify how to construct the windows near the end points.
This is done by passing an input argument of type :type:`gsl_movstat_end_t`:

.. type:: gsl_movstat_end_t

   This data type specifies how to construct windows near end points and can
   be selected from the following choices:

   .. macro:: GSL_MOVSTAT_END_PADZERO

      With this option, a full window of length :math:`K` will be constructed
      by inserting zeros into the window near the signal end points. Effectively,
      the input signal is modified to

      .. math:: \tilde{x} = \{ \underbrace{0, \dots, 0}_{H \textrm{ zeros}}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{0, \dots, 0}_{J \textrm{ zeros} } \}

      to ensure a well-defined window for all :math:`x_i`.

   .. macro:: GSL_MOVSTAT_END_PADVALUE

      With this option, a full window of length :math:`K` will be constructed
      by padding the window with the first and last sample in the input signal.
      Effectively, the input signal is modified to

      .. math:: \tilde{x} = \{ \underbrace{x_1, \dots, x_1}_{H}, x_1, x_2, \dots, x_{n-1}, x_n, \underbrace{x_n, \dots, x_n}_{J} \}

   .. macro:: GSL_MOVSTAT_END_TRUNCATE

      With this option, no padding is performed, and the windows are simply truncated
      as the end points are approached.

Moving Median
=============

The moving median calculates the median of the window :math:`W_i^{H,J}` for
each sample :math:`x_i`:

.. math:: y_i = \textrm{median} \left( W_i^{H,J} \right)

.. function:: gsl_movstat_median_workspace * gsl_movstat_median_alloc(const size_t K)

   This function allocates a workspace for computing a symmetric, centered moving median with a window
   length of :math:`K` samples. In this case, :math:`H = J = K/2`. The size of the workspace
   is :math:`O(4K)`.

.. function:: gsl_movstat_median_workspace * gsl_movstat_median_alloc2(const size_t H, const size_t J)

   This function allocates a workspace for computing a moving median using a window with :math:`H`
   samples prior to the current sample, and :math:`J` samples after the current sample. The
   total window size is :math:`K = H + J + 1`. The size of the workspace is :math:`O(4K)`.

.. function:: void * gsl_movstat_median_free(gsl_movstat_median_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_movstat_median(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_median_workspace * w)

   This function computes the moving median of the input vector :data:`x`, storing
   the output in :data:`y`. The parameter :data:`endtype` specifies how windows near
   the ends of the input should be handled.

Moving MAD
==========

The median absolute deviation (MAD) for the window :math:`W_i^{H,J}` is defined
to be the median of the absolute deviations from the window's median:

.. math:: MAD_i = \textrm{median} \left( \left| W_i^{H,J} - \textrm{median} \left( W_i^{H,J} \right) \right| \right)

In words, first the median of all samples in :math:`W_i^{H,J}` is computed. Then the median
is subtracted from all samples in the window to find the deviation of each sample
from the window median. The median of all absolute deviations is then the MAD.

.. function:: gsl_movstat_mad_workspace * gsl_movstat_mad_alloc(const size_t K)

   This function allocates a workspace for computing a symmetric, centered moving MAD with a window
   length of :math:`K` samples. In this case, :math:`H = J = K/2`. The size of the workspace
   is :math:`O(6K)`.

.. function:: gsl_movstat_mad_workspace * gsl_movstat_mad_alloc2(const size_t H, const size_t J)

   This function allocates a workspace for computing a moving MAD using a window with :math:`H`
   samples prior to the current sample, and :math:`J` samples after the current sample. The
   total window size is :math:`K = H + J + 1`. The size of the workspace is :math:`O(6K)`.

.. function:: void * gsl_movstat_mad_free(gsl_movstat_mad_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_movstat_mad(const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad, gsl_movstat_mad_workspace * w)

   This function computes the moving MAD of the input vector :data:`x` and stores the result
   in :data:`xmad`. The medians of each window :math:`W_i^{H,J}` are stored in :data:`xmedian`
   on output. The inputs :data:`x`, :data:`xmedian`, and :data:`xmad` must all be the same length.

Accumulators
============

Many of the moving statistics routines in this chapter are based on an accumulator design,
which track the desired statistic of a fixed-size window. Each time a new sample is
added to the accumulator (indicating the window is sliding forward), the oldest sample
is discarded and the relevant statistic is updated to incorporate the new sample. Most
users will likely want to use the routines described above, in which they can
input an entire vector of data, calculate the moving statistic, and obtain an output vector.
However the routines below are provided for specialized applications which may need
to implement a sliding window.

Median Accumulator
------------------

The median accumulator uses an efficient heap-based algorithm by Härdle and Steiger
to update the median of a window each time a new sample is added.

.. function:: gsl_movstat_medacc_workspace * gsl_movstat_medacc_alloc(const size_t K)

   This function allocates a workspace for tracking the median of a window of
   length :math:`K`. The size of the workspace is :math:`O(3K)`.

.. function:: void gsl_movstat_medacc_free(gsl_movstat_medacc_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_movstat_medacc_insert(const double x, gsl_movstat_medacc_workspace * w)

   This function inserts a single data sample :data:`x` into the accumulator
   and computes the median of all samples in the current window. If the window
   is full, the oldest sample is discarded.

.. function:: int gsl_movstat_medacc_reset(gsl_movstat_medacc_workspace * w)

   This function resets the accumulator to its initial state to begin working on
   a new window of data.

.. function:: double gsl_movstat_medacc_median(const gsl_movstat_medacc_workspace * w)

   This function returns the median of the current data window.

Min/Max Accumulator
-------------------

The min/max accumulator efficiently tracks the minimum and maximum values of the
current sliding window, using the algorithm by D. Lemire.

.. function:: gsl_movstat_minmaxacc_workspace * gsl_movstat_minmaxacc_alloc(const size_t K)

   This function allocates a workspace for tracking the minimum and maximum of a window of
   length :math:`K`. The size of the workspace is :math:`O(3K)`.

.. function:: void gsl_movstat_minmaxacc_free(gsl_movstat_minmaxacc_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_movstat_minmaxacc_insert(const double x, gsl_movstat_minmaxacc_workspace * w)

   This function inserts a single data sample :data:`x` into the accumulator
   and computes the minimum and maximum of all samples in the current window. If the window
   is full, the oldest sample is discarded.

.. function:: int gsl_movstat_minmaxacc_minmax(double * min, double * max, const gsl_movstat_minmaxacc_workspace * w)

   This function stores the minimum and maximum values of the current window in
   :data:`min` and :data:`max` respectively.

.. function:: double gsl_movstat_minmaxacc_min(const gsl_movstat_minmaxacc_workspace * w)

   This function returns the minimum value of the current window.

.. function:: double gsl_movstat_minmaxacc_max(const gsl_movstat_minmaxacc_workspace * w)

   This function returns the maximum value of the current window.

References and Further Reading
==============================

The following publications are relevant to the algorithms described
in this section,

* W. Härdle and W. Steiger, *Optimal Median Smoothing*, Appl. Statist. 44 (2), 1995.

* D. Lemire, *Streaming Maximum-Minimum Filter Using No More than Three Comparisons per Element*,
  Nordic Journal of Computing, 13 (4), 2006 (https://arxiv.org/abs/cs/0610046).
