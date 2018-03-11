*****************
Digital Filtering
*****************

Introduction
============

The filters discussed in this chapter are based on the following moving data
window which is centered on :math:`i`-th sample:

.. math:: W_i^H = \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

Here, :math:`H` is a non-negative integer called the *window half-length*, which
represents the number of samples before and after sample :math:`i`.
The total window length is :math:`K = 2 H + 1`.

Nonlinear Digital Filters
=========================

The nonlinear digital filters described below are based on the window median, which is given
by

.. math:: m_i = \textrm{median} \left\{ W_i^H \right\} = \textrm{median} \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

The median is considered robust to local outliers, unlike the mean.
Median filters can preserve sharp edges while at the same removing signal noise, and are used
in a wide range of applications.

Standard Median Filter
----------------------

The *standard median filter* (SMF) simply replaces the sample :math:`x_i` by the median
:math:`m_i` of the window :math:`W_i^H`: This filter has one tuning parameter given
by :math:`H`. The standard median filter is considered highly resistant to
local outliers and local noise in the data sequence :math:`\{x_i\}`.

.. function:: gsl_filter_median_workspace * gsl_filter_median_alloc(const size_t K)

   This function initializes a workspace for standard median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(3K)`.

.. function:: void gsl_filter_median_free(gsl_filter_median_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_median(const gsl_vector * x, gsl_vector * y, gsl_filter_median_workspace * w)

   This function applies a standard median filter to the input :data:`x`, storing the output in :data:`y`. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Recursive Median Filter
-----------------------

The recursive median filter (RMF) is a modification of the SMF to include previous filter outputs
in the window before computing the median. The filter's response is

.. math:: y_i = \textrm{median} \left( y_{i-H}, \dots, y_{i-1}, x_i, x_{i+1}, \dots, x_{i+H} \right)

Sometimes, the SMF must be applied several times in a row to achieve adequate smoothing (i.e. a cascade filter).
The RMF, on the other hand, converges to a *root sequence* in one pass,
and can sometimes provide a smoother result than several passes of the SMF. A root sequence is an input which is
left unchanged by the filter.  So there is no need to apply an RMF filter twice to an input vector.

.. function:: gsl_filter_rmedian_workspace * gsl_filter_rmedian_alloc(const size_t K)

   This function initializes a workspace for recursive median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(K)`.

.. function:: void gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_rmedian(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)

   This function applies a recursive median filter to the input :data:`x`, storing the output in :data:`y`. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Hampel Filter
-------------

The Hampel filter is closely related to the standard median filter, which is designed to detect
outliers based on the median and median absolute deviation (MAD) scale estimator.
This filter's response is

.. math:: y_i = \left\{
                  \begin{array}{cc}
                    x_i, & |x_i - m_i| \le t S_i \\
                    m_i, & |x_i - m_i| > t S_i
                  \end{array}
                \right.

where :math:`m_i` is the median value of the window :math:`W_i^H`, and :math:`S_i` is the
MAD scale estimate, defined by

.. math:: S_i = 1.4826 \times \textrm{median} \left\{ | W_i^H - m_i | \right\}

In other words, it takes the median of all the absolute deviations of each sample in the window :math:`W_i^H`
from its local window median :math:`m_i`. The factor :math:`1.4826` makes :math:`S_i` an unbiased estimate of
the standard deviation for Gaussian data. The MAD statistic is used instead of the window standard deviation
since it is much less sensitive to outliers. Finally, :math:`t` is a tunable parameter for deciding
how far a particular sample should be from the MAD scale estimate to identify it as an outlier. Identified
outliers are replaced by the local median, while unflagged data are left unchanged.

Note that when :math:`t = 0`, the Hampel filter is equivalent to the standard median filter. When
:math:`t \rightarrow \infty`, it becomes the identity filter. This means the Hampel filter can
be viewed as a "less aggressive" version of the standard median filter, becoming less aggressive as :math:`t` is
increased. Note that the Hampel filter modifies only samples identified as outliers, while the standard median
filter changes all samples to the local median, regardless of whether they are outliers. This fact, plus
the additional flexibility offered by the additional tuning parameter :math:`t` can make the Hampel filter
a better choice for some applications.

.. warning::

   While the MAD is much less sensitive to outliers than the standard deviation, it can suffer from an
   effect called *MAD implosion*. The standard deviation of a window :math:`W_i^H` will be zero
   if and only if all samples in the window are equal. However, it is possible for the MAD of a window
   to be zero even if all the samples in the window are not equal. For example, if :math:`K/2 + 1` or more
   of the :math:`K` samples in the window are equal to some value :math:`x^{*}`, then the window median will
   be equal to :math:`x^{*}`. Consequently, at least :math:`K/2 + 1` of the absolute deviations
   :math:`|x_j - x^{*}|` will be zero, and so the statistic :math:`S_i` will be zero. In such a case, the Hampel
   filter will act like the standard median filter regardless of the value of :math:`t`. Caution should also
   be exercised if dividing by :math:`S_i`.

Because of the possibility of MAD implosion, GSL offers a routine :func:`gsl_filter_hampel2` where
the user can input an additional parameter :data:`epsilon`. This parameter is used as a lower bound
on the :math:`S_i`. So for this function, the filter's response is

.. math:: y_i = \left\{
                  \begin{array}{cc}
                    x_i, & |x_i - m_i| \le t S_i \textrm{ or } S_i < \epsilon \\
                    m_i, & |x_i - m_i| > t S_i
                  \end{array}
                \right.

The function :func:`gsl_filter_hampel` sets :math:`\epsilon = 0`.

.. function:: gsl_filter_hampel_workspace * gsl_filter_hampel_alloc(const size_t K)

   This function initializes a workspace for Hampel filtering using a symmetric moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(6K)`.

.. function:: void gsl_filter_hampel_free(gsl_filter_hampel_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_hampel(const gsl_filter_end_t endtype, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_hampel_workspace * w)
.. function:: int gsl_filter_hampel2(const gsl_filter_end_t endtype, const double epsilon, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_hampel_workspace * w)

   These functions apply a Hampel filter to the input vector :data:`x`, storing the filtered output
   in :data:`y`.  The tuning parameter :math:`t` is provided in :data:`t`. The lower
   bound :math:`\epsilon` for the MAD scale estimates :math:`S_i` is provided in :data:`epsilon`.
   The window medians :math:`m_i` are stored in :data:`xmedian` and the :math:`S_i` are stored in :data:`xsigma` on output.
   The number of outliers detected is stored in :data:`noutlier` on output, while
   the locations of flagged outliers are stored in the boolean array :data:`ioutlier`. The input
   :data:`ioutlier` may be :code:`NULL` if not desired. It  is allowed to have :data:`x` = :data:`y` for an
   in-place filter.
