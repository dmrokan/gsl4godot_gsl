noinst_LTLIBRARIES = libgslmovstat.la

pkginclude_HEADERS = gsl_movstat.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslmovstat_la_SOURCES = \
  alloc.c                  \
	medacc.c                 \
	mmacc.c                  \
	movmad.c                 \
	movmean.c                \
	movmedian.c              \
	movminmax.c              \
	movsum.c                 \
	movSn.c                  \
	movQn.c                  \
	movqqr.c                 \
	movvariance.c            \
	mvacc.c                  \
	qnacc.c                  \
	snacc.c                  \
	sumacc.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmovstat.la ../statistics/libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la ../vector/libgslvector.la ../block/libgslblock.la


