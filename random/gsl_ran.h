/* $Id$ */
#ifndef _gsl_ran_RANDOM_H_
#define _gsl_ran_RANDOM_H_

double gsl_ran_uniform(void);
void gsl_ran_seed(int);

unsigned long gsl_ran_random(void);
double gsl_ran_max();

void *gsl_ran_getRandomState(void);
void gsl_ran_setRandomState(void *);

unsigned long gsl_ran_random_wstate(void *);
double gsl_ran_uniform_wstate(void *);
void gsl_ran_seed_wstate(void *, int);

/* Note: it's kind of ugly to have the gaussian random state defined
 * so explicitly, if we used void *'s we would lose some typechecking
 * ability, but it would look better.  maybe if we just called it
 * gsl_ran_gState, it wouldn't be so bad...  */
typedef struct {
    int ng;
    double g;
    void *randomState;
} gsl_ran_gaussianRandomState;

void gsl_ran_copyGaussState(gsl_ran_gaussianRandomState *,
                            gsl_ran_gaussianRandomState *);
gsl_ran_gaussianRandomState *gsl_ran_getGaussState(void);
void gsl_ran_setGaussState(gsl_ran_gaussianRandomState *);
double gsl_ran_gaussian_wstate(gsl_ran_gaussianRandomState *);
double gsl_ran_gaussian();

#endif /* _gsl_ran_RANDOM_H_ */
