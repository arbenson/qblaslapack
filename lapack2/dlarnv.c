#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlarnv_(int *idist, int *iseed, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlarnv(int *idist, int *iseed, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlarnv_(int *idist, int *iseed, int *n, 
#endif

	LONG DOUBLE *x)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARNV returns a vector of n random real numbers from a uniform or   
    normal distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  uniform (0,1)   
            = 2:  uniform (-1,1)   
            = 3:  normal (0,1)   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    N       (input) INTEGER   
            The number of random numbers to be generated.   

    X       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The generated random numbers.   

    Further Details   
    ===============   

    This routine calls the auxiliary routine DLARUV to generate random   
    real numbers from a uniform (0,1) distribution, in batches of up to   
    128 using vectorisable code. The Box-Muller method is used to   
    transform numbers from a uniform to a normal distribution.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1, i__2, i__3;
    /* Builtin functions */

    /* Local variables */
    static int i;
    static LONG DOUBLE u[128];
    static int il, iv;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaruv_(int *, int *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaruv(int *, int *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaruv_(int *, int *, LONG DOUBLE *);
#endif

    static int il2;


#define U(I) u[(I)]
#define X(I) x[(I)-1]
#define ISEED(I) iseed[(I)-1]


    i__1 = *n;
    for (iv = 1; iv <= *n; iv += 64) {
/* Computing MIN */
	i__2 = 64, i__3 = *n - iv + 1;
	il = MIN(i__2,i__3);
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

/*        Call DLARUV to generate IL2 numbers from a uniform (0,1)   
          distribution (IL2 <= LV) */


#ifdef PETSC_PREFIX_SUFFIX
	dlaruv_(&ISEED(1), &il2, u);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaruv(&ISEED(1), &il2, u);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaruv_(&ISEED(1), &il2, u);
#endif


	if (*idist == 1) {

/*           Copy generated numbers */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = U(i - 1);
/* L10: */
	    }
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribut
ion */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = U(i - 1) * 2. - 1.;
/* L20: */
	    }
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distributio
n */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = sqrt(log(U((i << 1) - 2)) * -2.) * cos(U((i <<
			 1) - 1) * 6.2831853071795864769252867663);
/* L30: */
	    }
	}
/* L40: */
    }
    return;

/*     End of DLARNV */

} /* dlarnv_ */

