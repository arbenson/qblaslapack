#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>
#include <stdio.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlamch_(char *cmach)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlamch(char *cmach)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlamch_(char *cmach)
#endif

{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMCH determines LONG DOUBLE precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DLAMCH:   
            = 'E' or 'e',   DLAMCH := eps   
            = 'S' or 's ,   DLAMCH := sfmin   
            = 'B' or 'b',   DLAMCH := base   
            = 'P' or 'p',   DLAMCH := eps*base   
            = 'N' or 'n',   DLAMCH := t   
            = 'R' or 'r',   DLAMCH := rnd   
            = 'M' or 'm',   DLAMCH := emin   
            = 'U' or 'u',   DLAMCH := rmin   
            = 'L' or 'l',   DLAMCH := emax   
            = 'O' or 'o',   DLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static long int first = 1;
    /* System generated locals */
    int i__1;
    LONG DOUBLE ret_val;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE base;
    static int beta;
    static LONG DOUBLE emin, prec, emax;
    static int imin, imax;
    static long int lrnd;
    static LONG DOUBLE rmin, rmax, t, rmach;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE small, sfmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlamc2_(int *, int *, long int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc2(int *, int *, long int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc2_(int *, int *, long int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static int it;
    static LONG DOUBLE rnd, eps;



    if (first) {
	first = 0;

#ifdef PETSC_PREFIX_SUFFIX
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc2(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
#endif

	base = (LONG DOUBLE) beta;
	t = (LONG DOUBLE) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow(base, (LONG DOUBLE)i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow(base, (LONG DOUBLE)i__1);
	}
	prec = eps * base;
	emin = (LONG DOUBLE) imin;
	emax = (LONG DOUBLE) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }

    if (lsame_(cmach, "E")) {
	rmach = eps;
    } else if (lsame_(cmach, "S")) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B")) {
	rmach = base;
    } else if (lsame_(cmach, "P")) {
	rmach = prec;
    } else if (lsame_(cmach, "N")) {
	rmach = t;
    } else if (lsame_(cmach, "R")) {
	rmach = rnd;
    } else if (lsame_(cmach, "M")) {
	rmach = emin;
    } else if (lsame_(cmach, "U")) {
	rmach = rmin;
    } else if (lsame_(cmach, "L")) {
	rmach = emax;
    } else if (lsame_(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlamc1_(int *beta, int *t, long int *rnd, long int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlamc1(int *beta, int *t, long int *rnd, long int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlamc1_(int *beta, int *t, long int *rnd, long int 
#endif

	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static long int first = 1;
    /* System generated locals */
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static long int lrnd;
    static LONG DOUBLE a, b, c, f;
    static int lbeta;
    static LONG DOUBLE savec;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static long int lieee1;
    static LONG DOUBLE t1, t2;
    static int lt;
    static LONG DOUBLE one, qtr;



    if (first) {
	first = 0;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive int m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&a, &one);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&a, &one);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&a, &one);
#endif

	    d__1 = -a;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&c, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&c, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&c, &d__1);
#endif

	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive int 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.;

#ifdef PETSC_PREFIX_SUFFIX
	c = dlamc3_(&a, &b);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c = qlamc3(&a, &b);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c = qlamc3_(&a, &b);
#endif


/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&a, &b);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&a, &b);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&a, &b);
#endif

	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	d__1 = -a;

#ifdef PETSC_PREFIX_SUFFIX
	c = dlamc3_(&c, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c = qlamc3(&c, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c = qlamc3_(&c, &d__1);
#endif

	lbeta = (int) (c + qtr);

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (LONG DOUBLE) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;

#ifdef PETSC_PREFIX_SUFFIX
	f = dlamc3_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	f = qlamc3(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	f = qlamc3_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	c = dlamc3_(&f, &a);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c = qlamc3(&f, &a);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c = qlamc3_(&f, &a);
#endif

	if (c == a) {
	    lrnd = 1;
	} else {
	    lrnd = 0;
	}
	d__1 = b / 2;
	d__2 = b / 100;

#ifdef PETSC_PREFIX_SUFFIX
	f = dlamc3_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	f = qlamc3(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	f = qlamc3_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	c = dlamc3_(&f, &a);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c = qlamc3(&f, &a);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c = qlamc3_(&f, &a);
#endif

	if (lrnd && c == a) {
	    lrnd = 0;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;

#ifdef PETSC_PREFIX_SUFFIX
	t1 = dlamc3_(&d__1, &a);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	t1 = qlamc3(&d__1, &a);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	t1 = qlamc3_(&d__1, &a);
#endif

	d__1 = b / 2;

#ifdef PETSC_PREFIX_SUFFIX
	t2 = dlamc3_(&d__1, &savec);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	t2 = qlamc3(&d__1, &savec);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	t2 = qlamc3_(&d__1, &savec);
#endif

	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  int part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive int 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&a, &one);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&a, &one);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&a, &one);
#endif

	    d__1 = -a;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&c, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&c, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&c, &d__1);
#endif

	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return;

/*     End of DLAMC1 */

} /* dlamc1_ */



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlamc2_(int *beta, int *t, long int *rnd, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlamc2(int *beta, int *t, long int *rnd, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlamc2_(int *beta, int *t, long int *rnd, 
#endif

	LONG DOUBLE *eps, int *emin, LONG DOUBLE *rmin, int *emax, 
	LONG DOUBLE *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) LONG DOUBLE PRECISION   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) LONG DOUBLE PRECISION   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) LONG DOUBLE PRECISION   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    
    /* Initialized data */
    static long int first = 1;
    static long int iwarn = 0;
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    /* Local variables */
    static long int ieee;
    static LONG DOUBLE half;
    static long int lrnd;
    static LONG DOUBLE leps, zero, a, b, c;
    static int i, lbeta;
    static LONG DOUBLE rbase;
    static int lemin, lemax, gnmin;
    static LONG DOUBLE small;
    static int gpmin;
    static LONG DOUBLE third, lrmin, lrmax, sixth;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlamc1_(int *, int *, long int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc1(int *, int *, long int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc1_(int *, int *, long int *, 
#endif

	    long int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static long int lieee1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlamc4_(int *, LONG DOUBLE *, int *), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc4(int *, LONG DOUBLE *, int *), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamc4_(int *, LONG DOUBLE *, int *), 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dlamc5_(int *, int *, int *, long int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlamc5(int *, int *, int *, long int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlamc5_(int *, int *, int *, long int *, int *, 
#endif

	    LONG DOUBLE *);
    static int lt, ngnmin, ngpmin;
    static LONG DOUBLE one, two;



    if (first) {
	first = 0;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/


#ifdef PETSC_PREFIX_SUFFIX
	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc1(&lbeta, &lt, &lrnd, &lieee1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc1_(&lbeta, &lt, &lrnd, &lieee1);
#endif


/*        Start to find EPS. */

	b = (LONG DOUBLE) lbeta;
	i__1 = -lt;
	a = pow(b, (LONG DOUBLE)i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;

#ifdef PETSC_PREFIX_SUFFIX
	sixth = dlamc3_(&b, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	sixth = qlamc3(&b, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	sixth = qlamc3_(&b, &d__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	third = dlamc3_(&sixth, &sixth);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	third = qlamc3(&sixth, &sixth);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	third = qlamc3_(&sixth, &sixth);
#endif

	d__1 = -half;

#ifdef PETSC_PREFIX_SUFFIX
	b = dlamc3_(&third, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	b = qlamc3(&third, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	b = qlamc3_(&third, &d__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	b = dlamc3_(&b, &sixth);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	b = qlamc3(&b, &sixth);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	b = qlamc3_(&b, &sixth);
#endif

	b = ABS(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&d__1, &d__2);
#endif

	    d__1 = -c;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&half, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&half, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&half, &d__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    b = dlamc3_(&half, &c);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    b = qlamc3(&half, &c);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    b = qlamc3_(&half, &c);
#endif

	    d__1 = -b;

#ifdef PETSC_PREFIX_SUFFIX
	    c = dlamc3_(&half, &d__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    c = qlamc3(&half, &d__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    c = qlamc3_(&half, &d__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    b = dlamc3_(&half, &c);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    b = qlamc3(&half, &c);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    b = qlamc3_(&half, &c);
#endif

	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    d__1 = small * rbase;

#ifdef PETSC_PREFIX_SUFFIX
	    small = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    small = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    small = qlamc3_(&d__1, &zero);
#endif

/* L20: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	a = dlamc3_(&one, &small);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	a = qlamc3(&one, &small);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	a = qlamc3_(&one, &small);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlamc4_(&ngpmin, &one, &lbeta);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc4(&ngpmin, &one, &lbeta);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc4_(&ngpmin, &one, &lbeta);
#endif

	d__1 = -one;

#ifdef PETSC_PREFIX_SUFFIX
	dlamc4_(&ngnmin, &d__1, &lbeta);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc4(&ngnmin, &d__1, &lbeta);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc4_(&ngnmin, &d__1, &lbeta);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlamc4_(&gpmin, &a, &lbeta);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc4(&gpmin, &a, &lbeta);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc4_(&gpmin, &a, &lbeta);
#endif

	d__1 = -a;

#ifdef PETSC_PREFIX_SUFFIX
	dlamc4_(&gnmin, &d__1, &lbeta);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc4(&gnmin, &d__1, &lbeta);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc4_(&gnmin, &d__1, &lbeta);
#endif

	ieee = 0;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = 1;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = MIN(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, ABS(i__1)) == 1) {
		lemin = MAX(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = MIN(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else if ((i__1 = ngpmin - ngnmin, ABS(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - MIN(ngpmin,ngnmin) == 3) {
		lemin = MAX(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = MIN(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else {
/* Computing MIN */
	    i__1 = MIN(ngpmin,ngnmin), i__1 = MIN(i__1,gpmin);
	    lemin = MIN(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = 1;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = 1;
	    fprintf(stdout,"\n\n WARNING. The value EMIN may be incorrect:- ");
	    fprintf(stdout,"EMIN = %8i\n",lemin);
	    fprintf(stdout,"If, after inspection, the value EMIN looks acceptable");
            fprintf(stdout,"please comment out \n the IF block as marked within the"); 
            fprintf(stdout,"code of routine DLAMC2, \n otherwise supply EMIN"); 
            fprintf(stdout,"explicitly.\n");
            fflush(stdout);
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine DLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    d__1 = lrmin * rbase;

#ifdef PETSC_PREFIX_SUFFIX
	    lrmin = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    lrmin = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    lrmin = qlamc3_(&d__1, &zero);
#endif

/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */


#ifdef PETSC_PREFIX_SUFFIX
	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamc5(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
#endif

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return;


/*     End of DLAMC2 */

} /* dlamc2_ */



#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlamc3_(LONG DOUBLE *a, LONG DOUBLE *b)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlamc3(LONG DOUBLE *a, LONG DOUBLE *b)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlamc3_(LONG DOUBLE *a, LONG DOUBLE *b)
#endif

{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) LONG DOUBLE PRECISION   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    LONG DOUBLE ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlamc4_(int *emin, LONG DOUBLE *start, int *base)
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlamc4(int *emin, LONG DOUBLE *start, int *base)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlamc4_(int *emin, LONG DOUBLE *start, int *base)
#endif

{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC4 is a service routine for DLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) LONG DOUBLE PRECISION   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Local variables */
    static LONG DOUBLE zero, a;
    static int i;
    static LONG DOUBLE rbase, b1, b2, c1, c2, d1, d2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE one;



    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;

#ifdef PETSC_PREFIX_SUFFIX
    b1 = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    b1 = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    b1 = qlamc3_(&d__1, &zero);
#endif

    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;

#ifdef PETSC_PREFIX_SUFFIX
	b1 = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	b1 = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	b1 = qlamc3_(&d__1, &zero);
#endif

	d__1 = b1 * *base;

#ifdef PETSC_PREFIX_SUFFIX
	c1 = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c1 = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c1 = qlamc3_(&d__1, &zero);
#endif

	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;

#ifdef PETSC_PREFIX_SUFFIX
	b2 = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	b2 = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	b2 = qlamc3_(&d__1, &zero);
#endif

	d__1 = b2 / rbase;

#ifdef PETSC_PREFIX_SUFFIX
	c2 = dlamc3_(&d__1, &zero);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	c2 = qlamc3(&d__1, &zero);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	c2 = qlamc3_(&d__1, &zero);
#endif

	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return;

/*     End of DLAMC4 */

} /* dlamc4_ */



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlamc5_(int *beta, int *p, int *emin, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlamc5(int *beta, int *p, int *emin, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlamc5_(int *beta, int *p, int *emin, 
#endif

	long int *ieee, int *emax, LONG DOUBLE *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + ABS(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A long int flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) LONG DOUBLE PRECISION   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       ABS(EMIN). We then assume that EMAX + ABS(EMIN) will sum   
       approximately to the bound that is closest to ABS(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static LONG DOUBLE c_b5 = 0.;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Local variables */
    static int lexp;
    static LONG DOUBLE oldy;
    static int uexp, i;
    static LONG DOUBLE y, z;
    static int nbits;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE recbas;
    static int exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1. / *beta;
    z = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.) {
	    oldy = y;
	}

#ifdef PETSC_PREFIX_SUFFIX
	y = dlamc3_(&y, &z);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	y = qlamc3(&y, &z);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	y = qlamc3_(&y, &z);
#endif

/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= *emax; ++i) {
	d__1 = y * *beta;

#ifdef PETSC_PREFIX_SUFFIX
	y = dlamc3_(&d__1, &c_b5);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	y = qlamc3(&d__1, &c_b5);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	y = qlamc3_(&d__1, &c_b5);
#endif

/* L30: */
    }

    *rmax = y;
    return;

/*     End of DLAMC5 */

} /* dlamc5_ */

