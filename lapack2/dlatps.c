#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlatps_(char *uplo, char *trans, char *diag, char *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlatps(char *uplo, char *trans, char *diag, char *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlatps_(char *uplo, char *trans, char *diag, char *
#endif

	normin, int *n, LONG DOUBLE *ap, LONG DOUBLE *x, LONG DOUBLE *scale, 
	LONG DOUBLE *cnorm, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DLATPS solves one of the triangular systems   

       A *x = s*b  or  A'*x = s*b   

    with scaling to prevent overflow, where A is an upper or lower   
    triangular matrix stored in packed form.  Here A' denotes the   
    transpose of A, x and b are n-element vectors, and s is a scaling   
    factor, usually less than or equal to 1, chosen so that the   
    components of x will be less than the overflow threshold.  If the   
    unscaled problem will not cause overflow, the Level 2 BLAS routine   
    DTPSV is called. If the matrix A is singular (A(j,j) = 0 for some j), 
  
    then s is set to 0 and a non-trivial solution to A*x = 0 is returned. 
  

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    TRANS   (input) CHARACTER*1   
            Specifies the operation applied to A.   
            = 'N':  Solve A * x = s*b  (No transpose)   
            = 'T':  Solve A'* x = s*b  (Transpose)   
            = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose) 
  

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    NORMIN  (input) CHARACTER*1   
            Specifies whether CNORM has been set or not.   
            = 'Y':  CNORM contains the column norms on entry   
            = 'N':  CNORM is not set on entry.  On exit, the norms will   
                    be computed and stored in CNORM.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the right hand side b of the triangular system.   
            On exit, X is overwritten by the solution vector x.   

    SCALE   (output) LONG DOUBLE PRECISION   
            The scaling factor s for the triangular system   
               A * x = s*b  or  A'* x = s*b.   
            If SCALE = 0, the matrix A is singular or badly scaled, and   
            the vector x is an exact or approximate solution to A*x = 0. 
  

    CNORM   (input or output) LONG DOUBLE PRECISION array, dimension (N)   

            If NORMIN = 'Y', CNORM is an input argument and CNORM(j)   
            contains the norm of the off-diagonal part of the j-th column 
  
            of A.  If TRANS = 'N', CNORM(j) must be greater than or equal 
  
            to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)   
            must be greater than or equal to the 1-norm.   

            If NORMIN = 'N', CNORM is an output argument and CNORM(j)   
            returns the 1-norm of the offdiagonal part of the j-th column 
  
            of A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -k, the k-th argument had an illegal value   

    Further Details   
    ======= =======   

    A rough bound on x is computed; if that is less than overflow, DTPSV 
  
    is called, otherwise, specific code is used which checks for possible 
  
    overflow or divide-by-zero at every operation.   

    A columnwise scheme is used for solving A*x = b.  The basic algorithm 
  
    if A is lower triangular is   

         x[1:n] := b[1:n]   
         for j = 1, ..., n   
              x(j) := x(j) / A(j,j)   
              x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]   
         end   

    Define bounds on the components of x after j iterations of the loop: 
  
       M(j) = bound on x[1:j]   
       G(j) = bound on x[j+1:n]   
    Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.   

    Then for iteration j+1 we have   
       M(j+1) <= G(j) / | A(j+1,j+1) |   
       G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |   
              <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )   

    where CNORM(j+1) is greater than or equal to the infinity-norm of   
    column j+1 of A, not counting the diagonal.  Hence   

       G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )   
                    1<=i<=j   
    and   

       |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) 
  
                                     1<=i< j   

    Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTPSV if the   
    reciprocal of the largest M(j), j=1,..,n, is larger than   
    MAX(underflow, 1/overflow).   

    The bound on x(j) is also used to determine when a step in the   
    columnwise method can be performed without fear of overflow.  If   
    the computed bound is greater than a large constant, x is scaled to   
    prevent overflow, but if the bound overflows, x is set to 0, x(j) to 
  
    1, and scale to 0, and a non-trivial solution to A*x = 0 is found.   

    Similarly, a row-wise scheme is used to solve A'*x = b.  The basic   
    algorithm for A upper triangular is   

         for j = 1, ..., n   
              x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)   
         end   

    We simultaneously compute two bounds   
         G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j   
         M(j) = bound on x(i), 1<=i<=j   

    The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we   
    add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.   
    Then the bound on x(j) is   

         M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |   

              <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )   
                        1<=i<=j   

    and we can safely call DTPSV if 1/M(n) and 1/G(n) are both greater   
    than MAX(underflow, 1/overflow).   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b36 = .5;
    
    /* System generated locals */
    int i__1, i__2, i__3;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int jinc, jlen;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE xbnd;
    static int imax;
    static LONG DOUBLE tmax, tjjs, xmax, grow, sumj;
    static int i, j;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);
    static LONG DOUBLE tscal, uscal;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dasum_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qasum(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qasum_(int *, LONG DOUBLE *, int *);
#endif

    static int jlast;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtpsv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv_(char *, char *, char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int ip;
    static LONG DOUBLE xj;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum;
    static long int notran;
    static int jfirst;
    static LONG DOUBLE smlnum;
    static long int nounit;
    static LONG DOUBLE rec, tjj;



#define CNORM(I) cnorm[(I)-1]
#define X(I) x[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");

/*     Test the input parameters. */

    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (! lsame_(normin, "Y") && ! lsame_(normin, "N"))
	     {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLATPS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Determine machine dependent parameters to control overflow. */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("Safe minimum") / dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("Safe minimum") / dlamch_("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("Safe minimum") / dlamch_("Precision");
#endif

    bignum = 1. / smlnum;
    *scale = 1.;

    if (lsame_(normin, "N")) {

/*        Compute the 1-norm of each column, not including the diagona
l. */

	if (upper) {

/*           A is upper triangular. */

	    ip = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		CNORM(j) = dasum_(&i__2, &AP(ip), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		CNORM(j) = qasum(&i__2, &AP(ip), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		CNORM(j) = qasum_(&i__2, &AP(ip), &c__1);
#endif

		ip += j;
/* L10: */
	    }
	} else {

/*           A is lower triangular. */

	    ip = 1;
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		CNORM(j) = dasum_(&i__2, &AP(ip + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		CNORM(j) = qasum(&i__2, &AP(ip + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		CNORM(j) = qasum_(&i__2, &AP(ip + 1), &c__1);
#endif

		ip = ip + *n - j + 1;
/* L20: */
	    }
	    CNORM(*n) = 0.;
	}
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is 
  
       greater than BIGNUM. */


#ifdef PETSC_PREFIX_SUFFIX
    imax = idamax_(n, &CNORM(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    imax = iqamax(n, &CNORM(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    imax = iqamax_(n, &CNORM(1), &c__1);
#endif

    tmax = CNORM(imax);
    if (tmax <= bignum) {
	tscal = 1.;
    } else {
	tscal = 1. / (smlnum * tmax);

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &tscal, &CNORM(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &tscal, &CNORM(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &tscal, &CNORM(1), &c__1);
#endif

    }

/*     Compute a bound on the computed solution vector to see if the   
       Level 2 BLAS routine DTPSV can be used. */


#ifdef PETSC_PREFIX_SUFFIX
    j = idamax_(n, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    j = iqamax(n, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    j = iqamax_(n, &X(1), &c__1);
#endif

    xmax = (d__1 = X(j), ABS(d__1));
    xbnd = xmax;
    if (notran) {

/*        Compute the growth in A * x = b. */

	if (upper) {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	} else {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L50;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, G(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / MAX(xbnd,smlnum);
	    xbnd = grow;
	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = *n;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              M(j) = G(j-1) / ABS(A(j,j)) */

		tjj = (d__1 = AP(ip), ABS(d__1));
/* Computing MIN */
		d__1 = xbnd, d__2 = MIN(1.,tjj) * grow;
		xbnd = MIN(d__1,d__2);
		if (tjj + CNORM(j) >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / ABS(A(j,
j)) ) */

		    grow *= tjj / (tjj + CNORM(j));
		} else {

/*                 G(j) could overflow, set GROW to 0. */

		    grow = 0.;
		}
		ip += jinc * jlen;
		--jlen;
/* L30: */
	    }
	    grow = xbnd;
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}.   

   Computing MIN */
	    d__1 = 1., d__2 = 1. / MAX(xbnd,smlnum);
	    grow = MIN(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

		grow *= 1. / (CNORM(j) + 1.);
/* L40: */
	    }
	}
L50:

	;
    } else {

/*        Compute the growth in A' * x = b. */

	if (upper) {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	} else {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L80;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, M(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / MAX(xbnd,smlnum);
	    xbnd = grow;
	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = 1;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = MAX( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) 
*/

		xj = CNORM(j) + 1.;
/* Computing MIN */
		d__1 = grow, d__2 = xbnd / xj;
		grow = MIN(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / ABS(A(j,j)) 
*/

		tjj = (d__1 = AP(ip), ABS(d__1));
		if (xj > tjj) {
		    xbnd *= tjj / xj;
		}
		++jlen;
		ip += jinc * jlen;
/* L60: */
	    }
	    grow = MIN(grow,xbnd);
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}.   

   Computing MIN */
	    d__1 = 1., d__2 = 1. / MAX(xbnd,smlnum);
	    grow = MIN(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

		xj = CNORM(j) + 1.;
		grow /= xj;
/* L70: */
	    }
	}
L80:
	;
    }

    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on
   
          elements of X is not too small. */


#ifdef PETSC_PREFIX_SUFFIX
	dtpsv_(uplo, trans, diag, n, &AP(1), &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtpsv(uplo, trans, diag, n, &AP(1), &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtpsv_(uplo, trans, diag, n, &AP(1), &X(1), &c__1);
#endif

    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal
 to   
             BIGNUM in absolute value. */

	    *scale = bignum / xmax;

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(n, scale, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(n, scale, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(n, scale, &X(1), &c__1);
#endif

	    xmax = bignum;
	}

	if (notran) {

/*           Solve A * x = b */

	    ip = jfirst * (jfirst + 1) / 2;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if nec
essary. */

		xj = (d__1 = X(j), ABS(d__1));
		if (nounit) {
		    tjjs = AP(ip) * tscal;
		} else {
		    tjjs = tscal;
		    if (tscal == 1.) {
			goto L100;
		    }
		}
		tjj = ABS(tjjs);
		if (tjj > smlnum) {

/*                    ABS(A(j,j)) > SMLNUM: */

		    if (tjj < 1.) {
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

			    rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    X(j) /= tjjs;
		    xj = (d__1 = X(j), ABS(d__1));
		} else if (tjj > 0.) {

/*                    0 < ABS(A(j,j)) <= SMLNUM: */

		    if (xj > tjj * bignum) {

/*                       Scale x by (1/ABS(x(j)))*ABS(
A(j,j))*BIGNUM   
                         to avoid overflow when dividi
ng by A(j,j). */

			rec = tjj * bignum / xj;
			if (CNORM(j) > 1.) {

/*                          Scale by 1/CNORM(j) to
 avoid overflow when   
                            multiplying x(j) times
 column j. */

			    rec /= CNORM(j);
			}

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(n, &rec, &X(1), &c__1);
#endif

			*scale *= rec;
			xmax *= rec;
		    }
		    X(j) /= tjjs;
		    xj = (d__1 = X(j), ABS(d__1));
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 
1, and   
                      scale = 0, and compute a solution to
 A*x = 0. */

		    i__3 = *n;
		    for (i = 1; i <= *n; ++i) {
			X(i) = 0.;
/* L90: */
		    }
		    X(j) = 1.;
		    xj = 1.;
		    *scale = 0.;
		    xmax = 0.;
		}
L100:

/*              Scale x if necessary to avoid overflow when ad
ding a   
                multiple of column j of A. */

		if (xj > 1.) {
		    rec = 1. / xj;
		    if (CNORM(j) > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*ABS(x(j))). */

			rec *= .5;

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(n, &rec, &X(1), &c__1);
#endif

			*scale *= rec;
		    }
		} else if (xj * CNORM(j) > bignum - xmax) {

/*                 Scale x by 1/2. */


#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(n, &c_b36, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(n, &c_b36, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(n, &c_b36, &X(1), &c__1);
#endif

		    *scale *= .5;
		}

		if (upper) {
		    if (j > 1) {

/*                    Compute the update   
                         x(1:j-1) := x(1:j-1) - x(j) *
 A(1:j-1,j) */

			i__3 = j - 1;
			d__1 = -X(j) * tscal;

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__3, &d__1, &AP(ip - j + 1), &c__1, &X(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__3, &d__1, &AP(ip - j + 1), &c__1, &X(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__3, &d__1, &AP(ip - j + 1), &c__1, &X(1), &
#endif

				c__1);
			i__3 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
			i = idamax_(&i__3, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			i = iqamax(&i__3, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			i = iqamax_(&i__3, &X(1), &c__1);
#endif

			xmax = (d__1 = X(i), ABS(d__1));
		    }
		    ip -= j;
		} else {
		    if (j < *n) {

/*                    Compute the update   
                         x(j+1:n) := x(j+1:n) - x(j) *
 A(j+1:n,j) */

			i__3 = *n - j;
			d__1 = -X(j) * tscal;

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__3, &d__1, &AP(ip + 1), &c__1, &X(j + 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__3, &d__1, &AP(ip + 1), &c__1, &X(j + 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__3, &d__1, &AP(ip + 1), &c__1, &X(j + 1), &
#endif

				c__1);
			i__3 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
			i = j + idamax_(&i__3, &X(j + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			i = j + iqamax(&i__3, &X(j + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			i = j + iqamax_(&i__3, &X(j + 1), &c__1);
#endif

			xmax = (d__1 = X(i), ABS(d__1));
		    }
		    ip = ip + *n - j + 1;
		}
/* L110: */
	    }

	} else {

/*           Solve A' * x = b */

	    ip = jfirst * (jfirst + 1) / 2;
	    jlen = 1;
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k).   
                                      k<>j */

		xj = (d__1 = X(j), ABS(d__1));
		uscal = tscal;
		rec = 1. / MAX(xmax,1.);
		if (CNORM(j) > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2
*XMAX). */

		    rec *= .5;
		    if (nounit) {
			tjjs = AP(ip) * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = ABS(tjjs);
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling
 x if A(j,j) > 1.   

   Computing MIN */
			d__1 = 1., d__2 = rec * tjj;
			rec = MIN(d__1,d__2);
			uscal /= tjjs;
		    }
		    if (rec < 1.) {

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(n, &rec, &X(1), &c__1);
#endif

			*scale *= rec;
			xmax *= rec;
		    }
		}

		sumj = 0.;
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot 
product is 1,   
                   call DDOT to perform the dot product. 
*/

		    if (upper) {
			i__3 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
			sumj = ddot_(&i__3, &AP(ip - j + 1), &c__1, &X(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			sumj = qdot(&i__3, &AP(ip - j + 1), &c__1, &X(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			sumj = qdot_(&i__3, &AP(ip - j + 1), &c__1, &X(1), &
#endif

				c__1);
		    } else if (j < *n) {
			i__3 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
			sumj = ddot_(&i__3, &AP(ip + 1), &c__1, &X(j + 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			sumj = qdot(&i__3, &AP(ip + 1), &c__1, &X(j + 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			sumj = qdot_(&i__3, &AP(ip + 1), &c__1, &X(j + 1), &
#endif

				c__1);
		    }
		} else {

/*                 Otherwise, use in-line code for the dot
 product. */

		    if (upper) {
			i__3 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    sumj += AP(ip - j + i) * uscal * X(i);
/* L120: */
			}
		    } else if (j < *n) {
			i__3 = *n - j;
			for (i = 1; i <= *n-j; ++i) {
			    sumj += AP(ip + i) * uscal * X(j + i);
/* L130: */
			}
		    }
		}

		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j
) if 1/A(j,j)   
                   was not used to scale the dotproduct. 
*/

		    X(j) -= sumj;
		    xj = (d__1 = X(j), ABS(d__1));
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), sc
aling if necessary. */

			tjjs = AP(ip) * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == 1.) {
			    goto L150;
			}
		    }
		    tjj = ABS(tjjs);
		    if (tjj > smlnum) {

/*                       ABS(A(j,j)) > SMLNUM: */

			if (tjj < 1.) {
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/ab
s(x(j)). */

				rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
				dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
				qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
				qscal_(n, &rec, &X(1), &c__1);
#endif

				*scale *= rec;
				xmax *= rec;
			    }
			}
			X(j) /= tjjs;
		    } else if (tjj > 0.) {

/*                       0 < ABS(A(j,j)) <= SMLNUM: */

			if (xj > tjj * bignum) {

/*                          Scale x by (1/ABS(x(j)
))*ABS(A(j,j))*BIGNUM. */

			    rec = tjj * bignum / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
			X(j) /= tjjs;
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, 
x(j) = 1, and   
                         scale = 0, and compute a solu
tion to A'*x = 0. */

			i__3 = *n;
			for (i = 1; i <= *n; ++i) {
			    X(i) = 0.;
/* L140: */
			}
			X(j) = 1.;
			*scale = 0.;
			xmax = 0.;
		    }
L150:
		    ;
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj i
f the dot   
                   product has already been divided by 1/A
(j,j). */

		    X(j) = X(j) / tjjs - sumj;
		}
/* Computing MAX */
		d__2 = xmax, d__3 = (d__1 = X(j), ABS(d__1));
		xmax = MAX(d__2,d__3);
		++jlen;
		ip += jinc * jlen;
/* L160: */
	    }
	}
	*scale /= tscal;
    }

/*     Scale the column norms by 1/TSCAL for return. */

    if (tscal != 1.) {
	d__1 = 1. / tscal;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__1, &CNORM(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__1, &CNORM(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__1, &CNORM(1), &c__1);
#endif

    }

    return;

/*     End of DLATPS */

} /* dlatps_ */

