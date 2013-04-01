#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaein_(long int *rightv, long int *noinit, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaein(long int *rightv, long int *noinit, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaein_(long int *rightv, long int *noinit, int *n, 
#endif

	LONG DOUBLE *h, int *ldh, LONG DOUBLE *wr, LONG DOUBLE *wi, 
	LONG DOUBLE *vr, LONG DOUBLE *vi, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *work, LONG DOUBLE *eps3, LONG DOUBLE *smlnum, LONG DOUBLE *
	bignum, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAEIN uses inverse iteration to find a right or left eigenvector   
    corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg   
    matrix H.   

    Arguments   
    =========   

    RIGHTV   (input) LOGICAL   
            = .TRUE. : compute right eigenvector;   
            = .FALSE.: compute left eigenvector.   

    NOINIT   (input) LOGICAL   
            = .TRUE. : no initial vector supplied in (VR,VI).   
            = .FALSE.: initial vector supplied in (VR,VI).   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    H       (input) LONG DOUBLE PRECISION array, dimension (LDH,N)   
            The upper Hessenberg matrix H.   

    LDH     (input) INTEGER   
            The leading dimension of the array H.  LDH >= MAX(1,N).   

    WR      (input) LONG DOUBLE PRECISION   
    WI      (input) LONG DOUBLE PRECISION   
            The real and imaginary parts of the eigenvalue of H whose   
            corresponding right or left eigenvector is to be computed.   

    VR      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
    VI      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain   
            a real starting vector for inverse iteration using the real   
            eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI   
            must contain the real and imaginary parts of a complex   
            starting vector for inverse iteration using the complex   
            eigenvalue (WR,WI); otherwise VR and VI need not be set.   
            On exit, if WI = 0.0 (real eigenvalue), VR contains the   
            computed real eigenvector; if WI.ne.0.0 (complex eigenvalue), 
  
            VR and VI contain the real and imaginary parts of the   
            computed complex eigenvector. The eigenvector is normalized   
            so that the component of largest magnitude has magnitude 1;   
            here the magnitude of a complex number (x,y) is taken to be   
            |x| + |y|.   
            VI is not referenced if WI = 0.0.   

    B       (workspace) LONG DOUBLE PRECISION array, dimension (LDB,N)   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= N+1.   

    WORK   (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    EPS3    (input) LONG DOUBLE PRECISION   
            A small machine-dependent value which is used to perturb   
            close eigenvalues, and to replace zero pivots.   

    SMLNUM  (input) LONG DOUBLE PRECISION   
            A machine-dependent value close to the underflow threshold.   

    BIGNUM  (input) LONG DOUBLE PRECISION   
            A machine-dependent value close to the overflow threshold.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            = 1:  inverse iteration did not converge; VR is set to the   
                  last iterate, and so is VI if WI.ne.0.0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static int ierr;
    static LONG DOUBLE temp, norm, vmax;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

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
    static LONG DOUBLE scale, w, x, y;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dasum_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qasum(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qasum_(int *, LONG DOUBLE *, int *);
#endif

    static char trans[1];
    static LONG DOUBLE vcrit;
    static int i1, i2, i3;
    static LONG DOUBLE rootn, vnorm, w1;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE ei, ej, absbii, absbjj, xi;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE xr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatrs_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrs(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrs_(char *, char *, char *, char *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *);
    static char normin[1];
    static LONG DOUBLE nrmsml, growto, rec;
    static int its;



#define VR(I) vr[(I)-1]
#define VI(I) vi[(I)-1]
#define WORK(I) work[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an   
       eigenvector. */

    rootn = sqrt((LONG DOUBLE) (*n));
    growto = .1 / rootn;
/* Computing MAX */
    d__1 = 1., d__2 = *eps3 * rootn;
    nrmsml = MAX(d__1,d__2) * *smlnum;

/*     Form B = H - (WR,WI)*I (except that the subdiagonal elements and   
       the imaginary parts of the diagonal elements are not stored). */

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    B(i,j) = H(i,j);
/* L10: */
	}
	B(j,j) = H(j,j) - *wr;
/* L20: */
    }

    if (*wi == 0.) {

/*        Real eigenvalue. */

	if (*noinit) {

/*           Set initial vector. */

	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		VR(i) = *eps3;
/* L30: */
	    }
	} else {

/*           Scale supplied initial vector. */


#ifdef PETSC_PREFIX_SUFFIX
	    vnorm = dnrm2_(n, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    vnorm = qnrm2(n, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    vnorm = qnrm2_(n, &VR(1), &c__1);
#endif

	    d__1 = *eps3 * rootn / MAX(vnorm,nrmsml);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(n, &d__1, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(n, &d__1, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(n, &d__1, &VR(1), &c__1);
#endif

	}

	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3. */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		ei = H(i+1,i);
		if ((d__1 = B(i,i), ABS(d__1)) < ABS(ei)) {

/*                 Interchange rows and eliminate. */

		    x = B(i,i) / ei;
		    B(i,i) = ei;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			temp = B(i+1,j);
			B(i+1,j) = B(i,j) - x * temp;
			B(i,j) = temp;
/* L40: */
		    }
		} else {

/*                 Eliminate without interchange. */

		    if (B(i,i) == 0.) {
			B(i,i) = *eps3;
		    }
		    x = ei / B(i,i);
		    if (x != 0.) {
			i__2 = *n;
			for (j = i + 1; j <= *n; ++j) {
			    B(i+1,j) -= x * B(i,j);
/* L50: */
			}
		    }
		}
/* L60: */
	    }
	    if (B(*n,*n) == 0.) {
		B(*n,*n) = *eps3;
	    }

	    *(unsigned char *)trans = 'N';

	} else {

/*           UL decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3. */

	    for (j = *n; j >= 2; --j) {
		ej = H(j,j-1);
		if ((d__1 = B(j,j), ABS(d__1)) < ABS(ej)) {

/*                 Interchange columns and eliminate. */

		    x = B(j,j) / ej;
		    B(j,j) = ej;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp = B(i,j-1);
			B(i,j-1) = B(i,j) - x * 
				temp;
			B(i,j) = temp;
/* L70: */
		    }
		} else {

/*                 Eliminate without interchange. */

		    if (B(j,j) == 0.) {
			B(j,j) = *eps3;
		    }
		    x = ej / B(j,j);
		    if (x != 0.) {
			i__1 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    B(i,j-1) -= x * B(i,j);
/* L80: */
			}
		    }
		}
/* L90: */
	    }
	    if (B(1,1) == 0.) {
		B(1,1) = *eps3;
	    }

	    *(unsigned char *)trans = 'T';

	}

	*(unsigned char *)normin = 'N';
	i__1 = *n;
	for (its = 1; its <= *n; ++its) {

/*           Solve U*x = scale*v for a right eigenvector   
               or U'*x = scale*v for a left eigenvector,   
             overwriting x on v. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatrs_("Upper", trans, "Nonunit", normin, n, &B(1,1), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrs("Upper", trans, "Nonunit", normin, n, &B(1,1), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrs_("Upper", trans, "Nonunit", normin, n, &B(1,1), ldb, &
#endif

		    VR(1), &scale, &WORK(1), &ierr);
	    *(unsigned char *)normin = 'Y';

/*           Test for sufficient growth in the norm of v. */


#ifdef PETSC_PREFIX_SUFFIX
	    vnorm = dasum_(n, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    vnorm = qasum(n, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    vnorm = qasum_(n, &VR(1), &c__1);
#endif

	    if (vnorm >= growto * scale) {
		goto L120;
	    }

/*           Choose new orthogonal starting vector and try again. 
*/

	    temp = *eps3 / (rootn + 1.);
	    VR(1) = *eps3;
	    i__2 = *n;
	    for (i = 2; i <= *n; ++i) {
		VR(i) = temp;
/* L100: */
	    }
	    VR(*n - its + 1) -= *eps3 * rootn;
/* L110: */
	}

/*        Failure to find eigenvector in N iterations. */

	*info = 1;

L120:

/*        Normalize eigenvector. */


#ifdef PETSC_PREFIX_SUFFIX
	i = idamax_(n, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	i = iqamax(n, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	i = iqamax_(n, &VR(1), &c__1);
#endif

	d__2 = 1. / (d__1 = VR(i), ABS(d__1));

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__2, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__2, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__2, &VR(1), &c__1);
#endif

    } else {

/*        Complex eigenvalue. */

	if (*noinit) {

/*           Set initial vector. */

	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		VR(i) = *eps3;
		VI(i) = 0.;
/* L130: */
	    }
	} else {

/*           Scale supplied initial vector. */


#ifdef PETSC_PREFIX_SUFFIX
	    d__1 = dnrm2_(n, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    d__1 = qnrm2(n, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    d__1 = qnrm2_(n, &VR(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    d__2 = dnrm2_(n, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    d__2 = qnrm2(n, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    d__2 = qnrm2_(n, &VI(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    norm = dlapy2_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    norm = qlapy2(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    norm = qlapy2_(&d__1, &d__2);
#endif

	    rec = *eps3 * rootn / MAX(norm,nrmsml);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(n, &rec, &VR(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(n, &rec, &VI(1), &c__1);
#endif

	}

	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3.   

             The imaginary part of the (i,j)-th element of U is st
ored in   
             B(j+1,i). */

	    B(2,1) = -(*wi);
	    i__1 = *n;
	    for (i = 2; i <= *n; ++i) {
		B(i+1,1) = 0.;
/* L140: */
	    }

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
		absbii = dlapy2_(&B(i,i), &B(i+1,i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		absbii = qlapy2(&B(i,i), &B(i+1,i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		absbii = qlapy2_(&B(i,i), &B(i+1,i));
#endif

		ei = H(i+1,i);
		if (absbii < ABS(ei)) {

/*                 Interchange rows and eliminate. */

		    xr = B(i,i) / ei;
		    xi = B(i+1,i) / ei;
		    B(i,i) = ei;
		    B(i+1,i) = 0.;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			temp = B(i+1,j);
			B(i+1,j) = B(i,j) - xr * temp;
			B(j+1,i+1) = B(j+1,i) - 
				xi * temp;
			B(i,j) = temp;
			B(j+1,i) = 0.;
/* L150: */
		    }
		    B(i+2,i) = -(*wi);
		    B(i+1,i+1) -= xi * *wi;
		    B(i+2,i+1) += xr * *wi;
		} else {

/*                 Eliminate without interchanging rows. 
*/

		    if (absbii == 0.) {
			B(i,i) = *eps3;
			B(i+1,i) = 0.;
			absbii = *eps3;
		    }
		    ei = ei / absbii / absbii;
		    xr = B(i,i) * ei;
		    xi = -B(i+1,i) * ei;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			B(i+1,j) = B(i+1,j) - xr * 
				B(i,j) + xi * B(j+1,i)
				;
			B(j+1,i+1) = -xr * B(j+1,i) - xi * B(i,j);
/* L160: */
		    }
		    B(i+2,i+1) -= *wi;
		}

/*              Compute 1-norm of offdiagonal elements of i-th
 row. */

		i__2 = *n - i;
		i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		WORK(i) = dasum_(&i__2, &B(i,i+1), ldb) + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		WORK(i) = qasum(&i__2, &B(i,i+1), ldb) + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		WORK(i) = qasum_(&i__2, &B(i,i+1), ldb) + 
#endif


#ifdef PETSC_PREFIX_SUFFIX
			dasum_(&i__3, &B(i+2,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qasum(&i__3, &B(i+2,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qasum_(&i__3, &B(i+2,i), &c__1);
#endif

/* L170: */
	    }
	    if (B(*n,*n) == 0. && B(*n+1,*n) == 0.) {
		B(*n,*n) = *eps3;
	    }
	    WORK(*n) = 0.;

	    i1 = *n;
	    i2 = 1;
	    i3 = -1;
	} else {

/*           UL decomposition with partial pivoting of conjg(B), 
  
             replacing zero pivots by EPS3.   

             The imaginary part of the (i,j)-th element of U is st
ored in   
             B(j+1,i). */

	    B(*n+1,*n) = *wi;
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		B(*n+1,j) = 0.;
/* L180: */
	    }

	    for (j = *n; j >= 2; --j) {
		ej = H(j,j-1);

#ifdef PETSC_PREFIX_SUFFIX
		absbjj = dlapy2_(&B(j,j), &B(j+1,j));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		absbjj = qlapy2(&B(j,j), &B(j+1,j));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		absbjj = qlapy2_(&B(j,j), &B(j+1,j));
#endif

		if (absbjj < ABS(ej)) {

/*                 Interchange columns and eliminate */

		    xr = B(j,j) / ej;
		    xi = B(j+1,j) / ej;
		    B(j,j) = ej;
		    B(j+1,j) = 0.;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp = B(i,j-1);
			B(i,j-1) = B(i,j) - xr * 
				temp;
			B(j,i) = B(j+1,i) - xi * temp;
			B(i,j) = temp;
			B(j+1,i) = 0.;
/* L190: */
		    }
		    B(j+1,j-1) = *wi;
		    B(j-1,j-1) += xi * *wi;
		    B(j,j-1) -= xr * *wi;
		} else {

/*                 Eliminate without interchange. */

		    if (absbjj == 0.) {
			B(j,j) = *eps3;
			B(j+1,j) = 0.;
			absbjj = *eps3;
		    }
		    ej = ej / absbjj / absbjj;
		    xr = B(j,j) * ej;
		    xi = -B(j+1,j) * ej;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			B(i,j-1) = B(i,j-1) - 
				xr * B(i,j) + xi * B(j+1,i);
			B(j,i) = -xr * B(j+1,i) - xi *
				 B(i,j);
/* L200: */
		    }
		    B(j,j-1) += *wi;
		}

/*              Compute 1-norm of offdiagonal elements of j-th
 column. */

		i__1 = j - 1;
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		WORK(j) = dasum_(&i__1, &B(1,j), &c__1) + dasum_(&
#endif
#ifdef Q_C_PREFIX_SUFFIX
		WORK(j) = qasum(&i__1, &B(1,j), &c__1) + dasum_(&
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		WORK(j) = qasum_(&i__1, &B(1,j), &c__1) + dasum_(&
#endif

			i__2, &B(j+1,1), ldb);
/* L210: */
	    }
	    if (B(1,1) == 0. && B(2,1) == 0.) {
		B(1,1) = *eps3;
	    }
	    WORK(1) = 0.;

	    i1 = 1;
	    i2 = *n;
	    i3 = 1;
	}

	i__1 = *n;
	for (its = 1; its <= *n; ++its) {
	    scale = 1.;
	    vmax = 1.;
	    vcrit = *bignum;

/*           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvect
or,   
               or U'*(xr,xi) = scale*(vr,vi) for a left eigenvecto
r,   
             overwriting (xr,xi) on (vr,vi). */

	    i__2 = i2;
	    i__3 = i3;
	    for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {

		if (WORK(i) > vcrit) {
		    rec = 1. / vmax;

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(n, &rec, &VR(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(n, &rec, &VI(1), &c__1);
#endif

		    scale *= rec;
		    vmax = 1.;
		    vcrit = *bignum;
		}

		xr = VR(i);
		xi = VI(i);
		if (*rightv) {
		    i__4 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			xr = xr - B(i,j) * VR(j) + B(j+1,i) * VI(j);
			xi = xi - B(i,j) * VI(j) - B(j+1,i) * VR(j);
/* L220: */
		    }
		} else {
		    i__4 = i - 1;
		    for (j = 1; j <= i-1; ++j) {
			xr = xr - B(j,i) * VR(j) + B(i+1,j) * VI(j);
			xi = xi - B(j,i) * VI(j) - B(i+1,j) * VR(j);
/* L230: */
		    }
		}

		w = (d__1 = B(i,i), ABS(d__1)) + (d__2 = B(i+1,i), ABS(d__2));
		if (w > *smlnum) {
		    if (w < 1.) {
			w1 = ABS(xr) + ABS(xi);
			if (w1 > w * *bignum) {
			    rec = 1. / w1;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &VR(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &VI(1), &c__1);
#endif

			    xr = VR(i);
			    xi = VI(i);
			    scale *= rec;
			    vmax *= rec;
			}
		    }

/*                 Divide by diagonal element of B. */


#ifdef PETSC_PREFIX_SUFFIX
		    dladiv_(&xr, &xi, &B(i,i), &B(i+1,i), &VR(i), &VI(i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qladiv(&xr, &xi, &B(i,i), &B(i+1,i), &VR(i), &VI(i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qladiv_(&xr, &xi, &B(i,i), &B(i+1,i), &VR(i), &VI(i));
#endif

/* Computing MAX */
		    d__3 = (d__1 = VR(i), ABS(d__1)) + (d__2 = VI(i), ABS(
			    d__2));
		    vmax = MAX(d__3,vmax);
		    vcrit = *bignum / vmax;
		} else {
		    i__4 = *n;
		    for (j = 1; j <= *n; ++j) {
			VR(j) = 0.;
			VI(j) = 0.;
/* L240: */
		    }
		    VR(i) = 1.;
		    VI(i) = 1.;
		    scale = 0.;
		    vmax = 1.;
		    vcrit = *bignum;
		}
/* L250: */
	    }

/*           Test for sufficient growth in the norm of (VR,VI). */


#ifdef PETSC_PREFIX_SUFFIX
	    vnorm = dasum_(n, &VR(1), &c__1) + dasum_(n, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    vnorm = qasum(n, &VR(1), &c__1) + dasum_(n, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    vnorm = qasum_(n, &VR(1), &c__1) + dasum_(n, &VI(1), &c__1);
#endif

	    if (vnorm >= growto * scale) {
		goto L280;
	    }

/*           Choose a new orthogonal starting vector and try again
. */

	    y = *eps3 / (rootn + 1.);
	    VR(1) = *eps3;
	    VI(1) = 0.;

	    i__3 = *n;
	    for (i = 2; i <= *n; ++i) {
		VR(i) = y;
		VI(i) = 0.;
/* L260: */
	    }
	    VR(*n - its + 1) -= *eps3 * rootn;
/* L270: */
	}

/*        Failure to find eigenvector in N iterations */

	*info = 1;

L280:

/*        Normalize eigenvector. */

	vnorm = 0.;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__3 = vnorm, d__4 = (d__1 = VR(i), ABS(d__1)) + (d__2 = VI(i), 
		    ABS(d__2));
	    vnorm = MAX(d__3,d__4);
/* L290: */
	}
	d__1 = 1. / vnorm;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__1, &VR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__1, &VR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__1, &VR(1), &c__1);
#endif

	d__1 = 1. / vnorm;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__1, &VI(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__1, &VI(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__1, &VI(1), &c__1);
#endif


    }

    return;

/*     End of DLAEIN */

} /* dlaein_ */

