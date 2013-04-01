#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dstein_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qstein(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qstein_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif

	int *m, LONG DOUBLE *w, int *iblock, int *isplit, 
	LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int *iwork, 
	int *ifail, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEIN computes the eigenvectors of a real symmetric tridiagonal   
    matrix T corresponding to specified eigenvalues, using inverse   
    iteration.   

    The maximum number of iterations allowed for each eigenvector is   
    specified by an internal parameter MAXITS (currently set to 5).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix T.   

    E       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The (n-1) subdiagonal elements of the tridiagonal matrix   
            T, in elements 1 to N-1.  E(N) need not be set.   

    M       (input) INTEGER   
            The number of eigenvectors to be found.  0 <= M <= N.   

    W       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The first M elements of W contain the eigenvalues for   
            which eigenvectors are to be computed.  The eigenvalues   
            should be grouped by split-off block and ordered from   
            smallest to largest within the block.  ( The output array   
            W from DSTEBZ with ORDER = 'B' is expected here. )   

    IBLOCK  (input) INTEGER array, dimension (N)   
            The submatrix indices associated with the corresponding   
            eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
            the first submatrix from the top, =2 if W(i) belongs to   
            the second submatrix, etc.  ( The output array IBLOCK   
            from DSTEBZ is expected here. )   

    ISPLIT  (input) INTEGER array, dimension (N)   
            The splitting points, at which T breaks up into submatrices. 
  
            The first submatrix consists of rows/columns 1 to   
            ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
            through ISPLIT( 2 ), etc.   
            ( The output array ISPLIT from DSTEBZ is expected here. )   

    Z       (output) LONG DOUBLE PRECISION array, dimension (LDZ, M)   
            The computed eigenvectors.  The eigenvector associated   
            with the eigenvalue W(i) is stored in the i-th column of   
            Z.  Any vector which fails to converge is set to its current 
  
            iterate after MAXITS iterations.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    IFAIL   (output) INTEGER array, dimension (M)   
            On normal exit, all elements of IFAIL are zero.   
            If one or more eigenvectors fail to converge after   
            MAXITS iterations, then their indices are stored in   
            array IFAIL.   

    INFO    (output) INTEGER   
            = 0: successful exit.   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, then i eigenvectors failed to converge   
                 in MAXITS iterations.  Their indices are stored in   
                 array IFAIL.   

    Internal Parameters   
    ===================   

    MAXITS  INTEGER, default = 5   
            The maximum number of iterations performed.   

    EXTRA   INTEGER, default = 2   
            The number of iterations performed after norm growth   
            criterion is satisfied, should be at least 1.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__2 = 2;
    static int c__1 = 1;
    static int c_n1 = -1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    /* Local variables */
    static int jblk, nblk;

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
    static int jmax;

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
    static int iseed[4], gpind, iinfo;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dasum_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qasum(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qasum_(int *, LONG DOUBLE *, int *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static int b1;

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
    static int P_j1;
    static LONG DOUBLE ortol;
    static int indrv1, indrv2, indrv3, indrv4, indrv5, bn;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlagtf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlagtf(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlagtf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *
	    , int *);
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


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dlagts_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlagts(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlagts_(
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int nrmchk;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarnv_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarnv(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarnv_(int *, int *, int *, 
#endif

	    LONG DOUBLE *);
    static int blksiz;
    static LONG DOUBLE onenrm, dtpcrt, pertol, scl, eps, sep, nrm, tol;
    static int its;
    static LONG DOUBLE xjm, ztr, eps1;



#define ISEED(I) iseed[(I)]
#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define W(I) w[(I)-1]
#define IBLOCK(I) iblock[(I)-1]
#define ISPLIT(I) isplit[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]
#define IFAIL(I) ifail[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;
    i__1 = *m;
    for (i = 1; i <= *m; ++i) {
	IFAIL(i) = 0;
/* L10: */
    }

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -4;
    } else if (*ldz < MAX(1,*n)) {
	*info = -9;
    } else {
	i__1 = *m;
	for (j = 2; j <= *m; ++j) {
	    if (IBLOCK(j) < IBLOCK(j - 1)) {
		*info = -6;
		goto L30;
	    }
	    if (IBLOCK(j) == IBLOCK(j - 1) && W(j) < W(j - 1)) {
		*info = -5;
		goto L30;
	    }
/* L20: */
	}
L30:
	;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEIN", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return;
    } else if (*n == 1) {
	Z(1,1) = 1.;
	return;
    }

/*     Get machine constants. */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Precision");
#endif


/*     Initialize seed for random number generator DLARNV. */

    for (i = 1; i <= 4; ++i) {
	ISEED(i - 1) = 1;
/* L40: */
    }

/*     Initialize pointers. */

    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

    P_j1 = 1;
    i__1 = IBLOCK(*m);
    for (nblk = 1; nblk <= IBLOCK(*m); ++nblk) {

/*        Find starting and ending indices of block nblk. */

	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = ISPLIT(nblk - 1) + 1;
	}
	bn = ISPLIT(nblk);
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    goto L60;
	}
	gpind = b1;

/*        Compute reorthogonalization criterion and stopping criterion
. */

	onenrm = (d__1 = D(b1), ABS(d__1)) + (d__2 = E(b1), ABS(d__2));
/* Computing MAX */
	d__3 = onenrm, d__4 = (d__1 = D(bn), ABS(d__1)) + (d__2 = E(bn - 1), 
		ABS(d__2));
	onenrm = MAX(d__3,d__4);
	i__2 = bn - 1;
	for (i = b1 + 1; i <= bn-1; ++i) {
/* Computing MAX */
	    d__4 = onenrm, d__5 = (d__1 = D(i), ABS(d__1)) + (d__2 = E(i - 1),
		     ABS(d__2)) + (d__3 = E(i), ABS(d__3));
	    onenrm = MAX(d__4,d__5);
/* L50: */
	}
	ortol = onenrm * .001;

	dtpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

L60:
	jblk = 0;
	i__2 = *m;
	for (j = P_j1; j <= *m; ++j) {
	    if (IBLOCK(j) != nblk) {
		P_j1 = j;
		goto L160;
	    }
	    ++jblk;
	    xj = W(j);

/*           Skip all the work if the block size is one. */

	    if (blksiz == 1) {
		WORK(indrv1 + 1) = 1.;
		goto L120;
	    }

/*           If eigenvalues j and j-1 are too close, add a relativ
ely   
             small perturbation. */

	    if (jblk > 1) {
		eps1 = (d__1 = eps * xj, ABS(d__1));
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

/*           Get random starting vector. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarnv_(&c__2, iseed, &blksiz, &WORK(indrv1 + 1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarnv(&c__2, iseed, &blksiz, &WORK(indrv1 + 1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarnv_(&c__2, iseed, &blksiz, &WORK(indrv1 + 1));
#endif


/*           Copy the matrix T so it won't be destroyed in factori
zation. */


#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&blksiz, &D(b1), &c__1, &WORK(indrv4 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&blksiz, &D(b1), &c__1, &WORK(indrv4 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&blksiz, &D(b1), &c__1, &WORK(indrv4 + 1), &c__1);
#endif

	    i__3 = blksiz - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__3, &E(b1), &c__1, &WORK(indrv2 + 2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__3, &E(b1), &c__1, &WORK(indrv2 + 2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__3, &E(b1), &c__1, &WORK(indrv2 + 2), &c__1);
#endif

	    i__3 = blksiz - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__3, &E(b1), &c__1, &WORK(indrv3 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__3, &E(b1), &c__1, &WORK(indrv3 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__3, &E(b1), &c__1, &WORK(indrv3 + 1), &c__1);
#endif


/*           Compute LU factors with partial pivoting  ( PT = LU )
 */

	    tol = 0.;

#ifdef PETSC_PREFIX_SUFFIX
	    dlagtf_(&blksiz, &WORK(indrv4 + 1), &xj, &WORK(indrv2 + 2), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlagtf(&blksiz, &WORK(indrv4 + 1), &xj, &WORK(indrv2 + 2), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlagtf_(&blksiz, &WORK(indrv4 + 1), &xj, &WORK(indrv2 + 2), &WORK(
#endif

		    indrv3 + 1), &tol, &WORK(indrv5 + 1), &IWORK(1), &iinfo);

/*           Update iteration count. */

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

/*           Normalize and scale the righthand side vector Pb.   

   Computing MAX */
	    d__2 = eps, d__3 = (d__1 = WORK(indrv4 + blksiz), ABS(d__1));

#ifdef PETSC_PREFIX_SUFFIX
	    scl = blksiz * onenrm * MAX(d__2,d__3) / dasum_(&blksiz, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    scl = blksiz * onenrm * MAX(d__2,d__3) / qasum(&blksiz, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    scl = blksiz * onenrm * MAX(d__2,d__3) / qasum_(&blksiz, &WORK(
#endif

		    indrv1 + 1), &c__1);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif


/*           Solve the system LU = Pb. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlagts_(&c_n1, &blksiz, &WORK(indrv4 + 1), &WORK(indrv2 + 2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlagts(&c_n1, &blksiz, &WORK(indrv4 + 1), &WORK(indrv2 + 2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlagts_(&c_n1, &blksiz, &WORK(indrv4 + 1), &WORK(indrv2 + 2), &
#endif

		    WORK(indrv3 + 1), &WORK(indrv5 + 1), &IWORK(1), &WORK(
		    indrv1 + 1), &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalu
es are   
             close enough. */

	    if (jblk == 1) {
		goto L90;
	    }
	    if ((d__1 = xj - xjm, ABS(d__1)) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i = gpind; i <= j-1; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
		    ztr = -ddot_(&blksiz, &WORK(indrv1 + 1), &c__1, &Z(b1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    ztr = -qdot(&blksiz, &WORK(indrv1 + 1), &c__1, &Z(b1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    ztr = -qdot_(&blksiz, &WORK(indrv1 + 1), &c__1, &Z(b1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&blksiz, &ztr, &Z(b1,i), &c__1, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&blksiz, &ztr, &Z(b1,i), &c__1, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&blksiz, &ztr, &Z(b1,i), &c__1, &WORK(
#endif

			    indrv1 + 1), &c__1);
/* L80: */
		}
	    }

/*           Check the infinity norm of the iterate. */

L90:

#ifdef PETSC_PREFIX_SUFFIX
	    jmax = idamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    jmax = iqamax(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    jmax = iqamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif

	    nrm = (d__1 = WORK(indrv1 + jmax), ABS(d__1));

/*           Continue for additional iterations after norm reaches
   
             stopping criterion. */

	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L110;

/*           If stopping criterion was not satisfied, update info 
and   
             store eigenvector number in array ifail. */

L100:
	    ++(*info);
	    IFAIL(*info) = j;

/*           Accept iterate as jth eigenvector. */

L110:

#ifdef PETSC_PREFIX_SUFFIX
	    scl = 1. / dnrm2_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    scl = 1. / qnrm2(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    scl = 1. / qnrm2_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    jmax = idamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    jmax = iqamax(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    jmax = iqamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
#endif

	    if (WORK(indrv1 + jmax) < 0.) {
		scl = -scl;
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
#endif

L120:
	    i__3 = *n;
	    for (i = 1; i <= *n; ++i) {
		Z(i,j) = 0.;
/* L130: */
	    }
	    i__3 = blksiz;
	    for (i = 1; i <= blksiz; ++i) {
		Z(b1+i-1,j) = WORK(indrv1 + i);
/* L140: */
	    }

/*           Save the shift to check eigenvalue spacing at next   
             iteration. */

	    xjm = xj;

/* L150: */
	}
L160:
	;
    }

    return;

/*     End of DSTEIN */

} /* dstein_ */

