#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dhsein_(char *side, char *eigsrc, char *initv, long int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qhsein(char *side, char *eigsrc, char *initv, long int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qhsein_(char *side, char *eigsrc, char *initv, long int *
#endif

	select, int *n, LONG DOUBLE *h, int *ldh, LONG DOUBLE *wr, 
	LONG DOUBLE *wi, LONG DOUBLE *vl, int *ldvl, LONG DOUBLE *vr, 
	int *ldvr, int *mm, int *m, LONG DOUBLE *work, int *
	ifaill, int *ifailr, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DHSEIN uses inverse iteration to find specified right and/or left   
    eigenvectors of a real upper Hessenberg matrix H.   

    The right eigenvector x and the left eigenvector y of the matrix H   
    corresponding to an eigenvalue w are defined by:   

                 H * x = w * x,     y**h * H = w * y**h   

    where y**h denotes the conjugate transpose of the vector y.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'R': compute right eigenvectors only;   
            = 'L': compute left eigenvectors only;   
            = 'B': compute both right and left eigenvectors.   

    EIGSRC  (input) CHARACTER*1   
            Specifies the source of eigenvalues supplied in (WR,WI):   
            = 'Q': the eigenvalues were found using DHSEQR; thus, if   
                   H has zero subdiagonal elements, and so is   
                   block-triangular, then the j-th eigenvalue can be   
                   assumed to be an eigenvalue of the block containing   
                   the j-th row/column.  This property allows DHSEIN to   
                   perform inverse iteration on just one diagonal block. 
  
            = 'N': no assumptions are made on the correspondence   
                   between eigenvalues and diagonal blocks.  In this   
                   case, DHSEIN must always perform inverse iteration   
                   using the whole matrix H.   

    INITV   (input) CHARACTER*1   
            = 'N': no initial vectors are supplied;   
            = 'U': user-supplied initial vectors are stored in the arrays 
  
                   VL and/or VR.   

    SELECT  (input/output) LOGICAL array, dimension (N)   
            Specifies the eigenvectors to be computed. To select the   
            real eigenvector corresponding to a real eigenvalue WR(j),   
            SELECT(j) must be set to .TRUE.. To select the complex   
            eigenvector corresponding to a complex eigenvalue   
            (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),   
            either SELECT(j) or SELECT(j+1) or both must be set to   
            .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is   
            .FALSE..   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    H       (input) LONG DOUBLE PRECISION array, dimension (LDH,N)   
            The upper Hessenberg matrix H.   

    LDH     (input) INTEGER   
            The leading dimension of the array H.  LDH >= MAX(1,N).   

    WR      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
    WI      (input) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the real and imaginary parts of the eigenvalues of 
  
            H; a complex conjugate pair of eigenvalues must be stored in 
  
            consecutive elements of WR and WI.   
            On exit, WR may have been altered since close eigenvalues   
            are perturbed slightly in searching for independent   
            eigenvectors.   

    VL      (input/output) LONG DOUBLE PRECISION array, dimension (LDVL,MM)   
            On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must   
            contain starting vectors for the inverse iteration for the   
            left eigenvectors; the starting vector for each eigenvector   
            must be in the same column(s) in which the eigenvector will   
            be stored.   
            On exit, if SIDE = 'L' or 'B', the left eigenvectors   
            specified by SELECT will be stored consecutively in the   
            columns of VL, in the same order as their eigenvalues. A   
            complex eigenvector corresponding to a complex eigenvalue is 
  
            stored in two consecutive columns, the first holding the real 
  
            part and the second the imaginary part.   
            If SIDE = 'R', VL is not referenced.   

    LDVL    (input) INTEGER   
            The leading dimension of the array VL.   
            LDVL >= MAX(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.   

    VR      (input/output) LONG DOUBLE PRECISION array, dimension (LDVR,MM)   
            On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must   
            contain starting vectors for the inverse iteration for the   
            right eigenvectors; the starting vector for each eigenvector 
  
            must be in the same column(s) in which the eigenvector will   
            be stored.   
            On exit, if SIDE = 'R' or 'B', the right eigenvectors   
            specified by SELECT will be stored consecutively in the   
            columns of VR, in the same order as their eigenvalues. A   
            complex eigenvector corresponding to a complex eigenvalue is 
  
            stored in two consecutive columns, the first holding the real 
  
            part and the second the imaginary part.   
            If SIDE = 'L', VR is not referenced.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= MAX(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.   

    MM      (input) INTEGER   
            The number of columns in the arrays VL and/or VR. MM >= M.   

    M       (output) INTEGER   
            The number of columns in the arrays VL and/or VR required to 
  
            store the eigenvectors; each selected real eigenvector   
            occupies one column and each selected complex eigenvector   
            occupies two columns.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension ((N+2)*N)   

    IFAILL  (output) INTEGER array, dimension (MM)   
            If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left   
            eigenvector in the i-th column of VL (corresponding to the   
            eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the   
            eigenvector converged satisfactorily. If the i-th and (i+1)th 
  
            columns of VL hold a complex eigenvector, then IFAILL(i) and 
  
            IFAILL(i+1) are set to the same value.   
            If SIDE = 'R', IFAILL is not referenced.   

    IFAILR  (output) INTEGER array, dimension (MM)   
            If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right   
            eigenvector in the i-th column of VR (corresponding to the   
            eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the   
            eigenvector converged satisfactorily. If the i-th and (i+1)th 
  
            columns of VR hold a complex eigenvector, then IFAILR(i) and 
  
            IFAILR(i+1) are set to the same value.   
            If SIDE = 'L', IFAILR is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, i is the number of eigenvectors which   
                  failed to converge; see IFAILL and IFAILR for further   
                  details.   

    Further Details   
    ===============   

    Each eigenvector is normalized so that the element of largest   
    magnitude has magnitude 1; here the magnitude of a complex number   
    (x,y) is taken to be |x|+|y|.   

    ===================================================================== 
  


       Decode and test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static long int c_false = 0;
    static long int c_true = 1;
    
    /* System generated locals */
    int  i__1, 
	    i__2;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static long int pair;
    static LONG DOUBLE unfl;
    static int i, k;
    extern long int lsame_(char *, char *);
    static int iinfo;
    static long int leftv, bothv;
    static LONG DOUBLE hnorm;
    static int kl;

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
    extern /* Subroutine */ void dlaein_(long int *, long int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaein(long int *, long int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaein_(long int *, long int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *
	    , LONG DOUBLE *, LONG DOUBLE *, int *);
    static int kr;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *);
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum;
    static long int noinit;
    static int ldwork;
    static long int rightv, fromqr;
    static LONG DOUBLE smlnum;
    static int kln, ksi;
    static LONG DOUBLE wki;
    static int ksr;
    static LONG DOUBLE ulp, wkr, eps3;



#define SELECT(I) select[(I)-1]
#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]
#define WORK(I) work[(I)-1]
#define IFAILL(I) ifaill[(I)-1]
#define IFAILR(I) ifailr[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    bothv = lsame_(side, "B");
    rightv = lsame_(side, "R") || bothv;
    leftv = lsame_(side, "L") || bothv;

    fromqr = lsame_(eigsrc, "Q");

    noinit = lsame_(initv, "N");

/*     Set M to the number of columns required to store the selected   
       eigenvectors, and standardize the array SELECT. */

    *m = 0;
    pair = 0;
    i__1 = *n;
    for (k = 1; k <= *n; ++k) {
	if (pair) {
	    pair = 0;
	    SELECT(k) = 0;
	} else {
	    if (WI(k) == 0.) {
		if (SELECT(k)) {
		    ++(*m);
		}
	    } else {
		pair = 1;
		if (SELECT(k) || SELECT(k + 1)) {
		    SELECT(k) = 1;
		    *m += 2;
		}
	    }
	}
/* L10: */
    }

    *info = 0;
    if (! rightv && ! leftv) {
	*info = -1;
    } else if (! fromqr && ! lsame_(eigsrc, "N")) {
	*info = -2;
    } else if (! noinit && ! lsame_(initv, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -5;
    } else if (*ldh < MAX(1,*n)) {
	*info = -7;
    } else if (*ldvl < 1 || (leftv && *ldvl < *n)) {
	*info = -11;
    } else if (*ldvr < 1 || (rightv && *ldvr < *n)) {
	*info = -13;
    } else if (*mm < *m) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DHSEIN", &i__1);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return;
    }

/*     Set machine-dependent constants. */


#ifdef PETSC_PREFIX_SUFFIX
    unfl = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    unfl = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    unfl = qlamch_("Safe minimum");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    ulp = dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    ulp = qlamch("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    ulp = qlamch_("Precision");
#endif

    smlnum = unfl * (*n / ulp);
    bignum = (1. - ulp) / smlnum;

    ldwork = *n + 1;

    kl = 1;
    kln = 0;
    if (fromqr) {
	kr = 0;
    } else {
	kr = *n;
    }
    ksr = 1;

    i__1 = *n;
    for (k = 1; k <= *n; ++k) {
	if (SELECT(k)) {

/*           Compute eigenvector(s) corresponding to W(K). */

	    if (fromqr) {

/*              If affiliation of eigenvalues is known, check 
whether   
                the matrix splits.   

                Determine KL and KR such that 1 <= KL <= K <= 
KR <= N   
                and H(KL,KL-1) and H(KR+1,KR) are zero (or KL 
= 1 or   
                KR = N).   

                Then inverse iteration can be performed with t
he   
                submatrix H(KL:N,KL:N) for a left eigenvector,
 and with   
                the submatrix H(1:KR,1:KR) for a right eigenve
ctor. */

		i__2 = kl + 1;
		for (i = k; i >= kl+1; --i) {
		    if (H(i,i-1) == 0.) {
			goto L30;
		    }
/* L20: */
		}
L30:
		kl = i;
		if (k > kr) {
		    i__2 = *n - 1;
		    for (i = k; i <= *n-1; ++i) {
			if (H(i+1,i) == 0.) {
			    goto L50;
			}
/* L40: */
		    }
L50:
		    kr = i;
		}
	    }

	    if (kl != kln) {
		kln = kl;

/*              Compute infinity-norm of submatrix H(KL:KR,KL:
KR) if it   
                has not ben computed before. */

		i__2 = kr - kl + 1;

#ifdef PETSC_PREFIX_SUFFIX
		hnorm = dlanhs_("I", &i__2, &H(kl,kl), ldh, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		hnorm = qlanhs("I", &i__2, &H(kl,kl), ldh, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		hnorm = qlanhs_("I", &i__2, &H(kl,kl), ldh, &WORK(
#endif

			1));
		if (hnorm > 0.) {
		    eps3 = hnorm * ulp;
		} else {
		    eps3 = smlnum;
		}
	    }

/*           Perturb eigenvalue if it is close to any previous   
             selected eigenvalues affiliated to the submatrix   
             H(KL:KR,KL:KR). Close roots are modified by EPS3. */

	    wkr = WR(k);
	    wki = WI(k);
L60:
	    i__2 = kl;
	    for (i = k - 1; i >= kl; --i) {
		if (SELECT(i) && (d__1 = WR(i) - wkr, ABS(d__1)) + (d__2 = WI(
			i) - wki, ABS(d__2)) < eps3) {
		    wkr += eps3;
		    goto L60;
		}
/* L70: */
	    }
	    WR(k) = wkr;

	    pair = wki != 0.;
	    if (pair) {
		ksi = ksr + 1;
	    } else {
		ksi = ksr;
	    }
	    if (leftv) {

/*              Compute left eigenvector. */

		i__2 = *n - kl + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlaein_(&c_false, &noinit, &i__2, &H(kl,kl), ldh, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaein(&c_false, &noinit, &i__2, &H(kl,kl), ldh, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaein_(&c_false, &noinit, &i__2, &H(kl,kl), ldh, &
#endif

			wkr, &wki, &VL(kl,ksr), &VL(kl,ksi), &WORK(1), &ldwork, &WORK(*n * *n + *n + 1), 
			&eps3, &smlnum, &bignum, &iinfo);
		if (iinfo > 0) {
		    if (pair) {
			*info += 2;
		    } else {
			++(*info);
		    }
		    IFAILL(ksr) = k;
		    IFAILL(ksi) = k;
		} else {
		    IFAILL(ksr) = 0;
		    IFAILL(ksi) = 0;
		}
		i__2 = kl - 1;
		for (i = 1; i <= kl-1; ++i) {
		    VL(i,ksr) = 0.;
/* L80: */
		}
		if (pair) {
		    i__2 = kl - 1;
		    for (i = 1; i <= kl-1; ++i) {
			VL(i,ksi) = 0.;
/* L90: */
		    }
		}
	    }
	    if (rightv) {

/*              Compute right eigenvector. */


#ifdef PETSC_PREFIX_SUFFIX
		dlaein_(&c_true, &noinit, &kr, &H(1,1), ldh, &wkr, &wki, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaein(&c_true, &noinit, &kr, &H(1,1), ldh, &wkr, &wki, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaein_(&c_true, &noinit, &kr, &H(1,1), ldh, &wkr, &wki, 
#endif

			&VR(1,ksr), &VR(1,ksi), &WORK(
			1), &ldwork, &WORK(*n * *n + *n + 1), &eps3, &smlnum, 
			&bignum, &iinfo);
		if (iinfo > 0) {
		    if (pair) {
			*info += 2;
		    } else {
			++(*info);
		    }
		    IFAILR(ksr) = k;
		    IFAILR(ksi) = k;
		} else {
		    IFAILR(ksr) = 0;
		    IFAILR(ksi) = 0;
		}
		i__2 = *n;
		for (i = kr + 1; i <= *n; ++i) {
		    VR(i,ksr) = 0.;
/* L100: */
		}
		if (pair) {
		    i__2 = *n;
		    for (i = kr + 1; i <= *n; ++i) {
			VR(i,ksi) = 0.;
/* L110: */
		    }
		}
	    }

	    if (pair) {
		ksr += 2;
	    } else {
		++ksr;
	    }
	}
/* L120: */
    }

    return;

/*     End of DHSEIN */

} /* dhsein_ */

