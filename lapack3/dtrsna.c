#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrsna_(char *job, char *howmny, long int *select, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrsna(char *job, char *howmny, long int *select, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrsna_(char *job, char *howmny, long int *select, 
#endif

	int *n, LONG DOUBLE *t, int *ldt, LONG DOUBLE *vl, int *
	ldvl, LONG DOUBLE *vr, int *ldvr, LONG DOUBLE *s, LONG DOUBLE *sep, 
	int *mm, int *m, LONG DOUBLE *work, int *ldwork, int *
	iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTRSNA estimates reciprocal condition numbers for specified   
    eigenvalues and/or right eigenvectors of a real upper   
    quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q   
    orthogonal).   

    T must be in Schur canonical form (as returned by DHSEQR), that is,   
    block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each   
    2-by-2 diagonal block has its diagonal elements equal and its   
    off-diagonal elements of opposite sign.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies whether condition numbers are required for   
            eigenvalues (S) or eigenvectors (SEP):   
            = 'E': for eigenvalues only (S);   
            = 'V': for eigenvectors only (SEP);   
            = 'B': for both eigenvalues and eigenvectors (S and SEP).   

    HOWMNY  (input) CHARACTER*1   
            = 'A': compute condition numbers for all eigenpairs;   
            = 'S': compute condition numbers for selected eigenpairs   
                   specified by the array SELECT.   

    SELECT  (input) LOGICAL array, dimension (N)   
            If HOWMNY = 'S', SELECT specifies the eigenpairs for which   
            condition numbers are required. To select condition numbers   
            for the eigenpair corresponding to a real eigenvalue w(j),   
            SELECT(j) must be set to .TRUE.. To select condition numbers 
  
            corresponding to a complex conjugate pair of eigenvalues w(j) 
  
            and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be 
  
            set to .TRUE..   
            If HOWMNY = 'A', SELECT is not referenced.   

    N       (input) INTEGER   
            The order of the matrix T. N >= 0.   

    T       (input) LONG DOUBLE PRECISION array, dimension (LDT,N)   
            The upper quasi-triangular matrix T, in Schur canonical form. 
  

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= MAX(1,N).   

    VL      (input) LONG DOUBLE PRECISION array, dimension (LDVL,M)   
            If JOB = 'E' or 'B', VL must contain left eigenvectors of T   
            (or of any Q*T*Q**T with Q orthogonal), corresponding to the 
  
            eigenpairs specified by HOWMNY and SELECT. The eigenvectors   
            must be stored in consecutive columns of VL, as returned by   
            DHSEIN or DTREVC.   
            If JOB = 'V', VL is not referenced.   

    LDVL    (input) INTEGER   
            The leading dimension of the array VL.   
            LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.   

    VR      (input) LONG DOUBLE PRECISION array, dimension (LDVR,M)   
            If JOB = 'E' or 'B', VR must contain right eigenvectors of T 
  
            (or of any Q*T*Q**T with Q orthogonal), corresponding to the 
  
            eigenpairs specified by HOWMNY and SELECT. The eigenvectors   
            must be stored in consecutive columns of VR, as returned by   
            DHSEIN or DTREVC.   
            If JOB = 'V', VR is not referenced.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.   

    S       (output) LONG DOUBLE PRECISION array, dimension (MM)   
            If JOB = 'E' or 'B', the reciprocal condition numbers of the 
  
            selected eigenvalues, stored in consecutive elements of the   
            array. For a complex conjugate pair of eigenvalues two   
            consecutive elements of S are set to the same value. Thus   
            S(j), SEP(j), and the j-th columns of VL and VR all   
            correspond to the same eigenpair (but not in general the   
            j-th eigenpair, unless all eigenpairs are selected).   
            If JOB = 'V', S is not referenced.   

    SEP     (output) LONG DOUBLE PRECISION array, dimension (MM)   
            If JOB = 'V' or 'B', the estimated reciprocal condition   
            numbers of the selected eigenvectors, stored in consecutive   
            elements of the array. For a complex eigenvector two   
            consecutive elements of SEP are set to the same value. If   
            the eigenvalues cannot be reordered to compute SEP(j), SEP(j) 
  
            is set to 0; this can only occur when the true value would be 
  
            very small anyway.   
            If JOB = 'E', SEP is not referenced.   

    MM      (input) INTEGER   
            The number of elements in the arrays S (if JOB = 'E' or 'B') 
  
             and/or SEP (if JOB = 'V' or 'B'). MM >= M.   

    M       (output) INTEGER   
            The number of elements of the arrays S and/or SEP actually   
            used to store the estimated condition numbers.   
            If HOWMNY = 'A', M is set to N.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LDWORK,N+1)   
            If JOB = 'E', WORK is not referenced.   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.   

    IWORK   (workspace) INTEGER array, dimension (N)   
            If JOB = 'E', IWORK is not referenced.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The reciprocal of the condition number of an eigenvalue lambda is   
    defined as   

            S(lambda) = |v'*u| / (norm(u)*norm(v))   

    where u and v are the right and left eigenvectors of T corresponding 
  
    to lambda; v' denotes the conjugate-transpose of v, and norm(u)   
    denotes the Euclidean norm. These reciprocal condition numbers always 
  
    lie between zero (very badly conditioned) and one (very well   
    conditioned). If n = 1, S(lambda) is defined to be 1.   

    An approximate error bound for a computed eigenvalue W(i) is given by 
  

                        EPS * norm(T) / S(i)   

    where EPS is the machine precision.   

    The reciprocal of the condition number of the right eigenvector u   
    corresponding to lambda is defined as follows. Suppose   

                T = ( lambda  c  )   
                    (   0    T22 )   

    Then the reciprocal condition number is   

            SEP( lambda, T22 ) = sigma-MIN( T22 - lambda*I )   

    where sigma-min denotes the smallest singular value. We approximate   
    the smallest singular value by the reciprocal of an estimate of the   
    one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is   
    defined to be ABS(T(1,1)).   

    An approximate error bound for a computed right eigenvector VR(i)   
    is given by   

                        EPS * norm(T) / SEP(i)   

    ===================================================================== 
  


       Decode and test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static long int c_true = 1;
    static long int c_false = 0;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static int kase;
    static LONG DOUBLE cond;

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
    static long int pair;
    static int ierr;
    static LONG DOUBLE dumm, prod;
    static int ifst;
    static LONG DOUBLE lnrm;
    static int ilst;
    static LONG DOUBLE rnrm;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE prod1, prod2;
    static int i, j, k;
    static LONG DOUBLE scale, delta;
    extern long int lsame_(char *, char *);
    static long int wants;
    static LONG DOUBLE dummy[1];
    static int n2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE cs;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int nn, ks;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, int *);
    static LONG DOUBLE sn, mu;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);
    static LONG DOUBLE bignum;
    static long int wantbh;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaqtr_(long int *, long int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqtr(long int *, long int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqtr_(long int *, long int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, int *), dtrexc_(char *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, int *), qtrexc(char *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, int *), qtrexc_(char *, int *
#endif

	    , LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *, 
	    int *, LONG DOUBLE *, int *);
    static long int somcon;
    static LONG DOUBLE smlnum;
    static long int wantsp;
    static LONG DOUBLE eps, est;



#define DUMMY(I) dummy[(I)]
#define SELECT(I) select[(I)-1]
#define S(I) s[(I)-1]
#define SEP(I) sep[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]
#define WORK(I,J) work[(I)-1 + ((J)-1)* ( *ldwork)]

    wantbh = lsame_(job, "B");
    wants = lsame_(job, "E") || wantbh;
    wantsp = lsame_(job, "V") || wantbh;

    somcon = lsame_(howmny, "S");

    *info = 0;
    if (! wants && ! wantsp) {
	*info = -1;
    } else if (! lsame_(howmny, "A") && ! somcon) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldt < MAX(1,*n)) {
	*info = -6;
    } else if (*ldvl < 1 || (wants && *ldvl < *n)) {
	*info = -8;
    } else if (*ldvr < 1 || (wants && *ldvr < *n)) {
	*info = -10;
    } else {

/*        Set M to the number of eigenpairs for which condition number
s   
          are required, and test MM. */

	if (somcon) {
	    *m = 0;
	    pair = 0;
	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {
		if (pair) {
		    pair = 0;
		} else {
		    if (k < *n) {
			if (T(k+1,k) == 0.) {
			    if (SELECT(k)) {
				++(*m);
			    }
			} else {
			    pair = 1;
			    if (SELECT(k) || SELECT(k + 1)) {
				*m += 2;
			    }
			}
		    } else {
			if (SELECT(*n)) {
			    ++(*m);
			}
		    }
		}
/* L10: */
	    }
	} else {
	    *m = *n;
	}

	if (*mm < *m) {
	    *info = -13;
	} else if (*ldwork < 1 || (wantsp && *ldwork < *n)) {
	    *info = -16;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTRSNA", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (somcon) {
	    if (! SELECT(1)) {
		return;
	    }
	}
	if (wants) {
	    S(1) = 1.;
	}
	if (wantsp) {
	    SEP(1) = (d__1 = T(1,1), ABS(d__1));
	}
	return;
    }

/*     Get machine constants */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("P");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("P");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("P");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S") / eps;
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S") / eps;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S") / eps;
#endif

    bignum = 1. / smlnum;

#ifdef PETSC_PREFIX_SUFFIX
    dlabad_(&smlnum, &bignum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlabad(&smlnum, &bignum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlabad_(&smlnum, &bignum);
#endif


    ks = 0;
    pair = 0;
    i__1 = *n;
    for (k = 1; k <= *n; ++k) {

/*        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block. */

	if (pair) {
	    pair = 0;
	    goto L60;
	} else {
	    if (k < *n) {
		pair = T(k+1,k) != 0.;
	    }
	}

/*        Determine whether condition numbers are required for the k-t
h   
          eigenpair. */

	if (somcon) {
	    if (pair) {
		if (! SELECT(k) && ! SELECT(k + 1)) {
		    goto L60;
		}
	    } else {
		if (! SELECT(k)) {
		    goto L60;
		}
	    }
	}

	++ks;

	if (wants) {

/*           Compute the reciprocal condition number of the k-th 
  
             eigenvalue. */

	    if (! pair) {

/*              Real eigenvalue. */


#ifdef PETSC_PREFIX_SUFFIX
		prod = ddot_(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		prod = qdot(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		prod = qdot_(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		rnrm = dnrm2_(n, &VR(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		rnrm = qnrm2(n, &VR(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		rnrm = qnrm2_(n, &VR(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		lnrm = dnrm2_(n, &VL(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		lnrm = qnrm2(n, &VL(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		lnrm = qnrm2_(n, &VL(1,ks), &c__1);
#endif

		S(ks) = ABS(prod) / (rnrm * lnrm);
	    } else {

/*              Complex eigenvalue. */


#ifdef PETSC_PREFIX_SUFFIX
		prod1 = ddot_(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		prod1 = qdot(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		prod1 = qdot_(n, &VR(1,ks), &c__1, &VL(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		prod1 += ddot_(n, &VR(1,ks+1), &c__1, &VL(1,ks+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		prod1 += qdot(n, &VR(1,ks+1), &c__1, &VL(1,ks+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		prod1 += qdot_(n, &VR(1,ks+1), &c__1, &VL(1,ks+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		prod2 = ddot_(n, &VL(1,ks), &c__1, &VR(1,ks+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		prod2 = qdot(n, &VL(1,ks), &c__1, &VR(1,ks+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		prod2 = qdot_(n, &VL(1,ks), &c__1, &VR(1,ks+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		prod2 -= ddot_(n, &VL(1,ks+1), &c__1, &VR(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		prod2 -= qdot(n, &VL(1,ks+1), &c__1, &VR(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		prod2 -= qdot_(n, &VL(1,ks+1), &c__1, &VR(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__1 = dnrm2_(n, &VR(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__1 = qnrm2(n, &VR(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__1 = qnrm2_(n, &VR(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__2 = dnrm2_(n, &VR(1,ks+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__2 = qnrm2(n, &VR(1,ks+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__2 = qnrm2_(n, &VR(1,ks+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		rnrm = dlapy2_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		rnrm = qlapy2(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		rnrm = qlapy2_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__1 = dnrm2_(n, &VL(1,ks), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__1 = qnrm2(n, &VL(1,ks), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__1 = qnrm2_(n, &VL(1,ks), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__2 = dnrm2_(n, &VL(1,ks+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__2 = qnrm2(n, &VL(1,ks+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__2 = qnrm2_(n, &VL(1,ks+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		lnrm = dlapy2_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		lnrm = qlapy2(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		lnrm = qlapy2_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		cond = dlapy2_(&prod1, &prod2) / (rnrm * lnrm);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		cond = qlapy2(&prod1, &prod2) / (rnrm * lnrm);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		cond = qlapy2_(&prod1, &prod2) / (rnrm * lnrm);
#endif

		S(ks) = cond;
		S(ks + 1) = cond;
	    }
	}

	if (wantsp) {

/*           Estimate the reciprocal condition number of the k-th 
  
             eigenvector.   

             Copy the matrix T to the array WORK and swap the diag
onal   
             block beginning at T(k,k) to the (1,1) position. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("Full", n, n, &T(1,1), ldt, &WORK(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("Full", n, n, &T(1,1), ldt, &WORK(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("Full", n, n, &T(1,1), ldt, &WORK(1,1), 
#endif

		    ldwork);
	    ifst = k;
	    ilst = 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtrexc_("No Q", n, &WORK(1,1), ldwork, dummy, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrexc("No Q", n, &WORK(1,1), ldwork, dummy, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrexc_("No Q", n, &WORK(1,1), ldwork, dummy, &c__1, &
#endif

		    ifst, &ilst, &WORK(1,*n+1), &ierr);

	    if (ierr == 1 || ierr == 2) {

/*              Could not swap because blocks not well separat
ed */

		scale = 1.;
		est = bignum;
	    } else {

/*              Reordering successful */

		if (WORK(2,1) == 0.) {

/*                 Form C = T22 - lambda*I in WORK(2:N,2:N
). */

		    i__2 = *n;
		    for (i = 2; i <= *n; ++i) {
			WORK(i,i) -= WORK(1,1);
/* L20: */
		    }
		    n2 = 1;
		    nn = *n - 1;
		} else {

/*                 Triangularize the 2 by 2 block by unita
ry   
                   transformation U = [  cs   i*ss ]   
                                      [ i*ss   cs  ].   
                   such that the (1,1) position of WORK is
 complex   
                   eigenvalue lambda with positive imagina
ry part. (2,2)   
                   position of WORK is the complex eigenva
lue lambda   
                   with negative imaginary  part. */

		    d__1 = WORK(1,2);d__2 = WORK(2,1);
                    mu = sqrt(( ABS(d__1))) 
			    * sqrt(( ABS(d__2)));

#ifdef PETSC_PREFIX_SUFFIX
		    delta = dlapy2_(&mu, &WORK(2,1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    delta = qlapy2(&mu, &WORK(2,1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    delta = qlapy2_(&mu, &WORK(2,1));
#endif

		    cs = mu / delta;
		    sn = -WORK(2,1) / delta;

/*                 Form   

                   C' = WORK(2:N,2:N) + i*[rwork(1) ..... 
rwork(n-1) ]   
                                          [   mu          
           ]   
                                          [         ..    
           ]   
                                          [             ..
           ]   
                                          [               
   mu      ]   
                   where C' is conjugate transpose of comp
lex matrix C,   
                   and RWORK is stored starting in the N+1
-st column of   
                   WORK. */

		    i__2 = *n;
		    for (j = 3; j <= *n; ++j) {
			WORK(2,j) = cs * WORK(2,j)
				;
			WORK(j,j) -= WORK(1,1);
/* L30: */
		    }
		    WORK(2,2) = 0.;

		    WORK(1,*n+1) = mu * 2.;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			WORK(i,*n+1) = sn * WORK(1,i+1);
/* L40: */
		    }
		    n2 = 2;
		    nn = (*n - 1) << 1;
		}

/*              Estimate norm(inv(C')) */

		est = 0.;
		kase = 0;
L50:

#ifdef PETSC_PREFIX_SUFFIX
		dlacon_(&nn, &WORK(1,*n+2), &WORK(1,*n+4), &IWORK(1), &est, &kase);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacon(&nn, &WORK(1,*n+2), &WORK(1,*n+4), &IWORK(1), &est, &kase);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacon_(&nn, &WORK(1,*n+2), &WORK(1,*n+4), &IWORK(1), &est, &kase);
#endif

		if (kase != 0) {
		    if (kase == 1) {
			if (n2 == 1) {

/*                       Real eigenvalue: solve C'
*x = scale*c. */

			    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dlaqtr_(&c_true, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &dumm, &scale, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlaqtr(&c_true, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &dumm, &scale, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlaqtr_(&c_true, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &dumm, &scale, 
#endif

				    &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
			} else {

/*                       Complex eigenvalue: solve
   
                         C'*(p+iq) = scale*(c+id) 
in real arithmetic. */

			    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dlaqtr_(&c_true, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlaqtr(&c_true, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlaqtr_(&c_true, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif

			}
		    } else {
			if (n2 == 1) {

/*                       Real eigenvalue: solve C*
x = scale*c. */

			    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dlaqtr_(&c_false, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlaqtr(&c_false, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlaqtr_(&c_false, &c_true, &i__2, &WORK(2,2), ldwork, dummy, &
#endif

				    dumm, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &
				    ierr);
			} else {

/*                       Complex eigenvalue: solve
   
                         C*(p+iq) = scale*(c+id) i
n real arithmetic. */

			    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dlaqtr_(&c_false, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlaqtr(&c_false, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlaqtr_(&c_false, &c_false, &i__2, &WORK(2,2), ldwork, &WORK(1,*n+1), &mu, &scale, &WORK(1,*n+4), &WORK(1,*n+6), &ierr);
#endif


			}
		    }

		    goto L50;
		}
	    }

	    SEP(ks) = scale / MAX(est,smlnum);
	    if (pair) {
		SEP(ks + 1) = SEP(ks);
	    }
	}

	if (pair) {
	    ++ks;
	}

L60:
	;
    }
    return;

/*     End of DTRSNA */

} /* dtrsna_ */

