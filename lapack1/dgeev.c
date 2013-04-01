#include <math.h>

#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgeev_(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgeev(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgeev_(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *wr, LONG DOUBLE *wi, LONG DOUBLE *vl, 
	int *ldvl, LONG DOUBLE *vr, int *ldvr, LONG DOUBLE *work, 
	int *lwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEEV computes for an N-by-N real nonsymmetric matrix A, the   
    eigenvalues and, optionally, the left and/or right eigenvectors.   

    The right eigenvector v(j) of A satisfies   
                     A * v(j) = lambda(j) * v(j)   
    where lambda(j) is its eigenvalue.   
    The left eigenvector u(j) of A satisfies   
                  u(j)**H * A = lambda(j) * u(j)**H   
    where u(j)**H denotes the conjugate transpose of u(j).   

    The computed eigenvectors are normalized to have Euclidean norm   
    equal to 1 and largest component real.   

    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N': left eigenvectors of A are not computed;   
            = 'V': left eigenvectors of A are computed.   

    JOBVR   (input) CHARACTER*1   
            = 'N': right eigenvectors of A are not computed;   
            = 'V': right eigenvectors of A are computed.   

    N       (input) INT   
            The order of the matrix A. N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten.   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= max(1,N).   

    WR      (output) LONG DOUBLE PRECISION array, dimension (N)   
    WI      (output) LONG DOUBLE PRECISION array, dimension (N)   
            WR and WI contain the real and imaginary parts,   
            respectively, of the computed eigenvalues.  Complex   
            conjugate pairs of eigenvalues appear consecutively   
            with the eigenvalue having the positive imaginary part   
            first.   

    VL      (output) LONG DOUBLE PRECISION array, dimension (LDVL,N)   
            If JOBVL = 'V', the left eigenvectors u(j) are stored one   
            after another in the columns of VL, in the same order   
            as their eigenvalues.   
            If JOBVL = 'N', VL is not referenced.   
            If the j-th eigenvalue is real, then u(j) = VL(:,j),   
            the j-th column of VL.   
            If the j-th and (j+1)-st eigenvalues form a complex   
            conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and   
            u(j+1) = VL(:,j) - i*VL(:,j+1).   

    LDVL    (input) INT   
            The leading dimension of the array VL.  LDVL >= 1; if   
            JOBVL = 'V', LDVL >= N.   

    VR      (output) LONG DOUBLE PRECISION array, dimension (LDVR,N)   
            If JOBVR = 'V', the right eigenvectors v(j) are stored one   
            after another in the columns of VR, in the same order   
            as their eigenvalues.   
            If JOBVR = 'N', VR is not referenced.   
            If the j-th eigenvalue is real, then v(j) = VR(:,j),   
            the j-th column of VR.   
            If the j-th and (j+1)-st eigenvalues form a complex   
            conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and   
            v(j+1) = VR(:,j) - i*VR(:,j+1).   

    LDVR    (input) INT   
            The leading dimension of the array VR.  LDVR >= 1; if   
            JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INT   
            The dimension of the array WORK.  LWORK >= max(1,3*N), and   
            if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good   
            performance, LWORK must generally be larger.   

    INFO    (output) INT   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the QR algorithm failed to compute all the 
  
                  eigenvalues, and no eigenvectors have been computed;   
                  elements i+1:N of WR and WI contain eigenvalues which   
                  have converged.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c__0 = 0;
    static int c__8 = 8;
    static int c_n1 = -1;
    static int c__4 = 4;
    
    /* System generated locals */
    int i__1,i__2, i__3, i__4;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */

    /* Local variables */
    static int ibal;
    static char side[1];
    static int maxb;
    static LONG DOUBLE anrm;
    static int ierr, itau;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *);
    static int iwrk, nout;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

    static int i, k;
    static LONG DOUBLE r;

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
    extern /* Subroutine */ void dlabad_(LONG DOUBLE *, LONG DOUBLE *), dgebak_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad(LONG DOUBLE *, LONG DOUBLE *), dgebak_(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad_(LONG DOUBLE *, LONG DOUBLE *), dgebak_(
#endif

	    char *, char *, int *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dgebal_(char *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgebal(char *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgebal_(char *, int *, LONG DOUBLE *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static LONG DOUBLE cs;
    static long int scalea;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE cscale;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlange_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlange(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlange_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgehrd_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgehrd(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgehrd_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    int *);
    static LONG DOUBLE sn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlascl_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, int *);

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
    extern /* Subroutine */ void dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *), xerbla_(char *, int *);
    static long int select[1];
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorghr_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorghr(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorghr_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dhseqr_(char *, char *, int *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qhseqr(char *, char *, int *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qhseqr_(char *, char *, int *, int *, int 
#endif

	    *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), dtrevc_(char *, char *, long int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qtrevc(char *, char *, long int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qtrevc_(char *, char *, long int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *, int *, LONG DOUBLE *, int *);
    static int minwrk, maxwrk;
    static long int wantvl;
    static LONG DOUBLE smlnum;
    static int hswork;
    static long int wantvr;
    static int ihi;
    static LONG DOUBLE scl;
    static int ilo;
    static LONG DOUBLE dum[1], eps;



#define DUM(I) dum[(I)]
#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    *info = 0;
    wantvl = lsame_(jobvl, "V");
    wantvr = lsame_(jobvr, "V");
    if (! wantvl && ! lsame_(jobvl, "N")) {
	*info = -1;
    } else if (! wantvr && ! lsame_(jobvr, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldvl < 1 || (wantvl && *ldvl < *n)) {
	*info = -9;
    } else if (*ldvr < 1 || (wantvr && *ldvr < *n)) {
	*info = -11;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.   
         HSWORK refers to the workspace preferred by DHSEQR, as   
         calculated below. HSWORK is computed assuming ILO=1 and IHI=N,   
         the worst case.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "DGEHRD", " ", n, &c__1, n, &
		c__0, 6L, 1L);
	if (! wantvl && ! wantvr) {
/* Computing MAX */
	    i__1 = 1, i__2 = *n * 3;
	    minwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "DHSEQR", "EN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = MAX(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "DHSEQR", "EN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = MIN(maxb,*n), i__2 = MAX(i__3,i__4);
	    k = MIN(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + 1, i__1 = MAX(i__1,i__2), i__2 = *n + 
		    hswork;
	    maxwrk = MAX(i__1,i__2);
	} else {
/* Computing MAX */
	    i__1 = 1, i__2 = *n << 2;
	    minwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "DOR"
		    "GHR", " ", n, &c__1, n, &c_n1, 6L, 1L);
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "DHSEQR", "SV", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = MAX(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "DHSEQR", "SV", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = MIN(maxb,*n), i__2 = MAX(i__3,i__4);
	    k = MIN(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + 1, i__1 = MAX(i__1,i__2), i__2 = *n + 
		    hswork;
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n << 2;
	    maxwrk = MAX(i__1,i__2);
	}
	WORK(1) = (LONG DOUBLE) maxwrk;
    }
    if (*lwork < minwrk) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEEV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
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
    smlnum = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S");
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

    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", n, n, &A(1,1), lda, dum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", n, n, &A(1,1), lda, dum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", n, n, &A(1,1), lda, dum);
#endif

    scalea = 0;
    if (anrm > 0. && anrm < smlnum) {
	scalea = 1;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = 1;
	cscale = bignum;
    }
    if (scalea) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &cscale, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &A(1,1), lda, &
#endif

		ierr);
    }

/*     Balance the matrix   
       (Workspace: need N) */

    ibal = 1;

#ifdef PETSC_PREFIX_SUFFIX
    dgebal_("B", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgebal("B", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgebal_("B", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif


/*     Reduce to upper Hessenberg form   
       (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = ibal + *n;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dgehrd_(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgehrd(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgehrd_(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
#endif

	     &ierr);

    if (wantvl) {

/*        Want left eigenvectors   
          Copy Householder vectors to VL */

	*(unsigned char *)side = 'L';

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("L", n, n, &A(1,1), lda, &VL(1,1), ldvl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("L", n, n, &A(1,1), lda, &VL(1,1), ldvl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("L", n, n, &A(1,1), lda, &VL(1,1), ldvl);
#endif


/*        Generate orthogonal matrix in VL   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorghr_(n, &ilo, &ihi, &VL(1,1), ldvl, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorghr(n, &ilo, &ihi, &VL(1,1), ldvl, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorghr_(n, &ilo, &ihi, &VL(1,1), ldvl, &WORK(itau), &WORK(iwrk),
#endif

		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dhseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qhseqr("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qhseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif

		VL(1,1), ldvl, &WORK(iwrk), &i__1, info);

	if (wantvr) {

/*           Want left and right eigenvectors   
             Copy Schur vectors to VR */

	    *(unsigned char *)side = 'B';

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("F", n, n, &VL(1,1), ldvl, &VR(1,1), ldvr)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("F", n, n, &VL(1,1), ldvl, &VR(1,1), ldvr)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("F", n, n, &VL(1,1), ldvl, &VR(1,1), ldvr)
#endif

		    ;
	}

    } else if (wantvr) {

/*        Want right eigenvectors   
          Copy Householder vectors to VR */

	*(unsigned char *)side = 'R';

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("L", n, n, &A(1,1), lda, &VR(1,1), ldvr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("L", n, n, &A(1,1), lda, &VR(1,1), ldvr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("L", n, n, &A(1,1), lda, &VR(1,1), ldvr);
#endif


/*        Generate orthogonal matrix in VR   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorghr_(n, &ilo, &ihi, &VR(1,1), ldvr, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorghr(n, &ilo, &ihi, &VR(1,1), ldvr, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorghr_(n, &ilo, &ihi, &VR(1,1), ldvr, &WORK(itau), &WORK(iwrk),
#endif

		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dhseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qhseqr("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qhseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif

		VR(1,1), ldvr, &WORK(iwrk), &i__1, info);

    } else {

/*        Compute eigenvalues only   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dhseqr_("E", "N", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qhseqr("E", "N", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qhseqr_("E", "N", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
#endif

		VR(1,1), ldvr, &WORK(iwrk), &i__1, info);
    }

/*     If INFO > 0 from DHSEQR, then quit */

    if (*info > 0) {
	goto L50;
    }

    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors   
          (Workspace: need 4*N) */


#ifdef PETSC_PREFIX_SUFFIX
	dtrevc_(side, "B", select, n, &A(1,1), lda, &VL(1,1), ldvl,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrevc(side, "B", select, n, &A(1,1), lda, &VL(1,1), ldvl,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrevc_(side, "B", select, n, &A(1,1), lda, &VL(1,1), ldvl,
#endif

		 &VR(1,1), ldvr, n, &nout, &WORK(iwrk), &ierr);
    }

    if (wantvl) {

/*        Undo balancing of left eigenvectors   
          (Workspace: need N) */


#ifdef PETSC_PREFIX_SUFFIX
	dgebak_("B", "L", n, &ilo, &ihi, &WORK(ibal), n, &VL(1,1), ldvl,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgebak("B", "L", n, &ilo, &ihi, &WORK(ibal), n, &VL(1,1), ldvl,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgebak_("B", "L", n, &ilo, &ihi, &WORK(ibal), n, &VL(1,1), ldvl,
#endif

		 &ierr);

/*        Normalize left eigenvectors and make largest component real 
*/

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WI(i) == 0.) {

#ifdef PETSC_PREFIX_SUFFIX
		scl = 1. / dnrm2_(n, &VL(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		scl = 1. / qnrm2(n, &VL(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		scl = 1. / qnrm2_(n, &VL(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VL(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VL(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VL(1,i), &c__1);
#endif

	    } else if (WI(i) > 0.) {

#ifdef PETSC_PREFIX_SUFFIX
		d__1 = dnrm2_(n, &VL(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__1 = qnrm2(n, &VL(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__1 = qnrm2_(n, &VL(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__2 = dnrm2_(n, &VL(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__2 = qnrm2(n, &VL(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__2 = qnrm2_(n, &VL(1,i+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		scl = 1. / dlapy2_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		scl = 1. / qlapy2(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		scl = 1. / qlapy2_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VL(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VL(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VL(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VL(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VL(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VL(1,i+1), &c__1);
#endif

		i__2 = *n;
		for (k = 1; k <= *n; ++k) {
/* Computing 2nd power */
		    d__1 = VL(k,i);
/* Computing 2nd power */
		    d__2 = VL(k,i+1);
		    WORK(iwrk + k - 1) = d__1 * d__1 + d__2 * d__2;
/* L10: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		k = idamax_(n, &WORK(iwrk), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		k = iqamax(n, &WORK(iwrk), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		k = iqamax_(n, &WORK(iwrk), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&VL(k,i), &VL(k,i+1), &cs,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&VL(k,i), &VL(k,i+1), &cs,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&VL(k,i), &VL(k,i+1), &cs,
#endif

			 &sn, &r);

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &VL(1,i), &c__1, &VL(1,i+1), &c__1, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &VL(1,i), &c__1, &VL(1,i+1), &c__1, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &VL(1,i), &c__1, &VL(1,i+1), &c__1, &cs, &sn);
#endif

		VL(k,i+1) = 0.;
	    }
/* L20: */
	}
    }

    if (wantvr) {

/*        Undo balancing of right eigenvectors   
          (Workspace: need N) */


#ifdef PETSC_PREFIX_SUFFIX
	dgebak_("B", "R", n, &ilo, &ihi, &WORK(ibal), n, &VR(1,1), ldvr,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgebak("B", "R", n, &ilo, &ihi, &WORK(ibal), n, &VR(1,1), ldvr,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgebak_("B", "R", n, &ilo, &ihi, &WORK(ibal), n, &VR(1,1), ldvr,
#endif

		 &ierr);

/*        Normalize right eigenvectors and make largest component real
 */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WI(i) == 0.) {

#ifdef PETSC_PREFIX_SUFFIX
		scl = 1. / dnrm2_(n, &VR(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		scl = 1. / qnrm2(n, &VR(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		scl = 1. / qnrm2_(n, &VR(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VR(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VR(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VR(1,i), &c__1);
#endif

	    } else if (WI(i) > 0.) {

#ifdef PETSC_PREFIX_SUFFIX
		d__1 = dnrm2_(n, &VR(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__1 = qnrm2(n, &VR(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__1 = qnrm2_(n, &VR(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		d__2 = dnrm2_(n, &VR(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		d__2 = qnrm2(n, &VR(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		d__2 = qnrm2_(n, &VR(1,i+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		scl = 1. / dlapy2_(&d__1, &d__2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		scl = 1. / qlapy2(&d__1, &d__2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		scl = 1. / qlapy2_(&d__1, &d__2);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VR(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VR(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VR(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dscal_(n, &scl, &VR(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(n, &scl, &VR(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(n, &scl, &VR(1,i+1), &c__1);
#endif

		i__2 = *n;
		for (k = 1; k <= *n; ++k) {
/* Computing 2nd power */
		    d__1 = VR(k,i);
/* Computing 2nd power */
		    d__2 = VR(k,i+1);
		    WORK(iwrk + k - 1) = d__1 * d__1 + d__2 * d__2;
/* L30: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		k = idamax_(n, &WORK(iwrk), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		k = iqamax(n, &WORK(iwrk), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		k = iqamax_(n, &WORK(iwrk), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&VR(k,i), &VR(k,i+1), &cs,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&VR(k,i), &VR(k,i+1), &cs,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&VR(k,i), &VR(k,i+1), &cs,
#endif

			 &sn, &r);

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &VR(1,i), &c__1, &VR(1,i+1), &c__1, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &VR(1,i), &c__1, &VR(1,i+1), &c__1, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &VR(1,i), &c__1, &VR(1,i+1), &c__1, &cs, &sn);
#endif

		VR(k,i+1) = 0.;
	    }
/* L40: */
	}
    }

/*     Undo scaling if necessary */

L50:
    if (scalea) {
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = MAX(i__3,1);

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(*info + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(*info + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(*info + 
#endif

		1), &i__2, &ierr);
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = MAX(i__3,1);

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(*info + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(*info + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(*info + 
#endif

		1), &i__2, &ierr);
	if (*info > 0) {
	    i__1 = ilo - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(1), 
#endif

		    n, &ierr);
	    i__1 = ilo - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(1), 
#endif

		    n, &ierr);
	}
    }

    WORK(1) = (LONG DOUBLE) maxwrk;
    return;

/*     End of DGEEV */

} /* dgeev_ */

