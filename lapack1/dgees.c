#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgees_(char *jobvs, char *sort, long int (*select)(LONG DOUBLE *,LONG DOUBLE*), int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgees(char *jobvs, char *sort, long int (*select)(LONG DOUBLE *,LONG DOUBLE*), int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgees_(char *jobvs, char *sort, long int (*select)(LONG DOUBLE *,LONG DOUBLE*), int *n, 
#endif

	LONG DOUBLE *a, int *lda, int *sdim, LONG DOUBLE *wr, 
	LONG DOUBLE *wi, LONG DOUBLE *vs, int *ldvs, LONG DOUBLE *work, 
	int *lwork, long int *bwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGEES computes for an N-by-N real nonsymmetric matrix A, the   
    eigenvalues, the real Schur form T, and, optionally, the matrix of   
    Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T). 
  

    Optionally, it also orders the eigenvalues on the diagonal of the   
    real Schur form so that selected eigenvalues are at the top left.   
    The leading columns of Z then form an orthonormal basis for the   
    invariant subspace corresponding to the selected eigenvalues.   

    A matrix is in real Schur form if it is upper quasi-triangular with   
    1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the   
    form   
            [  a  b  ]   
            [  c  a  ]   

    where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).   

    Arguments   
    =========   

    JOBVS   (input) CHARACTER*1   
            = 'N': Schur vectors are not computed;   
            = 'V': Schur vectors are computed.   

    SORT    (input) CHARACTER*1   
            Specifies whether or not to order the eigenvalues on the   
            diagonal of the Schur form.   
            = 'N': Eigenvalues are not ordered;   
            = 'S': Eigenvalues are ordered (see SELECT).   

    SELECT  (input) LOGICAL FUNCTION of two LONG DOUBLE PRECISION arguments   
            SELECT must be declared EXTERNAL in the calling subroutine.   
            If SORT = 'S', SELECT is used to select eigenvalues to sort   
            to the top left of the Schur form.   
            If SORT = 'N', SELECT is not referenced.   
            An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if   
            SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex 
  
            conjugate pair of eigenvalues is selected, then both complex 
  
            eigenvalues are selected.   
            Note that a selected complex eigenvalue may no longer   
            satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since   
            ordering may change the value of complex eigenvalues   
            (especially if the eigenvalue is ill-conditioned); in this   
            case INFO is set to N+2 (see INFO below).   

    N       (input) INTEGER   
            The order of the matrix A. N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten by its real Schur form T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    SDIM    (output) INTEGER   
            If SORT = 'N', SDIM = 0.   
            If SORT = 'S', SDIM = number of eigenvalues (after sorting)   
                           for which SELECT is true. (Complex conjugate   
                           pairs for which SELECT is true for either   
                           eigenvalue count as 2.)   

    WR      (output) LONG DOUBLE PRECISION array, dimension (N)   
    WI      (output) LONG DOUBLE PRECISION array, dimension (N)   
            WR and WI contain the real and imaginary parts,   
            respectively, of the computed eigenvalues in the same order   
            that they appear on the diagonal of the output Schur form T. 
  
            Complex conjugate pairs of eigenvalues will appear   
            consecutively with the eigenvalue having the positive   
            imaginary part first.   

    VS      (output) LONG DOUBLE PRECISION array, dimension (LDVS,N)   
            If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur 
  
            vectors.   
            If JOBVS = 'N', VS is not referenced.   

    LDVS    (input) INTEGER   
            The leading dimension of the array VS.  LDVS >= 1; if   
            JOBVS = 'V', LDVS >= N.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) contains the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= MAX(1,3*N).   
            For good performance, LWORK must generally be larger.   

    BWORK   (workspace) LOGICAL array, dimension (N)   
            Not referenced if SORT = 'N'.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value.   
            > 0: if INFO = i, and i is   
               <= N: the QR algorithm failed to compute all the   
                     eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI 
  
                     contain those eigenvalues which have converged; if   
                     JOBVS = 'V', VS contains the matrix which reduces A 
  
                     to its partially converged Schur form.   
               = N+1: the eigenvalues could not be reordered because some 
  
                     eigenvalues were too close to separate (the problem 
  
                     is very ill-conditioned);   
               = N+2: after reordering, roundoff changed values of some   
                     complex eigenvalues so that leading eigenvalues in   
                     the Schur form no longer satisfy SELECT=.TRUE.  This 
  
                     could also be caused by underflow due to scaling.   

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
    int i__1, i__2, i__3, i__4;
    /* Builtin functions */
    /* Local variables */
    static int ibal, maxb;
    static LONG DOUBLE anrm;
    static int idum[1], ierr, itau, iwrk, inxt, i, k;
    static LONG DOUBLE s;
    static int icond, ieval;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dswap_(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap_(int *, LONG DOUBLE *, int 
#endif

	    *, LONG DOUBLE *, int *);
    static long int cursl;
    static int i1, i2;

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
    static long int lst2sl, scalea;
    static int ip;
    static LONG DOUBLE cscale;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *), dlange_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlange_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlange_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *);

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

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlascl_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlascl(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlascl_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);
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
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), dtrsen_(char *, char *, long int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qtrsen(char *, char *, long int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qtrsen_(char *, char *, long int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *, int *, int *, int *);
    static long int lastsl;
    static int minwrk, maxwrk;
    static LONG DOUBLE smlnum;
    static int hswork;
    static long int wantst, wantvs;
    static int ihi, ilo;
    static LONG DOUBLE dum[1], eps, sep;



#define IDUM(I) idum[(I)]
#define DUM(I) dum[(I)]
#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]
#define WORK(I) work[(I)-1]
#define BWORK(I) bwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define VS(I,J) vs[(I)-1 + ((J)-1)* ( *ldvs)]

    *info = 0;
    wantvs = lsame_(jobvs, "V");
    wantst = lsame_(sort, "S");
    if (! wantvs && ! lsame_(jobvs, "N")) {
	*info = -1;
    } else if (! wantst && ! lsame_(sort, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < MAX(1,*n)) {
	*info = -6;
    } else if (*ldvs < 1 || (wantvs && *ldvs < *n)) {
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
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3;
	minwrk = MAX(i__1,i__2);
	if (! wantvs) {
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "DHSEQR", "SN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = MAX(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "DHSEQR", "SN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = MIN(maxb,*n), i__2 = MAX(i__3,i__4);
	    k = MIN(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + hswork, i__1 = MAX(i__1,i__2);
	    maxwrk = MAX(i__1,1);
	} else {
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "DOR"
		    "GHR", " ", n, &c__1, n, &c_n1, 6L, 1L);
	    maxwrk = MAX(i__1,i__2);
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
	    i__1 = maxwrk, i__2 = *n + hswork, i__1 = MAX(i__1,i__2);
	    maxwrk = MAX(i__1,1);
	}
	WORK(1) = (LONG DOUBLE) maxwrk;
    }
    if (*lwork < minwrk) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEES ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	*sdim = 0;
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

/*     Permute the matrix to make it more nearly triangular   
       (Workspace: need N) */

    ibal = 1;

#ifdef PETSC_PREFIX_SUFFIX
    dgebal_("P", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgebal("P", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgebal_("P", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);
#endif


/*     Reduce to upper Hessenberg form   
       (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = *n + ibal;
    iwrk = *n + itau;
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

    if (wantvs) {

/*        Copy Householder vectors to VS */


#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("L", n, n, &A(1,1), lda, &VS(1,1), ldvs);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("L", n, n, &A(1,1), lda, &VS(1,1), ldvs);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("L", n, n, &A(1,1), lda, &VS(1,1), ldvs);
#endif


/*        Generate orthogonal matrix in VS   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorghr_(n, &ilo, &ihi, &VS(1,1), ldvs, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorghr(n, &ilo, &ihi, &VS(1,1), ldvs, &WORK(itau), &WORK(iwrk),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorghr_(n, &ilo, &ihi, &VS(1,1), ldvs, &WORK(itau), &WORK(iwrk),
#endif

		 &i__1, &ierr);
    }

    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired 
  
       (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

    iwrk = itau;
    i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dhseqr_("S", jobvs, n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &VS(1,1), ldvs, &WORK(iwrk), &i__1, &ieval);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qhseqr("S", jobvs, n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &VS(1,1), ldvs, &WORK(iwrk), &i__1, &ieval);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qhseqr_("S", jobvs, n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &VS(1,1), ldvs, &WORK(iwrk), &i__1, &ieval);
#endif

    if (ieval > 0) {
	*info = ieval;
    }

/*     Sort eigenvalues if desired */

    if (wantst && *info == 0) {
	if (scalea) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WR(1), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WR(1), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WR(1), n, &
#endif

		    ierr);

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WI(1), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WI(1), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WI(1), n, &
#endif

		    ierr);
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    BWORK(i) = (*select)(&WR(i), &WI(i));
/* L10: */
	}

/*        Reorder eigenvalues and transform Schur vectors   
          (Workspace: none needed) */

	i__1 = *lwork - iwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dtrsen_("N", jobvs, &BWORK(1), n, &A(1,1), lda, &VS(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrsen("N", jobvs, &BWORK(1), n, &A(1,1), lda, &VS(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrsen_("N", jobvs, &BWORK(1), n, &A(1,1), lda, &VS(1,1), 
#endif

		ldvs, &WR(1), &WI(1), sdim, &s, &sep, &WORK(iwrk), &i__1, 
		idum, &c__1, &icond);
	if (icond > 0) {
	    *info = *n + icond;
	}
    }

    if (wantvs) {

/*        Undo balancing   
          (Workspace: need N) */


#ifdef PETSC_PREFIX_SUFFIX
	dgebak_("P", "R", n, &ilo, &ihi, &WORK(ibal), n, &VS(1,1), ldvs,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgebak("P", "R", n, &ilo, &ihi, &WORK(ibal), n, &VS(1,1), ldvs,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgebak_("P", "R", n, &ilo, &ihi, &WORK(ibal), n, &VS(1,1), ldvs,
#endif

		 &ierr);
    }

    if (scalea) {

/*        Undo scaling for the Schur form of A */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("H", &c__0, &c__0, &cscale, &anrm, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &A(1,1), lda, &
#endif

		ierr);
	i__1 = *lda + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &A(1,1), &i__1, &WR(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &A(1,1), &i__1, &WR(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &A(1,1), &i__1, &WR(1), &c__1);
#endif

	if (cscale == smlnum) {

/*           If scaling back towards underflow, adjust WI if an   
             offdiagonal element of a 2-by-2 block in the Schur fo
rm   
             underflows. */

	    if (ieval > 0) {
		i1 = ieval + 1;
		i2 = ihi - 1;
		i__1 = ilo - 1;
/* Computing MAX */
		i__3 = ilo - 1;
		i__2 = MAX(i__3,1);

#ifdef PETSC_PREFIX_SUFFIX
		dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(
#endif

			1), &i__2, &ierr);
	    } else if (wantst) {
		i1 = 1;
		i2 = *n - 1;
	    } else {
		i1 = ilo;
		i2 = ihi - 1;
	    }
	    inxt = i1 - 1;
	    i__1 = i2;
	    for (i = i1; i <= i2; ++i) {
		if (i < inxt) {
		    goto L20;
		}
		if (WI(i) == 0.) {
		    inxt = i + 1;
		} else {
		    if (A(i+1,i) == 0.) {
			WI(i) = 0.;
			WI(i + 1) = 0.;
		    } else if (A(i+1,i) != 0. && A(i,i+1) == 0.) {
			WI(i) = 0.;
			WI(i + 1) = 0.;
			if (i > 1) {
			    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dswap_(&i__2, &A(1,i), &c__1, &A(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qswap(&i__2, &A(1,i), &c__1, &A(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qswap_(&i__2, &A(1,i), &c__1, &A(1,i+1), &c__1);
#endif

			}
			if (*n > i + 1) {
			    i__2 = *n - i - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dswap_(&i__2, &A(i,i+2), lda, &A(i+1,i+2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qswap(&i__2, &A(i,i+2), lda, &A(i+1,i+2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qswap_(&i__2, &A(i,i+2), lda, &A(i+1,i+2), lda);
#endif

			}

#ifdef PETSC_PREFIX_SUFFIX
			dswap_(n, &VS(1,i), &c__1, &VS(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qswap(n, &VS(1,i), &c__1, &VS(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qswap_(n, &VS(1,i), &c__1, &VS(1,i+1), &c__1);
#endif

			A(i,i+1) = A(i+1,i);
			A(i+1,i) = 0.;
		    }
		    inxt = i + 2;
		}
L20:
		;
	    }
	}

/*        Undo scaling for the imaginary part of the eigenvalues */

	i__1 = *n - ieval;
/* Computing MAX */
	i__3 = *n - ieval;
	i__2 = MAX(i__3,1);

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(ieval + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(ieval + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(ieval + 
#endif

		1), &i__2, &ierr);
    }

    if (wantst && *info == 0) {

/*        Check if reordering successful */

	lastsl = 1;
	lst2sl = 1;
	*sdim = 0;
	ip = 0;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    cursl = (*select)(&WR(i), &WI(i));
	    if (WI(i) == 0.) {
		if (cursl) {
		    ++(*sdim);
		}
		ip = 0;
		if (cursl && ! lastsl) {
		    *info = *n + 2;
		}
	    } else {
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

		    cursl = cursl || lastsl;
		    lastsl = cursl;
		    if (cursl) {
			*sdim += 2;
		    }
		    ip = -1;
		    if (cursl && ! lst2sl) {
			*info = *n + 2;
		    }
		} else {

/*                 First eigenvalue of conjugate pair */

		    ip = 1;
		}
	    }
	    lst2sl = lastsl;
	    lastsl = cursl;
/* L30: */
	}
    }

    WORK(1) = (LONG DOUBLE) maxwrk;
    return;

/*     End of DGEES */

} /* dgees_ */

