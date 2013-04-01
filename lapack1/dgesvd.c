#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgesvd(char *jobu, char *jobvt, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgesvd_(char *jobu, char *jobvt, int *m, int *n, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *s, LONG DOUBLE *u, int *
	ldu, LONG DOUBLE *vt, int *ldvt, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGESVD computes the singular value decomposition (SVD) of a real   
    M-by-N matrix A, optionally computing the left and/or right singular 
  
    vectors. The SVD is written   

         A = U * SIGMA * transpose(V)   

    where SIGMA is an M-by-N matrix which is zero except for its   
    MIN(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and   
    V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA   
    are the singular values of A; they are real and non-negative, and   
    are returned in descending order.  The first MIN(m,n) columns of   
    U and V are the left and right singular vectors of A.   

    Note that the routine returns V**T, not V.   

    Arguments   
    =========   

    JOBU    (input) CHARACTER*1   
            Specifies options for computing all or part of the matrix U: 
  
            = 'A':  all M columns of U are returned in array U:   
            = 'S':  the first MIN(m,n) columns of U (the left singular   
                    vectors) are returned in the array U;   
            = 'O':  the first MIN(m,n) columns of U (the left singular   
                    vectors) are overwritten on the array A;   
            = 'N':  no columns of U (no left singular vectors) are   
                    computed.   

    JOBVT   (input) CHARACTER*1   
            Specifies options for computing all or part of the matrix   
            V**T:   
            = 'A':  all N rows of V**T are returned in the array VT;   
            = 'S':  the first MIN(m,n) rows of V**T (the right singular   
                    vectors) are returned in the array VT;   
            = 'O':  the first MIN(m,n) rows of V**T (the right singular   
                    vectors) are overwritten on the array A;   
            = 'N':  no rows of V**T (no right singular vectors) are   
                    computed.   

            JOBVT and JOBU cannot both be 'O'.   

    M       (input) INT   
            The number of rows of the input matrix A.  M >= 0.   

    N       (input) INT   
            The number of columns of the input matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit,   
            if JOBU = 'O',  A is overwritten with the first MIN(m,n)   
                            columns of U (the left singular vectors,   
                            stored columnwise);   
            if JOBVT = 'O', A is overwritten with the first MIN(m,n)   
                            rows of V**T (the right singular vectors,   
                            stored rowwise);   
            if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A   
                            are destroyed.   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= max(1,M).   

    S       (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The singular values of A, sorted so that S(i) >= S(i+1).   

    U       (output) LONG DOUBLE PRECISION array, dimension (LDU,UCOL)   
            (LDU,M) if JOBU = 'A' or (LDU,MIN(M,N)) if JOBU = 'S'.   
            If JOBU = 'A', U contains the M-by-M orthogonal matrix U;   
            if JOBU = 'S', U contains the first MIN(m,n) columns of U   
            (the left singular vectors, stored columnwise);   
            if JOBU = 'N' or 'O', U is not referenced.   

    LDU     (input) INT   
            The leading dimension of the array U.  LDU >= 1; if   
            JOBU = 'S' or 'A', LDU >= M.   

    VT      (output) LONG DOUBLE PRECISION array, dimension (LDVT,N)   
            If JOBVT = 'A', VT contains the N-by-N orthogonal matrix   
            V**T;   
            if JOBVT = 'S', VT contains the first MIN(m,n) rows of   
            V**T (the right singular vectors, stored rowwise);   
            if JOBVT = 'N' or 'O', VT is not referenced.   

    LDVT    (input) INT   
            The leading dimension of the array VT.  LDVT >= 1; if   
            JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= MIN(M,N).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK;   
            if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged   
            superdiagonal elements of an upper bidiagonal matrix B   
            whose diagonal is in S (not necessarily sorted). B   
            satisfies A = U * B * VT, so it has the same singular values 
  
            as A, and singular vectors related by U and VT.   

    LWORK   (input) INT   
            The dimension of the array WORK. LWORK >= 1.   
            LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N)-4).   
            For good performance, LWORK should generally be larger.   

    INFO    (output) INT   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if DBDSQR did not converge, INFO specifies how many   
                  superdiagonals of an intermediate bidiagonal form B   
                  did not converge to zero. See the description of WORK   
                  above for details.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__6 = 6;
    static int c__0 = 0;
    static int c__1 = 1;
    static int c_n1 = -1;
    static LONG DOUBLE c_b416 = 0.;
    static LONG DOUBLE c_b438 = 1.;
    
    /* System generated locals */
    int  i__2, i__3, i__4;
    char ch__1[3];
    /* Builtin functions */
    /* Local variables */
    static int iscl;
    static LONG DOUBLE anrm;
    static int ierr, itau, ncvt, nrvt, i;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemm_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    extern long int lsame_(char *, char *);
    static int chunk, minmn, wrkbl, itaup, itauq, mnthr, iwork;
    static long int wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    static int ie;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgebrd_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebrd(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebrd_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, int *, int *);

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
    static int ir, bdspac, iu;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgelqf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelqf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelqf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, int *),

#ifdef PETSC_PREFIX_SUFFIX
	     dgeqrf_(int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     qgeqrf(int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     qgeqrf_(int *, int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), dlacpy_(char *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlacpy(char *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlacpy_(char *,
#endif

	     int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dbdsqr_(char *, int *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qbdsqr(char *, int *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qbdsqr_(char *, int *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), dorgbr_(char *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qorgbr(char *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qorgbr_(char *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    int *);
    static LONG DOUBLE bignum;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormbr_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormbr(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormbr_(char *, char *, char *, int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), dorglq_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qorglq(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qorglq_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dorgqr_(int *, int *, int *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qorgqr(int *, int *, int *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qorgqr_(int *, int *, int *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static int ldwrkr, minwrk, ldwrku, maxwrk;
    static LONG DOUBLE smlnum;
    static long int wntuas, wntvas;
    static int blk, ncu;
    static LONG DOUBLE dum[1], eps;
    static int nru;



#define DUM(I) dum[(I)]
#define S(I) s[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define VT(I,J) vt[(I)-1 + ((J)-1)* ( *ldvt)]

    *info = 0;
    minmn = MIN(*m,*n);
/* Writing concatenation */
    /*
    i__1[0] = 1, a__1[0] = jobu;
    i__1[1] = 1, a__1[1] = jobvt;
    s_cat(ch__1, a__1, i__1, &c__2, 2L); */
    ch__1[0] = *jobu; ch__1[1] = *jobvt; ch__1[2] = 0;
    mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, 6L, 2L);
    wntua = lsame_(jobu, "A");
    wntus = lsame_(jobu, "S");
    wntuas = wntua || wntus;
    wntuo = lsame_(jobu, "O");
    wntun = lsame_(jobu, "N");
    wntva = lsame_(jobvt, "A");
    wntvs = lsame_(jobvt, "S");
    wntvas = wntva || wntvs;
    wntvo = lsame_(jobvt, "O");
    wntvn = lsame_(jobvt, "N");
    minwrk = 1;

    if (! (wntua || wntus || wntuo || wntun)) {
	*info = -1;
    } else if (! (wntva || wntvs || wntvo || wntvn) || (wntvo && wntuo)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < MAX(1,*m)) {
	*info = -6;
    } else if (*ldu < 1 || (wntuas && *ldu < *m)) {
	*info = -9;
    } else if (*ldvt < 1 || (wntva && *ldvt < *n) || (wntvs && *ldvt < minmn)) {
	*info = -11;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.) */

    if (*info == 0 && *lwork >= 1 && *m > 0 && *n > 0) {
	if (*m >= *n) {

/*           Compute space needed for DBDSQR   

   Computing MAX */
	    i__2 = *n * 3, i__3 = *n * 5 - 4;
	    bdspac = MAX(i__2,i__3);
	    if (*m >= mnthr) {
		if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N')
 */

		    maxwrk = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		    if (wntvo || wntvas) {
/* Computing MAX */
			i__2 = maxwrk, i__3 = *n * 3 + (*n - 1) * ilaenv_(&
				c__1, "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
			maxwrk = MAX(i__2,i__3);
		    }
		    maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
		    i__2 = *n << 2;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O',
 JOBVT='N') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
/* Computing MAX */
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
		    maxwrk = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O',
 JOBVT='S' or   
                   'A') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
/* Computing MAX */
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
		    maxwrk = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S',
 JOBVT='N') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S',
 JOBVT='O') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S',
 JOBVT='S' or   
                   'A') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A',
 JOBVT='N') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *m * ilaenv_(&c__1, "DORGQR", 
			    " ", m, m, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A',
 JOBVT='O') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *m * ilaenv_(&c__1, "DORGQR", 
			    " ", m, m, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A',
 JOBVT='S' or   
                   'A') */

		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n + *m * ilaenv_(&c__1, "DORGQR", 
			    " ", m, m, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORGBR"
			    , "Q", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
		    i__2 = *n * 3 + *m;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		}
	    } else {

/*              Path 10 (M at least N, but not much larger) */

		maxwrk = *n * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m,
			 n, &c_n1, &c_n1, 6L, 1L);
		if (wntus || wntuo) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *n * 3 + *n * ilaenv_(&c__1, "DORG"
			    "BR", "Q", m, n, n, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		if (wntua) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *n * 3 + *m * ilaenv_(&c__1, "DORG"
			    "BR", "Q", m, m, n, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		if (! wntvn) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", n, n, n, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
		i__2 = *n * 3 + *m;
		minwrk = MAX(i__2,bdspac);
		maxwrk = MAX(maxwrk,minwrk);
	    }
	} else {

/*           Compute space needed for DBDSQR   

   Computing MAX */
	    i__2 = *m * 3, i__3 = *m * 5 - 4;
	    bdspac = MAX(i__2,i__3);
	    if (*n >= mnthr) {
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N'
) */

		    maxwrk = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		    if (wntuo || wntuas) {
/* Computing MAX */
			i__2 = maxwrk, i__3 = *m * 3 + *m * ilaenv_(&c__1, 
				"DORGBR", "Q", m, m, m, &c_n1, 6L, 1L);
			maxwrk = MAX(i__2,i__3);
		    }
		    maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
		    i__2 = *m << 2;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N',
 JOBVT='O') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
/* Computing MAX */
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
		    maxwrk = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' 
or 'A',   
                   JOBVT='O') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR"
			    , "Q", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
/* Computing MAX */
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
		    maxwrk = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N',
 JOBVT='S') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O',
 JOBVT='S') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR"
			    , "Q", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' 
or 'A',   
                   JOBVT='S') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR"
			    , "Q", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N',
 JOBVT='A') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *n * ilaenv_(&c__1, "DORGLQ", 
			    " ", n, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O',
 JOBVT='A') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *n * ilaenv_(&c__1, "DORGLQ", 
			    " ", n, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR"
			    , "Q", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' 
or 'A',   
                   JOBVT='A') */

		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m + *n * ilaenv_(&c__1, "DORGLQ", 
			    " ", n, n, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR"
			    , "Q", m, m, m, &c_n1, 6L, 1L);
		    wrkbl = MAX(i__2,i__3);
		    wrkbl = MAX(wrkbl,bdspac);
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
		    i__2 = *m * 3 + *n;
		    minwrk = MAX(i__2,bdspac);
		    maxwrk = MAX(maxwrk,minwrk);
		}
	    } else {

/*              Path 10t(N greater than M, but not much larger
) */

		maxwrk = *m * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m,
			 n, &c_n1, &c_n1, 6L, 1L);
		if (wntvs || wntvo) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *m * 3 + *m * ilaenv_(&c__1, "DORG"
			    "BR", "P", m, n, m, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		if (wntva) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *m * 3 + *n * ilaenv_(&c__1, "DORG"
			    "BR", "P", n, n, m, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		if (! wntun) {
/* Computing MAX */
		    i__2 = maxwrk, i__3 = *m * 3 + (*m - 1) * ilaenv_(&c__1, 
			    "DORGBR", "Q", m, m, m, &c_n1, 6L, 1L);
		    maxwrk = MAX(i__2,i__3);
		}
		maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
		i__2 = *m * 3 + *n;
		minwrk = MAX(i__2,bdspac);
		maxwrk = MAX(maxwrk,minwrk);
	    }
	}
	WORK(1) = (LONG DOUBLE) maxwrk;
    }

    if (*lwork < minwrk) {
	*info = -13;
    }
    if (*info != 0) {
	i__2 = -(*info);
	xerbla_("DGESVD", &i__2);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	if (*lwork >= 1) {
	    WORK(1) = 1.;
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
    smlnum = sqrt(dlamch_("S")) / eps;
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = sqrt(qlamch("S")) / eps;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = sqrt(qlamch_("S")) / eps;
#endif

    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", m, n, &A(1,1), lda, dum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", m, n, &A(1,1), lda, dum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", m, n, &A(1,1), lda, dum);
#endif

    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, &
#endif

		ierr);
    } else if (anrm > bignum) {
	iscl = 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, &
#endif

		ierr);
    }

    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficientl
y   
          more rows than columns, first reduce using the QR   
          decomposition (if sufficient workspace available) */

	if (*m >= mnthr) {

	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N')   
                No left singular vectors to be computed */

		itau = 1;
		iwork = itau + *n;

/*              Compute A=Q*R   
                (Workspace: need 2*N, prefer N+N*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif

			i__2, &ierr);

/*              Zero out below R */

		i__2 = *n - 1;
		i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), 
#endif

			lda);
		ie = 1;
		itauq = ie + *n;
		itaup = itauq + *n;
		iwork = itaup + *n;

/*              Bidiagonalize R in A   
                (Workspace: need 4*N, prefer 3*N+2*N*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgebrd(n, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif

			itauq), &WORK(itaup), &WORK(iwork), &i__2, &ierr);
		ncvt = 0;
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, gene
rate P'.   
                   (Workspace: need 4*N-1, prefer 3*N+(N-1
)*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", n, n, n, &A(1,1), lda, &WORK(itaup), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &
#endif

			    WORK(iwork), &i__2, &ierr);
		    ncvt = *n;
		}
		iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing rig
ht   
                singular vectors of A in A if desired   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("U", n, &ncvt, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(iwork), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("U", n, &ncvt, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(iwork), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("U", n, &ncvt, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(iwork), 
#endif

			info);

/*              If right singular vectors desired in VT, copy 
them there */

		if (wntvas) {

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("F", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("F", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("F", n, n, &A(1,1), lda, &VT(1,1), 
#endif

			    ldvt);
		}

	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT=
'N')   
                N left singular vectors to be overwritten on A
 and   
                no right singular vectors to be computed   

   Computing MAX */
		i__2 = *n << 2;
		if (*lwork >= *n * *n + MAX(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorit
hm */

		    ir = 1;
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *lda * *n + *n;
		    if (*lwork >= MAX(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) i
s LDA by N */

			ldwrku = *lda;
			ldwrkr = *lda;
		    } else /* if(complicated condition) */ {
/* Computing MAX */
			i__2 = wrkbl, i__3 = *lda * *n + *n;
			if (*lwork >= MAX(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(I
R) is N by N */

			    ldwrku = *lda;
			    ldwrkr = *n;
			} else {

/*                    WORK(IU) is LDWRKU by N, WOR
K(IR) is N by N */

			    ldwrku = (*lwork - *n * *n - *n) / *n;
			    ldwrkr = *n;
			}
		    }
		    itau = ir + ldwrkr * *n;
		    iwork = itau + *n;

/*                 Compute A=Q*R   
                   (Workspace: need N*N+2*N, prefer N*N+N+
N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below i
t */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("U", n, n, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif

		    i__2 = *n - 1;
		    i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 1)
#endif

			    , &ldwrkr);

/*                 Generate Q in A   
                   (Workspace: need N*N+2*N, prefer N*N+N+
N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__2, &ierr);
		    ie = itau;
		    itauq = ie + *n;
		    itaup = itauq + *n;
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR)   
                   (Workspace: need N*N+4*N, prefer N*N+3*
N+2*N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R
   
                   (Workspace: need N*N+4*N, prefer N*N+3*
N+N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif

			    WORK(iwork), &i__2, &ierr);
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of R in WORK(IR)   
                   (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), dum, &
#endif

			    c__1, &WORK(ir), &ldwrkr, dum, &c__1, &WORK(iwork)
			    , info);
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vector
s of R in   
                   WORK(IR), storing result in WORK(IU) an
d copying to A   
                   (Workspace: need N*N+2*N, prefer N*N+M*
N+N) */

		    i__2 = *m;
		    i__3 = ldwrku;
		    for (i = 1; ldwrku < 0 ? i >= *m : i <= *m; i += ldwrku) {
/* Computing MIN */
			i__4 = *m - i + 1;
			chunk = MIN(i__4,ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif

				, lda, &WORK(ir), &ldwrkr, &c_b416, &WORK(iu),
				 &ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif

/* L10: */
		    }

		} else {

/*                 Insufficient workspace for a fast algor
ithm */

		    ie = 1;
		    itauq = ie + *n;
		    itaup = itauq + *n;
		    iwork = itaup + *n;

/*                 Bidiagonalize A   
                   (Workspace: need 3*N+M, prefer 3*N+(M+N
)*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A
   
                   (Workspace: need 4*N, prefer 3*N+N*NB) 
*/

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &
#endif

			    WORK(iwork), &i__3, &ierr);
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of A in A   
                   (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), dum, &
#endif

			    c__1, &A(1,1), lda, dum, &c__1, &WORK(iwork),
			     info);

		}

	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT=
'S' or 'A')   
                N left singular vectors to be overwritten on A
 and   
                N right singular vectors to be computed in VT 
  

   Computing MAX */
		i__3 = *n << 2;
		if (*lwork >= *n * *n + MAX(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorit
hm */

		    ir = 1;
/* Computing MAX */
		    i__3 = wrkbl, i__2 = *lda * *n + *n;
		    if (*lwork >= MAX(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR
) is LDA by N */

			ldwrku = *lda;
			ldwrkr = *lda;
		    } else /* if(complicated condition) */ {
/* Computing MAX */
			i__3 = wrkbl, i__2 = *lda * *n + *n;
			if (*lwork >= MAX(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WOR
K(IR) is N by N */

			    ldwrku = *lda;
			    ldwrkr = *n;
			} else {

/*                    WORK(IU) is LDWRKU by N and 
WORK(IR) is N by N */

			    ldwrku = (*lwork - *n * *n - *n) / *n;
			    ldwrkr = *n;
			}
		    }
		    itau = ir + ldwrkr * *n;
		    iwork = itau + *n;

/*                 Compute A=Q*R   
                   (Workspace: need N*N+2*N, prefer N*N+N+
N*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif

			    ldvt);
		    i__3 = *n - 1;
		    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("L", &i__3, &i__2, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("L", &i__3, &i__2, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("L", &i__3, &i__2, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif


/*                 Generate Q in A   
                   (Workspace: need N*N+2*N, prefer N*N+N+
N*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__3, &ierr);
		    ie = itau;
		    itauq = ie + *n;
		    itaup = itauq + *n;
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result t
o WORK(IR)   
                   (Workspace: need N*N+4*N, prefer N*N+3*
N+2*N*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif

			    WORK(itauq), &WORK(itaup), &WORK(iwork), &i__3, &
			    ierr);

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("L", n, n, &VT(1,1), ldvt, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("L", n, n, &VT(1,1), ldvt, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("L", n, n, &VT(1,1), ldvt, &WORK(ir), &
#endif

			    ldwrkr);

/*                 Generate left vectors bidiagonalizing R
 in WORK(IR)   
                   (Workspace: need N*N+4*N, prefer N*N+3*
N+N*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq), &
#endif

			    WORK(iwork), &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing 
R in VT   
                   (Workspace: need N*N+4*N-1, prefer N*N+
3*N+(N-1)*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif

			    &WORK(iwork), &i__3, &ierr);
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of R in WORK(IR) and c
omputing right   
                   singular vectors of R in VT   
                   (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(ir), &ldwrkr, dum, &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(ir), &ldwrkr, dum, &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(ir), &ldwrkr, dum, &c__1, 
#endif

			    &WORK(iwork), info);
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vector
s of R in   
                   WORK(IR), storing result in WORK(IU) an
d copying to A   
                   (Workspace: need N*N+2*N, prefer N*N+M*
N+N) */

		    i__3 = *m;
		    i__2 = ldwrku;
		    for (i = 1; ldwrku < 0 ? i >= *m : i <= *m; i += ldwrku) {
/* Computing MIN */
			i__4 = *m - i + 1;
			chunk = MIN(i__4,ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", &chunk, n, n, &c_b438, &A(i,1)
#endif

				, lda, &WORK(ir), &ldwrkr, &c_b416, &WORK(iu),
				 &ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", &chunk, n, &WORK(iu), &ldwrku, &A(i,1), lda);
#endif

/* L20: */
		    }

		} else {

/*                 Insufficient workspace for a fast algor
ithm */

		    itau = 1;
		    iwork = itau + *n;

/*                 Compute A=Q*R   
                   (Workspace: need 2*N, prefer N+N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif

			    ldvt);
		    i__2 = *n - 1;
		    i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif


/*                 Generate Q in A   
                   (Workspace: need 2*N, prefer N+N*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__2, &ierr);
		    ie = itau;
		    itauq = ie + *n;
		    itaup = itauq + *n;
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT   
                   (Workspace: need 4*N, prefer 3*N+2*N*NB
) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), &
#endif

			    WORK(itauq), &WORK(itaup), &WORK(iwork), &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiago
nalizing R   
                   (Workspace: need 3*N+M, prefer 3*N+M*NB
) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qormbr("Q", "R", "N", m, n, n, &VT(1,1), ldvt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, &
#endif

			    WORK(itauq), &A(1,1), lda, &WORK(iwork), &
			    i__2, &ierr);

/*                 Generate right vectors bidiagonalizing 
R in VT   
                   (Workspace: need 4*N-1, prefer 3*N+(N-1
)*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), 
#endif

			    &WORK(iwork), &i__2, &ierr);
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of A in A and computin
g right   
                   singular vectors of A in VT   
                   (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif

			    WORK(iwork), info);

		}

	    } else if (wntus) {

		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S',
 JOBVT='N')   
                   N left singular vectors to be computed 
in U and   
                   no right singular vectors to be compute
d   

   Computing MAX */
		    i__2 = *n << 2;
		    if (*lwork >= *n * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			ir = 1;
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

			    ldwrkr = *lda;
			} else {

/*                       WORK(IR) is N by N */

			    ldwrkr = *n;
			}
			itau = ir + ldwrkr * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif

				ldwrkr);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif

				+ 1), &ldwrkr);

/*                    Generate Q in A   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Generate left vectors bidiagonal
izing R in WORK(IR)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IR
)   
                      (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif

				dum, &c__1, &WORK(ir), &ldwrkr, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply Q in A by left singular
 vectors of R in   
                      WORK(IR), storing result in U   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif

				&WORK(ir), &ldwrkr, &c_b416, &U(1,1), 
				ldu);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Zero out below R in A */

			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif


/*                    Bidiagonalize R in A   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors 
bidiagonalizing R   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif

				WORK(itauq), &U(1,1), ldu, &WORK(iwork), 
				&i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif

				dum, &c__1, &U(1,1), ldu, dum, &c__1, &
				WORK(iwork), info);

		    }

		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S',
 JOBVT='O')   
                   N left singular vectors to be computed 
in U and   
                   N right singular vectors to be overwrit
ten on A   

   Computing MAX */
		    i__2 = *n << 2;
		    if (*lwork >= (*n << 1) * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and 
WORK(IR) is LDA by N */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *lda;
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and 
WORK(IR) is N by N */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *n;
			} else {

/*                       WORK(IU) is N by N and WO
RK(IR) is N by N */

			    ldwrku = *n;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *n;
			}
			itau = ir + ldwrkr * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R   
                      (Workspace: need 2*N*N+2*N, pref
er 2*N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ 1), &ldwrku);

/*                    Generate Q in A   
                      (Workspace: need 2*N*N+2*N, pref
er 2*N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), cop
ying result to   
                      WORK(IR)   
                      (Workspace: need 2*N*N+4*N,   
                                  prefer 2*N*N+3*N+2*N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif

				ldwrkr);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IU)   
                      (Workspace: need 2*N*N+4*N, pref
er 2*N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in WORK(IR)   
                      (Workspace: need 2*N*N+4*N-1,   
                                  prefer 2*N*N+3*N+(N-
1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IU
) and computing   
                      right singular vectors of R in W
ORK(IR)   
                      (Workspace: need 2*N*N+BDSPAC) 
*/


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				ir), &ldwrkr, &WORK(iu), &ldwrku, dum, &c__1, 
				&WORK(iwork), info);

/*                    Multiply Q in A by left singular
 vectors of R in   
                      WORK(IU), storing result in U   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif

				&WORK(iu), &ldwrku, &c_b416, &U(1,1), 
				ldu);

/*                    Copy right singular vectors of R
 to A   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif

				lda);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Zero out below R in A */

			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif


/*                    Bidiagonalize R in A   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors 
bidiagonalizing R   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif

				WORK(itauq), &U(1,1), ldu, &WORK(iwork), 
				&i__2, &ierr);

/*                    Generate right vectors bidiagona
lizing R in A   
                      (Workspace: need 4*N-1, prefer 3
*N+(N-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in A   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif

				 &WORK(iwork), info);

		    }

		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S',
 JOBVT='S'   
                           or 'A')   
                   N left singular vectors to be computed 
in U and   
                   N right singular vectors to be computed
 in VT   

   Computing MAX */
		    i__2 = *n << 2;
		    if (*lwork >= *n * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

			    ldwrku = *lda;
			} else {

/*                       WORK(IU) is N by N */

			    ldwrku = *n;
			}
			itau = iu + ldwrku * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ 1), &ldwrku);

/*                    Generate Q in A   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), cop
ying result to VT   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif

				 ldvt);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IU)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in VT   
                      (Workspace: need N*N+4*N-1,   
                                  prefer N*N+3*N+(N-1)
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif

				itaup), &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IU
) and computing   
                      right singular vectors of R in V
T   
                      (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif

				c__1, &WORK(iwork), info);

/*                    Multiply Q in A by left singular
 vectors of R in   
                      WORK(IU), storing result in U   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &A(1,1), lda, 
#endif

				&WORK(iu), &ldwrku, &c_b416, &U(1,1), 
				ldu);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, n, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy R to VT, zeroing out below 
it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif

			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif

				&WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagon
alizing vectors   
                      in VT   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif

				&WORK(itauq), &U(1,1), ldu, &WORK(iwork),
				 &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in VT   
                      (Workspace: need 4*N-1, prefer 3
*N+(N-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif

				itaup), &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		}

	    } else if (wntua) {

		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A',
 JOBVT='N')   
                   M left singular vectors to be computed 
in U and   
                   no right singular vectors to be compute
d   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= *n * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			ir = 1;
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

			    ldwrkr = *lda;
			} else {

/*                       WORK(IR) is N by N */

			    ldwrkr = *n;
			}
			itau = ir + ldwrkr * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Copy R to WORK(IR), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(ir), &
#endif

				ldwrkr);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif

				+ 1), &ldwrkr);

/*                    Generate Q in U   
                      (Workspace: need N*N+N+M, prefer
 N*N+N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IR)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IR
)   
                      (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, &c__0, n, &c__0, &S(1), &WORK(ie), 
#endif

				dum, &c__1, &WORK(ir), &ldwrkr, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply Q in U by left singular
 vectors of R in   
                      WORK(IR), storing result in A   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif

				&WORK(ir), &ldwrkr, &c_b416, &A(1,1), 
				lda);

/*                    Copy left singular vectors of A 
from A to U */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need N+M, prefer N+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Zero out below R in A */

			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif


/*                    Bidiagonalize R in A   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagon
alizing vectors   
                      in A   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif

				WORK(itauq), &U(1,1), ldu, &WORK(iwork), 
				&i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, &c__0, m, &c__0, &S(1), &WORK(ie), 
#endif

				dum, &c__1, &U(1,1), ldu, dum, &c__1, &
				WORK(iwork), info);

		    }

		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A',
 JOBVT='O')   
                   M left singular vectors to be computed 
in U and   
                   N right singular vectors to be overwrit
ten on A   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= (*n << 1) * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and 
WORK(IR) is LDA by N */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *lda;
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and 
WORK(IR) is N by N */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *n;
			} else {

/*                       WORK(IU) is N by N and WO
RK(IR) is N by N */

			    ldwrku = *n;
			    ir = iu + ldwrku * *n;
			    ldwrkr = *n;
			}
			itau = ir + ldwrkr * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N*N+2*N, pref
er 2*N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need 2*N*N+N+M, pref
er 2*N*N+N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ 1), &ldwrku);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), cop
ying result to   
                      WORK(IR)   
                      (Workspace: need 2*N*N+4*N,   
                                  prefer 2*N*N+3*N+2*N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &WORK(iu), &ldwrku, &WORK(ir), &
#endif

				ldwrkr);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IU)   
                      (Workspace: need 2*N*N+4*N, pref
er 2*N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in WORK(IR)   
                      (Workspace: need 2*N*N+4*N-1,   
                                  prefer 2*N*N+3*N+(N-
1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IU
) and computing   
                      right singular vectors of R in W
ORK(IR)   
                      (Workspace: need 2*N*N+BDSPAC) 
*/


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				ir), &ldwrkr, &WORK(iu), &ldwrku, dum, &c__1, 
				&WORK(iwork), info);

/*                    Multiply Q in U by left singular
 vectors of R in   
                      WORK(IU), storing result in A   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif

				&WORK(iu), &ldwrku, &c_b416, &A(1,1), 
				lda);

/*                    Copy left singular vectors of A 
from A to U */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Copy right singular vectors of R
 from WORK(IR) to A */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", n, n, &WORK(ir), &ldwrkr, &A(1,1), 
#endif

				lda);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need N+M, prefer N+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Zero out below R in A */

			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &A(2,1), lda);
#endif


/*                    Bidiagonalize R in A   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagon
alizing vectors   
                      in A   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &A(1,1), lda, &
#endif

				WORK(itauq), &U(1,1), ldu, &WORK(iwork), 
				&i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in A   
                      (Workspace: need 4*N-1, prefer 3
*N+(N-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in A   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1,
#endif

				 &WORK(iwork), info);

		    }

		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A',
 JOBVT='S'   
                           or 'A')   
                   M left singular vectors to be computed 
in U and   
                   N right singular vectors to be computed
 in VT   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= *n * *n + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

			    ldwrku = *lda;
			} else {

/*                       WORK(IU) is N by N */

			    ldwrku = *n;
			}
			itau = iu + ldwrku * *n;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need N*N+2*N, prefer
 N*N+N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need N*N+N+M, prefer
 N*N+N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ 1), &ldwrku);
			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), cop
ying result to VT   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &WORK(iu), &ldwrku, &VT(1,1),
#endif

				 ldvt);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IU)   
                      (Workspace: need N*N+4*N, prefer
 N*N+3*N+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", n, n, n, &WORK(iu), &ldwrku, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in VT   
                      (Workspace: need N*N+4*N-1,   
                                  prefer N*N+3*N+(N-1)
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif

				itaup), &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of R in WORK(IU
) and computing   
                      right singular vectors of R in V
T   
                      (Workspace: need N*N+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, n, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &WORK(iu), &ldwrku, dum, &
#endif

				c__1, &WORK(iwork), info);

/*                    Multiply Q in U by left singular
 vectors of R in   
                      WORK(IU), storing result in A   
                      (Workspace: need N*N) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, n, &c_b438, &U(1,1), ldu, 
#endif

				&WORK(iu), &ldwrku, &c_b416, &A(1,1), 
				lda);

/*                    Copy left singular vectors of A 
from A to U */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to
 U   
                      (Workspace: need 2*N, prefer N+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, n, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), 
#endif

				ldu);

/*                    Generate Q in U   
                      (Workspace: need N+M, prefer N+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgqr(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgqr_(m, m, n, &U(1,1), ldu, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out
 below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);
			i__2 = *n - 1;
			i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("L", &i__2, &i__3, &c_b416, &c_b416, &VT(2,1), ldvt);
#endif

			ie = itau;
			itauq = ie + *n;
			itaup = itauq + *n;
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT   
                      (Workspace: need 4*N, prefer 3*N
+2*N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(n, n, &VT(1,1), ldvt, &S(1), &WORK(ie), 
#endif

				&WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagon
alizing vectors   
                      in VT   
                      (Workspace: need 3*N+M, prefer 3
*N+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("Q", "R", "N", m, n, n, &VT(1,1), ldvt, 
#endif

				&WORK(itauq), &U(1,1), ldu, &WORK(iwork),
				 &i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in VT   
                      (Workspace: need 4*N-1, prefer 3
*N+(N-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(
#endif

				itaup), &WORK(iwork), &i__2, &ierr);
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", n, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		}

	    }

	} else {

/*           M .LT. MNTHR   

             Path 10 (M at least N, but not much larger)   
             Reduce to bidiagonal form without QR decomposition */

	    ie = 1;
	    itauq = ie + *n;
	    itaup = itauq + *n;
	    iwork = itaup + *n;

/*           Bidiagonalize A   
             (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgebrd(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif

		    WORK(itaup), &WORK(iwork), &i__2, &ierr);
	    if (wntuas) {

/*              If left singular vectors desired in U, copy re
sult to U   
                and generate left bidiagonalizing vectors in U
   
                (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB) 
*/


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("L", m, n, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("L", m, n, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("L", m, n, &A(1,1), lda, &U(1,1), ldu);
#endif

		if (wntus) {
		    ncu = *n;
		}
		if (wntua) {
		    ncu = *m;
		}
		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("Q", m, &ncu, n, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("Q", m, &ncu, n, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("Q", m, &ncu, n, &U(1,1), ldu, &WORK(itauq), &
#endif

			WORK(iwork), &i__2, &ierr);
	    }
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy 
result to   
                VT and generate right bidiagonalizing vectors 
in VT   
                (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) 
*/


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("U", n, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("U", n, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("P", n, n, n, &VT(1,1), ldvt, &WORK(itaup), &
#endif

			WORK(iwork), &i__2, &ierr);
	    }
	    if (wntuo) {

/*              If left singular vectors desired in A, generat
e left   
                bidiagonalizing vectors in A   
                (Workspace: need 4*N, prefer 3*N+N*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("Q", m, n, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif

			iwork), &i__2, &ierr);
	    }
	    if (wntvo) {

/*              If right singular vectors desired in A, genera
te right   
                bidiagonalizing vectors in A   
                (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) 
*/

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(
#endif

			iwork), &i__2, &ierr);
	    }
	    iwork = ie + *n;
	    if (wntuas || wntuo) {
		nru = *m;
	    }
	    if (wntun) {
		nru = 0;
	    }
	    if (wntvas || wntvo) {
		ncvt = *n;
	    }
	    if (wntvn) {
		ncvt = 0;
	    }
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in U and computing right
 singular   
                vectors in VT   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif

			WORK(iwork), info);
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in U and computing right
 singular   
                vectors in A   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif

			iwork), info);
	    } else {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in A and computing right
 singular   
                vectors in VT   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("U", n, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif

			WORK(iwork), info);
	    }

	}

    } else {

/*        A has more columns than rows. If A has sufficiently more   
          columns than rows, first reduce using the LQ decomposition (
if   
          sufficient workspace available) */

	if (*n >= mnthr) {

	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N')   
                No right singular vectors to be computed */

		itau = 1;
		iwork = itau + *m;

/*              Compute A=L*Q   
                (Workspace: need 2*M, prefer M+M*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &
#endif

			i__2, &ierr);

/*              Zero out above L */

		i__2 = *m - 1;
		i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif

		ie = 1;
		itauq = ie + *m;
		itaup = itauq + *m;
		iwork = itaup + *m;

/*              Bidiagonalize L in A   
                (Workspace: need 4*M, prefer 3*M+2*M*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgebrd(m, m, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif

			itauq), &WORK(itaup), &WORK(iwork), &i__2, &ierr);
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, gener
ate Q   
                   (Workspace: need 4*M, prefer 3*M+M*NB) 
*/

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", m, m, m, &A(1,1), lda, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq), &
#endif

			    WORK(iwork), &i__2, &ierr);
		}
		iwork = ie + *m;
		nru = 0;
		if (wntuo || wntuas) {
		    nru = *m;
		}

/*              Perform bidiagonal QR iteration, computing lef
t singular   
                vectors of A in A if desired   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("U", m, &c__0, &nru, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("U", m, &c__0, &nru, &c__0, &S(1), &WORK(ie), dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("U", m, &c__0, &nru, &c__0, &S(1), &WORK(ie), dum, &
#endif

			c__1, &A(1,1), lda, dum, &c__1, &WORK(iwork), 
			info);

/*              If left singular vectors desired in U, copy th
em there */

		if (wntuas) {

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("F", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("F", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("F", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif

		}

	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT=
'O')   
                M right singular vectors to be overwritten on 
A and   
                no left singular vectors to be computed   

   Computing MAX */
		i__2 = *m << 2;
		if (*lwork >= *m * *m + MAX(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorit
hm */

		    ir = 1;
/* Computing MAX */
		    i__2 = wrkbl, i__3 = *lda * *n + *m;
		    if (*lwork >= MAX(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR
) is LDA by M */

			ldwrku = *lda;
			chunk = *n;
			ldwrkr = *lda;
		    } else /* if(complicated condition) */ {
/* Computing MAX */
			i__2 = wrkbl, i__3 = *lda * *n + *m;
			if (*lwork >= MAX(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WOR
K(IR) is M by M */

			    ldwrku = *lda;
			    chunk = *n;
			    ldwrkr = *m;
			} else {

/*                    WORK(IU) is M by CHUNK and W
ORK(IR) is M by M */

			    ldwrku = *m;
			    chunk = (*lwork - *m * *m - *m) / *m;
			    ldwrkr = *m;
			}
		    }
		    itau = ir + ldwrkr * *m;
		    iwork = itau + *m;

/*                 Compute A=L*Q   
                   (Workspace: need M*M+2*M, prefer M*M+M+
M*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above i
t */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("L", m, m, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &ldwrkr);
#endif

		    i__2 = *m - 1;
		    i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir + 
#endif

			    ldwrkr), &ldwrkr);

/*                 Generate Q in A   
                   (Workspace: need M*M+2*M, prefer M*M+M+
M*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__2, &ierr);
		    ie = itau;
		    itauq = ie + *m;
		    itaup = itauq + *m;
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR)   
                   (Workspace: need M*M+4*M, prefer M*M+3*
M+2*M*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing 
L   
                   (Workspace: need M*M+4*M-1, prefer M*M+
3*M+(M-1)*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif

			    WORK(iwork), &i__2, &ierr);
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, comput
ing right   
                   singular vectors of L in WORK(IR)   
                   (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &WORK(
#endif

			    ir), &ldwrkr, dum, &c__1, dum, &c__1, &WORK(iwork)
			    , info);
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in
 WORK(IR) by Q   
                   in A, storing result in WORK(IU) and co
pying to A   
                   (Workspace: need M*M+2*M, prefer M*M+M*
N+M) */

		    i__2 = *n;
		    i__3 = chunk;
		    for (i = 1; chunk < 0 ? i >= *n : i <= *n; i += chunk) {
/* Computing MIN */
			i__4 = *n - i + 1;
			blk = MIN(i__4,chunk);

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif

				ldwrkr, &A(1,i), lda, &c_b416, &
				WORK(iu), &ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif

/* L30: */
		    }

		} else {

/*                 Insufficient workspace for a fast algor
ithm */

		    ie = 1;
		    itauq = ie + *m;
		    itaup = itauq + *m;
		    iwork = itaup + *m;

/*                 Bidiagonalize A   
                   (Workspace: need 3*M+N, prefer 3*M+(M+N
)*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing 
A   
                   (Workspace: need 4*M, prefer 3*M+M*NB) 
*/

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", m, n, m, &A(1,1), lda, &WORK(itaup), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &
#endif

			    WORK(iwork), &i__3, &ierr);
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, comput
ing right   
                   singular vectors of A in A   
                   (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("L", m, n, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("L", m, n, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("L", m, n, &c__0, &c__0, &S(1), &WORK(ie), &A(1,1), lda, dum, &c__1, dum, &c__1, &WORK(
#endif

			    iwork), info);

		}

	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A',
 JOBVT='O')   
                M right singular vectors to be overwritten on 
A and   
                M left singular vectors to be computed in U   

   Computing MAX */
		i__3 = *m << 2;
		if (*lwork >= *m * *m + MAX(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorit
hm */

		    ir = 1;
/* Computing MAX */
		    i__3 = wrkbl, i__2 = *lda * *n + *m;
		    if (*lwork >= MAX(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR
) is LDA by M */

			ldwrku = *lda;
			chunk = *n;
			ldwrkr = *lda;
		    } else /* if(complicated condition) */ {
/* Computing MAX */
			i__3 = wrkbl, i__2 = *lda * *n + *m;
			if (*lwork >= MAX(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WOR
K(IR) is M by M */

			    ldwrku = *lda;
			    chunk = *n;
			    ldwrkr = *m;
			} else {

/*                    WORK(IU) is M by CHUNK and W
ORK(IR) is M by M */

			    ldwrku = *m;
			    chunk = (*lwork - *m * *m - *m) / *m;
			    ldwrkr = *m;
			}
		    }
		    itau = ir + ldwrkr * *m;
		    iwork = itau + *m;

/*                 Compute A=L*Q   
                   (Workspace: need M*M+2*M, prefer M*M+M+
M*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif

		    i__3 = *m - 1;
		    i__2 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("U", &i__3, &i__2, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("U", &i__3, &i__2, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("U", &i__3, &i__2, &c_b416, &c_b416, &U(1,2), ldu);
#endif


/*                 Generate Q in A   
                   (Workspace: need M*M+2*M, prefer M*M+M+
M*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__3, &ierr);
		    ie = itau;
		    itauq = ie + *m;
		    itaup = itauq + *m;
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to
 WORK(IR)   
                   (Workspace: need M*M+4*M, prefer M*M+3*
M+2*M*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__3, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("U", m, m, &U(1,1), ldu, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("U", m, m, &U(1,1), ldu, &WORK(ir), &ldwrkr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("U", m, m, &U(1,1), ldu, &WORK(ir), &ldwrkr);
#endif


/*                 Generate right vectors bidiagonalizing 
L in WORK(IR)   
                   (Workspace: need M*M+4*M-1, prefer M*M+
3*M+(M-1)*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup), &
#endif

			    WORK(iwork), &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing L
 in U   
                   (Workspace: need M*M+4*M, prefer M*M+3*
M+M*NB) */

		    i__3 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif

			    WORK(iwork), &i__3, &ierr);
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of L in U, and computi
ng right   
                   singular vectors of L in WORK(IR)   
                   (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(ir), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(ir), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(ir), 
#endif

			    &ldwrkr, &U(1,1), ldu, dum, &c__1, &WORK(
			    iwork), info);
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in
 WORK(IR) by Q   
                   in A, storing result in WORK(IU) and co
pying to A   
                   (Workspace: need M*M+2*M, prefer M*M+M*
N+M)) */

		    i__3 = *n;
		    i__2 = chunk;
		    for (i = 1; chunk < 0 ? i >= *n : i <= *n; i += chunk) {
/* Computing MIN */
			i__4 = *n - i + 1;
			blk = MIN(i__4,chunk);

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, &blk, m, &c_b438, &WORK(ir), &
#endif

				ldwrkr, &A(1,i), lda, &c_b416, &
				WORK(iu), &ldwrku);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, &blk, &WORK(iu), &ldwrku, &A(1,i), lda);
#endif

/* L40: */
		    }

		} else {

/*                 Insufficient workspace for a fast algor
ithm */

		    itau = 1;
		    iwork = itau + *m;

/*                 Compute A=L*Q   
                   (Workspace: need 2*M, prefer M+M*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork)
#endif

			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */


#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif

		    i__2 = *m - 1;
		    i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif


/*                 Generate Q in A   
                   (Workspace: need 2*M, prefer M+M*NB) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &WORK(
#endif

			    iwork), &i__2, &ierr);
		    ie = itau;
		    itauq = ie + *m;
		    itaup = itauq + *m;
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U   
                   (Workspace: need 4*M, prefer 3*M+2*M*NB
) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgebrd(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &WORK(
#endif

			    itauq), &WORK(itaup), &WORK(iwork), &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing 
L by Q in A   
                   (Workspace: need 3*M+N, prefer 3*M+N*NB
) */

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qormbr("P", "L", "T", m, n, m, &U(1,1), ldu, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &WORK(
#endif

			    itaup), &A(1,1), lda, &WORK(iwork), &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing L
 in U   
                   (Workspace: need 4*M, prefer 3*M+M*NB) 
*/

		    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq), &
#endif

			    WORK(iwork), &i__2, &ierr);
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, comput
ing left   
                   singular vectors of A in U and computin
g right   
                   singular vectors of A in A   
                   (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		    dbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qbdsqr("U", m, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &
#endif

			    WORK(iwork), info);

		}

	    } else if (wntvs) {

		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N',
 JOBVT='S')   
                   M right singular vectors to be computed
 in VT and   
                   no left singular vectors to be computed
   

   Computing MAX */
		    i__2 = *m << 2;
		    if (*lwork >= *m * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			ir = 1;
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

			    ldwrkr = *lda;
			} else {

/*                       WORK(IR) is M by M */

			    ldwrkr = *m;
			}
			itau = ir + ldwrkr * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out 
above it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif

				ldwrkr);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif

				+ ldwrkr), &ldwrkr);

/*                    Generate Q in A   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR)   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Generate right vectors bidiagona
lizing L in   
                      WORK(IR)   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+(M-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing right   
                      singular vectors of L in WORK(IR
)   
                      (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif

				WORK(ir), &ldwrkr, dum, &c__1, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IR) by   
                      Q in A, storing result in VT   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif

				 &A(1,1), lda, &c_b416, &VT(1,1), 
				ldvt);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy result to VT */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Zero out above L in A */

			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif


/*                    Bidiagonalize L in A   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right vectors bidiagona
lizing L by Q in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif

				VT(1,1), ldvt, dum, &c__1, dum, &c__1, &
				WORK(iwork), info);

		    }

		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O',
 JOBVT='S')   
                   M right singular vectors to be computed
 in VT and   
                   M left singular vectors to be overwritt
en on A   

   Computing MAX */
		    i__2 = *m << 2;
		    if (*lwork >= (*m << 1) * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and 
WORK(IR) is LDA by M */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *lda;
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and 
WORK(IR) is M by M */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *m;
			} else {

/*                       WORK(IU) is M by M and WO
RK(IR) is M by M */

			    ldwrku = *m;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *m;
			}
			itau = ir + ldwrkr * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q   
                      (Workspace: need 2*M*M+2*M, pref
er 2*M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out 
below it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ ldwrku), &ldwrku);

/*                    Generate Q in A   
                      (Workspace: need 2*M*M+2*M, pref
er 2*M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), cop
ying result to   
                      WORK(IR)   
                      (Workspace: need 2*M*M+4*M,   
                                  prefer 2*M*M+3*M+2*M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif

				ldwrkr);

/*                    Generate right bidiagonalizing v
ectors in WORK(IU)   
                      (Workspace: need 2*M*M+4*M-1,   
                                  prefer 2*M*M+3*M+(M-
1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IR)   
                      (Workspace: need 2*M*M+4*M, pref
er 2*M*M+3*M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of L in WORK(IR
) and computing   
                      right singular vectors of L in W
ORK(IU)   
                      (Workspace: need 2*M*M+BDSPAC) 
*/


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				iu), &ldwrku, &WORK(ir), &ldwrkr, dum, &c__1, 
				&WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IU) by   
                      Q in A, storing result in VT   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif

				 &A(1,1), lda, &c_b416, &VT(1,1), 
				ldvt);

/*                    Copy left singular vectors of L 
to A   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif

				lda);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Zero out above L in A */

			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif


/*                    Bidiagonalize L in A   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right vectors bidiagona
lizing L by Q in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors of L in A   
                      (Workspace: need 4*M, prefer 3*M
+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 compute left   
                      singular vectors of A in A and c
ompute right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' 
or 'A',   
                           JOBVT='S')   
                   M right singular vectors to be computed
 in VT and   
                   M left singular vectors to be computed 
in U   

   Computing MAX */
		    i__2 = *m << 2;
		    if (*lwork >= *m * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

			    ldwrku = *lda;
			} else {

/*                       WORK(IU) is LDA by M */

			    ldwrku = *m;
			}
			itau = iu + ldwrku * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out 
above it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ ldwrku), &ldwrku);

/*                    Generate Q in A   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &A(1,1), lda, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), cop
ying result to U   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif

				ldu);

/*                    Generate right bidiagonalizing v
ectors in WORK(IU)   
                      (Workspace: need M*M+4*M-1,   
                                  prefer M*M+3*M+(M-1)
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in U   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of L in U and c
omputing right   
                      singular vectors of L in WORK(IU
)   
                      (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				iu), &ldwrku, &U(1,1), ldu, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IU) by   
                      Q in A, storing result in VT   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif

				 &A(1,1), lda, &c_b416, &VT(1,1), 
				ldvt);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(m, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy L to U, zeroing out above i
t */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &U(1,1), 
#endif

				ldu);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif

			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in U   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing v
ectors in U by Q   
                      in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in U   
                      (Workspace: need 4*M, prefer 3*M
+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		}

	    } else if (wntva) {

		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N',
 JOBVT='A')   
                   N right singular vectors to be computed
 in VT and   
                   no left singular vectors to be computed
   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= *m * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			ir = 1;
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

			    ldwrkr = *lda;
			} else {

/*                       WORK(IR) is M by M */

			    ldwrkr = *m;
			}
			itau = ir + ldwrkr * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Copy L to WORK(IR), zeroing out 
above it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(ir), &
#endif

				ldwrkr);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(ir 
#endif

				+ ldwrkr), &ldwrkr);

/*                    Generate Q in VT   
                      (Workspace: need M*M+M+N, prefer
 M*M+M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR)   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(ir), &ldwrkr, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Generate right bidiagonalizing v
ectors in WORK(IR)   
                      (Workspace: need M*M+4*M-1,   
                                  prefer M*M+3*M+(M-1)
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(ir), &ldwrkr, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing right   
                      singular vectors of L in WORK(IR
)   
                      (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, &c__0, &c__0, &S(1), &WORK(ie), &
#endif

				WORK(ir), &ldwrkr, dum, &c__1, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IR) by   
                      Q in VT, storing result in A   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(ir), &ldwrkr,
#endif

				 &VT(1,1), ldvt, &c_b416, &A(1,1), 
				lda);

/*                    Copy right singular vectors of A
 from A to VT */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need M+N, prefer M+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Zero out above L in A */

			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif


/*                    Bidiagonalize L in A   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing v
ectors in A by Q   
                      in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, &c__0, &c__0, &S(1), &WORK(ie), &
#endif

				VT(1,1), ldvt, dum, &c__1, dum, &c__1, &
				WORK(iwork), info);

		    }

		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O',
 JOBVT='A')   
                   N right singular vectors to be computed
 in VT and   
                   M left singular vectors to be overwritt
en on A   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= (*m << 1) * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and 
WORK(IR) is LDA by M */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *lda;
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and 
WORK(IR) is M by M */

			    ldwrku = *lda;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *m;
			} else {

/*                       WORK(IU) is M by M and WO
RK(IR) is M by M */

			    ldwrku = *m;
			    ir = iu + ldwrku * *m;
			    ldwrkr = *m;
			}
			itau = ir + ldwrkr * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M*M+2*M, pref
er 2*M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need 2*M*M+M+N, pref
er 2*M*M+M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out 
above it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ ldwrku), &ldwrku);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), cop
ying result to   
                      WORK(IR)   
                      (Workspace: need 2*M*M+4*M,   
                                  prefer 2*M*M+3*M+2*M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &WORK(iu), &ldwrku, &WORK(ir), &
#endif

				ldwrkr);

/*                    Generate right bidiagonalizing v
ectors in WORK(IU)   
                      (Workspace: need 2*M*M+4*M-1,   
                                  prefer 2*M*M+3*M+(M-
1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in WORK(IR)   
                      (Workspace: need 2*M*M+4*M, pref
er 2*M*M+3*M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &WORK(ir), &ldwrkr, &WORK(itauq)
#endif

				, &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of L in WORK(IR
) and computing   
                      right singular vectors of L in W
ORK(IU)   
                      (Workspace: need 2*M*M+BDSPAC) 
*/


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				iu), &ldwrku, &WORK(ir), &ldwrkr, dum, &c__1, 
				&WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IU) by   
                      Q in VT, storing result in A   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif

				 &VT(1,1), ldvt, &c_b416, &A(1,1), 
				lda);

/*                    Copy right singular vectors of A
 from A to VT */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Copy left singular vectors of A 
from WORK(IR) to A */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, m, &WORK(ir), &ldwrkr, &A(1,1), 
#endif

				lda);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need M+N, prefer M+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Zero out above L in A */

			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &A(1,2), lda);
#endif


/*                    Bidiagonalize L in A   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &A(1,1), lda, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing v
ectors in A by Q   
                      in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &A(1,1), lda, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in A   
                      (Workspace: need 4*M, prefer 3*M
+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &A(1,1), lda, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in A and c
omputing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' 
or 'A',   
                           JOBVT='A')   
                   N right singular vectors to be computed
 in VT and   
                   M left singular vectors to be computed 
in U   

   Computing MAX */
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = MAX(i__2,i__3);
		    if (*lwork >= *m * *m + MAX(i__2,bdspac)) {

/*                    Sufficient workspace for a fast 
algorithm */

			iu = 1;
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

			    ldwrku = *lda;
			} else {

/*                       WORK(IU) is M by M */

			    ldwrku = *m;
			}
			itau = iu + ldwrku * *m;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need M*M+2*M, prefer
 M*M+M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need M*M+M+N, prefer
 M*M+M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out 
above it */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &WORK(iu), &
#endif

				ldwrku);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &WORK(iu 
#endif

				+ ldwrku), &ldwrku);
			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), cop
ying result to U   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &WORK(iu), &ldwrku, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &WORK(iu), &ldwrku, &U(1,1), 
#endif

				ldu);

/*                    Generate right bidiagonalizing v
ectors in WORK(IU)   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+(M-1)*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("P", m, m, m, &WORK(iu), &ldwrku, &WORK(itaup)
#endif

				, &WORK(iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in U   
                      (Workspace: need M*M+4*M, prefer
 M*M+3*M+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of L in U and c
omputing right   
                      singular vectors of L in WORK(IU
)   
                      (Workspace: need M*M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, m, m, &c__0, &S(1), &WORK(ie), &WORK(
#endif

				iu), &ldwrku, &U(1,1), ldu, dum, &c__1, &
				WORK(iwork), info);

/*                    Multiply right singular vectors 
of L in WORK(IU) by   
                      Q in VT, storing result in A   
                      (Workspace: need M*M) */


#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("N", "N", m, n, m, &c_b438, &WORK(iu), &ldwrku,
#endif

				 &VT(1,1), ldvt, &c_b416, &A(1,1), 
				lda);

/*                    Copy right singular vectors of A
 from A to VT */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("F", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

		    } else {

/*                    Insufficient workspace for a fas
t algorithm */

			itau = 1;
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to
 VT   
                      (Workspace: need 2*M, prefer M+M
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(
#endif

				iwork), &i__2, &ierr);

#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), 
#endif

				ldvt);

/*                    Generate Q in VT   
                      (Workspace: need M+N, prefer M+N
*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorglq(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorglq_(n, n, m, &VT(1,1), ldvt, &WORK(itau), &
#endif

				WORK(iwork), &i__2, &ierr);

/*                    Copy L to U, zeroing out above i
t */


#ifdef PETSC_PREFIX_SUFFIX
			dlacpy_("L", m, m, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlacpy("L", m, m, &A(1,1), lda, &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlacpy_("L", m, m, &A(1,1), lda, &U(1,1), 
#endif

				ldu);
			i__2 = *m - 1;
			i__3 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlaset("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlaset_("U", &i__2, &i__3, &c_b416, &c_b416, &U(1,2), ldu);
#endif

			ie = itau;
			itauq = ie + *m;
			itaup = itauq + *m;
			iwork = itaup + *m;

/*                    Bidiagonalize L in U   
                      (Workspace: need 4*M, prefer 3*M
+2*M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgebrd(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgebrd_(m, m, &U(1,1), ldu, &S(1), &WORK(ie), &
#endif

				WORK(itauq), &WORK(itaup), &WORK(iwork), &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing v
ectors in U by Q   
                      in VT   
                      (Workspace: need 3*M+N, prefer 3
*M+N*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qormbr("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qormbr_("P", "L", "T", m, n, m, &U(1,1), ldu, &
#endif

				WORK(itaup), &VT(1,1), ldvt, &WORK(
				iwork), &i__2, &ierr);

/*                    Generate left bidiagonalizing ve
ctors in U   
                      (Workspace: need 4*M, prefer 3*M
+M*NB) */

			i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qorgbr("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qorgbr_("Q", m, m, m, &U(1,1), ldu, &WORK(itauq),
#endif

				 &WORK(iwork), &i__2, &ierr);
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration,
 computing left   
                      singular vectors of A in U and c
omputing right   
                      singular vectors of A in VT   
                      (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
			dbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qbdsqr("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qbdsqr_("U", m, n, m, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &
#endif

				c__1, &WORK(iwork), info);

		    }

		}

	    }

	} else {

/*           N .LT. MNTHR   

             Path 10t(N greater than M, but not much larger)   
             Reduce to bidiagonal form without LQ decomposition */

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    iwork = itaup + *m;

/*           Bidiagonalize A   
             (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgebrd(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif

		    WORK(itaup), &WORK(iwork), &i__2, &ierr);
	    if (wntuas) {

/*              If left singular vectors desired in U, copy re
sult to U   
                and generate left bidiagonalizing vectors in U
   
                (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) 
*/


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("L", m, m, &A(1,1), lda, &U(1,1), ldu);
#endif

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("Q", m, m, n, &U(1,1), ldu, &WORK(itauq), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("Q", m, m, n, &U(1,1), ldu, &WORK(itauq), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("Q", m, m, n, &U(1,1), ldu, &WORK(itauq), &WORK(
#endif

			iwork), &i__2, &ierr);
	    }
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy 
result to   
                VT and generate right bidiagonalizing vectors 
in VT   
                (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB)
 */


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("U", m, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("U", m, n, &A(1,1), lda, &VT(1,1), ldvt);
#endif

		if (wntva) {
		    nrvt = *n;
		}
		if (wntvs) {
		    nrvt = *m;
		}
		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("P", &nrvt, n, m, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("P", &nrvt, n, m, &VT(1,1), ldvt, &WORK(itaup), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("P", &nrvt, n, m, &VT(1,1), ldvt, &WORK(itaup), 
#endif

			&WORK(iwork), &i__2, &ierr);
	    }
	    if (wntuo) {

/*              If left singular vectors desired in A, generat
e left   
                bidiagonalizing vectors in A   
                (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) 
*/

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("Q", m, m, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("Q", m, m, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("Q", m, m, n, &A(1,1), lda, &WORK(itauq), &WORK(
#endif

			iwork), &i__2, &ierr);
	    }
	    if (wntvo) {

/*              If right singular vectors desired in A, genera
te right   
                bidiagonalizing vectors in A   
                (Workspace: need 4*M, prefer 3*M+M*NB) */

		i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qorgbr("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif

			iwork), &i__2, &ierr);
	    }
	    iwork = ie + *m;
	    if (wntuas || wntuo) {
		nru = *m;
	    }
	    if (wntun) {
		nru = 0;
	    }
	    if (wntvas || wntvo) {
		ncvt = *n;
	    }
	    if (wntvn) {
		ncvt = 0;
	    }
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in U and computing right
 singular   
                vectors in VT   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &U(1,1), ldu, dum, &c__1, &
#endif

			WORK(iwork), info);
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in U and computing right
 singular   
                vectors in A   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &A(1,1), lda, &U(1,1), ldu, dum, &c__1, &WORK(
#endif

			iwork), info);
	    } else {

/*              Perform bidiagonal QR iteration, if desired, c
omputing   
                left singular vectors in A and computing right
 singular   
                vectors in VT   
                (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qbdsqr("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qbdsqr_("L", m, &ncvt, &nru, &c__0, &S(1), &WORK(ie), &VT(1,1), ldvt, &A(1,1), lda, dum, &c__1, &
#endif

			WORK(iwork), info);
	    }

	}

    }

/*     If DBDSQR failed to converge, copy unconverged superdiagonals   
       to WORK( 2:MINMN ) */

    if (*info != 0) {
	if (ie > 2) {
	    i__2 = minmn - 1;
	    for (i = 1; i <= minmn-1; ++i) {
		WORK(i + 1) = WORK(i + ie - 1);
/* L50: */
	    }
	}
	if (ie < 2) {
	    for (i = minmn - 1; i >= 1; --i) {
		WORK(i + 1) = WORK(i + ie - 1);
/* L60: */
	    }
	}
    }

/*     Undo scaling if necessary */

    if (iscl == 1) {
	if (anrm > bignum) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif

		    minmn, &ierr);
	}
	if (*info != 0 && anrm > bignum) {
	    i__2 = minmn - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &WORK(2),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &WORK(2),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &WORK(2),
#endif

		     &minmn, &ierr);
	}
	if (anrm < smlnum) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif

		    minmn, &ierr);
	}
	if (*info != 0 && anrm < smlnum) {
	    i__2 = minmn - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &WORK(2),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &WORK(2),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &WORK(2),
#endif

		     &minmn, &ierr);
	}
    }

/*     Return optimal workspace in WORK(1) */

    WORK(1) = (LONG DOUBLE) maxwrk;

    return;

/*     End of DGESVD */

} /* dgesvd_ */

