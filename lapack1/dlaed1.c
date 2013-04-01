#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed1_(int *n, LONG DOUBLE *d, LONG DOUBLE *q, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed1(int *n, LONG DOUBLE *d, LONG DOUBLE *q, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed1_(int *n, LONG DOUBLE *d, LONG DOUBLE *q, 
#endif

	int *ldq, int *indxq, LONG DOUBLE *rho, int *cutpnt, 
	LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED1 computes the updated eigensystem of a diagonal   
    matrix after modification by a rank-one symmetric matrix.  This   
    routine is used only for the eigenproblem which requires all   
    eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles 
  
    the case in which eigenvalues only or eigenvalues and eigenvectors   
    of a full symmetric matrix (which was reduced to tridiagonal form)   
    are desired.   

      T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out) 
  

       where Z = Q'u, u is a vector of length N with ones in the   
       CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.   

       The eigenvectors of the original matrix are stored in Q, and the   
       eigenvalues are in D.  The algorithm consists of three stages:   

          The first stage consists of deflating the size of the problem   
          when there are multiple eigenvalues or if there is a zero in   
          the Z vector.  For each such occurence the dimension of the   
          secular equation problem is reduced by one.  This stage is   
          performed by the routine DLAED2.   

          The second stage consists of calculating the updated   
          eigenvalues. This is done by finding the roots of the secular   
          equation via the routine DLAED4 (as called by SLAED3).   
          This routine also calculates the eigenvectors of the current   
          problem.   

          The final stage consists of computing the updated eigenvectors 
  
          directly using the updated eigenvalues.  The eigenvectors for   
          the current problem are multiplied with the eigenvectors from   
          the overall problem.   

    Arguments   
    =========   

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
           On entry, the eigenvalues of the rank-1-perturbed matrix.   
           On exit, the eigenvalues of the repaired matrix.   

    Q      (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
           On entry, the eigenvectors of the rank-1-perturbed matrix.   
           On exit, the eigenvectors of the repaired tridiagonal matrix. 
  

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  LDQ >= MAX(1,N).   

    INDXQ  (input/output) INTEGER array, dimension (N)   
           On entry, the permutation which separately sorts the two   
           subproblems in D into ascending order.   
           On exit, the permutation which will reintegrate the   
           subproblems back into sorted order,   
           i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.   

    RHO    (input) LONG DOUBLE PRECISION   
           The subdiagonal entry used to create the rank-1 modification. 
  

    CUTPNT (input) INTEGER   
           The location of the last eigenvalue in the leading sub-matrix. 
  
           MIN(1,N) <= CUTPNT <= N.   

    WORK   (workspace) LONG DOUBLE PRECISION array, dimension (3*N+2*N**2)   

    IWORK  (workspace) INTEGER array, dimension (4*N)   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, an eigenvalue did not converge   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int indx, i, k, indxc;

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
    static int indxp;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaed2_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed2(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed2_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *, 
	    LONG DOUBLE *, int *, int *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlaed3_(int *, int *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaed3(int *, int *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaed3_(int *, int *, int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, int *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, int *);
    static int n1, n2, idlmda, is, iw, iz;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlamrg_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamrg(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamrg_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, int *), xerbla_(char *, int *);
    static int coltyp, iq2, ldq2, zpp1;



#define D(I) d[(I)-1]
#define INDXQ(I) indxq[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*ldq < MAX(1,*n)) {
	*info = -4;
    } else if (MIN(1,*n) > *cutpnt || *n < *cutpnt) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED1", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     The following values are for bookkeeping purposes only.  They are 
  
       int pointers which indicate the portion of the workspace   
       used by a particular array in DLAED2 and SLAED3. */

    ldq2 = *n;

    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    is = iq2 + *n * ldq2;

    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;


/*     Form the z-vector which consists of the last row of Q_1 and the   
       first row of Q_2. */


#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(cutpnt, &Q(*cutpnt,1), ldq, &WORK(iz), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(cutpnt, &Q(*cutpnt,1), ldq, &WORK(iz), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(cutpnt, &Q(*cutpnt,1), ldq, &WORK(iz), &c__1);
#endif

    zpp1 = *cutpnt + 1;
    i__1 = *n - *cutpnt;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&i__1, &Q(zpp1,zpp1), ldq, &WORK(iz + *cutpnt), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&i__1, &Q(zpp1,zpp1), ldq, &WORK(iz + *cutpnt), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&i__1, &Q(zpp1,zpp1), ldq, &WORK(iz + *cutpnt), &c__1);
#endif


/*     Deflate eigenvalues. */


#ifdef PETSC_PREFIX_SUFFIX
    dlaed2_(&k, n, &D(1), &Q(1,1), ldq, &INDXQ(1), rho, cutpnt, &WORK(iz)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlaed2(&k, n, &D(1), &Q(1,1), ldq, &INDXQ(1), rho, cutpnt, &WORK(iz)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlaed2_(&k, n, &D(1), &Q(1,1), ldq, &INDXQ(1), rho, cutpnt, &WORK(iz)
#endif

	    , &WORK(idlmda), &WORK(iq2), &ldq2, &IWORK(indxc), &WORK(iw), &
	    IWORK(indxp), &IWORK(indx), &IWORK(coltyp), info);
    if (*info != 0) {
	goto L20;
    }

/*     Solve Secular Equation. */

    if (k != 0) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaed3_(&k, &c__1, &k, n, &D(1), &Q(1,1), ldq, rho, cutpnt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaed3(&k, &c__1, &k, n, &D(1), &Q(1,1), ldq, rho, cutpnt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaed3_(&k, &c__1, &k, n, &D(1), &Q(1,1), ldq, rho, cutpnt, &
#endif

		WORK(idlmda), &WORK(iq2), &ldq2, &IWORK(indxc), &IWORK(coltyp)
		, &WORK(iw), &WORK(is), &k, info);
	if (*info != 0) {
	    goto L20;
	}

/*     Prepare the INDXQ sorting permutation. */

	n1 = k;
	n2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
	dlamrg_(&n1, &n2, &D(1), &c__1, &c_n1, &INDXQ(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlamrg(&n1, &n2, &D(1), &c__1, &c_n1, &INDXQ(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlamrg_(&n1, &n2, &D(1), &c__1, &c_n1, &INDXQ(1));
#endif

    } else {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    INDXQ(i) = i;
/* L10: */
	}
    }

L20:
    return;

/*     End of DLAED1 */

} /* dlaed1_ */

