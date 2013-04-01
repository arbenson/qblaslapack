#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrexc_(char *compq, int *n, LONG DOUBLE *t, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrexc(char *compq, int *n, LONG DOUBLE *t, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrexc_(char *compq, int *n, LONG DOUBLE *t, int *
#endif

	ldt, LONG DOUBLE *q, int *ldq, int *ifst, int *ilst, 
	LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTREXC reorders the real Schur factorization of a real matrix   
    A = Q*T*Q**T, so that the diagonal block of T with row index IFST is 
  
    moved to row ILST.   

    The real Schur form T is reordered by an orthogonal similarity   
    transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors 
  
    is updated by postmultiplying it with Z.   

    T must be in Schur canonical form (as returned by DHSEQR), that is,   
    block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each   
    2-by-2 diagonal block has its diagonal elements equal and its   
    off-diagonal elements of opposite sign.   

    Arguments   
    =========   

    COMPQ   (input) CHARACTER*1   
            = 'V':  update the matrix Q of Schur vectors;   
            = 'N':  do not update Q.   

    N       (input) INTEGER   
            The order of the matrix T. N >= 0.   

    T       (input/output) LONG DOUBLE PRECISION array, dimension (LDT,N)   
            On entry, the upper quasi-triangular matrix T, in Schur   
            Schur canonical form.   
            On exit, the reordered upper quasi-triangular matrix, again   
            in Schur canonical form.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= MAX(1,N).   

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            On entry, if COMPQ = 'V', the matrix Q of Schur vectors.   
            On exit, if COMPQ = 'V', Q has been postmultiplied by the   
            orthogonal transformation matrix Z which reorders T.   
            If COMPQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= MAX(1,N).   

    IFST    (input/output) INTEGER   
    ILST    (input/output) INTEGER   
            Specify the reordering of the diagonal blocks of T.   
            The block with row index IFST is moved to row ILST, by a   
            sequence of transpositions between adjacent blocks.   
            On exit, if IFST pointed on entry to the second row of a   
            2-by-2 block, it is changed to point to the first row; ILST   
            always points to the first row of the block in its final   
            position (which may differ from its input value by +1 or -1). 
  
            1 <= IFST <= N; 1 <= ILST <= N.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            = 1:  two adjacent blocks were too close to swap (the problem 
  
                  is very ill-conditioned); T may have been partially   
                  reordered, and ILST points to the first row of the   
                  current position of the block being moved.   

    ===================================================================== 
  


       Decode and test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c__2 = 2;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int here;
    extern long int lsame_(char *, char *);
    static long int wantq;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaexc_(long int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaexc(long int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaexc_(long int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *, int *, int *, int 
	    *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    static int nbnext, nbf, nbl;



#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;
    wantq = lsame_(compq, "V");
    if (! wantq && ! lsame_(compq, "N")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldt < MAX(1,*n)) {
	*info = -4;
    } else if (*ldq < 1 || (wantq && *ldq < MAX(1,*n))) {
	*info = -6;
    } else if (*ifst < 1 || *ifst > *n) {
	*info = -7;
    } else if (*ilst < 1 || *ilst > *n) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTREXC", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return;
    }

/*     Determine the first row of specified block   
       and find out it is 1 by 1 or 2 by 2. */

    if (*ifst > 1) {
	if (T(*ifst,*ifst-1) != 0.) {
	    --(*ifst);
	}
    }
    nbf = 1;
    if (*ifst < *n) {
	if (T(*ifst+1,*ifst) != 0.) {
	    nbf = 2;
	}
    }

/*     Determine the first row of the final block   
       and find out it is 1 by 1 or 2 by 2. */

    if (*ilst > 1) {
	if (T(*ilst,*ilst-1) != 0.) {
	    --(*ilst);
	}
    }
    nbl = 1;
    if (*ilst < *n) {
	if (T(*ilst+1,*ilst) != 0.) {
	    nbl = 2;
	}
    }

    if (*ifst == *ilst) {
	return;
    }

    if (*ifst < *ilst) {

/*        Update ILST */

	if (nbf == 2 && nbl == 1) {
	    --(*ilst);
	}
	if (nbf == 1 && nbl == 2) {
	    ++(*ilst);
	}

	here = *ifst;

L10:

/*        Swap block with next one below */

	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

	    nbnext = 1;
	    if (here + nbf + 1 <= *n) {
		if (T(here+nbf+1,here+nbf) != 0.) {
		    nbnext = 2;
		}
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &here, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &here, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &here, &
#endif

		    nbf, &nbnext, &WORK(1), info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    here += nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

	    if (nbf == 2) {
		if (T(here+1,here) == 0.) {
		    nbf = 3;
		}
	    }

	} else {

/*           Current block consists of two 1 by 1 blocks each of w
hich   
             must be swapped individually */

	    nbnext = 1;
	    if (here + 3 <= *n) {
		if (T(here+3,here+2) != 0.) {
		    nbnext = 2;
		}
	    }
	    i__1 = here + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif

		    c__1, &nbnext, &WORK(1), info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible 
*/


#ifdef PETSC_PREFIX_SUFFIX
		dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			here, &c__1, &nbnext, &WORK(1), info);
		++here;
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

		if (T(here+2,here+1) == 0.) {
		    nbnext = 1;
		}
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */


#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    here, &c__1, &nbnext, &WORK(1), info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    here += 2;
		} else {

/*                 2 by 2 Block did split */


#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    here, &c__1, &c__1, &WORK(1), info);
		    i__1 = here + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    i__1, &c__1, &c__1, &WORK(1), info);
		    here += 2;
		}
	    }
	}
	if (here < *ilst) {
	    goto L10;
	}

    } else {

	here = *ifst;
L20:

/*        Swap block with next one above */

	if (nbf == 1 || nbf == 2) {

/*           Current block either 1 by 1 or 2 by 2 */

	    nbnext = 1;
	    if (here >= 3) {
		if (T(here-1,here-2) != 0.) {
		    nbnext = 2;
		}
	    }
	    i__1 = here - nbnext;

#ifdef PETSC_PREFIX_SUFFIX
	    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif

		    nbnext, &nbf, &WORK(1), info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    here -= nbnext;

/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

	    if (nbf == 2) {
		if (T(here+1,here) == 0.) {
		    nbf = 3;
		}
	    }

	} else {

/*           Current block consists of two 1 by 1 blocks each of w
hich   
             must be swapped individually */

	    nbnext = 1;
	    if (here >= 3) {
		if (T(here-1,here-2) != 0.) {
		    nbnext = 2;
		}
	    }
	    i__1 = here - nbnext;

#ifdef PETSC_PREFIX_SUFFIX
	    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &i__1, &
#endif

		    nbnext, &c__1, &WORK(1), info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    if (nbnext == 1) {

/*              Swap two 1 by 1 blocks, no problems possible 
*/


#ifdef PETSC_PREFIX_SUFFIX
		dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			here, &nbnext, &c__1, &WORK(1), info);
		--here;
	    } else {

/*              Recompute NBNEXT in case 2 by 2 split */

		if (T(here,here-1) == 0.) {
		    nbnext = 1;
		}
		if (nbnext == 2) {

/*                 2 by 2 Block did not split */

		    i__1 = here - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    i__1, &c__2, &c__1, &WORK(1), info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    here += -2;
		} else {

/*                 2 by 2 Block did split */


#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    here, &c__1, &c__1, &WORK(1), info);
		    i__1 = here - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaexc(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaexc_(&wantq, n, &T(1,1), ldt, &Q(1,1), ldq, &
#endif

			    i__1, &c__1, &c__1, &WORK(1), info);
		    here += -2;
		}
	    }
	}
	if (here > *ilst) {
	    goto L20;
	}
    }
    *ilst = here;

    return;

/*     End of DTREXC */

} /* dtrexc_ */

