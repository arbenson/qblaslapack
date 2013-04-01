#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaeda_(int *n, int *tlvls, int *curlvl, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaeda(int *n, int *tlvls, int *curlvl, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaeda_(int *n, int *tlvls, int *curlvl, 
#endif

	int *curpbm, int *prmptr, int *perm, int *givptr, 
	int *givcol, LONG DOUBLE *givnum, LONG DOUBLE *q, int *qptr, 
	LONG DOUBLE *z, LONG DOUBLE *ztemp, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAEDA computes the Z vector corresponding to the merge step in the   
    CURLVLth step of the merge process with TLVLS steps for the CURPBMth 
  
    problem.   

    Arguments   
    =========   

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    TLVLS  (input) INTEGER   
           The total number of merging levels in the overall divide and   
           conquer tree.   

    CURLVL (input) INTEGER   
           The current level in the overall merge routine,   
           0 <= curlvl <= tlvls.   

    CURPBM (input) INTEGER   
           The current problem in the current level in the overall   
           merge routine (counting from upper left to lower right).   

    PRMPTR (input) INTEGER array, dimension (N lg N)   
           Contains a list of pointers which indicate where in PERM a   
           level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)   
           indicates the size of the permutation and incidentally the   
           size of the full, non-deflated problem.   

    PERM   (input) INTEGER array, dimension (N lg N)   
           Contains the permutations (from deflation and sorting) to be   
           applied to each eigenblock.   

    GIVPTR (input) INTEGER array, dimension (N lg N)   
           Contains a list of pointers which indicate where in GIVCOL a   
           level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i) 
  
           indicates the number of Givens rotations.   

    GIVCOL (input) INTEGER array, dimension (2, N lg N)   
           Each pair of numbers indicates a pair of columns to take place 
  
           in a Givens rotation.   

    GIVNUM (input) LONG DOUBLE PRECISION array, dimension (2, N lg N)   
           Each number indicates the S value to be used in the   
           corresponding Givens rotation.   

    Q      (input) LONG DOUBLE PRECISION array, dimension (N**2)   
           Contains the square eigenblocks from previous levels, the   
           starting positions for blocks are given by QPTR.   

    QPTR   (input) INTEGER array, dimension (N+2)   
           Contains a list of pointers which indicate where in Q an   
           eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates   
           the size of the block.   

    Z      (output) LONG DOUBLE PRECISION array, dimension (N)   
           On output this vector contains the updating vector (the last   
           row of the first sub-eigenvector matrix and the first row of   
           the second sub-eigenvector matrix).   

    ZTEMP  (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static int c__2 = 2;
    static int c__1 = 1;
    static LONG DOUBLE c_b24 = 1.;
    static LONG DOUBLE c_b26 = 0.;
    
    /* System generated locals */
    int i__1, i__2, i__3;
    /* Builtin functions */
    /* Local variables */

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
    static int curr, bsiz1, bsiz2, psiz1, psiz2, i, k, zptr1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dcopy_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy_(int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), xerbla_(char *,
	     int *);
    static int mid, ptr;


    --ztemp;
    --z;
    --qptr;
    --q;
    givnum -= 3;
    givcol -= 3;
    --givptr;
    --perm;
    --prmptr;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -1;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAEDA", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Determine location of first number in second half. */

    mid = *n / 2 + 1;

/*     Gather last/first rows of appropriate eigenblocks into center of Z 
*/

    ptr = 1;

/*     Determine location of lowest level subproblem in the full storage 
  
       scheme */

    i__1 = *curlvl - 1;
    curr = ptr + *curpbm * (int)pow((LONG DOUBLE)c__2,(LONG DOUBLE) *curlvl)+(int)pow((LONG DOUBLE)c__2, (LONG DOUBLE)i__1)-1;

/*     Determine size of these matrices.  We add HALF to the value of   
       the SQRT in case the machine underestimates one of these square   
       roots. */

    bsiz1 = (int) (sqrt((LONG DOUBLE) (qptr[curr + 1] - qptr[curr])) + .5);
    bsiz2 = (int) (sqrt((LONG DOUBLE) (qptr[curr + 2] - qptr[curr + 1])) + 
	    .5);
    i__1 = mid - bsiz1 - 1;
    for (k = 1; k <= i__1; ++k) {
	z[k] = 0.;
/* L10: */
    }

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z[mid - bsiz1], &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z[mid - bsiz1], &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z[mid - bsiz1], &c__1)
#endif

	    ;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z[mid], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z[mid], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z[mid], &c__1);
#endif

    i__1 = *n;
    for (k = mid + bsiz2; k <= i__1; ++k) {
	z[k] = 0.;
/* L20: */
    }

/*     Loop thru remaining levels 1 -> CURLVL applying the Givens   
       rotations and permutation and then multiplying the center matrices 
  
       against the current Z. */

    ptr = (int)pow((LONG DOUBLE)c__2, (LONG DOUBLE)*tlvls) + 1;
    i__1 = *curlvl - 1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *curlvl - k;
	i__3 = *curlvl - k - 1;
	curr = ptr + *curpbm * (int)pow((LONG DOUBLE)c__2, (LONG DOUBLE)i__2) + (int)pow((LONG DOUBLE)c__2, (LONG DOUBLE)i__3) - 
		1;
	psiz1 = prmptr[curr + 1] - prmptr[curr];
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
	zptr1 = mid - psiz1;

/*       Apply Givens at CURR and CURR+1 */

	i__2 = givptr[curr + 1] - 1;
	for (i = givptr[curr]; i <= i__2; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&c__1, &z[zptr1 + givcol[(i << 1) + 1] - 1], &c__1, &z[
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&c__1, &z[zptr1 + givcol[(i << 1) + 1] - 1], &c__1, &z[
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&c__1, &z[zptr1 + givcol[(i << 1) + 1] - 1], &c__1, &z[
#endif

		    zptr1 + givcol[(i << 1) + 2] - 1], &c__1, &givnum[(i << 1)
		     + 1], &givnum[(i << 1) + 2]);
/* L30: */
	}
	i__2 = givptr[curr + 2] - 1;
	for (i = givptr[curr + 1]; i <= i__2; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&c__1, &z[mid - 1 + givcol[(i << 1) + 1]], &c__1, &z[mid - 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&c__1, &z[mid - 1 + givcol[(i << 1) + 1]], &c__1, &z[mid - 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&c__1, &z[mid - 1 + givcol[(i << 1) + 1]], &c__1, &z[mid - 
#endif

		    1 + givcol[(i << 1) + 2]], &c__1, &givnum[(i << 1) + 1], &
		    givnum[(i << 1) + 2]);
/* L40: */
	}
	psiz1 = prmptr[curr + 1] - prmptr[curr];
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
	i__2 = psiz1 - 1;
	for (i = 0; i <= i__2; ++i) {
	    ztemp[i + 1] = z[zptr1 + perm[prmptr[curr] + i] - 1];
/* L50: */
	}
	i__2 = psiz2 - 1;
	for (i = 0; i <= i__2; ++i) {
	    ztemp[psiz1 + i + 1] = z[mid + perm[prmptr[curr + 1] + i] - 1];
/* L60: */
	}

/*        Multiply Blocks at CURR and CURR+1   

          Determine size of these matrices.  We add HALF to the value 
of   
          the SQRT in case the machine underestimates one of these   
          square roots. */

	bsiz1 = (int) (sqrt((LONG DOUBLE) (qptr[curr + 1] - qptr[curr])) + 
		.5);
	bsiz2 = (int) (sqrt((LONG DOUBLE) (qptr[curr + 2] - qptr[curr + 1])
		) + .5);
	if (bsiz1 > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, &
#endif

		    ztemp[1], &c__1, &c_b26, &z[zptr1], &c__1);
	}
	i__2 = psiz1 - bsiz1;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(&i__2, &ztemp[bsiz1 + 1], &c__1, &z[zptr1 + bsiz1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(&i__2, &ztemp[bsiz1 + 1], &c__1, &z[zptr1 + bsiz1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(&i__2, &ztemp[bsiz1 + 1], &c__1, &z[zptr1 + bsiz1], &c__1);
#endif

	if (bsiz2 > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, &
#endif

		    ztemp[psiz1 + 1], &c__1, &c_b26, &z[mid], &c__1);
	}
	i__2 = psiz2 - bsiz2;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z[mid + bsiz2], &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z[mid + bsiz2], &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z[mid + bsiz2], &
#endif

		c__1);

	i__2 = *tlvls - k;
	ptr += (int)pow((LONG DOUBLE)c__2, (LONG DOUBLE)i__2);
/* L70: */
    }

    return;

/*     End of DLAEDA */

} /* dlaeda_ */

