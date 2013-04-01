#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsbtrd_(char *vect, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsbtrd(char *vect, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsbtrd_(char *vect, char *uplo, int *n, int *kd, 
#endif

	LONG DOUBLE *ab, int *ldab, LONG DOUBLE *d, LONG DOUBLE *e, 
	LONG DOUBLE *q, int *ldq, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSBTRD reduces a real symmetric band matrix A to symmetric   
    tridiagonal form T by an orthogonal similarity transformation:   
    Q**T * A * Q = T.   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            = 'N':  do not form Q;   
            = 'V':  form Q;   
            = 'U':  update a matrix X, by forming X*Q.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  
            On exit, the diagonal elements of AB are overwritten by the   
            diagonal elements of the tridiagonal matrix T; if KD > 0, the 
  
            elements on the first superdiagonal (if UPLO = 'U') or the   
            first subdiagonal (if UPLO = 'L') are overwritten by the   
            off-diagonal elements of T; the rest of AB is overwritten by 
  
            values generated during the reduction.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    D       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T.   

    E       (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. 
  

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            On entry, if VECT = 'U', then Q must contain an N-by-N   
            matrix X; if VECT = 'N' or 'V', then Q need not be set.   

            On exit:   
            if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;   
            if VECT = 'U', Q contains the product X*Q;   
            if VECT = 'N', the array Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b9 = 0.;
    static LONG DOUBLE c_b10 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int inca;
    static LONG DOUBLE temp;

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
    static int i, j, k, l;
    extern long int lsame_(char *, char *);
    static long int initq, wantq, upper;
    static int P_j1, j2;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlar2v_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlar2v(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlar2v_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int nr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *), xerbla_(char *, int *), dlargv_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *), xerbla_(char *, int *), qlargv(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *), xerbla_(char *, int *), qlargv_(
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dlartv_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlartv(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlartv_(int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *);
    static int kd1, kdn, nrt;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    initq = lsame_(vect, "V");
    wantq = initq || lsame_(vect, "U");
    upper = lsame_(uplo, "U");
    kd1 = *kd + 1;
    *info = 0;
    if (! wantq && ! lsame_(vect, "N")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*kd < 0) {
	*info = -4;
    } else if (*ldab < kd1) {
	*info = -6;
    } else if (*ldq < MAX(1,*n) && wantq) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSBTRD", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Initialize Q to the unit matrix, if needed */

    if (initq) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b9, &c_b10, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b9, &c_b10, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b9, &c_b10, &Q(1,1), ldq);
#endif

    }

/*     Wherever possible, plane rotations are generated and applied in   
       vector operations of length NR over the index set P_J1:J2:KD1.   

       The cosines and sines of the plane rotations are stored in the   
       arrays D and WORK. */

    inca = kd1 * *ldab;
/* Computing MIN */
    i__1 = *n - 1;
    kdn = MIN(i__1,*kd);
    if (upper) {

	if (*kd > 1) {

/*           Reduce to tridiagonal form, working with upper triang
le */

	    nr = 0;
	    P_j1 = kdn + 2;
	    j2 = 1;

	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {

/*              Reduce i-th row of matrix to tridiagonal form 
*/

		for (k = kdn + 1; k >= 2; --k) {
		    P_j1 += kdn;
		    j2 += kdn;

		    if (nr > 0) {

/*                    generate plane rotations to anni
hilate nonzero   
                      elements which have been created
 outside the band */


#ifdef PETSC_PREFIX_SUFFIX
			dlargv_(&nr, &AB(1,P_j1-1), &inca, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlargv(&nr, &AB(1,P_j1-1), &inca, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlargv_(&nr, &AB(1,P_j1-1), &inca, &
#endif

				WORK(P_j1), &kd1, &D(P_j1), &kd1);

/*                    apply rotations from the right 
*/

			i__2 = *kd - 1;
			for (l = 1; l <= *kd-1; ++l) {

#ifdef PETSC_PREFIX_SUFFIX
			    dlartv_(&nr, &AB(l+1,P_j1-1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartv(&nr, &AB(l+1,P_j1-1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartv_(&nr, &AB(l+1,P_j1-1), &
#endif

				    inca, &AB(l,P_j1), &inca, &D(P_j1)
				    , &WORK(P_j1), &kd1);
/* L10: */
			}
		    }

		    if (k > 2) {
			if (k <= *n - i + 1) {

/*                       generate plane rotation t
o annihilate a(i,i+k-1)   
                         within the band */


#ifdef PETSC_PREFIX_SUFFIX
			    dlartg_(&AB(*kd-k+3,i+k-2), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartg(&AB(*kd-k+3,i+k-2), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartg_(&AB(*kd-k+3,i+k-2), 
#endif

				    &AB(*kd-k+2,i+k-1), 
				    &D(i + k - 1), &WORK(i + k - 1), &temp);
			    AB(*kd-k+3,i+k-2) = temp;

/*                       apply rotation from the r
ight */

			    i__2 = k - 3;

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(&i__2, &AB(*kd-k+4,i+k-2), &c__1, &AB(*kd-k+3,i+k-1), &c__1, &D(i + k - 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(&i__2, &AB(*kd-k+4,i+k-2), &c__1, &AB(*kd-k+3,i+k-1), &c__1, &D(i + k - 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(&i__2, &AB(*kd-k+4,i+k-2), &c__1, &AB(*kd-k+3,i+k-1), &c__1, &D(i + k - 1), &
#endif

				    WORK(i + k - 1));
			}
			++nr;
			P_j1 = P_j1 - kdn - 1;
		    }

/*                 apply plane rotations from both sides t
o diagonal   
                   blocks */

		    if (nr > 0) {

#ifdef PETSC_PREFIX_SUFFIX
			dlar2v_(&nr, &AB(kd1,P_j1-1), &AB(kd1,P_j1), &AB(*kd,P_j1), &inca,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlar2v(&nr, &AB(kd1,P_j1-1), &AB(kd1,P_j1), &AB(*kd,P_j1), &inca,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlar2v_(&nr, &AB(kd1,P_j1-1), &AB(kd1,P_j1), &AB(*kd,P_j1), &inca,
#endif

				 &D(P_j1), &WORK(P_j1), &kd1);
		    }

/*                 apply plane rotations from the left */

		    i__2 = *kd - 1;
		    for (l = 1; l <= *kd-1; ++l) {
			if (j2 + l > *n) {
			    nrt = nr - 1;
			} else {
			    nrt = nr;
			}
			if (nrt > 0) {

#ifdef PETSC_PREFIX_SUFFIX
			    dlartv_(&nrt, &AB(*kd-l,P_j1+l), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartv(&nrt, &AB(*kd-l,P_j1+l), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartv_(&nrt, &AB(*kd-l,P_j1+l), &
#endif

				    inca, &AB(*kd-l+1,P_j1+l), &inca, &D(P_j1), &WORK(P_j1), &kd1);
			}
/* L20: */
		    }

		    if (wantq) {

/*                    accumulate product of plane rota
tions in Q */

			i__2 = j2;
			i__3 = kd1;
			for (j = P_j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += 
				i__3) {

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif

/* L30: */
			}
		    }

		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bou
nds of the matrix */

			--nr;
			j2 = j2 - kdn - 1;
		    }

		    i__3 = j2;
		    i__2 = kd1;
		    for (j = P_j1; kd1 < 0 ? j >= j2 : j <= j2; j += kd1) 
			    {

/*                    create nonzero element a(j-1,j+k
d) outside the band   
                      and store it in WORK */

			WORK(j + *kd) = WORK(j) * AB(1,j+*kd);
			AB(1,j+*kd) = D(j) * AB(1,j+*kd);
/* L40: */
		    }
/* L50: */
		}
/* L60: */
	    }
	}

	if (*kd > 0) {

/*           copy off-diagonal elements to E */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = AB(*kd,i+1);
/* L70: */
	    }
	} else {

/*           set E to zero if original matrix was diagonal */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = 0.;
/* L80: */
	    }
	}

/*        copy diagonal elements to D */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    D(i) = AB(kd1,i);
/* L90: */
	}

    } else {

	if (*kd > 1) {

/*           Reduce to tridiagonal form, working with lower triang
le */

	    nr = 0;
	    P_j1 = kdn + 2;
	    j2 = 1;

	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {

/*              Reduce i-th column of matrix to tridiagonal fo
rm */

		for (k = kdn + 1; k >= 2; --k) {
		    P_j1 += kdn;
		    j2 += kdn;

		    if (nr > 0) {

/*                    generate plane rotations to anni
hilate nonzero   
                      elements which have been created
 outside the band */


#ifdef PETSC_PREFIX_SUFFIX
			dlargv_(&nr, &AB(kd1,P_j1-kd1), &inca, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlargv(&nr, &AB(kd1,P_j1-kd1), &inca, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlargv_(&nr, &AB(kd1,P_j1-kd1), &inca, &
#endif

				WORK(P_j1), &kd1, &D(P_j1), &kd1);

/*                    apply plane rotations from one s
ide */

			i__2 = *kd - 1;
			for (l = 1; l <= *kd-1; ++l) {

#ifdef PETSC_PREFIX_SUFFIX
			    dlartv_(&nr, &AB(kd1-l,P_j1-kd1+l), &inca, &AB(kd1-l+1,P_j1-kd1+l), &inca, &D(P_j1), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartv(&nr, &AB(kd1-l,P_j1-kd1+l), &inca, &AB(kd1-l+1,P_j1-kd1+l), &inca, &D(P_j1), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartv_(&nr, &AB(kd1-l,P_j1-kd1+l), &inca, &AB(kd1-l+1,P_j1-kd1+l), &inca, &D(P_j1), &WORK(
#endif

				    P_j1), &kd1);
/* L100: */
			}
		    }

		    if (k > 2) {
			if (k <= *n - i + 1) {

/*                       generate plane rotation t
o annihilate a(i+k-1,i)   
                         within the band */


#ifdef PETSC_PREFIX_SUFFIX
			    dlartg_(&AB(k-1,i), &AB(k,i), &D(i + k - 1), &WORK(i + k - 1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartg(&AB(k-1,i), &AB(k,i), &D(i + k - 1), &WORK(i + k - 1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartg_(&AB(k-1,i), &AB(k,i), &D(i + k - 1), &WORK(i + k - 1),
#endif

				     &temp);
			    AB(k-1,i) = temp;

/*                       apply rotation from the l
eft */

			    i__2 = k - 3;
			    i__3 = *ldab - 1;
			    i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(&i__2, &AB(k-2,i+1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(&i__2, &AB(k-2,i+1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(&i__2, &AB(k-2,i+1), &
#endif

				    i__3, &AB(k-1,i+1), &
				    i__4, &D(i + k - 1), &WORK(i + k - 1));
			}
			++nr;
			P_j1 = P_j1 - kdn - 1;
		    }

/*                 apply plane rotations from both sides t
o diagonal   
                   blocks */

		    if (nr > 0) {

#ifdef PETSC_PREFIX_SUFFIX
			dlar2v_(&nr, &AB(1,P_j1-1), &AB(1,P_j1), &AB(2,P_j1-1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlar2v(&nr, &AB(1,P_j1-1), &AB(1,P_j1), &AB(2,P_j1-1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlar2v_(&nr, &AB(1,P_j1-1), &AB(1,P_j1), &AB(2,P_j1-1), &
#endif

				inca, &D(P_j1), &WORK(P_j1), &kd1);
		    }

/*                 apply plane rotations from the right */

		    i__2 = *kd - 1;
		    for (l = 1; l <= *kd-1; ++l) {
			if (j2 + l > *n) {
			    nrt = nr - 1;
			} else {
			    nrt = nr;
			}
			if (nrt > 0) {

#ifdef PETSC_PREFIX_SUFFIX
			    dlartv_(&nrt, &AB(l+2,P_j1-1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qlartv(&nrt, &AB(l+2,P_j1-1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qlartv_(&nrt, &AB(l+2,P_j1-1), &
#endif

				    inca, &AB(l+1,P_j1), &inca, &
				    D(P_j1), &WORK(P_j1), &kd1);
			}
/* L110: */
		    }

		    if (wantq) {

/*                    accumulate product of plane rota
tions in Q */

			i__2 = j2;
			i__3 = kd1;
			for (j = P_j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += 
				i__3) {

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
#endif

/* L120: */
			}
		    }

		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bou
nds of the matrix */

			--nr;
			j2 = j2 - kdn - 1;
		    }

		    i__3 = j2;
		    i__2 = kd1;
		    for (j = P_j1; kd1 < 0 ? j >= j2 : j <= j2; j += kd1) 
			    {

/*                    create nonzero element a(j+kd,j-
1) outside the   
                      band and store it in WORK */

			WORK(j + *kd) = WORK(j) * AB(kd1,j);
			AB(kd1,j) = D(j) * AB(kd1,j);
/* L130: */
		    }
/* L140: */
		}
/* L150: */
	    }
	}

	if (*kd > 0) {

/*           copy off-diagonal elements to E */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = AB(2,i);
/* L160: */
	    }
	} else {

/*           set E to zero if original matrix was diagonal */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = 0.;
/* L170: */
	    }
	}

/*        copy diagonal elements to D */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    D(i) = AB(1,i);
/* L180: */
	}
    }

    return;

/*     End of DSBTRD */

} /* dsbtrd_ */

