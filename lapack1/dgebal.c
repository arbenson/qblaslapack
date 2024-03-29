#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgebal_(char *job, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgebal(char *job, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgebal_(char *job, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *ilo, int *ihi, LONG DOUBLE *scale, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEBAL balances a general real matrix A.  This involves, first,   
    permuting A by a similarity transformation to isolate eigenvalues   
    in the first 1 to ILO-1 and last IHI+1 to N elements on the   
    diagonal; and second, applying a diagonal similarity transformation   
    to rows and columns ILO to IHI to make the rows and columns as   
    close in norm as possible.  Both steps are optional.   

    Balancing may reduce the 1-norm of the matrix, and improve the   
    accuracy of the computed eigenvalues and/or eigenvectors.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the operations to be performed on A:   
            = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0   
                    for i = 1,...,N;   
            = 'P':  permute only;   
            = 'S':  scale only;   
            = 'B':  both permute and scale.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the input matrix A.   
            On exit,  A is overwritten by the balanced matrix.   
            If JOB = 'N', A is not referenced.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    ILO     (output) INTEGER   
    IHI     (output) INTEGER   
            ILO and IHI are set to ints such that on exit   
            A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.   
            If JOB = 'N' or 'S', ILO = 1 and IHI = N.   

    SCALE   (output) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied to   
            A.  If P(j) is the index of the row and column interchanged   
            with row and column j and D(j) is the scaling factor   
            applied to row and column j, then   
            SCALE(j) = P(j)    for j = 1,...,ILO-1   
                     = D(j)    for j = ILO,...,IHI   
                     = P(j)    for j = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The permutations consist of row and column interchanges which put   
    the matrix in the form   

               ( T1   X   Y  )   
       P A P = (  0   B   Z  )   
               (  0   0   T2 )   

    where T1 and T2 are upper triangular matrices whose eigenvalues lie   
    along the diagonal.  The column indices ILO and IHI mark the starting 
  
    and ending columns of the submatrix B. Balancing consists of applying 
  
    a diagonal similarity transformation inv(D) * B * D to make the   
    1-norms of each row of B and its corresponding column nearly equal.   
    The output matrix is   

       ( T1     X*D          Y    )   
       (  0  inv(D)*B*D  inv(D)*Z ).   
       (  0      0           T2   )   

    Information about the permutations P and the diagonal matrix D is   
    returned in the vector SCALE.   

    This subroutine is based on the EISPACK routine BALANC.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static int iexc;
    static LONG DOUBLE c, f, g;
    static int i, j, k, l, m;
    static LONG DOUBLE r, s;

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
    extern /* Subroutine */ void dswap_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static LONG DOUBLE sfmin1, sfmin2, sfmax1, sfmax2, ca, ra;

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
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int noconv;
    static int ica, ira;



#define SCALE(I) scale[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(
	    job, "S") && ! lsame_(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEBAL", &i__1);
	return;
    }

    k = 1;
    l = *n;

    if (*n == 0) {
	goto L210;
    }

    if (lsame_(job, "N")) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    SCALE(i) = 1.;
/* L10: */
	}
	goto L210;
    }

    if (lsame_(job, "S")) {
	goto L120;
    }

/*     Permutation to isolate eigenvalues if possible */

    goto L50;

/*     Row and column exchange. */

L20:
    SCALE(m) = (LONG DOUBLE) j;
    if (j == m) {
	goto L30;
    }


#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif

    i__1 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&i__1, &A(j,k), lda, &A(m,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&i__1, &A(j,k), lda, &A(m,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&i__1, &A(j,k), lda, &A(m,k), lda);
#endif


L30:
    switch (iexc) {
	case 1:  goto L40;
	case 2:  goto L80;
    }

/*     Search for rows isolating an eigenvalue and push them down. */

L40:
    if (l == 1) {
	goto L210;
    }
    --l;

L50:
    for (j = l; j >= 1; --j) {

	i__1 = l;
	for (i = 1; i <= l; ++i) {
	    if (i == j) {
		goto L60;
	    }
	    if (A(j,i) != 0.) {
		goto L70;
	    }
L60:
	    ;
	}

	m = l;
	iexc = 1;
	goto L20;
L70:
	;
    }

    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

L80:
    ++k;

L90:
    i__1 = l;
    for (j = k; j <= l; ++j) {

	i__2 = l;
	for (i = k; i <= l; ++i) {
	    if (i == j) {
		goto L100;
	    }
	    if (A(i,j) != 0.) {
		goto L110;
	    }
L100:
	    ;
	}

	m = k;
	iexc = 2;
	goto L20;
L110:
	;
    }

L120:
    i__1 = l;
    for (i = k; i <= l; ++i) {
	SCALE(i) = 1.;
/* L130: */
    }

    if (lsame_(job, "P")) {
	goto L210;
    }

/*     Balance the submatrix in rows K to L.   

       Iterative loop for norm reduction */


#ifdef PETSC_PREFIX_SUFFIX
    sfmin1 = dlamch_("S") / dlamch_("P");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sfmin1 = qlamch("S") / dlamch_("P");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sfmin1 = qlamch_("S") / dlamch_("P");
#endif

    sfmax1 = 1. / sfmin1;
    sfmin2 = sfmin1 * 10.;
    sfmax2 = 1. / sfmin2;
L140:
    noconv = 0;

    i__1 = l;
    for (i = k; i <= l; ++i) {
	c = 0.;
	r = 0.;

	i__2 = l;
	for (j = k; j <= l; ++j) {
	    if (j == i) {
		goto L150;
	    }
	    c += (d__1 = A(j,i), ABS(d__1));
	    r += (d__1 = A(i,j), ABS(d__1));
L150:
	    ;
	}

#ifdef PETSC_PREFIX_SUFFIX
	ica = idamax_(&l, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	ica = iqamax(&l, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	ica = iqamax_(&l, &A(1,i), &c__1);
#endif

	ca = (d__1 = A(ica,i), ABS(d__1));
	i__2 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
	ira = idamax_(&i__2, &A(i,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	ira = iqamax(&i__2, &A(i,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	ira = iqamax_(&i__2, &A(i,k), lda);
#endif

	ra = (d__1 = A(i,ira+k-1), ABS(d__1));

/*        Guard against zero C or R due to underflow. */

	if (c == 0. || r == 0.) {
	    goto L200;
	}
	g = r / 10.;
	f = 1.;
	s = c + r;
L160:
/* Computing MAX */
	d__1 = MAX(f,c);
/* Computing MIN */
	d__2 = MIN(r,g);
	if (c >= g || MAX(d__1,ca) >= sfmax2 || MIN(d__2,ra) <= sfmin2) {
	    goto L170;
	}
	f *= 10.;
	c *= 10.;
	ca *= 10.;
	r /= 10.;
	g /= 10.;
	ra /= 10.;
	goto L160;

L170:
	g = c / 10.;
L180:
/* Computing MIN */
	d__1 = MIN(f,c), d__1 = MIN(d__1,g);
	if (g < r || MAX(r,ra) >= sfmax2 || MIN(d__1,ca) <= sfmin2) {
	    goto L190;
	}
	f /= 10.;
	c /= 10.;
	g /= 10.;
	ca /= 10.;
	r *= 10.;
	ra *= 10.;
	goto L180;

/*        Now balance. */

L190:
	if (c + r >= s * .95) {
	    goto L200;
	}
	if (f < 1. && SCALE(i) < 1.) {
	    if (f * SCALE(i) <= sfmin1) {
		goto L200;
	    }
	}
	if (f > 1. && SCALE(i) > 1.) {
	    if (SCALE(i) >= sfmax1 / f) {
		goto L200;
	    }
	}
	g = 1. / f;
	SCALE(i) *= f;
	noconv = 1;

	i__2 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__2, &g, &A(i,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__2, &g, &A(i,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__2, &g, &A(i,k), lda);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&l, &f, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&l, &f, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&l, &f, &A(1,i), &c__1);
#endif


L200:
	;
    }

    if (noconv) {
	goto L140;
    }

L210:
    *ilo = k;
    *ihi = l;

    return;

/*     End of DGEBAL */

} /* dgebal_ */

