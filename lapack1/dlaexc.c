#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaexc_(long int *wantq, int *n, LONG DOUBLE *t, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaexc(long int *wantq, int *n, LONG DOUBLE *t, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaexc_(long int *wantq, int *n, LONG DOUBLE *t, 
#endif

	int *ldt, LONG DOUBLE *q, int *ldq, int *P_j1, int *n1, 
	int *n2, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in 
  
    an upper quasi-triangular matrix T by an orthogonal similarity   
    transformation.   

    T must be in Schur canonical form, that is, block upper triangular   
    with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block   
    has its diagonal elemnts equal and its off-diagonal elements of   
    opposite sign.   

    Arguments   
    =========   

    WANTQ   (input) LOGICAL   
            = .TRUE. : accumulate the transformation in the matrix Q;   
            = .FALSE.: do not accumulate the transformation.   

    N       (input) INTEGER   
            The order of the matrix T. N >= 0.   

    T       (input/output) LONG DOUBLE PRECISION array, dimension (LDT,N)   
            On entry, the upper quasi-triangular matrix T, in Schur   
            canonical form.   
            On exit, the updated matrix T, again in Schur canonical form. 
  

    LDT     (input)  INTEGER   
            The leading dimension of the array T. LDT >= MAX(1,N).   

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            On entry, if WANTQ is .TRUE., the orthogonal matrix Q.   
            On exit, if WANTQ is .TRUE., the updated matrix Q.   
            If WANTQ is .FALSE., Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.   

    P_J1      (input) INTEGER   
            The index of the first row of the first block T11.   

    N1      (input) INTEGER   
            The order of the first block T11. N1 = 0, 1 or 2.   

    N2      (input) INTEGER   
            The order of the second block T22. N2 = 0, 1 or 2.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            = 1: the transformed matrix T would be too far from Schur   
                 form; the blocks are not swapped and T and Q are   
                 unchanged.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c__4 = 4;
    static long int c_false = 0;
    static int c_n1 = -1;
    static int c__2 = 2;
    static int c__3 = 3;
    
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int ierr;
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
    static LONG DOUBLE d[16]	/* was [4][4] */;
    static int k;
    static LONG DOUBLE u[3], scale, x[4]	/* was [2][2] */, dnorm;
    static int j2, j3, j4;
    static LONG DOUBLE xnorm, u1[3], u2[3];

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlanv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlanv2(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlanv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), dlasy2_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlasy2(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlasy2_(
#endif

	    long int *, long int *, int *, int *, int *, LONG DOUBLE 
	    *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static int nd;
    static LONG DOUBLE cs, t11, t22;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE t33;

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
    extern /* Subroutine */ void dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *);
    static LONG DOUBLE sn;

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


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *), dlarfx_(char *, int *, int *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfx(char *, int *, int *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfx_(char *, int *, int *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static LONG DOUBLE thresh, smlnum, wi1, wi2, wr1, wr2, eps, tau, tau1, 
	    tau2;



#define D(I) d[(I)]
#define WAS(I) was[(I)]
#define U(I) u[(I)]
#define U1(I) u1[(I)]
#define U2(I) u2[(I)]
#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;

/*     Quick return if possible */

    if (*n == 0 || *n1 == 0 || *n2 == 0) {
	return;
    }
    if (*P_j1 + *n1 > *n) {
	return;
    }

    j2 = *P_j1 + 1;
    j3 = *P_j1 + 2;
    j4 = *P_j1 + 3;

    if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

	t11 = T(*P_j1,*P_j1);
	t22 = T(j2,j2);

/*        Determine the transformation to perform the interchange. */

	d__1 = t22 - t11;

#ifdef PETSC_PREFIX_SUFFIX
	dlartg_(&T(*P_j1,j2), &d__1, &cs, &sn, &temp);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlartg(&T(*P_j1,j2), &d__1, &cs, &sn, &temp);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlartg_(&T(*P_j1,j2), &d__1, &cs, &sn, &temp);
#endif


/*        Apply transformation to the matrix T. */

	if (j3 <= *n) {
	    i__1 = *n - *P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__1, &T(*P_j1,j3), ldt, &T(j2,j3), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__1, &T(*P_j1,j3), ldt, &T(j2,j3), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__1, &T(*P_j1,j3), ldt, &T(j2,j3), 
#endif

		    ldt, &cs, &sn);
	}
	i__1 = *P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	drot_(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qrot(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qrot_(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &c__1, 
#endif

		&cs, &sn);

	T(*P_j1,*P_j1) = t22;
	T(j2,j2) = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */


#ifdef PETSC_PREFIX_SUFFIX
	    drot_(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &c__1, 
#endif

		    &cs, &sn);
	}

    } else {

/*        Swapping involves at least one 2-by-2 block.   

          Copy the diagonal block of order N1+N2 to the local array D 
  
          and compute its norm. */

	nd = *n1 + *n2;

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("Full", &nd, &nd, &T(*P_j1,*P_j1), ldt, d, &c__4);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("Full", &nd, &nd, &T(*P_j1,*P_j1), ldt, d, &c__4);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("Full", &nd, &nd, &T(*P_j1,*P_j1), ldt, d, &c__4);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dnorm = dlange_("Max", &nd, &nd, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	dnorm = qlange("Max", &nd, &nd, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	dnorm = qlange_("Max", &nd, &nd, d, &c__4, &WORK(1));
#endif


/*        Compute machine-dependent threshold for test for accepting 
  
          swap. */


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

/* Computing MAX */
	d__1 = eps * 10. * dnorm;
	thresh = MAX(d__1,smlnum);

/*        Solve T11*X - X*T22 = scale*T12 for X. */


#ifdef PETSC_PREFIX_SUFFIX
	dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d, &c__4, &D(*n1 + 1 + ((*
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlasy2(&c_false, &c_false, &c_n1, n1, n2, d, &c__4, &D(*n1 + 1 + ((*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlasy2_(&c_false, &c_false, &c_n1, n1, n2, d, &c__4, &D(*n1 + 1 + ((*
#endif

		n1 + 1) << 2) - 5), &c__4, &D((((*n1 + 1) << 2)) - 4), &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

	k = *n1 + *n1 + *n2 - 3;
	switch (k) {
	    case 1:  goto L10;
	    case 2:  goto L20;
	    case 3:  goto L30;
	}

L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that:   

          ( scale, X11, X12 ) H = ( 0, 0, * ) */

	U(0) = scale;
	U(1) = x[0];
	U(2) = x[2];

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&c__3, &U(2), u, &c__1, &tau);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&c__3, &U(2), u, &c__1, &tau);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&c__3, &U(2), u, &c__1, &tau);
#endif

	U(2) = 1.;
	t11 = T(*P_j1,*P_j1);

/*        Perform swap provisionally on diagonal block in D. */


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif


/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = ABS(D(2)), d__3 = ABS(D(6)), d__2 = MAX(d__2,d__3), d__3 = (
		d__1 = D(10) - t11, ABS(d__1));
	if (MAX(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *P_j1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &i__1, u, &tau, &T(*P_j1,*P_j1), ldt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &i__1, u, &tau, &T(*P_j1,*P_j1), ldt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &i__1, u, &tau, &T(*P_j1,*P_j1), ldt, &
#endif

		WORK(1));

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &j2, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &j2, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &j2, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif


	T(j3,*P_j1) = 0.;
	T(j3,j2) = 0.;
	T(j3,j3) = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfx_("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfx("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfx_("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif

		    1));
	}
	goto L40;

L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that:   

          H (  -X11 ) = ( * )   
            (  -X21 ) = ( 0 )   
            ( scale ) = ( 0 ) */

	U(0) = -x[0];
	U(1) = -x[1];
	U(2) = scale;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&c__3, u, &U(1), &c__1, &tau);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&c__3, u, &U(1), &c__1, &tau);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&c__3, u, &U(1), &c__1, &tau);
#endif

	U(0) = 1.;
	t33 = T(j3,j3);

/*        Perform swap provisionally on diagonal block in D. */


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
#endif


/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = ABS(D(1)), d__3 = ABS(D(2)), d__2 = MAX(d__2,d__3), d__3 = (
		d__1 = D(0) - t33, ABS(d__1));
	if (MAX(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &j3, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &j3, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &j3, &c__3, u, &tau, &T(1,*P_j1), ldt, &WORK(1));
#endif

	i__1 = *n - *P_j1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &i__1, u, &tau, &T(*P_j1,j2), ldt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &i__1, u, &tau, &T(*P_j1,j2), ldt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &i__1, u, &tau, &T(*P_j1,j2), ldt, &WORK(
#endif

		1));

	T(*P_j1,*P_j1) = t33;
	T(j2,*P_j1) = 0.;
	T(j3,*P_j1) = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfx_("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfx("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfx_("R", n, &c__3, u, &tau, &Q(1,*P_j1), ldq, &WORK(
#endif

		    1));
	}
	goto L40;

L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2)
 so   
          that:   

          H(2) H(1) (  -X11  -X12 ) = (  *  * )   
                    (  -X21  -X22 )   (  0  * )   
                    ( scale    0  )   (  0  0 )   
                    (    0  scale )   (  0  0 ) */

	U1(0) = -x[0];
	U1(1) = -x[1];
	U1(2) = scale;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&c__3, u1, &U1(1), &c__1, &tau1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&c__3, u1, &U1(1), &c__1, &tau1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&c__3, u1, &U1(1), &c__1, &tau1);
#endif

	U1(0) = 1.;

	temp = -tau1 * (x[2] + U1(1) * x[3]);
	U2(0) = -temp * U1(1) - x[3];
	U2(1) = -temp * U1(2);
	U2(2) = scale;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&c__3, u2, &U2(1), &c__1, &tau2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&c__3, u2, &U2(1), &c__1, &tau2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&c__3, u2, &U2(1), &c__1, &tau2);
#endif

	U2(0) = 1.;

/*        Perform swap provisionally on diagonal block in D. */


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &c__4, u1, &tau1, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &c__4, u1, &tau1, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &c__4, u1, &tau1, d, &c__4, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &c__4, &c__3, u1, &tau1, d, &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &c__4, &c__3, u1, &tau1, d, &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &c__4, &c__3, u1, &tau1, d, &c__4, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &c__4, u2, &tau2, &D(1), &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &c__4, u2, &tau2, &D(1), &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &c__4, u2, &tau2, &D(1), &c__4, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &c__4, &c__3, u2, &tau2, &D(4), &c__4, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &c__4, &c__3, u2, &tau2, &D(4), &c__4, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &c__4, &c__3, u2, &tau2, &D(4), &c__4, &WORK(1));
#endif


/*        Test whether to reject swap.   

   Computing MAX */
	d__1 = ABS(D(2)), d__2 = ABS(D(6)), d__1 = MAX(d__1,d__2), d__2 = ABS(
		D(3)), d__1 = MAX(d__1,d__2), d__2 = ABS(D(7));
	if (MAX(d__1,d__2) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *P_j1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &i__1, u1, &tau1, &T(*P_j1,*P_j1), ldt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &i__1, u1, &tau1, &T(*P_j1,*P_j1), ldt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &i__1, u1, &tau1, &T(*P_j1,*P_j1), ldt, &
#endif

		WORK(1));

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &j4, &c__3, u1, &tau1, &T(1,*P_j1), ldt, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &j4, &c__3, u1, &tau1, &T(1,*P_j1), ldt, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &j4, &c__3, u1, &tau1, &T(1,*P_j1), ldt, &WORK(
#endif

		1));
	i__1 = *n - *P_j1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("L", &c__3, &i__1, u2, &tau2, &T(j2,*P_j1), ldt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("L", &c__3, &i__1, u2, &tau2, &T(j2,*P_j1), ldt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("L", &c__3, &i__1, u2, &tau2, &T(j2,*P_j1), ldt, &
#endif

		WORK(1));

#ifdef PETSC_PREFIX_SUFFIX
	dlarfx_("R", &j4, &c__3, u2, &tau2, &T(1,j2), ldt, &WORK(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfx("R", &j4, &c__3, u2, &tau2, &T(1,j2), ldt, &WORK(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfx_("R", &j4, &c__3, u2, &tau2, &T(1,j2), ldt, &WORK(1)
#endif

		);

	T(j3,*P_j1) = 0.;
	T(j3,j2) = 0.;
	T(j4,*P_j1) = 0.;
	T(j4,j2) = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfx_("R", n, &c__3, u1, &tau1, &Q(1,*P_j1), ldq, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfx("R", n, &c__3, u1, &tau1, &Q(1,*P_j1), ldq, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfx_("R", n, &c__3, u1, &tau1, &Q(1,*P_j1), ldq, &
#endif

		    WORK(1));

#ifdef PETSC_PREFIX_SUFFIX
	    dlarfx_("R", n, &c__3, u2, &tau2, &Q(1,j2), ldq, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfx("R", n, &c__3, u2, &tau2, &Q(1,j2), ldq, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfx_("R", n, &c__3, u2, &tau2, &Q(1,j2), ldq, &WORK(
#endif

		    1));
	}

L40:

	if (*n2 == 2) {

/*           Standardize new 2-by-2 block T11 */


#ifdef PETSC_PREFIX_SUFFIX
	    dlanv2_(&T(*P_j1,*P_j1), &T(*P_j1,j2), &T(j2,*P_j1), &T(j2,j2), &wr1, &wi1, &wr2, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlanv2(&T(*P_j1,*P_j1), &T(*P_j1,j2), &T(j2,*P_j1), &T(j2,j2), &wr1, &wi1, &wr2, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlanv2_(&T(*P_j1,*P_j1), &T(*P_j1,j2), &T(j2,*P_j1), &T(j2,j2), &wr1, &wi1, &wr2, &
#endif

		    wi2, &cs, &sn);
	    i__1 = *n - *P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__1, &T(*P_j1,*P_j1+2), ldt, &T(j2,*P_j1+2), ldt, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__1, &T(*P_j1,*P_j1+2), ldt, &T(j2,*P_j1+2), ldt, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__1, &T(*P_j1,*P_j1+2), ldt, &T(j2,*P_j1+2), ldt, &cs, &sn);
#endif

	    i__1 = *P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__1, &T(1,*P_j1), &c__1, &T(1,j2), &
#endif

		    c__1, &cs, &sn);
	    if (*wantq) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Q(1,*P_j1), &c__1, &Q(1,j2), &
#endif

			c__1, &cs, &sn);
	    }
	}

	if (*n1 == 2) {

/*           Standardize new 2-by-2 block T22 */

	    j3 = *P_j1 + *n2;
	    j4 = j3 + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlanv2_(&T(j3,j3), &T(j3,j4), &T(j4,j3), &T(j4,j4), &wr1, &wi1, &wr2, &wi2, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlanv2(&T(j3,j3), &T(j3,j4), &T(j4,j3), &T(j4,j4), &wr1, &wi1, &wr2, &wi2, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlanv2_(&T(j3,j3), &T(j3,j4), &T(j4,j3), &T(j4,j4), &wr1, &wi1, &wr2, &wi2, &
#endif

		    cs, &sn);
	    if (j3 + 2 <= *n) {
		i__1 = *n - j3 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		drot_(&i__1, &T(j3,j3+2), ldt, &T(j4,j3+2), ldt, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(&i__1, &T(j3,j3+2), ldt, &T(j4,j3+2), ldt, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(&i__1, &T(j3,j3+2), ldt, &T(j4,j3+2), ldt, &cs, &sn);
#endif

	    }
	    i__1 = j3 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__1, &T(1,j3), &c__1, &T(1,j4), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__1, &T(1,j3), &c__1, &T(1,j4), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__1, &T(1,j3), &c__1, &T(1,j4), &
#endif

		    c__1, &cs, &sn);
	    if (*wantq) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Q(1,j3), &c__1, &Q(1,j4), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Q(1,j3), &c__1, &Q(1,j4), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Q(1,j3), &c__1, &Q(1,j4), &
#endif

			c__1, &cs, &sn);
	    }
	}

    }
    return;

/*     Exit with INFO = 1 if swap was rejected. */

L50:
    *info = 1;
    return;

/*     End of DLAEXC */

} /* dlaexc_ */

