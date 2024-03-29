#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlatzm_(char *side, int *m, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlatzm(char *side, int *m, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlatzm_(char *side, int *m, int *n, LONG DOUBLE *
#endif

	v, int *incv, LONG DOUBLE *tau, LONG DOUBLE *c1, LONG DOUBLE *c2, 
	int *ldc, LONG DOUBLE *work)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLATZM applies a Householder matrix generated by DTZRQF to a matrix. 
  

    Let P = I - tau*u*u',   u = ( 1 ),   
                                ( v )   
    where v is an (m-1) vector if SIDE = 'L', or a (n-1) vector if   
    SIDE = 'R'.   

    If SIDE equals 'L', let   
           C = [ C1 ] 1   
               [ C2 ] m-1   
                 n   
    Then C is overwritten by P*C.   

    If SIDE equals 'R', let   
           C = [ C1, C2 ] m   
                  1  n-1   
    Then C is overwritten by C*P.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form P * C   
            = 'R': form C * P   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) LONG DOUBLE PRECISION array, dimension   
                    (1 + (M-1)*ABS(INCV)) if SIDE = 'L'   
                    (1 + (N-1)*ABS(INCV)) if SIDE = 'R'   
            The vector v in the representation of P. V is not used   
            if TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0   

    TAU     (input) LONG DOUBLE PRECISION   
            The value tau in the representation of P.   

    C1      (input/output) LONG DOUBLE PRECISION array, dimension   
                           (LDC,N) if SIDE = 'L'   
                           (M,1)   if SIDE = 'R'   
            On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1   
            if SIDE = 'R'.   

            On exit, the first row of P*C if SIDE = 'L', or the first   
            column of C*P if SIDE = 'R'.   

    C2      (input/output) LONG DOUBLE PRECISION array, dimension   
                           (LDC, N)   if SIDE = 'L'   
                           (LDC, N-1) if SIDE = 'R'   
            On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the   
            m x (n - 1) matrix C2 if SIDE = 'R'.   

            On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P 
  
            if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the arrays C1 and C2. LDC >= (1,M). 
  

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
                        (N) if SIDE = 'L'   
                        (M) if SIDE = 'R'   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b5 = 1.;
    
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dger_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qger(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qger_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    extern long int lsame_(char *, char *);

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


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), daxpy_(int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy(int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy_(int 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *)
	    ;



#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C2(I,J) c2[(I)-1 + ((J)-1)* ( *ldc)]
#define C1(I,J) c1[(I)-1 + ((J)-1)* ( *ldc)]

    if (MIN(*m,*n) == 0 || *tau == 0.) {
	return;
    }

    if (lsame_(side, "L")) {

/*        w := C1 + v' * C2 */


#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &C1(1,1), ldc, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &C1(1,1), ldc, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &C1(1,1), ldc, &WORK(1), &c__1);
#endif

	i__1 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("Transpose", &i__1, n, &c_b5, &C2(1,1), ldc, &V(1), incv,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("Transpose", &i__1, n, &c_b5, &C2(1,1), ldc, &V(1), incv,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("Transpose", &i__1, n, &c_b5, &C2(1,1), ldc, &V(1), incv,
#endif

		 &c_b5, &WORK(1), &c__1);

/*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w'   
          [ C2 ]    [ C2 ]        [ v ] */

	d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	daxpy_(n, &d__1, &WORK(1), &c__1, &C1(1,1), ldc);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qaxpy(n, &d__1, &WORK(1), &c__1, &C1(1,1), ldc);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qaxpy_(n, &d__1, &WORK(1), &c__1, &C1(1,1), ldc);
#endif

	i__1 = *m - 1;
	d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	dger_(&i__1, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C2(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qger(&i__1, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C2(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qger_(&i__1, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C2(1,1), 
#endif

		ldc);

    } else if (lsame_(side, "R")) {

/*        w := C1 + C2 * v */


#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(m, &C1(1,1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(m, &C1(1,1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(m, &C1(1,1), &c__1, &WORK(1), &c__1);
#endif

	i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("No transpose", m, &i__1, &c_b5, &C2(1,1), ldc, &V(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("No transpose", m, &i__1, &c_b5, &C2(1,1), ldc, &V(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("No transpose", m, &i__1, &c_b5, &C2(1,1), ldc, &V(1), 
#endif

		incv, &c_b5, &WORK(1), &c__1);

/*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v'] */

	d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	daxpy_(m, &d__1, &WORK(1), &c__1, &C1(1,1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qaxpy(m, &d__1, &WORK(1), &c__1, &C1(1,1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qaxpy_(m, &d__1, &WORK(1), &c__1, &C1(1,1), &c__1);
#endif

	i__1 = *n - 1;
	d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	dger_(m, &i__1, &d__1, &WORK(1), &c__1, &V(1), incv, &C2(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qger(m, &i__1, &d__1, &WORK(1), &c__1, &V(1), incv, &C2(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qger_(m, &i__1, &d__1, &WORK(1), &c__1, &V(1), incv, &C2(1,1), 
#endif

		ldc);
    }

    return;

/*     End of DLATZM */

} /* dlatzm_ */

