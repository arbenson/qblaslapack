#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlahqr_(long int *wantt, long int *wantz, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlahqr(long int *wantt, long int *wantz, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlahqr_(long int *wantt, long int *wantz, int *n, 
#endif

	int *ilo, int *ihi, LONG DOUBLE *h, int *ldh, LONG DOUBLE *
	wr, LONG DOUBLE *wi, int *iloz, int *ihiz, LONG DOUBLE *z, 
	int *ldz, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAHQR is an auxiliary routine called by DHSEQR to update the   
    eigenvalues and Schur decomposition already computed by DHSEQR, by   
    dealing with the Hessenberg submatrix in rows and columns ILO to IHI. 
  

    Arguments   
    =========   

    WANTT   (input) LOGICAL   
            = .TRUE. : the full Schur form T is required;   
            = .FALSE.: only eigenvalues are required.   

    WANTZ   (input) LOGICAL   
            = .TRUE. : the matrix of Schur vectors Z is required;   
            = .FALSE.: Schur vectors are not required.   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that H is already upper quasi-triangular in   
            rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless   
            ILO = 1). DLAHQR works primarily with the Hessenberg   
            submatrix in rows and columns ILO to IHI, but applies   
            transformations to all of H if WANTT is .TRUE..   
            1 <= ILO <= MAX(1,IHI); IHI <= N.   

    H       (input/output) LONG DOUBLE PRECISION array, dimension (LDH,N)   
            On entry, the upper Hessenberg matrix H.   
            On exit, if WANTT is .TRUE., H is upper quasi-triangular in   
            rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in 
  
            standard form. If WANTT is .FALSE., the contents of H are   
            unspecified on exit.   

    LDH     (input) INTEGER   
            The leading dimension of the array H. LDH >= MAX(1,N).   

    WR      (output) LONG DOUBLE PRECISION array, dimension (N)   
    WI      (output) LONG DOUBLE PRECISION array, dimension (N)   
            The real and imaginary parts, respectively, of the computed   
            eigenvalues ILO to IHI are stored in the corresponding   
            elements of WR and WI. If two eigenvalues are computed as a   
            complex conjugate pair, they are stored in consecutive   
            elements of WR and WI, say the i-th and (i+1)th, with   
            WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the   
            eigenvalues are stored in the same order as on the diagonal   
            of the Schur form returned in H, with WR(i) = H(i,i), and, if 
  
            H(i:i+1,i:i+1) is a 2-by-2 diagonal block,   
            WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).   

    ILOZ    (input) INTEGER   
    IHIZ    (input) INTEGER   
            Specify the rows of Z to which transformations must be   
            applied if WANTZ is .TRUE..   
            1 <= ILOZ <= ILO; IHI <= IHIZ <= N.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ,N)   
            If WANTZ is .TRUE., on entry Z must contain the current   
            matrix Z of transformations accumulated by DHSEQR, and on   
            exit Z has been updated; transformations are applied only to 
  
            the submatrix Z(ILOZ:IHIZ,ILO:IHI).   
            If WANTZ is .FALSE., Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z. LDZ >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            > 0: DLAHQR failed to compute all the eigenvalues ILO to IHI 
  
                 in a total of 30*(IHI-ILO+1) iterations; if INFO = i,   
                 elements i+1:ihi of WR and WI contain those eigenvalues 
  
                 which have been successfully computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static LONG DOUBLE h43h34, unfl, ovfl;

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
    static LONG DOUBLE work[1];
    static int i, j, k, l, m;
    static LONG DOUBLE s, v[3];

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
    static int i1, i2;
    static LONG DOUBLE t1, t2, t3, v1, v2, v3;

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
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), dlabad_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlabad(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlabad_(
#endif

	    LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE h00, h10, h11, h12, h21, h22, h33, h44;
    static int nh;
    static LONG DOUBLE cs;

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
    extern /* Subroutine */ void dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *);
    static int nr;
    static LONG DOUBLE sn;
    static int nz;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *);
    static LONG DOUBLE smlnum, h33s, h44s;
    static int itn, its;
    static LONG DOUBLE ulp, sum, tst1;



#define WORK(I) work[(I)]
#define V(I) v[(I)]
#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }
    if (*ilo == *ihi) {
	WR(*ilo) = H(*ilo,*ilo);
	WI(*ilo) = 0.;
	return;
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion.   
       If norm(H) <= sqrt(OVFL), overflow should not occur. */


#ifdef PETSC_PREFIX_SUFFIX
    unfl = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    unfl = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    unfl = qlamch_("Safe minimum");
#endif

    ovfl = 1. / unfl;

#ifdef PETSC_PREFIX_SUFFIX
    dlabad_(&unfl, &ovfl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlabad(&unfl, &ovfl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlabad_(&unfl, &ovfl);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    ulp = dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    ulp = qlamch("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    ulp = qlamch_("Precision");
#endif

    smlnum = unfl * (nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H   
       to which transformations must be applied. If eigenvalues only are 
  
       being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITN is the total number of QR iterations allowed. */

    itn = nh * 30;

/*     The main loop begins here. I is the loop index and decreases from 
  
       IHI to ILO in steps of 1 or 2. Each iteration of the loop works   
       with the active submatrix in rows and columns L to I.   
       Eigenvalues I+1 to IHI have already converged. Either L = ILO or   
       H(L,L-1) is negligible so that the matrix splits. */

    i = *ihi;
L10:
    l = *ilo;
    if (i < *ilo) {
	goto L150;
    }

/*     Perform QR iterations on rows and columns ILO to I until a   
       submatrix of order 1 or 2 splits off at the bottom because a   
       subdiagonal element has become negligible. */

    i__1 = itn;
    for (its = 0; its <= itn; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i; k >= l+1; --k) {
	    tst1 = (d__1 = H(k-1,k-1), ABS(d__1)) + (d__2 = 
		    H(k,k), ABS(d__2));
	    if (tst1 == 0.) {
		i__3 = i - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
		tst1 = dlanhs_("1", &i__3, &H(l,l), ldh, work);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		tst1 = qlanhs("1", &i__3, &H(l,l), ldh, work);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		tst1 = qlanhs_("1", &i__3, &H(l,l), ldh, work);
#endif

	    }
/* Computing MAX */
	    d__2 = ulp * tst1;
	    if ((d__1 = H(k,k-1), ABS(d__1)) <= MAX(d__2,
		    smlnum)) {
		goto L30;
	    }
/* L20: */
	}
L30:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    H(l,l-1) = 0.;
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off.
 */

	if (l >= i - 1) {
	    goto L140;
	}

/*        Now the active submatrix is in rows and columns L to I. If 
  
          eigenvalues only are being computed, only the active submatr
ix   
          need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i;
	}

	if (its == 10 || its == 20) {

/*           Exceptional shift. */

	    s = (d__1 = H(i,i-1), ABS(d__1)) + (d__2 = H(i-1,i-2), ABS(d__2));
	    h44 = s * .75;
	    h33 = h44;
	    h43h34 = s * -.4375 * s;
	} else {

/*           Prepare to use Wilkinson's LONG DOUBLE shift */

	    h44 = H(i,i);
	    h33 = H(i-1,i-1);
	    h43h34 = H(i,i-1) * H(i-1,i);
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l;
	for (m = i - 2; m >= l; --m) {

/*           Determine the effect of starting the LONG DOUBLE-shift QR 
  
             iteration at row M, and see if this would make H(M,M-
1)   
             negligible. */

	    h11 = H(m,m);
	    h22 = H(m+1,m+1);
	    h21 = H(m+1,m);
	    h12 = H(m,m+1);
	    h44s = h44 - h11;
	    h33s = h33 - h11;
	    v1 = (h33s * h44s - h43h34) / h21 + h12;
	    v2 = h22 - h11 - h33s - h44s;
	    v3 = H(m+2,m+1);
	    s = ABS(v1) + ABS(v2) + ABS(v3);
	    v1 /= s;
	    v2 /= s;
	    v3 /= s;
	    V(0) = v1;
	    V(1) = v2;
	    V(2) = v3;
	    if (m == l) {
		goto L50;
	    }
	    h00 = H(m-1,m-1);
	    h10 = H(m,m-1);
	    tst1 = ABS(v1) * (ABS(h00) + ABS(h11) + ABS(h22));
	    if (ABS(h10) * (ABS(v2) + ABS(v3)) <= ulp * tst1) {
		goto L50;
	    }
/* L40: */
	}
L50:

/*        LONG DOUBLE-shift QR step */

	i__2 = i - 1;
	for (k = m; k <= i-1; ++k) {

/*           The first iteration of this loop determines a reflect
ion G   
             from the vector V and applies it from left and right 
to H,   
             thus creating a nonzero bulge below the subdiagonal. 
  

             Each subsequent iteration determines a reflection G t
o   
             restore the Hessenberg form in the (K-1)th column, an
d thus   
             chases the bulge one step toward the bottom of the ac
tive   
             submatrix. NR is the order of G.   

   Computing MIN */
	    i__3 = 3, i__4 = i - k + 1;
	    nr = MIN(i__3,i__4);
	    if (k > m) {

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&nr, &H(k,k-1), &c__1, v, &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&nr, &H(k,k-1), &c__1, v, &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&nr, &H(k,k-1), &c__1, v, &c__1);
#endif

	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlarfg_(&nr, v, &V(1), &c__1, &t1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfg(&nr, v, &V(1), &c__1, &t1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfg_(&nr, v, &V(1), &c__1, &t1);
#endif

	    if (k > m) {
		H(k,k-1) = V(0);
		H(k+1,k-1) = 0.;
		if (k < i - 1) {
		    H(k+2,k-1) = 0.;
		}
	    } else if (m > l) {
		H(k,k-1) = -H(k,k-1);
	    }
	    v2 = V(1);
	    t2 = t1 * v2;
	    if (nr == 3) {
		v3 = V(2);
		t3 = t1 * v3;

/*              Apply G from the left to transform the rows of
 the matrix   
                in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i2; ++j) {
		    sum = H(k,j) + v2 * H(k+1,j) + v3 
			    * H(k+2,j);
		    H(k,j) -= sum * t1;
		    H(k+1,j) -= sum * t2;
		    H(k+2,j) -= sum * t3;
/* L60: */
		}

/*              Apply G from the right to transform the column
s of the   
                matrix in rows I1 to MIN(K+3,I).   

   Computing MIN */
		i__4 = k + 3;
		i__3 = MIN(i__4,i);
		for (j = i1; j <= MIN(k+3,i); ++j) {
		    sum = H(j,k) + v2 * H(j,k+1) + 
			    v3 * H(j,k+2);
		    H(j,k) -= sum * t1;
		    H(j,k+1) -= sum * t2;
		    H(j,k+2) -= sum * t3;
/* L70: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matri
x Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= *ihiz; ++j) {
			sum = Z(j,k) + v2 * Z(j,k+1)
				 + v3 * Z(j,k+2);
			Z(j,k) -= sum * t1;
			Z(j,k+1) -= sum * t2;
			Z(j,k+2) -= sum * t3;
/* L80: */
		    }
		}
	    } else if (nr == 2) {

/*              Apply G from the left to transform the rows of
 the matrix   
                in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i2; ++j) {
		    sum = H(k,j) + v2 * H(k+1,j);
		    H(k,j) -= sum * t1;
		    H(k+1,j) -= sum * t2;
/* L90: */
		}

/*              Apply G from the right to transform the column
s of the   
                matrix in rows I1 to MIN(K+3,I). */

		i__3 = i;
		for (j = i1; j <= i; ++j) {
		    sum = H(j,k) + v2 * H(j,k+1);
		    H(j,k) -= sum * t1;
		    H(j,k+1) -= sum * t2;
/* L100: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matri
x Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= *ihiz; ++j) {
			sum = Z(j,k) + v2 * Z(j,k+1)
				;
			Z(j,k) -= sum * t1;
			Z(j,k+1) -= sum * t2;
/* L110: */
		    }
		}
	    }
/* L120: */
	}

/* L130: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i;
    return;

L140:

    if (l == i) {

/*        H(I,I-1) is negligible: one eigenvalue has converged. */

	WR(i) = H(i,i);
	WI(i) = 0.;
    } else if (l == i - 1) {

/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converg
ed.   

          Transform the 2-by-2 submatrix to standard Schur form,   
          and compute and store the eigenvalues. */


#ifdef PETSC_PREFIX_SUFFIX
	dlanv2_(&H(i-1,i-1), &H(i-1,i), &H(i,i-1), &H(i,i), &WR(i - 1), &WI(i - 1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlanv2(&H(i-1,i-1), &H(i-1,i), &H(i,i-1), &H(i,i), &WR(i - 1), &WI(i - 1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlanv2_(&H(i-1,i-1), &H(i-1,i), &H(i,i-1), &H(i,i), &WR(i - 1), &WI(i - 1), 
#endif

		&WR(i), &WI(i), &cs, &sn);

	if (*wantt) {

/*           Apply the transformation to the rest of H. */

	    if (i2 > i) {
		i__1 = i2 - i;

#ifdef PETSC_PREFIX_SUFFIX
		drot_(&i__1, &H(i-1,i+1), ldh, &H(i,i+1), ldh, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(&i__1, &H(i-1,i+1), ldh, &H(i,i+1), ldh, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(&i__1, &H(i-1,i+1), ldh, &H(i,i+1), ldh, &cs, &sn);
#endif

	    }
	    i__1 = i - i1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__1, &H(i1,i-1), &c__1, &H(i1,i)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__1, &H(i1,i-1), &c__1, &H(i1,i)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__1, &H(i1,i-1), &c__1, &H(i1,i)
#endif

		    , &c__1, &cs, &sn);
	}
	if (*wantz) {

/*           Apply the transformation to Z. */


#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&nz, &Z(*iloz,i-1), &c__1, &Z(*iloz,i), &c__1, &cs, &sn);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&nz, &Z(*iloz,i-1), &c__1, &Z(*iloz,i), &c__1, &cs, &sn);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&nz, &Z(*iloz,i-1), &c__1, &Z(*iloz,i), &c__1, &cs, &sn);
#endif

	}
    }

/*     Decrement number of remaining iterations, and return to start of   
       the main loop with new value of I. */

    itn -= its;
    i = l - 1;
    goto L10;

L150:
    return;

/*     End of DLAHQR */

} /* dlahqr_ */

