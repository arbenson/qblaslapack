#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaebz_(int *ijob, int *nitmax, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaebz(int *ijob, int *nitmax, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaebz_(int *ijob, int *nitmax, int *n, 
#endif

	int *mmax, int *minp, int *nbmin, LONG DOUBLE *abstol, 
	LONG DOUBLE *reltol, LONG DOUBLE *pivmin, LONG DOUBLE *d, LONG DOUBLE *e, 
	LONG DOUBLE *e2, int *nval, LONG DOUBLE *ab, LONG DOUBLE *c, int 
	*mout, int *nab, LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAEBZ contains the iteration loops which compute and use the   
    function N(w), which is the count of eigenvalues of a symmetric   
    tridiagonal matrix T less than or equal to its argument  w.  It   
    performs a choice of two types of loops:   

    IJOB=1, followed by   
    IJOB=2: It takes as input a list of intervals and returns a list of   
            sufficiently small intervals whose union contains the same   
            eigenvalues as the union of the original intervals.   
            The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.   
            The output interval (AB(j,1),AB(j,2)] will contain   
            eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.   

    IJOB=3: It performs a binary search in each input interval   
            (AB(j,1),AB(j,2)] for a point  w(j)  such that   
            N(w(j))=NVAL(j), and uses  C(j)  as the starting point of   
            the search.  If such a w(j) is found, then on output   
            AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output 
  
            (AB(j,1),AB(j,2)] will be a small interval containing the   
            point where N(w) jumps through NVAL(j), unless that point   
            lies outside the initial interval.   

    Note that the intervals are in all cases half-open intervals,   
    i.e., of the form  (a,b] , which includes  b  but not  a .   

    To avoid underflow, the matrix should be scaled so that its largest   
    element is no greater than  overflow**(1/2) * underflow**(1/4)   
    in absolute value.  To assure the most accurate computation   
    of small eigenvalues, the matrix should be scaled to be   
    not much smaller than that, either.   

    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal   
    Matrix", Report CS41, Computer Science Dept., Stanford   
    University, July 21, 1966   

    Note: the arguments are, in general, *not* checked for unreasonable   
    values.   

    Arguments   
    =========   

    IJOB    (input) INTEGER   
            Specifies what is to be done:   
            = 1:  Compute NAB for the initial intervals.   
            = 2:  Perform bisection iteration to find eigenvalues of T.   
            = 3:  Perform bisection iteration to invert N(w), i.e.,   
                  to find a point which has a specified number of   
                  eigenvalues of T to its left.   
            Other values will cause DLAEBZ to return with INFO=-1.   

    NITMAX  (input) INTEGER   
            The maximum number of "levels" of bisection to be   
            performed, i.e., an interval of width W will not be made   
            smaller than 2^(-NITMAX) * W.  If not all intervals   
            have converged after NITMAX iterations, then INFO is set   
            to the number of non-converged intervals.   

    N       (input) INTEGER   
            The dimension n of the tridiagonal matrix T.  It must be at   
            least 1.   

    MMAX    (input) INTEGER   
            The maximum number of intervals.  If more than MMAX intervals 
  
            are generated, then DLAEBZ will quit with INFO=MMAX+1.   

    MINP    (input) INTEGER   
            The initial number of intervals.  It may not be greater than 
  
            MMAX.   

    NBMIN   (input) INTEGER   
            The smallest number of intervals that should be processed   
            using a vector loop.  If zero, then only the scalar loop   
            will be used.   

    ABSTOL  (input) LONG DOUBLE PRECISION   
            The minimum (absolute) width of an interval.  When an   
            interval is narrower than ABSTOL, or than RELTOL times the   
            larger (in magnitude) endpoint, then it is considered to be   
            sufficiently small, i.e., converged.  This must be at least   
            zero.   

    RELTOL  (input) LONG DOUBLE PRECISION   
            The minimum relative width of an interval.  When an interval 
  
            is narrower than ABSTOL, or than RELTOL times the larger (in 
  
            magnitude) endpoint, then it is considered to be   
            sufficiently small, i.e., converged.  Note: this should   
            always be at least radix*machine epsilon.   

    PIVMIN  (input) LONG DOUBLE PRECISION   
            The minimum absolute value of a "pivot" in the Sturm   
            sequence loop.  This *must* be at least  max |e(j)**2| *   
            safe_min  and at least safe_min, where safe_min is at least   
            the smallest number that can divide one without overflow.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T.   

    E       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The offdiagonal elements of the tridiagonal matrix T in   
            positions 1 through N-1.  E(N) is arbitrary.   

    E2      (input) LONG DOUBLE PRECISION array, dimension (N)   
            The squares of the offdiagonal elements of the tridiagonal   
            matrix T.  E2(N) is ignored.   

    NVAL    (input/output) INTEGER array, dimension (MINP)   
            If IJOB=1 or 2, not referenced.   
            If IJOB=3, the desired values of N(w).  The elements of NVAL 
  
            will be reordered to correspond with the intervals in AB.   
            Thus, NVAL(j) on output will not, in general be the same as   
            NVAL(j) on input, but it will correspond with the interval   
            (AB(j,1),AB(j,2)] on output.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (MMAX,2)   
            The endpoints of the intervals.  AB(j,1) is  a(j), the left   
            endpoint of the j-th interval, and AB(j,2) is b(j), the   
            right endpoint of the j-th interval.  The input intervals   
            will, in general, be modified, split, and reordered by the   
            calculation.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (MMAX)   
            If IJOB=1, ignored.   
            If IJOB=2, workspace.   
            If IJOB=3, then on input C(j) should be initialized to the   
            first search point in the binary search.   

    MOUT    (output) INTEGER   
            If IJOB=1, the number of eigenvalues in the intervals.   
            If IJOB=2 or 3, the number of intervals output.   
            If IJOB=3, MOUT will equal MINP.   

    NAB     (input/output) INTEGER array, dimension (MMAX,2)   
            If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)). 
  
            If IJOB=2, then on input, NAB(i,j) should be set.  It must   
               satisfy the condition:   
               N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),   
               which means that in interval i only eigenvalues   
               NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,   
               NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with   
               IJOB=1.   
               On output, NAB(i,j) will contain   
               MAX(na(k),MIN(nb(k),N(AB(i,j)))), where k is the index of 
  
               the input interval that the output interval   
               (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the   
               the input values of NAB(k,1) and NAB(k,2).   
            If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),   
               unless N(w) > NVAL(i) for all search points  w , in which 
  
               case NAB(i,1) will not be modified, i.e., the output   
               value will be the same as the input value (modulo   
               reorderings -- see NVAL and AB), or unless N(w) < NVAL(i) 
  
               for all search points  w , in which case NAB(i,2) will   
               not be modified.  Normally, NAB should be set to some   
               distinctive value(s) before DLAEBZ is called.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (MMAX)   
            Workspace.   

    IWORK   (workspace) INTEGER array, dimension (MMAX)   
            Workspace.   

    INFO    (output) INTEGER   
            = 0:       All intervals converged.   
            = 1--MMAX: The last INFO intervals did not converge.   
            = MMAX+1:  More than MMAX intervals were generated.   

    Further Details   
    ===============   

        This routine is intended to be called only by other LAPACK   
    routines, thus the interface is less user-friendly.  It is intended   
    for two purposes:   

    (a) finding eigenvalues.  In this case, DLAEBZ should have one or   
        more initial intervals set up in AB, and DLAEBZ should be called 
  
        with IJOB=1.  This sets up NAB, and also counts the eigenvalues. 
  
        Intervals with no eigenvalues would usually be thrown out at   
        this point.  Also, if not all the eigenvalues in an interval i   
        are desired, NAB(i,1) can be increased or NAB(i,2) decreased.   
        For example, set NAB(i,1)=NAB(i,2)-1 to get the largest   
        eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX   
        no smaller than the value of MOUT returned by the call with   
        IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1   
        through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the   
        tolerance specified by ABSTOL and RELTOL.   

    (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l). 
  
        In this case, start with a Gershgorin interval  (a,b).  Set up   
        AB to contain 2 search intervals, both initially (a,b).  One   
        NVAL element should contain  f-1  and the other should contain  l 
  
        , while C should contain a and b, resp.  NAB(i,1) should be -1   
        and NAB(i,2) should be N+1, to flag an error if the desired   
        interval does not lie in (a,b).  DLAEBZ is then called with   
        IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals -- 
  
        j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while   
        if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r   
        >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and   
        N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and   
        w(l-r)=...=w(l+k) are handled similarly.   

    ===================================================================== 
  


       Check for Errors   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3, i__4,  i__5, i__6;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Local variables */
    static int itmp1, itmp2, j, kfnew, klnew, kf, ji, kl, jp, jit;
    static LONG DOUBLE tmp1, tmp2;


#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define E2(I) e2[(I)-1]
#define NVAL(I) nval[(I)-1]
#define C(I) c[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define NAB(I,J) nab[(I)-1 + ((J)-1)* ( *mmax)]
#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *mmax)]

    *info = 0;
    if (*ijob < 1 || *ijob > 3) {
	*info = -1;
	return;
    }

/*     Initialize NAB */

    if (*ijob == 1) {

/*        Compute the number of eigenvalues in the initial intervals. 
*/

	*mout = 0;
	i__1 = *minp;
	for (ji = 1; ji <= *minp; ++ji) {
	    for (jp = 1; jp <= 2; ++jp) {
		tmp1 = D(1) - AB(ji,jp);
		if (ABS(tmp1) < *pivmin) {
		    tmp1 = -(*pivmin);
		}
		NAB(ji,jp) = 0;
		if (tmp1 <= 0.) {
		    NAB(ji,jp) = 1;
		}

		i__2 = *n;
		for (j = 2; j <= *n; ++j) {
		    tmp1 = D(j) - E2(j - 1) / tmp1 - AB(ji,jp);
		    if (ABS(tmp1) < *pivmin) {
			tmp1 = -(*pivmin);
		    }
		    if (tmp1 <= 0.) {
			++NAB(ji,jp);
		    }
/* L10: */
		}
/* L20: */
	    }
	    *mout = *mout + NAB(ji,2) - NAB(ji,1);
/* L30: */
	}
	return;
    }

/*     Initialize for loop   

       KF and KL have the following meaning:   
          Intervals 1,...,KF-1 have converged.   
          Intervals KF,...,KL  still need to be refined. */

    kf = 1;
    kl = *minp;

/*     If IJOB=2, initialize C.   
       If IJOB=3, use the user-supplied starting point. */

    if (*ijob == 2) {
	i__1 = *minp;
	for (ji = 1; ji <= *minp; ++ji) {
	    C(ji) = (AB(ji,1) + AB(ji,2)) * .5;
/* L40: */
	}
    }

/*     Iteration loop */

    i__1 = *nitmax;
    for (jit = 1; jit <= *nitmax; ++jit) {

/*        Loop over intervals */

	if (kl - kf + 1 >= *nbmin && *nbmin > 0) {

/*           Begin of Parallel Version of the loop */

	    i__2 = kl;
	    for (ji = kf; ji <= kl; ++ji) {

/*              Compute N(c), the number of eigenvalues less t
han c */

		WORK(ji) = D(1) - C(ji);
		IWORK(ji) = 0;
		if (WORK(ji) <= *pivmin) {
		    IWORK(ji) = 1;
/* Computing MIN */
		    d__1 = WORK(ji), d__2 = -(*pivmin);
		    WORK(ji) = MIN(d__1,d__2);
		}

		i__3 = *n;
		for (j = 2; j <= *n; ++j) {
		    WORK(ji) = D(j) - E2(j - 1) / WORK(ji) - C(ji);
		    if (WORK(ji) <= *pivmin) {
			++IWORK(ji);
/* Computing MIN */
			d__1 = WORK(ji), d__2 = -(*pivmin);
			WORK(ji) = MIN(d__1,d__2);
		    }
/* L50: */
		}
/* L60: */
	    }

	    if (*ijob <= 2) {

/*              IJOB=2: Choose all intervals containing eigenv
alues. */

		klnew = kl;
		i__2 = kl;
		for (ji = kf; ji <= kl; ++ji) {

/*                 Insure that N(w) is monotone   

   Computing MIN   
   Computing MAX */
		    i__5 = NAB(ji,1), i__6 = IWORK(ji);
		    i__3 = NAB(ji,2), i__4 = MAX(i__5,i__6);
		    IWORK(ji) = MIN(i__3,i__4);

/*                 Update the Queue -- add intervals if bo
th halves   
                   contain eigenvalues. */

		    if (IWORK(ji) == NAB(ji,2)) {

/*                    No eigenvalue in the upper inter
val:   
                      just use the lower interval. */

			AB(ji,2) = C(ji);

		    } else if (IWORK(ji) == NAB(ji,1)) {

/*                    No eigenvalue in the lower inter
val:   
                      just use the upper interval. */

			AB(ji,1) = C(ji);
		    } else {
			++klnew;
			if (klnew <= *mmax) {

/*                       Eigenvalue in both interv
als -- add upper to   
                         queue. */

			    AB(klnew,2) = AB(ji,2);
			    NAB(klnew,2) = NAB(ji,2);
			    AB(klnew,1) = C(ji);
			    NAB(klnew,1) = IWORK(ji);
			    AB(ji,2) = C(ji);
			    NAB(ji,2) = IWORK(ji);
			} else {
			    *info = *mmax + 1;
			}
		    }
/* L70: */
		}
		if (*info != 0) {
		    return;
		}
		kl = klnew;
	    } else {

/*              IJOB=3: Binary search.  Keep only the interval
 containing   
                        w   s.t. N(w) = NVAL */

		i__2 = kl;
		for (ji = kf; ji <= kl; ++ji) {
		    if (IWORK(ji) <= NVAL(ji)) {
			AB(ji,1) = C(ji);
			NAB(ji,1) = IWORK(ji);
		    }
		    if (IWORK(ji) >= NVAL(ji)) {
			AB(ji,2) = C(ji);
			NAB(ji,2) = IWORK(ji);
		    }
/* L80: */
		}
	    }

	} else {

/*           End of Parallel Version of the loop   

             Begin of Serial Version of the loop */

	    klnew = kl;
	    i__2 = kl;
	    for (ji = kf; ji <= kl; ++ji) {

/*              Compute N(w), the number of eigenvalues less t
han w */

		tmp1 = C(ji);
		tmp2 = D(1) - tmp1;
		itmp1 = 0;
		if (tmp2 <= *pivmin) {
		    itmp1 = 1;
/* Computing MIN */
		    d__1 = tmp2, d__2 = -(*pivmin);
		    tmp2 = MIN(d__1,d__2);
		}

/*              A series of compiler directives to defeat vect
orization   
                for the next loop   

   $PL$ CMCHAR=' '   
   DIR$          NEXTSCALAR   
   $DIR          SCALAR   
   DIR$          NEXT SCALAR   
   VD$L          NOVECTOR   
   DEC$          NOVECTOR   
   VD$           NOVECTOR   
   VDIR          NOVECTOR   
   VOCL          LOOP,SCALAR   
   IBM           PREFER SCALAR   
   $PL$ CMCHAR='*' */

		i__3 = *n;
		for (j = 2; j <= *n; ++j) {
		    tmp2 = D(j) - E2(j - 1) / tmp2 - tmp1;
		    if (tmp2 <= *pivmin) {
			++itmp1;
/* Computing MIN */
			d__1 = tmp2, d__2 = -(*pivmin);
			tmp2 = MIN(d__1,d__2);
		    }
/* L90: */
		}

		if (*ijob <= 2) {

/*                 IJOB=2: Choose all intervals containing
 eigenvalues.   

                   Insure that N(w) is monotone   

   Computing MIN   
   Computing MAX */
		    i__5 = NAB(ji,1);
		    i__3 = NAB(ji,2), i__4 = MAX(i__5,itmp1);
		    itmp1 = MIN(i__3,i__4);

/*                 Update the Queue -- add intervals if bo
th halves   
                   contain eigenvalues. */

		    if (itmp1 == NAB(ji,2)) {

/*                    No eigenvalue in the upper inter
val:   
                      just use the lower interval. */

			AB(ji,2) = tmp1;

		    } else if (itmp1 == NAB(ji,1)) {

/*                    No eigenvalue in the lower inter
val:   
                      just use the upper interval. */

			AB(ji,1) = tmp1;
		    } else if (klnew < *mmax) {

/*                    Eigenvalue in both intervals -- 
add upper to queue. */

			++klnew;
			AB(klnew,2) = AB(ji,2);
			NAB(klnew,2) = NAB(ji,2);
			AB(klnew,1) = tmp1;
			NAB(klnew,1) = itmp1;
			AB(ji,2) = tmp1;
			NAB(ji,2) = itmp1;
		    } else {
			*info = *mmax + 1;
			return;
		    }
		} else {

/*                 IJOB=3: Binary search.  Keep only the i
nterval   
                           containing  w  s.t. N(w) = NVAL
 */

		    if (itmp1 <= NVAL(ji)) {
			AB(ji,1) = tmp1;
			NAB(ji,1) = itmp1;
		    }
		    if (itmp1 >= NVAL(ji)) {
			AB(ji,2) = tmp1;
			NAB(ji,2) = itmp1;
		    }
		}
/* L100: */
	    }
	    kl = klnew;

/*           End of Serial Version of the loop */

	}

/*        Check for convergence */

	kfnew = kf;
	i__2 = kl;
	for (ji = kf; ji <= kl; ++ji) {
	    tmp1 = (d__1 = AB(ji,2) - AB(ji,1), ABS(
		    d__1));
/* Computing MAX */
	    d__3 = (d__1 = AB(ji,2), ABS(d__1)), d__4 = (d__2 =
		     AB(ji,1), ABS(d__2));
	    tmp2 = MAX(d__3,d__4);
/* Computing MAX */
	    d__1 = MAX(*abstol,*pivmin), d__2 = *reltol * tmp2;
	    if (tmp1 < MAX(d__1,d__2) || NAB(ji,1) >= NAB(ji,2)) {

/*              Converged -- Swap with position KFNEW,   
                             then increment KFNEW */

		if (ji > kfnew) {
		    tmp1 = AB(ji,1);
		    tmp2 = AB(ji,2);
		    itmp1 = NAB(ji,1);
		    itmp2 = NAB(ji,2);
		    AB(ji,1) = AB(kfnew,1);
		    AB(ji,2) = AB(kfnew,2);
		    NAB(ji,1) = NAB(kfnew,1);
		    NAB(ji,2) = NAB(kfnew,2);
		    AB(kfnew,1) = tmp1;
		    AB(kfnew,2) = tmp2;
		    NAB(kfnew,1) = itmp1;
		    NAB(kfnew,2) = itmp2;
		    if (*ijob == 3) {
			itmp1 = NVAL(ji);
			NVAL(ji) = NVAL(kfnew);
			NVAL(kfnew) = itmp1;
		    }
		}
		++kfnew;
	    }
/* L110: */
	}
	kf = kfnew;

/*        Choose Midpoints */

	i__2 = kl;
	for (ji = kf; ji <= kl; ++ji) {
	    C(ji) = (AB(ji,1) + AB(ji,2)) * .5;
/* L120: */
	}

/*        If no more intervals to refine, quit. */

	if (kf > kl) {
	    goto L140;
	}
/* L130: */
    }

/*     Converged */

L140:
/* Computing MAX */
    i__1 = kl + 1 - kf;
    *info = MAX(i__1,0);
    *mout = kl;

    return;

/*     End of DLAEBZ */

} /* dlaebz_ */

