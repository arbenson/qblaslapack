# $Id: makefile,v 1.11 2001/04/03 16:09:03 balay Exp $ 

SOURCEC  = \
 dpocon.c dpoequ.c dporfs.c dposv.c dposvx.c dpotf2.c dpotrf.c \
 dpotri.c dpotrs.c dppcon.c dppequ.c dpprfs.c dppsv.c dppsvx.c \
 dpptrf.c dpptri.c dpptrs.c dptcon.c dpteqr.c dptrfs.c dptsv.c \
 dptsvx.c dpttrf.c dpttrs.c drscl.c dsbev.c dsbevd.c dsbevx.c \
 dsbgst.c dsbgv.c dsbtrd.c dspcon.c dspev.c dspevd.c dspevx.c \
 dspgst.c dspgv.c dsprfs.c dspsv.c dspsvx.c dsptrd.c dsptrf.c \
 dsptri.c dsptrs.c dstebz.c dstedc.c dstein.c dsteqr.c dsterf.c \
 dstev.c dstevd.c dstevx.c dsycon.c dsyev.c dsyevd.c dsyevx.c \
 dsygs2.c dsygst.c dsygv.c dsyrfs.c dsysv.c dsysvx.c dsytd2.c \
 dsytf2.c dsytrd.c dsytrf.c dsytri.c dsytrs.c dtbcon.c dtbrfs.c \
 dtbtrs.c dtgevc.c dtgsja.c dtpcon.c dtprfs.c dtptri.c dtptrs.c \
 dtrcon.c dtrevc.c dtrexc.c dtrrfs.c dtrsen.c dtrsna.c dtrsyl.c \
 dtrti2.c dtrtri.c dtrtrs.c dtzrqf.c ilaenv.c
OBJSC    = \
 dpocon.o dpoequ.o dporfs.o dposv.o dposvx.o dpotf2.o dpotrf.o \
 dpotri.o dpotrs.o dppcon.o dppequ.o dpprfs.o dppsv.o dppsvx.o \
 dpptrf.o dpptri.o dpptrs.o dptcon.o dpteqr.o dptrfs.o dptsv.o \
 dptsvx.o dpttrf.o dpttrs.o drscl.o dsbev.o dsbevd.o dsbevx.o \
 dsbgst.o dsbgv.o dsbtrd.o dspcon.o dspev.o dspevd.o dspevx.o \
 dspgst.o dspgv.o dsprfs.o dspsv.o dspsvx.o dsptrd.o dsptrf.o \
 dsptri.o dsptrs.o dstebz.o dstedc.o dstein.o dsteqr.o dsterf.o \
 dstev.o dstevd.o dstevx.o dsycon.o dsyev.o dsyevd.o dsyevx.o \
 dsygs2.o dsygst.o dsygv.o dsyrfs.o dsysv.o dsysvx.o dsytd2.o \
 dsytf2.o dsytrd.o dsytrf.o dsytri.o dsytrs.o dtbcon.o dtbrfs.o \
 dtbtrs.o dtgevc.o dtgsja.o dtpcon.o dtprfs.o dtptri.o dtptrs.o \
 dtrcon.o dtrevc.o dtrexc.o dtrrfs.o dtrsen.o dtrsna.o dtrsyl.o \
 dtrti2.o dtrtri.o dtrtrs.o dtzrqf.o ilaenv.o

lib:
	$(CC) -c $(COPTFLAGS) $(SOURCEC)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSC)
	$(RM) $(OBJSC)
