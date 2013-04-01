# $Id: makefile,v 1.12 2001/04/03 16:08:59 balay Exp $ 

SOURCEC  = \
 dasum.c   dgbmv.c  drot.c   dspr.c   dsyr.c    dtbsv.c  dtrsm.c \
 daxpy.c   dgemm.c  drotg.c  dspr2.c  dsyr2.c   dtpmv.c  dtrsv.c \
 dcabs1.c  dgemv.c  dsbmv.c  dswap.c  dsyr2k.c  dtpsv.c  idamax.c \
 dcopy.c   dger.c   dscal.c  dsymm.c  dsyrk.c   dtrmm.c  lsame.c \
 ddot.c    dnrm2.c  dspmv.c  dsymv.c  dtbmv.c   dtrmv.c  xerbla.c f2cb.c
OBJSC    = \
 dasum.o   dgbmv.o  drot.o   dspr.o   dsyr.o    dtbsv.o  dtrsm.o \
 daxpy.o   dgemm.o  drotg.o  dspr2.o  dsyr2.o   dtpmv.o  dtrsv.o \
 dcabs1.o  dgemv.o  dsbmv.o  dswap.o  dsyr2k.o  dtpsv.o  idamax.o \
 dcopy.o   dger.o   dscal.o  dsymm.o  dsyrk.o   dtrmm.o  lsame.o \
 ddot.o    dnrm2.o  dspmv.o  dsymv.o  dtbmv.o   dtrmv.o  xerbla.o f2cb.o

lib:
	$(CC) -c $(COPTFLAGS) $(SOURCEC)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSC)
	$(RM) $(OBJSC)
