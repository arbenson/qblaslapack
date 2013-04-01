# $Id: makefile,v 1.1 2001/04/03 00:34:33 balay Exp $ 

SOURCEC  = \
dbdsqr.c  dgelqf.c  dhseqr.c  dlacpy.c  dlapy2.c  dlaset.c  dorghr.c \
dgehrd.c  dgetf2.c  dlabad.c  dlange.c  dlartg.c  dlaswp.c  dormbr.c
OBJSC    = \
dbdsqr.o  dgelqf.o  dhseqr.o  dlacpy.o  dlapy2.o  dlaset.o  dorghr.o \
dgehrd.o  dgetf2.o  dlabad.o  dlange.o  dlartg.o  dlaswp.o  dormbr.o

lib:
	$(CC) -c $(COPTFLAGS) $(SOURCEC)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSC)
	$(RM) $(OBJSC)
