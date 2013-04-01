# $Id: makefile,v 1.10 2001/04/03 16:09:00 balay Exp $ 

SOURCEC  = \
 ddisna.c dgbbrd.c dgbcon.c dgbequ.c dgbrfs.c dgbsv.c dgbsvx.c dgbtf2.c \
 dgbtrf.c dgbtrs.c dgebak.c dgebal.c dgebd2.c dgebrd.c dgecon.c dgeequ.c \
 dgees.c dgeesx.c dgeev.c dgeevx.c dgegs.c dgegv.c dgehd2.c dgelq2.c \
 dgels.c dgelss.c dgelsx.c dgeql2.c dgeqlf.c dgeqpf.c dgeqr2.c dgeqrf.c \
 dgerfs.c dgerq2.c dgerqf.c dgesv.c dgesvd.c dgesvx.c dgetrf.c dgetri.c \
 dgetrs.c dggbak.c dggbal.c dggglm.c dgghrd.c dgglse.c dggqrf.c dggrqf.c \
 dggsvd.c dggsvp.c dgtcon.c dgtrfs.c dgtsv.c dgtsvx.c dgttrf.c dgttrs.c \
 dhgeqz.c dhsein.c dlabrd.c dlacon.c dladiv.c dlae2.c dlaebz.c dlaed0.c \
 dlaed1.c dlaed2.c dlaed3.c dlaed4.c dlaed5.c dlaed6.c dlaed7.c dlaed8.c \
 dlaed9.c dlaeda.c dlaein.c dlaev2.c dlaexc.c dlag2.c dlags2.c dlagtf.c \
 dlagtm.c dlagts.c dlahqr.c dlahrd.c dlaic1.c 
OBJSC    = \
 ddisna.o dgbbrd.o dgbcon.o dgbequ.o dgbrfs.o dgbsv.o dgbsvx.o dgbtf2.o \
 dgbtrf.o dgbtrs.o dgebak.o dgebal.o dgebd2.o dgebrd.o dgecon.o dgeequ.o \
 dgees.o dgeesx.o dgeev.o dgeevx.o dgegs.o dgegv.o dgehd2.o dgelq2.o \
 dgels.o dgelss.o dgelsx.o dgeql2.o dgeqlf.o dgeqpf.o dgeqr2.o dgeqrf.o \
 dgerfs.o dgerq2.o dgerqf.o dgesv.o dgesvd.o dgesvx.o dgetrf.o dgetri.o \
 dgetrs.o dggbak.o dggbal.o dggglm.o dgghrd.o dgglse.o dggqrf.o dggrqf.o \
 dggsvd.o dggsvp.o dgtcon.o dgtrfs.o dgtsv.o dgtsvx.o dgttrf.o dgttrs.o \
 dhgeqz.o dhsein.o dlabrd.o dlacon.o dladiv.o dlae2.o dlaebz.o dlaed0.o \
 dlaed1.o dlaed2.o dlaed3.o dlaed4.o dlaed5.o dlaed6.o dlaed7.o dlaed8.o \
 dlaed9.o dlaeda.o dlaein.o dlaev2.o dlaexc.o dlag2.o dlags2.o dlagtf.o \
 dlagtm.o dlagts.o dlahqr.o dlahrd.o dlaic1.o 

lib:
	$(CC) -c $(COPTFLAGS) $(SOURCEC)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSC)
	$(RM) $(OBJSC)
