# $Id: makefile,v 1.10 2001/04/03 16:09:02 balay Exp $ 

SOURCEC  = \
 dlaln2.c dlamch.c dlamrg.c dlangb.c dlangt.c dlanhs.c dlansb.c \
 dlansp.c dlanst.c dlansy.c dlantb.c dlantp.c dlantr.c dlanv2.c \
 dlapll.c dlapmt.c dlapy3.c dlaqgb.c dlaqge.c dlaqsb.c dlaqsp.c \
 dlaqsy.c dlaqtr.c dlar2v.c dlarf.c dlarfb.c dlarfg.c dlarft.c \
 dlarfx.c dlargv.c dlarnv.c dlartv.c dlaruv.c dlas2.c dlascl.c \
 dlasq1.c dlasq2.c dlasq3.c dlasq4.c dlasr.c dlasrt.c dlassq.c \
 dlasv2.c dlasy2.c dlasyf.c dlatbs.c dlatps.c dlatrd.c dlatrs.c \
 dlatzm.c dlauu2.c dlauum.c dopgtr.c dopmtr.c dorg2l.c dorg2r.c \
 dorgbr.c dorgl2.c dorglq.c dorgql.c dorgqr.c dorgr2.c dorgrq.c \
 dorgtr.c dorm2l.c dorm2r.c dormhr.c dorml2.c dormlq.c dormql.c \
 dormqr.c dormr2.c dormrq.c dormtr.c dpbcon.c dpbequ.c dpbrfs.c \
 dpbstf.c dpbsv.c dpbsvx.c dpbtf2.c dpbtrf.c dpbtrs.c 
OBJSC    = \
 dlaln2.o dlamch.o dlamrg.o dlangb.o dlangt.o dlanhs.o dlansb.o \
 dlansp.o dlanst.o dlansy.o dlantb.o dlantp.o dlantr.o dlanv2.o \
 dlapll.o dlapmt.o dlapy3.o dlaqgb.o dlaqge.o dlaqsb.o dlaqsp.o \
 dlaqsy.o dlaqtr.o dlar2v.o dlarf.o dlarfb.o dlarfg.o dlarft.o \
 dlarfx.o dlargv.o dlarnv.o dlartv.o dlaruv.o dlas2.o dlascl.o \
 dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasr.o dlasrt.o dlassq.o \
 dlasv2.o dlasy2.o dlasyf.o dlatbs.o dlatps.o dlatrd.o dlatrs.o \
 dlatzm.o dlauu2.o dlauum.o dopgtr.o dopmtr.o dorg2l.o dorg2r.o \
 dorgbr.o dorgl2.o dorglq.o dorgql.o dorgqr.o dorgr2.o dorgrq.o \
 dorgtr.o dorm2l.o dorm2r.o dormhr.o dorml2.o dormlq.o dormql.o \
 dormqr.o dormr2.o dormrq.o dormtr.o dpbcon.o dpbequ.o dpbrfs.o \
 dpbstf.o dpbsv.o dpbsvx.o dpbtf2.o dpbtrf.o dpbtrs.o 

SOURCEC_SPECIAL = dlamch.c
OBJSC_SPECIAL   = dlamch.o

lib:
	$(CC) -c $(CNOOPT) $(SOURCEC_SPECIAL)
	$(CC) -c $(COPTFLAGS) $(SOURCEC)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSC_SPECIAL) $(OBJSC)
	$(RM) $(OBJSC_SPECIAL) $(OBJSC)
