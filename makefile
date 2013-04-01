# $Id: makefile,v 1.12 2001/11/29 19:38:51 balay Exp balay $ 

########################################################################################
# Please read the readme file before editing the makefile
########################################################################################

all: blas_lib lapack_lib

########################################################################################
# Specify options to compile and create libraries
########################################################################################
# NAME_TYPE is one of:
#
#   PETSC_PREFIX_SUFFIX
#   Q_NORMAL_PREFIX_SUFFIX
#   Q_C_PREFIX_SUFFIX
#
# which creates function names of the form:
#
#   daxpy_()
#   qaxpy_()
#   qaxpy()
#
# respectively.
NAME_TYPE  = Q_C_PREFIX_SUFFIX
COPTFLAGS  = -O -DDOUBLE=double -DLONG="" -D$(NAME_TYPE)=1
CNOOPT     = -O0 -DDOUBLE=double -DLONG="" -D$(NAME_TYPE)=1
CC         = cc
OMAKE      = make
RM         = /bin/rm
AR         = ar
AR_FLAGS   = cr
LIB_SUFFIX = a
RANLIB     = ranlib
########################################################################################
# By default, pick up the options from the PETSc configuration files
########################################################################################
BOPT             = O
BLASLAPACK_TYPE  = F2CBLASLAPACK

########################################################################################
# compile the source files and create the blas and lapack libs
########################################################################################

BLAS_LIB_NAME       = libqblas.$(LIB_SUFFIX)
LAPACK_LIB_NAME     = libqlapack.$(LIB_SUFFIX)
MAKE_OPTIONS        =  CC="$(CC)" COPTFLAGS="$(COPTFLAGS)" CNOOPT="$(CNOOPT)" AR="$(AR)" AR_FLAGS="$(AR_FLAGS)" RM="$(RM)"
MAKE_OPTIONS_BLAS   = $(MAKE_OPTIONS) LIBNAME="$(BLAS_LIB_NAME)"
MAKE_OPTIONS_LAPACK = $(MAKE_OPTIONS) LIBNAME="$(LAPACK_LIB_NAME)"

blas_lib:
	-@cd blas1;   $(OMAKE) lib $(MAKE_OPTIONS_BLAS)
	-@$(RANLIB) $(BLAS_LIB_NAME)

lapack_lib:
	-@cd lapack1; $(OMAKE) lib $(MAKE_OPTIONS_LAPACK)
	-@cd lapack2; $(OMAKE) lib $(MAKE_OPTIONS_LAPACK)
	-@cd lapack3; $(OMAKE) lib $(MAKE_OPTIONS_LAPACK)
	-@cd lapack4; $(OMAKE) lib $(MAKE_OPTIONS_LAPACK)
	-@$(RANLIB) $(LAPACK_LIB_NAME)

check_blaslapackdir:
	@if [ ! -f  lapack4/dbdsqr.c ]; then \
	echo "make invoked from incorrect location!"; \
	echo "Aborting build"; \
	false; fi

cleanblaslapck: check_blaslapackdir
	$(RM) */*.o

cleanlib: check_blaslapackdir
	$(RM) ./*.a ./*.lib



