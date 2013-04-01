/* dummy routine to identify f2cblaslapack library */


#ifdef PETSC_PREFIX_SUFFIX
int f2cblaslapack_id_(int *n)
#endif
#ifdef Q_C_PREFIX_SUFFIX
int f2cblaslapack_iq(int *n)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
int f2cblaslapack_iq_(int *n)
#endif

{
  *n = 1;
  return 0;
}
