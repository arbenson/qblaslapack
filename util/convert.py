#!/usr/bin/env python
"""
This script translates function names of the form:

   dfunction_name_(

to the form:
   
   BLAS_LAPACK_PREFIXfunction_nameBLAS_LAPACK_SUFFIX(

It does so for all files matching *.c in the current working directory and
stores the files in *.c.bak.

The point of this procedure is to be able to compile a BLAS/LAPACK library
with different prefixes and suffixes, in particular, for quad-precision blas.
For example, with BLAS_LAPACK_PREFIX and BLAS_LAPACK_SUFFIX defined to be
'q' and '' respectively, GEMM would change from dgemm_() to qgemm().

The following command can be used to change *.c.bak files to *.c:
for i in *.bak; do mv "$i" "${i%.bak}"; done
"""

import re
import os

func_header = re.compile('d\w*_\(')
src_files = filter(lambda x: x[-2:] == '.c', os.listdir('.'))

for src_file in src_files:
  with open(src_file, 'r') as f_src:
    with open(src_file + '.bak', 'w') as f_bak:
      for line in f_src:
        match = re.search(func_header, line)
        if match != None:
          f_bak.write('\n')
          f_bak.write('#ifdef PETSC_PREFIX_SUFFIX\n')
          f_bak.write(line)
          f_bak.write('#endif\n')
          f_bak.write('#ifdef Q_C_PREFIX_SUFFIX\n')
          f_bak.write('%s%s%s(%s' % (line[0:match.start()],
                                     'q',
                                     line[match.start()+1:match.end()-2],
                                     line[match.end():]))
          f_bak.write('#endif\n')
          f_bak.write('#ifdef Q_NORMAL_PREFIX_SUFFIX\n')
          f_bak.write('%s%s%s' % (line[0:match.start()],
                                  'q',
                                  line[match.start()+1:]))
          f_bak.write('#endif\n')
          f_bak.write('\n')
        else:
          f_bak.write(line)
