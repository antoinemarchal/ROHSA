c     Block length of the OpenMP loops in blas.f and lbfgsb.f. A sum over
c     a longer vector is computed block by block and the block sums are
c     added in block order, so the result does not depend on the number
c     of threads. Vectors of at most nblk elements are summed exactly as
c     in the serial code.
      integer nblk
      parameter (nblk = 65536)
