      subroutine der2y(q1,q2)
      use nvtx

      include 'par.inc'
      dimension q1(ny), q2(ny) 
      call nvtxStartRange('der2y',21)
      call drfftf(ny, q1, wsavey)

      do i = 2, ny-1, 2
         z1    = grady * float(i) / 2.0
         z2    = - z1 * z1 / float(ny)
         q2(i)   = q1(i)    * z2 
         q2(i+1) = q1(i+1)  * z2 
      enddo

      q2(1)  = 0.0
      q2(ny) = 0.0

      call drfftb(ny, q2, wsavey)
      call nvtxEndRange()
      return
      end
