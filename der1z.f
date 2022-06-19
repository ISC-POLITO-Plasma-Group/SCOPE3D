      subroutine der1z(q1,q2)

!	first derivative by FFT; author: F. Califano. Dec 2000

      include 'par.inc'
      dimension q1(nz), q2(nz) 

      call drfftf(nz, q1, wsavez)

      do i = 2, nz-1, 2
         zi    = gradz * float(i) / float(2 * nz)
         q2(i)   = - q1(i+1) * zi 
         q2(i+1) =   q1(i)  * zi 
      enddo

      q2(1)  = 0.0
      q2(nz) = 0.0

      call drfftb(nz, q2, wsavez)
      
      return
      end
