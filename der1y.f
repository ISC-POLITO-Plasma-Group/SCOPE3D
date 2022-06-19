       subroutine der1y(q1,q2)

!	first derivative by FFT; author: F. Califano. Dec 2000
       use nvtx
	include 'par.inc'
       dimension q1(ny), q2(ny) 
       call nvtxStartRange('der1y',18)
       call drfftf(ny, q1, wsavey)

       do i = 2, ny-1, 2
	zi    = grady * float(i) / float(2 * ny)
        q2(i)   = - q1(i+1) * zi 
        q2(i+1) =   q1(i)  * zi 
	enddo

       q2(1)  = 0.0
       q2(ny) = 0.0

       call drfftb(ny, q2, wsavey)
       call nvtxEndRange()
	return
	end
