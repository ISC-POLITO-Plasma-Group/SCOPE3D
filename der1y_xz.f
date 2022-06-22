       subroutine der1y_xz(aux_t)

!	first derivative by FFT; author: F. Califano. Dec 2000
       use nvtx
	include 'par.inc'
       dimension q1(ny), q2(ny), aux_t(nxl,ny,nzl) 
       call nvtxStartRange('der1y_xz',18)

       
        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 q1(iy) = aux_t(ix,iy,iz)
              enddo
              call drfftf(ny, q1, wsavey)
             do iy = 2, ny-1, 2
	        zi    = grady * float(iy) / float(2 * ny)
                q2(iy)   = - q1(iy+1) * zi 
                q2(iy+1) =   q1(iy)  * zi 
	     enddo
             q2(1)  = 0.0
             q2(ny) = 0.0

             call drfftb(ny, q2, wsavey)
             do iy = 1,ny
                 aux_t(ix,iy,iz) = q2(iy)
             enddo
           enddo
        enddo

       call nvtxEndRange()
	return
	end
