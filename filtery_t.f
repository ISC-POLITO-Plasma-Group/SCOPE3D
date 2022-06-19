          subroutine filtery_t(g)

          use nvtx
          include 'par.inc'


          dimension filtroy(ny)
          real*8 :: g(nx,nyl,nzl), gt(nxl,ny,nzl)
          real*8 :: gt1(nxl,ny,nzl)
c
c  filtro in y
c  
        call nvtxStartRange('filtery_t',15)
        call trasponi_yx(g, gt, 1)

! ...   Filtering along y

	do k = 1,nzl          
          do i = 1, nxl
            do j = 1, ny
              filtroy(j) = gt(i,j,k)
            enddo
         
            call drfftf(ny, filtroy, wsavey) 
        
          do m1 = 2, ny - 1, 2

            m2 = m1 + 1
            m  = m1 / 2

            filtroy(m1) = filtroy(m1) * work4(m)
            filtroy(m2) = filtroy(m2) * work4(m)

          enddo

          filtroy(ny) = 0.0

          filtroy = filtroy / ny

          call drfftb(ny, filtroy, wsavey)

            do j = 1, ny
              gt1(i,j,k) = filtroy(j)
            enddo
          enddo
        enddo

        call trasponi_yx(g, gt1, -1)
       call nvtxEndRange()
       return
       end 
