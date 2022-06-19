          subroutine filterz_t(g)
          use nvtx
          include 'par.inc'

          dimension filtroz(nz)
          real*8 :: g(nx,nyl,nzl), gt(nxl2,nyl,nz)
	  real*8 :: gt1(nxl2,nyl,nz)
c
c  filtro in z
c
        call nvtxStartRange('filterz_t',16)
	call trasponi_zx(g, gt, 1)

! ...   Filtering along z

	do j = 1,nyl
          do i = 1, nxl2
            do k = 1, nz
	      filtroz(k) = gt(i,j,k)
            enddo
          
            call drfftf(nz, filtroz, wsavez)

            do m1 = 2, nz - 1, 2
             
              m2 = m1 + 1
              m  = m1 / 2

              filtroz(m1) = filtroz(m1) * work5(m)
              filtroz(m2) = filtroz(m2) * work5(m)
 
            enddo
          
            filtroz(nz) = 0.0

            filtroz = filtroz / nz

            call drfftb(nz, filtroz, wsavez)


            do k = 1, nz
              gt1(i,j,k) = filtroz(k)
            enddo
	  enddo
	enddo

	call trasponi_zx(g, gt1, -1)
        call nvtxEndRange()
        return
        end 
