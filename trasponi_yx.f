       subroutine trasponi_yx(a, at, idir)

        use nvtx
        include 'par.inc'
        integer :: idir
        real*8 :: a(nx,nyl,nzl), at(nxl,ny,nzl)
        real*8 :: buffer(nyl*nzl*nxl)

!        write(6,*) ' DEBUG trasponi_yx ', idir
        call nvtxStartRange('trasponiyx',19)
        if( idir >= 0 ) then

          ibufsize = nzl * nyl * nxl ! size of the block
          write(*,*),'ibufsize=',ibufsize

          DO IP = 1, NPROCY
            ISOUR = MOD(INT(MPIME/NPROCZ)-IP+NPROCY,NPROCY)  !  ISOUR -> MPIME - 1 ... NPROC
            IDEST = MOD(INT(MPIME/NPROCZ)+IP      ,NPROCY)  !  IDEST -> MPIME + 1 ... NPROC
            ibuf = 1                           !  ISOUR, IDEST, MPIME start from 0
            inx_start = idest * nxl + 1
            inx_end   = inx_start + nxl - 1
            do i = inx_start , inx_end
              do j = 1, nyl 
                do k = 1 , nzl
                  buffer(ibuf) = a(i,j,k)
                  ibuf = ibuf + 1
                end do
              end do
            end do

	    idest1 = idest*NPROCZ+MOD(MPIME,NPROCZ) 
!+ nprocy*INT(MPIME/nprocy)
            isour1 = isour*NPROCZ+MOD(MPIME,NPROCZ) 
!+ nprocy*INT(MPIME/nprocy)

            call sendrecv(buffer, ibufsize, idest1, isour1,ip)
            ibuf = 1
            iny_start = isour * nyl + 1
            iny_end   = iny_start + nyl - 1
            do i = 1 , nxl
              do j = iny_start , iny_end
                do k = 1,nzl
                  at(i,j,k) = buffer(ibuf)
                  ibuf = ibuf + 1
                end do
              end do
            end do
          end do
	
        else

          ibufsize = nzl * nyl * nxl ! size of the block
          DO IP = 1, NPROCY
            ISOUR = MOD(INT(MPIME/NPROCZ)-IP+NPROCY,NPROCY)  !  ISOUR -> MPIME - 1 ... NPROC
            IDEST = MOD(INT(MPIME/NPROCZ)+IP      ,NPROCY)  !  IDEST -> MPIME + 1 ... NPROC
            ibuf = 1
            iny_start = isour * nyl + 1
            iny_end   = iny_start + nyl - 1
            do i = 1 , nxl
              do j = iny_start , iny_end
                do k = 1, nzl 
                  buffer(ibuf) = at(i,j,k)
                  ibuf = ibuf + 1
                end do
              end do
            end do

	    idest1 = idest*NPROCZ+MOD(MPIME,NPROCZ)
!+ nprocy*INT(MPIME/nprocy)
            isour1 = isour*NPROCZ+MOD(MPIME,NPROCZ)
!+ nprocy*INT(MPIME/nprocy)

            call sendrecv(buffer, ibufsize, isour1, idest1, ip)
            ibuf = 1                           !  ISOUR, IDEST, MPIME start from 0
            inx_start = idest * nxl + 1
            inx_end   = inx_start + nxl - 1
            do i = inx_start , inx_end
              do j = 1,nyl 
                do k = 1 , nzl
                  a(i,j,k) = buffer(ibuf)
                  ibuf = ibuf + 1
                end do
              end do
            end do
          end do

        end if
       call nvtxEndRange()
      return
      end


