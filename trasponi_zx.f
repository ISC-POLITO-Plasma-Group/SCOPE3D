       subroutine trasponi_zx(a, at, idir)
        use nvtx       
        include 'par.inc'
        
        integer :: idir
        real*8 :: a(nx,nyl,nzl), at(nxl2,nyl,nz)
        real*8 :: buffer(nyl*nzl*nxl2)

!        write(6,*) ' DEBUG trasponi_zx ', idir
        call nvtxStartRange('trasponi_zx',19)
        if( idir >= 0 ) then

          ibufsize = nzl * nyl * nxl2 ! size of the block
          DO IP = 1, NPROCZ
            ISOUR = MOD(MPIME-IP+NPROCZ,NPROCZ)  !  ISOUR -> MPIME - 1 ... NPROC
            IDEST = MOD(MPIME+IP      ,NPROCZ)  !  IDEST -> MPIME + 1 ... NPROC
            ibuf = 1                           !  ISOUR, IDEST, MPIME start from 0
            inx_start = idest * nxl2 + 1
            inx_end   = inx_start + nxl2 - 1
            do i = inx_start , inx_end
              do j = 1, nyl 
                do k = 1 , nzl
                  buffer(ibuf) = a(i,j,k)
                  ibuf = ibuf + 1
                end do
              end do
            end do

	    idest1 = idest+nprocz*INT(MPIME/nprocz)
            isour1 = isour+nprocz*INT(MPIME/nprocz)

            call sendrecv(buffer, ibufsize, idest1, isour1,ip)
            ibuf = 1
            inz_start = isour * nzl + 1
            inz_end   = inz_start + nzl - 1
            do i = 1 , nxl2
              do j = 1,nyl 
                do k = inz_start , inz_end
                  at(i,j,k) = buffer(ibuf)
                  ibuf = ibuf + 1
                end do
              end do
            end do
          end do
	
        else

          ibufsize = nzl * nyl * nxl2 ! size of the block
          DO IP = 1, NPROCZ
            ISOUR = MOD(MPIME-IP+NPROCZ,NPROCZ)  !  ISOUR -> MPIME - 1 ... NPROC
            IDEST = MOD(MPIME+IP      ,NPROCZ)  !  IDEST -> MPIME + 1 ... NPROC
            ibuf = 1
            inz_start = isour * nzl + 1
            inz_end   = inz_start + nzl - 1
            do i = 1 , nxl2
              do j = 1,nyl 
                do k = inz_start , inz_end
                  buffer(ibuf) = at(i,j,k)
                  ibuf = ibuf + 1
                end do
              end do
            end do

	    idest1 = idest+nprocz*INT(MPIME/nprocz)
            isour1 = isour+nprocz*INT(MPIME/nprocz)

            call sendrecv(buffer, ibufsize, isour1, idest1, ip)
            ibuf = 1                           !  ISOUR, IDEST, MPIME start from 0
            inx_start = idest * nxl2 + 1
            inx_end   = inx_start + nxl2 - 1
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


