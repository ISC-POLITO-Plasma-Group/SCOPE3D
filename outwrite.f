       subroutine outwrite(nunit,field)

       include 'par.inc'
       real*8 :: field(nx, nyl, nzl), field_rec(nx, nyl, nzl)
       integer nunit

       lx = 1 
       ly = 1
       lz = 1 
       isiz = nx*nyl*nzl


!       print*,time

!       IF( mpime == root ) then
!         write(nunit,201) time, nx/lx, ny/ly, nz/lz
!       END IF

       IF( mpime == root ) then
         write(nunit) time, nx/lx, ny/ly, nz/lz
       END IF
       
       DO ip = 1, nproc

!         write(6,*) ' DEBUG outwrite ',ip 
         it1 = ip

         IF( (ip-1) /= root ) THEN
           IF(mpime == (ip-1)) THEN
             CALL send_real( field, isiz, root, it1 )
           END IF
           IF ( mpime == root ) THEN
             CALL recv_real( field_rec, isiz, (ip-1), it1)
           END IF
         ELSE
           IF(mpime == root ) THEN
             field_rec = field
           END IF
         END IF

         CALL SYNCRONIZE()

!         IF(mpime == root) THEN
!           write(nunit,109) (((field_rec(ix,iy,iz),
!     &	        ix = 1, nx, lx),
!     &     iy = 1, nyl, ly),iz = 1, nzl, lz)  
!         END IF

!************UNFORMATTED************************
         IF(mpime == root) THEN
           write(nunit) (((field_rec(ix,iy,iz),
     &     ix = 1, nx, lx),
     &     iy = 1, nyl, ly),iz = 1, nzl, lz)  
         END IF
!***********************************************

       END DO


201    format(1x, 1e12.5, 1x, 3i4)
109    format(1x, 1(1x, 1e22.9))
       return
       end
