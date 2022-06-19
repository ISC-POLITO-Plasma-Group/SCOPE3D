      SUBROUTINE gpm_write(hp, hm,iunit)
        include 'par.inc'
        integer :: iunit
        real*8 :: hp(nx, nyl, nzl)
        real*8 :: hm(nx, nyl, nzl)

        call syncronize()
        call system_clock(it1)
        IF ( mpime == root ) then
          open(unit=iunit,status='unknown',
     &      file='DATA_1',form='unformatted')
        endif

        call write_out(hp, iunit)

        write(*,*) 'SCRITTO hp'

        call write_out(hm, iunit)

        write(*,*) 'SCRITTO hm'

        call write_out(phi, iunit)

        write(*,*) 'SCRITTO phi'

        call write_out(psi, iunit)

        write(*,*) 'SCRITTO psi'
      


        IF( mpime == root ) then
          write(iunit) oldtime, ioutt, ioutf

          write(*,*) 'SCRITTO oldtime, ioutt, ioutf'

          CLOSE(iunit)
        END IF

        call syncronize()
        call system_clock(it2,count_rate=iclk)
        write(6,100) dble(it2-it1)/dble(iclk)
 100    format('  restart file written in ',F8.2,' sec. ')
        RETURN
      END SUBROUTINE



      SUBROUTINE gpm_read(hp, hm, iunit)
        include 'par.inc'
        integer :: iunit
        real*8 :: hp(nx, nyl, nzl)
        real*8 :: hm(nx, nyl, nzl)
        call syncronize()
        call system_clock(it1)
        IF ( mpime == root ) then
          open(unit=iunit, status='old',
     &      file='DATA_1',form='unformatted')
        endif
        call read_in(hp, iunit)

        write(*,*) 'LETTO hp'

!        hp=hp*0.01

        call read_in(hm, iunit)

        write(*,*) 'LETTO hm'

!        hm=hm*0.01

        call read_in(phi, iunit)

        write(*,*) 'LETTO phi'

!        phi = phi*0.01

        call read_in(psi, iunit)

        write(*,*) 'LETTO psi'

!        psi = psi*0.01


        IF( mpime == root ) then
          read(iunit) oldtime, ioutt, ioutf

          write(*,*) 'LETTO oldtime, ioutt, ioutf'

          close(iunit)
        END IF
        call bcast_real(oldtime, 1, root)
        call bcast_integer(ioutt, 1, root)
        call bcast_integer(ioutf, 1, root)
        call syncronize()
        call system_clock(it2,count_rate=iclk)
        write(6,100) dble(it2-it1)/dble(iclk)
 100    format('  restart file read in ',F8.2,' sec. ')
        RETURN
      END SUBROUTINE


      SUBROUTINE write_out(a, iunit)
        include 'par.inc'
        integer :: ip, k
        real*8 :: a(nx, nyl, nzl) 
        real*8 :: a_rec(nx, nyl, nzl) 
        DO ip = 1, nproc
          IF( (ip-1) .NE. root ) THEN
            IF(mpime .EQ. (ip-1)) THEN
              CALL send_real( a, size(a), root, ip )
            END IF
            IF ( mpime .EQ. root) THEN
              CALL recv_real( a_rec, size(a_rec),(ip-1), ip)
            END IF
          ELSE
            IF(mpime .EQ. root) THEN
              a_rec = a
            END IF
          END IF
          CALL SYNCRONIZE()
          IF(mpime .EQ. root) THEN
            do k = 1, nzl
              write(iunit) (a_rec(i,1,k), i = 1, nx*nyl)  
            end do
          END IF
        END DO
        RETURN
      END SUBROUTINE

      SUBROUTINE read_in(a, iunit)
        include 'par.inc'
        integer :: ip, k
        real*8 :: a(nx, nyl, nzl)
        real*8 :: a_rec(nx, nyl, nzl)
        DO ip = 1, nproc
          IF(mpime == root) THEN
            do k = 1, nzl
              read(iunit) (a_rec(i,1,k), i = 1, nx*nyl)  
            end do
          END IF
          IF( (ip-1) /= root ) THEN
            IF ( mpime == root) THEN
              CALL send_real( a_rec, size(a_rec), (ip-1), ip)
            END IF
            IF(mpime == (ip-1)) THEN
              CALL recv_real( a, size(a), root, ip )
            END IF
          ELSE
            IF(mpime .EQ. root) THEN
              a = a_rec
            END IF
          END IF
          CALL SYNCRONIZE()
        END DO
        RETURN
      END SUBROUTINE

