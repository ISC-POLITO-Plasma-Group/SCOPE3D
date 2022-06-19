

! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

      subroutine sendrecv(array,size,dest,sour,ip)
        implicit none
        INCLUDE 'mpif.h'
        REAL*8 array(*)
        integer size, dest, sour, ip
        integer istatus(MPI_STATUS_SIZE), ierr
! ... 
        call mpi_sendrecv_replace(array, size, MPI_DOUBLE_PRECISION, 
     &    dest,ip,sour,ip,MPI_COMM_WORLD, ISTATUS, IERR)

        RETURN
      END SUBROUTINE

      subroutine fddsendrecv(fddin,fddout,dest,sour,ip)
        implicit none
        INCLUDE 'mpif.h'
        REAL*8 fddin(2), fddout(2)
        integer dest, sour, ip, ierr
        integer istatus(MPI_STATUS_SIZE)
        fddout(1) = fddin(1) 
        fddout(2) = fddin(2) 
        call mpi_sendrecv_replace(fddout, 2, MPI_DOUBLE_PRECISION,
     &     dest, ip, sour, ip, MPI_COMM_WORLD, ISTATUS, IERR)
        RETURN
      END SUBROUTINE


      subroutine send_real(array, size, dest, ip)
        implicit none
        INCLUDE 'mpif.h'
        REAL*8 array(*)
        integer size, dest, ip
        integer ierr
! ...
        call mpi_send(array, size, MPI_DOUBLE_PRECISION,
     &    dest, ip, MPI_COMM_WORLD, IERR)
        RETURN
      END SUBROUTINE

      subroutine recv_real(array, size, sour, ip)
        implicit none
        INCLUDE 'mpif.h'
        REAL*8 array(*)
        integer size, sour, ip, nrec, sourid
        integer istatus(MPI_STATUS_SIZE), ierr
! ... 
        IF( sour >= 0 ) THEN
          sourid = sour
        ELSE
          sourid = MPI_ANY_SOURCE
        END IF
        call mpi_recv(array, size, MPI_DOUBLE_PRECISION, 
     &    sourid, ip, MPI_COMM_WORLD, ISTATUS, IERR)
        CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrec, ierr)
! ...
        RETURN
      END SUBROUTINE

