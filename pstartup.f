         subroutine pstartup
         use openacc
         include 'par.inc'
         integer globalindex_y,globalindex_z
         integer localdim

c-----------------------------------------------------------------
c   nx+1=number of grid points in the x-direction
c   ny+1=number of grid points in the y-direction 
c
c   nvx+1=number of grid points in the vx-direction
c   nvy+1=number of grid points in the vy-direction
c
c----------------------------------------------------------------
        
        
        call parallel_startup(nproc, mpime, root, group)
         ndev= acc_get_num_devices(acc_device_nvidia)
         idev=mod(mpime,ndev)
         call acc_set_device_num(idev,acc_device_nvidia)
         write(*,*) 'hello from', idev, 'of', ndev
         ridev= acc_get_device_num(acc_device_nvidia)
         write(*,*) 'hello from', ridev, 'of', ndev

         nprocz = NINT(SQRT(REAL(nproc)))
         nprocz = 1
!         nprocz = 16
         nprocy = nproc / nprocz
         write(*,*) 'NPY, NPZ',nprocy, nprocz

        if( nproc .gt. nprocx ) then
          write(6,*) '*** too many processors ***'
          write(6,*) '***   increase nprocx   ***'
          call hangup
          stop
        end if
        if( MOD( nz, nprocz) /= 0  ) then
          write(6,*) '*** nproc not a divisor of nz ***'
          write(6,*) '*** nz = ', nz, ' nproc = ', nproc
          call hangup
          stop
        end if
        if( MOD( ny, nprocy) /= 0  ) then
          write(6,*) '*** nproc not a divisor of ny ***'
          write(6,*) '*** ny = ', ny, ' nproc = ', nproc
          call hangup
          stop
        end if
        nyl = ny / nprocy
        nxl = nx / nprocy
        nyg = ny
        iylg = globalindex_y(1,nyg,nprocy,nprocz,mpime)

        write(*,*) 'IYLG', iylg

        nzl = nz / nprocz
        nxl2 = nx / nprocz
        nzg = nz
        izlg = globalindex_z(1,nzg,nprocy,nprocz,mpime)
        

        write(*, 1000) mpime, nproc, nprocy,nprocz
        write(*, 1001) mpime, nzl, nyl, izlg, iylg, nzg, nyg

1000   FORMAT(' mpime, nproc ',6I4) 
1001   FORMAT(' mpime, nzl, nyl, izlg, iylg, nzg, nyg',7I4) 

       end

