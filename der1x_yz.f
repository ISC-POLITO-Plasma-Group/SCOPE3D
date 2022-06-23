! --------------------------------------------------- !
! Derivate compatte 3x3 con derivata al bordo         !
!     calcolata con formula compatta punti interni    !
! Usata per derivare G+ e G-                          !
! --------------------------------------------------- !

        subroutine der1x_yz(F0,f1,use_gpu_here,stream)

c  This routine calculates the first derivative in x direction (CFD 3 points)
        use nvtx
        use gpu_func
        use cudafor
        use cusparse
        use openacc
        include 'par.inc'   

        dimension F0(nx,nyl,nzl),f1(nx,nyl,nzl)

        ! gpu variables
        integer(4) ::  ierr
        logical :: use_gpu_here 
        integer(kind=cuda_stream_kind) :: stream
        integer :: strid
        

        call nvtxStartRange('der1x_y',17)

        ! Set the streamid equiv to the one created for cusparse
        call acc_set_cuda_stream(strid, stream)

!$acc parallel loop collapse(2) async(strid)       
        do iz=1,nzl
          do iy = 1, nyl
!$acc loop seq        
            do ix = 2, nx-1
               
               f1(ix,iy,iz) = aa_1_G(ix) * F0(ix-1,iy,iz) + bb_1_G(ix) 
     &                     * F0(ix,iy,iz) 
     &                     + cc_1_G(ix) * F0(ix+1,iy,iz) 
                
            enddo
        
   
          f1(1,iy,iz) = aa_1_G(1) * F0(1,iy,iz) + bb_1_G(1) 
     &         * F0(2,iy,iz) 
     &         + cc_1_G(1) * F0(3,iy,iz) + d_1_G(1) * F0(4,iy,iz)
          f1(nx,iy,iz) = aa_1_G(nx) * F0(nx,iy,iz) + bb_1_G(nx) 
     &         * F0(nx-1,iy,iz) 
     &         + cc_1_G(nx) * F0(nx-2,iy,iz) + d_1_G(nx) 
     &         * F0(nx-3,iy,iz)

        
!        write(*,*) 'OCCHIO ai COEFFICIENTI in DER1X'
!        write(*,*)  f1(1)
      
         enddo
        enddo
       
!        do i=1,nx
!        write(63,*) aux_alfa_1(i),aux_gamma_1(i),aux_beta_1(i)
!        enddo

        if (use_gpu_here) then 
#ifdef _OPENACC                

          ! Set the stream for the cusparse        
          ierr = ierr + cusparseSetStream(handle, stream)

!$acc host_data use_device(aux_alfa_2_G,aux_gamma_1_G,aux_beta_2_G,
!$acc& f1)

          if ( .not. allocated(pbuffer) ) then
            ! Compute buffersize for the cusparse and allocate
            ierr=ierr + cusparseDgtsv2_nopivot_bufferSize(handle,
     &      nx,nyl*nzl,aux_alfa_2_G,
     &      aux_gamma_1_G,aux_beta_2_G,f1,LDB,
     &      pBufferSizeInBytes)
            !write(*,*) 'buffersize=', pBufferSizeInBytes
            allocate(pBuffer(pBufferSizeInBytes))
          endif

          ! Do cusparse, nopivot algorithm
          ierr= ierr + cusparseDgtsv2_nopivot(handle, nx, nyl*nzl, 
     &    aux_alfa_2_G,
     &    aux_gamma_1_G,aux_beta_2_G, f1, LDB,
     &    pBuffer) 

!$acc end host_data

#endif
        else

!$acc update host(f1) 
          do iz=1,nzl              
              CALL  DGTTRS(TRANS,nx,nyl,aux_alfa_1_G,aux_gamma_1_G,
     &        aux_beta_1_G,ww_1_G,ipv_1_G,f1(:,:,iz),LDB,INFO)
          enddo

          if (info > 0 .or. info < 0) then
              write(*,*) 'Problemi soluzione, info:', info
              stop
          endif

        endif

        call nvtxEndRange()
      return

      end
