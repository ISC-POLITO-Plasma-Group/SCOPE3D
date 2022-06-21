! --------------------------------------------------- !
! Derivate compatte 3x3 con derivata al bordo         !
!     calcolata con formula compatta punti interni    !
! Usata per derivare G+ e G-                          !
! --------------------------------------------------- !

        subroutine der1x_y(F0,aux_f1,use_gpu_here)

c  This routine calculates the first derivative in x direction (CFD 3 points)
        use nvtx
        use gpu_func
        use cudafor
        use cusparse
        include 'par.inc'   

        dimension aux_f1(nx,nyl),F0(nx,nyl),f1(nx,nyl)
        integer(4)::  ierr
        integer(8):: pBufferSizeInBytes
        character(1), device, allocatable:: pBuffer(:)       
!        logical use_gpu
        logical use_gpu_here 

!       use_gpu_here  = .false.
!       if(present(use_gpu)) use_gpu_here = use_gpu 
!       costruzione termine noto
        

        call nvtxStartRange('der1x_free',17)
        do iy = 1, nyl
            do ix = 2, nx-1
               f1(ix,iy) = aa_1_G(ix) * F0(ix-1,iy) + bb_1_G(ix) 
     &                     * F0(ix,iy) 
     &                     + cc_1_G(ix) * F0(ix+1,iy)  
            enddo
        

   
          f1(1,iy) = aa_1_G(1) * F0(1,iy) + bb_1_G(1) * F0(2,iy) 
     &         + cc_1_G(1) * F0(3,iy) + d_1_G(1) * F0(4,iy)
          f1(nx,iy) = aa_1_G(nx) * F0(nx,iy) + bb_1_G(nx) * F0(nx-1,iy) 
     &         + cc_1_G(nx) * F0(nx-2,iy) + d_1_G(nx) * F0(nx-3,iy)

        
!        write(*,*) 'OCCHIO ai COEFFICIENTI in DER1X'
!        write(*,*)  f1(1)
      

*******************************************
       
           do ix = 1,nx         
              aux_f1(ix,iy) = f1(ix,iy)
           enddo
       enddo
        if (use_gpu_here) then 
!$acc enter data copyin(aux_f1) 
c*******************************

!        do i=1,nx
!        write(63,*) aux_alfa_1(i),aux_gamma_1(i),aux_beta_1(i) 
!        enddo

!$acc host_data use_device(aux_alfa_2_G,aux_gamma_1_G,aux_beta_2_G,
!$acc& aux_f1)

       ierr=ierr + cusparseDgtsv2_bufferSize(handle,nx,nyl,
     & aux_alfa_2_G,
     & aux_gamma_1_G,aux_beta_2_G,aux_f1,LDB,
     & pBufferSizeInBytes)

       write(*,*) 'buffersize=', pBufferSizeInBytes
       allocate(pBuffer(pBufferSizeInBytes))

       ierr= ierr + cusparseDgtsv2(handle, nx, nyl, aux_alfa_2_G,
     & aux_gamma_1_G,aux_beta_2_G, aux_f1, LDB,
     & pBuffer) 

!$acc end host_data
!$acc exit data copyout(aux_f1)
       else 
       CALL  DGTTRS(TRANS,nx,nyl,aux_alfa_1_G,aux_gamma_1_G,
     &        aux_beta_1_G,ww_1_G,ipv_1_G,aux_f1,LDB,INFO)

        

        if (info > 0 .or. info < 0) then
          write(*,*) 'Problemi soluzione, info:', info
        stop
        endif
      endif 

!!!!!!!!!!!    CPU END !!!!!!!!!!!

        call nvtxEndRange()

	return
        end
