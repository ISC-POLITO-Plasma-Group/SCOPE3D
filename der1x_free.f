! --------------------------------------------------- !
! Derivate compatte 3x3 con derivata al bordo         !
!     calcolata con formula compatta punti interni    !
! Usata per derivare G+ e G-                          !
! --------------------------------------------------- !

        subroutine der1x_free(F0,aux_f1)

c  This routine calculates the first derivative in x direction (CFD 3 points)
        use nvtx
        include 'par.inc'   

        dimension aux_f1(nx),F0(nx),f1(nx)

!       costruzione termine noto
        call nvtxStartRange('der1x_free',17)
        do ix = 2, nx-1
          f1(ix) = aa_1_G(ix) * F0(ix-1) + bb_1_G(ix) * F0(ix) 
     &               + cc_1_G(ix) * F0(ix+1)  
        enddo


          f1(1) = aa_1_G(1) * F0(1) + bb_1_G(1) * F0(2) 
     &         + cc_1_G(1) * F0(3) + d_1_G(1) * F0(4)
          f1(nx) = aa_1_G(nx) * F0(nx) + bb_1_G(nx) * F0(nx-1) 
     &         + cc_1_G(nx) * F0(nx-2) + d_1_G(nx) * F0(nx-3)

        
!        write(*,*) 'OCCHIO ai COEFFICIENTI in DER1X'
!        write(*,*)  f1(1)
      

*******************************************

        do ix = 1,nx         
          aux_f1(ix) = f1(ix)
        enddo

c*******************************

!        do i=1,nx
!        write(63,*) aux_alfa_1(i),aux_gamma_1(i),aux_beta_1(i) 
!        enddo


        CALL  DGTTRS(TRANS,nx,NRHS,aux_alfa_1_G,aux_gamma_1_G,
     &        aux_beta_1_G,ww_1_G,ipv_1_G,aux_f1,LDB,INFO)

        

        if (info > 0 .or. info < 0) then
          write(*,*) 'Problemi soluzione, info:', info
        stop
        endif

        call nvtxEndRange()
	return
        end
