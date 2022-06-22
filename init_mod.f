        subroutine init_mod
        use gpu_func
        include 'par.inc'

        dimension aa(nx-1), bb(nx), cc(nx-1), ww(nx-2)
        dimension aa_my(nx-1,nyl), bb_my(nx,nyl), cc_my(nx-1,nyl)
        dimension aux_aa(nx-1), aux_bb(nx), aux_cc(nx-1)
        dimension aa_my_1(nx1,nyl), bb_my_1(nx,nyl), cc_my_1(nx1,nyl)
        dimension aux_aa_1(nx1), aux_bb_1(nx), aux_cc_1(nx1)
        dimension aux_alfa_2(nx1), aux_beta_2(nx1), aux_gamma_2(nx)
        dimension ww_2(nx-2)

        tlength = dabs(t0) * 2.0d0
        xnx      = t0/gridbeta - gridalfa * dtanh(t0/gridgamma)
        x0       = - xnx               
        xlength  = 2.0d0 * xnx 
        dtx       = tlength / (nx-1)
        do i = 1, nx
         x(i) = t0+(i-1)*dtx
         x(i) = (x(i))/gridbeta -  gridalfa 
     &        * dtanh((x(i))/gridgamma)         
        enddo
        write(*,*) 'grid par:',tlength,xnx,dtx

!************** EQUISPACED GRID ************
!	t0 = -pi
!	dx_eq = dabs(t0) * 2.0d0 / (nx-1)
!	do ix = 1, nx
!	  x(ix) = t0 + (ix - 1) * dx_eq
!	enddo	
!*******************************************
!
        do i = 1,nx-1                                       
           dx(i) = x(i+1) - x(i)  
        enddo
        dx(nx)   = dx(nx-1)
!
        dy = yl / float(ny)
        dz = zl / float(nz)
!	
        do i = 1,nx
           dxy(i) = dx(i) * dy
           dxz(i) = dx(i) * dz
           dxyz(i) = dx(i)* dy * dz
        enddo
        dyz = dz * dy
!
        de2 = de * de
        rhos2 = rhos * rhos
!
        do iy = 1, ny
           y(iy) = -yl / 2.0 + (iy - 1) * dy
        enddo
!
        do iz = 1, nz
           z(iz) = -zl / 2.0 + (iz - 1) * dz
        enddo

        do iy = 1, ny/2+1
           yk(iy) = 2.*0d0*pi*(iy-1)/yl
        enddo

        do iz = 1, nz/2+1
           zk(iz) = 2.0d0*pi*(iz-1)/zl
        enddo
!
        do iy = 2, ny-1,2
           ky(iy) = iy / 2
           ky(iy+1) = ky(iy)
        enddo

        ky(1) = 0
        ky(ny)= 0

        do iz = 2, nz-1,2
           kz(iz) = iz / 2
           kz(iz+1) = kz(iz)
        enddo

        kz(1) = 0
        kz(nz)= 0


        call drffti(ny, wsavey)
        call drffti(nz, wsavez)


c   filter parameters
!
       zeta = (3.d0 - 2.d0*omega)/10.d0
       alpha = (2.d0 + 3.d0*omega)/4.d0
       beta1 = (6.d0 + 7.d0*omega)/8.d0
       gamma = (6.d0 + omega)/20.d0
       delta = (2.d0 - 3.d0*omega)/40.d0

       do kk = 1, ny/2
!         w = 2.d0 * pi * grady * dfloat(kk) / dfloat(ny)
          w = 2.d0 * pi * dfloat(kk) / dfloat(ny)
          work1 = alpha + beta1*dcos(w) +
     &          gamma*dcos(2.d0*w) + delta*dcos(3.d0*w)
          work2 = 1.d0+2.d0*omega*dcos(w)+2.d0*zeta*dcos(2.d0*w)
          work4(kk) = work1/work2
       enddo

       do kk = 1, nz/2
!         w = 2.d0 * pi * gradz * dfloat(kk) / dfloat(nz)
          w = 2.d0 * pi * dfloat(kk) / dfloat(nz)
          work1 = alpha + beta1*dcos(w) +
     &          gamma*dcos(2.d0*w) + delta*dcos(3.d0*w)
          work2 = 1.d0+2.d0*omega*dcos(w)+2.d0*zeta*dcos(2.d0*w)
          work5(kk) = work1/work2
       enddo


***** Setup solver tridiagonale LU  FILTRO in X hp, hm, psi, phi ***

       info_t = 0
       TRANS = 'N'
       nrhs  = 1
       ndb_t   = nx

!    R.H.S. coefficients of the tridiag. filter

        a_filt =   (5.0d0 + 6.0d0 * omegax) / 8.0d0
        b_filt =   (1.0d0 + 2.0d0 * omegax) / 4.0d0
        c_filt = - (1.0d0 - 2.0d0 * omegax) / 16.0d0

!    R.H.S boundary conditions coefficients

        a11_filt =  1.0d0-(1.0d0 - alfa1) / 16.0d0
        a12_filt =  alfa1 + (1.0d0 - alfa1) / 4.0d0
        a13_filt = -(6.0d0 / 16.0d0) * (1.0d0 - alfa1)
        a14_filt = (1.0d0 - alfa1) / 4.0d0
        a15_filt = -(1.0d0 - alfa1) / 16.0d0

        a21_filt =   1.0d0/ 16.0d0
        a22_filt =   3.0d0 / 4.0d0
        a23_filt =   3.0d0 / 8.0d0
        a24_filt = - 1.0d0 / 4.0d0
        a25_filt =   1.0d0/ 16.0d0

!    Diagonale filt_d, sotto diagonale filt_m e super diagonale filt_p
!    (filt_2 e' un working array)

        filt_d_t = 1.0d0

        filt_p_t = omegax
        filt_m_t = omegax

        filt_m_t(1) = 0.0d0
        filt_p_t(1)    = alfa1
        filt_p_t(2)    = 0.0d0
        filt_m_t(nx-1) = alfa1
        filt_m_t(nx-2) = 0.0d0
        filt_p_t(nx-1) = 0.0d0

        filt_2_t = 0.0d0

!     Fattorizzazione LU

         CALL DGTTRF(nx, filt_m_t, filt_d_t,
     &      filt_p_t, filt_2_t, jpv_t, info_t)

         if (info_t > 0 .or. info_t < 0) then
          write(*,*) 'Problemi fattorizzazione LU, info:', info
          stop
         endif


***************************************************************


!
*****   Setup solver Poisson tridiagonale LU  ***

       info = 0
       TRANS = 'N'
       nrhs  = 1
       ldb   = nx

*************************
! calcolo coefficienti interni

	write(*,*) '!!!!!!!!!IYLG!!!!!!!!!',iylg,mpime
        do ix = 2,nx-1

          dp1 = dx(ix)
          dm1 = dx(ix-1)
          dp2 = dp1 * dp1
          dm2 = dm1 * dm1
          dp3 = dp1 * dp2
          dm3 = dm1 * dm2

          temp1 = (dp3 + dm3) - (dp1 - dm1)
     &     * (5.0d0 * dp3 + 6.0d0 * dp2 * dm1 - dm3)
     &     / (2.0d0 * dp1 + dm1)
          temp2 = (5.0d0 * dm3 + 6.0d0 * dm2 * dp1 - dp3) +
     &     (dp1 + 2.0d0 * dm1) *
     &     (5.0d0 * dp3 + 6.0d0 * dp2 * dm1 - dm3) /
     &    (dm1 + 2.0d0 * dp1)
        
          alfa(ix) = temp1 / temp2
          beta(ix) = (alfa(ix) * (dp1 + 2.0d0 * dm1) +
     &    (dp1 - dm1)) / (dm1 + 2.0d0 * dp1)
          cc(ix) = 2.0d0*(1.0d0 + alfa(ix) + beta(ix)) /
     &    (dp1 + dm1) / dp1

          aa(ix-1) =   cc(ix) * dp1 / dm1
          bb(ix)   = - aa(ix-1) - cc(ix)

!! ATTENZIONE! nyl deve essere maggiore di 3 

        do my = 2,nyl-2,2
          myg = my + iylg - 1
          hy = 2.0d0 * pi * myg/yl
          hy1 = 2.0d0 * pi * (myg+2)/yl

          aa_my(ix-1,my) = aa(ix-1) - hy*hy*alfa(ix)/4.0d0
          aa_my(ix-1,my+1) = aa_my(ix-1,my)
          bb_my(ix,my) = bb(ix) - hy*hy/4.0d0
          bb_my(ix,my+1) = bb_my(ix,my)
          cc_my(ix,my) = cc(ix) - hy*hy*beta(ix)/4.0d0
          cc_my(ix,my+1) = cc_my(ix,my)
          aa_my(ix-1,my+2) = aa(ix-1) - hy1*hy1*alfa(ix)/4.0d0
          bb_my(ix,my+2) = bb(ix) - hy1*hy1/4.0d0
          cc_my(ix,my+2) = cc(ix) - hy1*hy1*beta(ix)/4.0d0
        enddo

        if(INT(mpime/nprocz) == 0) then
         aa_my(ix-1,1) = aa(ix-1)
         bb_my(ix,1) = bb(ix)
         cc_my(ix,1) = cc(ix)
        endif
        if(INT(mpime/nprocz).ne.0) then
!          hy2 = 2.0d0 * pi * (myg-2)/yl
          hy2 = hy-2.0d0 * pi * (nyl-2)/yl
          aa_my(ix-1,1) = aa(ix-1) - hy2*hy2*alfa(ix)/4.0d0
          bb_my(ix,1) = bb(ix) - hy2*hy2/4.0d0
          cc_my(ix,1) = cc(ix) - hy2*hy2*beta(ix)/4.0d0
        endif

        enddo


        alfa(1) = 0.0d0
        beta(1) = 0.0d0
! calcolo coefficienti al bordo (k = 0)

********** Von Neumann bc (onde ky,kz)*********

!        if (kbc_p.eq.1) then
!          bb(1) = 1.0d0
!          cc(1) = -bb(1)

!          bb(nx)   = 1.0d0
!          aa(nx-1) = -bb(nx)
!        endif

***********************************

********* dirichelet bc ***********

!        if (kbc_p.eq.0) then
          bb(1) = 1.0d0
          cc(1) = 0.0d0

          bb(nx)   = 1.0d0
          aa(nx-1) = 0.0d0
!        endif

***********************************
***********************************

        do my = 2,nyl-2,2

          bb_my(1,my) = bb(1)
          bb_my(1,my+1) = bb_my(1,my)
          cc_my(1,my) = cc(1)
          cc_my(1,my+1) = cc_my(1,my)
          aa_my(nx-1,my) = aa(nx-1)
          aa_my(nx-1,my+1) = aa_my(nx-1,my)
          bb_my(nx,my) = bb(nx)
          bb_my(nx,my+1) = bb_my(nx,my)

	  bb_my(1,my+2) = bb(1)
	  cc_my(1,my+2) = cc(1)
	  aa_my(nx-1,my+2) = aa(nx-1)
	  bb_my(nx,my+2) = bb(nx) 

        enddo
! calcolo coefficienti al bordo per (k.eq.0)

        bb_my(1,1) = 1.0d0
        cc_my(1,1) = 0.0d0

        aa_my(nx-1,1) = 0.0d0
        bb_my(nx,1) = 1.0d0


!     Fattorizzazione LU


        do my = 1,nyl

         do ix = 1,nx-1
           aux_aa(ix) = aa_my(ix,my)
           aux_bb(ix) = bb_my(ix,my)
           aux_cc(ix) = cc_my(ix,my)
        enddo

         aux_bb(nx) = bb_my(nx,my)

         CALL DGTTRF(nx, aux_aa, aux_bb, aux_cc, ww, ipv, info)

         do ix = 1,nx-1
          fact_aa(ix,my) = aux_aa(ix)
          fact_bb(ix,my) = aux_bb(ix)
          fact_cc(ix,my) = aux_cc(ix)
          ipv_my(ix,my) = ipv(ix)
         enddo

         do ix = 1,nx-2
          fact_ww(ix,my) = ww(ix)
         enddo

         fact_bb(nx,my) = aux_bb(nx)
         ipv_my(nx,my) = ipv(nx)


         if (info > 0 .or. info < 0) then
         write(*,*) 'Problemi LU - init_poiss - , info:', info,mpime,my
            stop
         endif
        enddo
	call SYNCRONIZE()
!
!******Setup solver HElmoltz tridiagonale LU*********

       info = 0
       TRANS = 'N'
       nrhs  = 1
       ldb   = nx

        if (de.ne.0.0) then  

        do ix = 2,nx-1
        do my = 2,nyl-2,2
	  myg = my + iylg - 1
          hy = 2.0d0 * pi * myg/yl
          hy1 = 2.0d0 * pi * (myg+2)/yl

          aa_my_1(ix-1,my) = aa(ix-1) 
     &     - (hy*hy/4.0d0+1.0d0/de2)*alfa(ix)
          aa_my_1(ix-1,my+1) = aa_my_1(ix-1,my)
          bb_my_1(ix,my) = bb(ix) -(hy*hy/4.0d0+1.0d0/de2)
          bb_my_1(ix,my+1) = bb_my_1(ix,my)
          cc_my_1(ix,my) = cc(ix) - (hy*hy/4.0d0+1.0d0/de2)*beta(ix)
          cc_my_1(ix,my+1) = cc_my_1(ix,my)
          aa_my_1(ix-1,my+2) = aa(ix-1)
     &          -(hy1*hy1/4.0d0+1.0d0/de2)*alfa(ix)
          bb_my_1(ix,my+2) = bb(ix) 
     &          - (hy1*hy1/4.0d0+1.0d0/de2)
          cc_my_1(ix,my+2) = cc(ix)
     &          -(hy1*hy1/4.0d0+1.0d0/de2)*beta(ix)
        enddo

        if(INT(mpime/nprocz) == 0) then
         aa_my_1(ix-1,1) = aa(ix-1)-(1.0d0/de2)*alfa(ix)
         bb_my_1(ix,1) = bb(ix)-(1.0d0/de2)
         cc_my_1(ix,1) = cc(ix)-(1.0d0/de2)*beta(ix)
        endif
        if(INT(mpime/nprocz).ne.0) then
!          hy2 = 2.0d0 * pi * (myg-2)/yl
          hy2 = hy-2.0d0 * pi * (nyl-2)/yl
          aa_my_1(ix-1,1) = aa(ix-1)
     &      - (hy2*hy2/4.0d0+1.0d0/de2)*alfa(ix)
          bb_my_1(ix,1) = bb(ix) - (hy2*hy2/4.0d0+1.0d0/de2)
          cc_my_1(ix,1) = cc(ix)
     &      - (hy2*hy2/4.0d0+1.0d0/de2)*beta(ix)
        endif

	enddo

        alfa(1) = 0.0d0
        beta(1) = 0.0d0

! calcolo coefficienti al bordo (k = 0)

********** Von Neumann bc (onde ky,kz)*********

!        if (kbc_p.eq.1) then
!          bb(1) = 1.0d0
!          cc(1) = -bb(1)

!          bb(nx)   = 1.0d0
!          aa(nx-1) = -bb(nx)
!        endif

***********************************

********* dirichelet bc ***********

!        if (kbc_p.eq.0) then
          bb(1) = 1.0d0
          cc(1) = 0.0d0

          bb(nx)   = 1.0d0
          aa(nx-1) = 0.0d0
!        endif

***********************************

        do my = 2,nyl-2,2

          bb_my_1(1,my) = bb(1)
          bb_my_1(1,my+1) = bb_my_1(1,my)
          cc_my_1(1,my) = cc(1)
          cc_my_1(1,my+1) = cc_my_1(1,my)
          aa_my_1(nx-1,my) = aa(nx-1)
          aa_my_1(nx-1,my+1) = aa_my_1(nx-1,my)
          bb_my_1(nx,my) = bb(nx)
          bb_my_1(nx,my+1) = bb_my_1(nx,my)

	  bb_my_1(1,my+2) = bb(1)
	  cc_my_1(1,my+2) = cc(1)
	  aa_my_1(nx-1,my+2) = aa(nx-1)
	  bb_my_1(nx,my+2) = bb(nx) 

        enddo

! calcolo coefficienti al bordo per (k.eq.0)

        bb_my_1(1,1) = 1.0d0
        cc_my_1(1,1) = 0.0d0

        aa_my_1(nx-1,1) = 0.0d0
        bb_my_1(nx,1) = 1.0d0

!     Fattorizzazione LU

        do my = 1,nyl

         do ix = 1,nx-1
           aux_aa_1(ix) = aa_my_1(ix,my)
           aux_bb_1(ix) = bb_my_1(ix,my)
           aux_cc_1(ix) = cc_my_1(ix,my)
        enddo

         aux_bb_1(nx) = bb_my_1(nx,my)

       CALL DGTTRF(nx, aux_aa_1, aux_bb_1, aux_cc_1, ww, ipv, info)

         do ix = 1,nx-1
          fact_aa_1(ix,my) = aux_aa_1(ix)
          fact_bb_1(ix,my) = aux_bb_1(ix)
          fact_cc_1(ix,my) = aux_cc_1(ix)
          ipv_my_1(ix,my) = ipv(ix)
         enddo

         do ix = 1,nx-2
          fact_ww_1(ix,my) = ww(ix)
         enddo

         fact_bb_1(nx,my) = aux_bb_1(nx)
         ipv_my_1(nx,my) = ipv(nx)


         if (info > 0 .or. info < 0) then
         write(*,*) 'Problemi LU - init_helm - , info:', info,mpime,my
            stop
         endif
        enddo
        
        endif

	call SYNCRONIZE()

!***************************************************************

***** Setup derivata x tridiagonale LU   ***

!
!     calcolo coefficienti interni 
       
      do ix = 2,nx-1
         
         dp1 = dx(ix)
         dm1 = dx(ix-1)
         dp2 = dp1 * dp1 
         dm2 = dm1 * dm1
         dp3 = dp1 * dp2
         dm3 = dm1 * dm2
         
         temp1 = - dp1 * dm2 * (dp1 + dm1 / 2.0d0) / 3.0d0 

         temp2 = ((dm1/2.0d0 + dp1/6.0d0) * (dp2 - dm2/6.0d0) - 
     &     dp1 * (2.0d0 * dm1 / 3.0d0 + dp1 / 4.0d0) * 
     &        (dm1 / 3.0d0 + dp1)) * dp2

         temp3 = dp2 * (dp1 / 6.0d0 + dm1 / 2.0d0)
         temp4 = (dm1 / 3.0d0 + dp1) * (dm2 / 2.0d0)

         cc_1_G(ix) = temp1 / temp2
         aa_1_G(ix) = (cc_1_G(ix) * temp3 - dp1 * dm1) / temp4
         beta_1_G(ix) = ((dm1 + 0.5d0 * dp1) * dp1 * cc_1_G(ix) - 
     &        0.5d0 * dm2 * aa_1_G(ix) -  dm1) / (dp1 + dm1)
         alfa_1_G(ix-1) = - dm1 * aa_1_G(ix) + dp1 * cc_1_G(ix) - 
     &        beta_1_G(ix) - 1.0d0
         
         bb_1_G(ix) = - aa_1_G(ix) - cc_1_G(ix)
          
      enddo

      do ix = 2,nx-1
         d_1_G(ix) = 0.0d0
      enddo

!!!   write(*,*) ' Derivata al bordo IMPLICITA '

      beta_1_G(1) = 3.0d0
      aa_1_G(1) = -17.0d0 / 6.0d0 / dx(1)
      bb_1_G(1) = 3.0d0 / 2.0d0 / dx(1)
      cc_1_G(1) = 3.0d0 / 2.0d0 / dx(1)
      d_1_G(1) = -1.0d0 / 6.0d0 / dx(1)


      alfa_1_G(nx-1) = 3.0d0
      aa_1_G(nx) = 17.0d0 / 6.0d0 / dx(nx-1)
      bb_1_G(nx) = -3.0d0 / 2.0d0 / dx(nx-1)
      cc_1_G(nx) = -3.0d0 / 2.0d0 / dx(nx-1)
      d_1_G(nx) = 1.0d0 / 6.0d0/ dx(nx-1)
      
      write(*,*) ' *************************** '
      write(*,*) ' Derivata al bordo ESPLICITA '
      write(*,*) ' *************************** '
      
      beta_1_G(1) = 0.0d0
      aa_1_G(1) = - 3.0d0 / 2.0d0 / dx(1)
      bb_1_G(1) = 2.0d0 / dx(1)
      cc_1_G(1) = - 1.0d0 / 2.0d0 / dx(1)
      d_1_G(1)  = 0.0d0
      
      alfa_1_G(nx-1) = 0.0d0
      aa_1_G(nx)     = 3.0d0 / 2.0d0 / dx(nx-1)
      bb_1_G(nx)     = - 2.0d0 / dx(nx-1)
      cc_1_G(nx)     = 1.0d0 / 2.0d0 / dx(nx-1)
      d_1_G(nx)      = 0.0d0
        

      write(*,*) ' *************************** '
      write(*,*) ' Derivata al bordo ESPLICITA neg '
      write(*,*) ' *************************** '
      
      beta_1_G(1) = 0.0d0
      aa_1_G(1) = -(2.0d0*dx(1)+dx(2))/
     &     (dx(1)*(dx(1)+dx(2)))
      bb_1_G(1) = (dx(1)+dx(2))/
     &     (dx(1)*dx(2))
      cc_1_G(1) = -dx(1)/
     &     (dx(2)*(dx(1)+dx(2)))
      d_1_G(1)  = 0.0d0
      
      alfa_1_G(nx-1) = 0.0d0
      aa_1_G(nx)     = (2.0d0*dx(nx-1)+dx(nx-2))/
     &     (dx(nx-1)*(dx(nx-1)+dx(nx-2)))
      bb_1_G(nx)     = -(dx(nx-1)+dx(nx-2))/
     &     (dx(nx-1)*dx(nx-2))
      cc_1_G(nx)     =  dx(nx-1)/
     &     (dx(nx-2)*(dx(nx-1)+dx(nx-2)))
      d_1_G(nx)      = 0.0d0

!     Fattorizzazione LU
      
      do ix = 1,nx-1
         aux_alfa_1_G(ix) = alfa_1_G(ix)
         aux_gamma_1_G(ix) = 1.0d0
         aux_beta_1_G(ix) = beta_1_G(ix)
      enddo
      aux_gamma_1_G(nx) = 1.0d0

!$acc enter data copyin(aa_1_G,bb_1_G,cc_1_G,d_1_G)      
      if (use_gpu) then

      aux_beta_2_G(1:nx-1)=aux_beta_1_G(1:nx-1)
         aux_beta_2_G(nx)=0.0
      aux_alfa_2_G(2:nx)=aux_alfa_1_G(1:nx-1)
         aux_alfa_2_G(1) =0.0
!$acc enter data copyin(aux_beta_2_G, aux_alfa_2_G, aux_gamma_1_G)      
      endif

      CALL DGTTRF(nx, aux_alfa_1_G, aux_gamma_1_G, aux_beta_1_G, 
     &     ww_1_G,ipv_1_G, info)

      
      if (info > 0 .or. info < 0) then
         write(*,*) 'Problemi LU - init_der1x - , info:', info
         stop
      endif

***************************************************

***** Setup derivata seconda x tridiagonale LU per Psi  ***
! calcolo coefficienti interni 
!
      info = 0
      TRANS = 'N'
      nrhs  = 1
      ldb   = nx

      do ix = 2,nx-1
         dp1 = dx(ix)
         dm1 = dx(ix-1)
         dp2 = dp1 * dp1
         dm2 = dm1 * dm1
         dp3 = dp1 * dp2
         dm3 = dm1 * dm2
         temp1 = (dp3 + dm3) - (dp1 - dm1)
     &        * (5.0d0 * dp3 + 6.0d0 * dp2 * dm1 - dm3)
     &        / (2.0d0 * dp1 + dm1)
         temp2 = (5.0d0 * dm3 + 6.0d0 * dm2 * dp1 - dp3) +
     &        (dp1 + 2.0d0 * dm1) *
     &        (5.0d0 * dp3 + 6.0d0 * dp2 * dm1 - dm3) /
     &        (dm1 + 2.0d0 * dp1)
         alfa_2(ix-1) = temp1 / temp2
         beta_2(ix) = (alfa_2(ix-1) * (dp1 + 2.0d0 * dm1) +
     &        (dp1 - dm1)) / (dm1 + 2.0d0 * dp1)
         cc_2(ix) = 2.0d0*(1.0d0 + alfa_2(ix-1) + beta_2(ix)) /
     &        (dp1 + dm1) / dp1
         aa_2(ix) =   cc_2(ix) * dp1 / dm1
         bb_2(ix)   = - aa_2(ix) - cc_2(ix)
      enddo
      beta_2(1) = 0.0d0
      alfa_2(nx-1) = 0.0d0

      do ix = 1,nx-1
         aux_alfa_2(ix) = alfa_2(ix)
         aux_beta_2(ix) = beta_2(ix)
         aux_gamma_2(ix) = 1.0d0
      enddo
      aux_gamma_2(nx) = 1.0d0

      CALL DGTTRF(nx,  aux_alfa_2,  aux_gamma_2,  aux_beta_2,
     &     ww_2, ipv2, info)
      print*,info
      fact_alfa_2 =  aux_alfa_2
      fact_beta_2 =  aux_beta_2
      fact_gamma_2 =  aux_gamma_2
      fact_ww_2 = ww_2
      if (info > 0 .or. info < 0) then
         write(*,*) 'Problemi LU - init_der_II - , info:', info
         stop
      endif

      dp  = dx(1)
      dp1 = dx(1) + dx(2)
      dp2 = dx(1) + dx(2) + dx(3)
      ddsecder1 = 2.0d0 * (dp + dp1) / dp2 / (dp2 - dp1)
     &     / (dp - dp2)
      ccsecder1 = (2.0d0 - ddsecder1 * dp2 * (dp2 - dp)) / dp1
     &     / (dp1 - dp)
      bbsecder1 = -(ccsecder1 * dp1 + ddsecder1 * dp2) / dp
      aasecder1 = -ccsecder1 - ddsecder1 -bbsecder1
      dm1 = dx(nx-1)
      dm2 = dx(nx-1) + dx(nx-2)
      dm3 = dx(nx-1) + dx(nx-2) + dx(nx-3)
      ddsecdernx = -2.0d0 * (dm1 + dm2) / dm3
     &     / (dm1 - dm3) / (dm2-dm3)
      ccsecdernx = (2.0d0 + ddsecdernx * dm3 *  (dm1 -dm3))
     &     / dm2 /(dm2 - dm1)
      bbsecdernx = -(ccsecdernx * dm2 + ddsecdernx * dm3) / dm1
      aasecdernx = -ccsecdernx - ddsecdernx -bbsecdernx

!**************************************************
	return

	end
