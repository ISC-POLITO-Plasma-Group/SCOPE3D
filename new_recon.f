      program recon
      use nvtx
      include 'par.inc'
           
      real*8, dimension (:,:,:), allocatable :: 
     &     flxm1, flxm2, flxm3     
     &     ,flxp1, flxp2, flxp3
     &     ,rhshpxy, rhshpz,rhshmxy, rhshmz
     &     ,hp, hm, psi_temp,flxp_av
     &     ,psiy_temp,hpy_temp
!     &     ,uu_temp 
      double precision t1, t2, usrtime

!-----------------------------------------------------------------
!   nx = number of grid points in the x-direction
!   ny = number of grid points in the y-direction
!   nz = number of grid points in the z-direction

!   xl = box dimension in the x-direction
!   yl=  box dimension in the y-direction
!   zl=  box dimension in the z-direction
!----------------------------------------------------------------
!
!  Read data from file='in.com'
           
     call nvtxStartRange('new_recon',5) 
      CALL pstartup()

!     print *, ' lecture des donnees ....'
      
      IF( mpime == root ) then
         open(unit=7,status='old',file='in.com')
         read (7,*)        
         read (7,*) dt, nstep, noutt, noutf
         read (7,*)
         read (7,*) xl, yl, zl
         read (7,*)
         read (7,*) t0,gridalfa,gridbeta,gridgamma
         read (7,*)
         read (7,*) de, rhos, beta_e, tau, eta, mm, pso, psoeq, 
     &              phieq, eq_l, asym
         read (7,*)   
         read (7,*) istart, omegax, omega, alfa1
         read (7,*)   
         read (7,*) kbc_p, kbc_h
         read (7,*)   
         read (7,*) front_l, front_r, delta_ps
         read (7,*)   
         read (7,*) irestart_zero
         close (unit = 7)

        write(*,*) 'PHIEQ: ',phieq

      end if
      
      call bcast_real(dt,1,root)
      call bcast_integer(nstep,1,root)
      call bcast_integer(noutt,1,root)
      call bcast_integer(noutf,1,root)
      call bcast_real(xl,1,root)
      call bcast_real(yl,1,root)
      call bcast_real(zl,1,root)
      call bcast_real(t0,1,root)
      call bcast_real(gridalfa,1,root)
      call bcast_real(gridbeta,1,root)
      call bcast_real(gridgamma,1,root)
      call bcast_real(de,1,root)
      call bcast_real(rhos,1,root)
      call bcast_real(tau,1,root)
      call bcast_real(eta,1,root)
      call bcast_real(beta_e,1,root)
      call bcast_integer(mm,1,root)
      call bcast_real(pso,1,root)
      call bcast_real(psoeq,1,root)
      call bcast_real(phieq,1,root)
      call bcast_real(eq_l,1,root)
      call bcast_real(asym,1,root)
      call bcast_integer(istart,1,root)
      call bcast_integer(irestart_zero,1,root)
      call bcast_integer(kbc_p,1,root)
      call bcast_integer(kbc_h,1,root)
      call bcast_real(front_l,1,root)
      call bcast_real(front_r,1,root)
      call bcast_real(delta_ps,1,root)
      call bcast_real(omegax,1,root)
      call bcast_real(omega,1,root)
      call bcast_real(alfa1,1,root)
      
!   xl    : dimensions de la boite en unites de 2*pi
!   yl    : dimensions de la boite en unites de 2*pi
!   zl    : dimensions de la boite en unites de 2*pi

! .. allocate memory
      allocate(psi(nx,nyl,nzl))
      allocate(psi_temp(nx,nyl,nzl))
      allocate(psiy_temp(nx,nyl,nzl))
      allocate(hpy_temp(nx,nyl,nzl))
      allocate(phi(nx,nyl,nzl))
      allocate(cur(nx,nyl,nzl))
      allocate(cur_old(nx,nyl,nzl))
      allocate(uu(nx,nyl,nzl))
      allocate(dpsi_dt(nx,nyl,nzl))
      allocate(flxp_av(nx,nyl,nzl))
      allocate(flxm1(nx,nyl,nzl))
      allocate(flxm2(nx,nyl,nzl))
      allocate(flxm3(nx,nyl,nzl))
      allocate(flxp1(nx,nyl,nzl))
      allocate(flxp2(nx,nyl,nzl))
      allocate(flxp3(nx,nyl,nzl))



      allocate(hp(nx,nyl,nzl))
      allocate(hm(nx,nyl,nzl))
      allocate(rhshpxy(nx,nyl,nzl))
      allocate(rhshmxy(nx,nyl,nzl))
      allocate(rhshpz(nx,nyl,nzl))
      allocate(rhshmz(nx,nyl,nzl))
      allocate(fact_aa(nx1,nyl))
      allocate(fact_bb(nx,nyl))
      allocate(fact_cc(nx1,nyl))
      allocate(fact_ww(nx-2,nyl))
      allocate(ipv_my(nx,nyl))
      allocate(fact_aa_1(nx1,nyl))
      allocate(fact_bb_1(nx,nyl))
      allocate(fact_cc_1(nx1,nyl))
      allocate(fact_ww_1(nx-2,nyl))
      allocate(ipv_my_1(nx,nyl))
      allocate(hpy(nx,nyl,nzl))
      allocate(hmy(nx,nyl,nzl))
      allocate(psiy(nx,nyl,nzl))
      allocate(phiy(nx,nyl,nzl))
      allocate(e_para(nx,nyl,nzl))
      allocate(abs_b(nx,nyl,nzl))

!      allocate(d2_uu(nx,nyl,nzl))
!      allocate(uu_temp(nx,nyl,nzl))

  
      if(istart.EQ.1) then
         
         if( mpime == root ) then
            open(unit=19,status='unknown',file='Energie.dat')
 903        read(19,110,end=904)
            goto 903
 904        backspace 19
            open(unit = 29, status='unknown', file='int_rhs.dat')
 1025       read(29,129,end=1026)
            goto 1025
 1026       backspace 29
         end if
         
      endif

      pi = dacos(-1.0d0)
      
      grady = 1.0d0 / yl
      gradz = 1.0d0 / zl
      
      xl = 2.0d0 * xl
      yl = 2.0d0 * pi * yl
      zl = 2.0d0 * pi * zl
      
      nsortie=0
      
      ab1   = dt * 5.d0 / 12.d0
      ab2   = - dt * 4.d0 / 3.d0
      ab3   = dt * 23.d0 / 12.d0 
      
      
      call init_mod

      if( mpime == root ) then
         open(unit = 18, status='unknown',
     &           file='grid.dat')

      endif

      call out_grid

      if( mpime == root ) then
         close (unit=18)
      endif


!	call condinit_mod_sh_wave_new(hp,hm)
!      call condinit_mod_dh(hp,hm)
!       call condinit_mod_sh(hp,hm)

       call condinit_kh(hp,hm)

!!!!!!!!!! INITIAL CONDITIONS MR!!!!!!!!!
!       call condinit_mod_ran(hp,hm)
!       call condinit_mod_ran_w3D(hp,hm,hp_pascal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call condinit(hp,hm,hp_pascal)

!       call condinit_mod_ran_dh(hp,hm)
!       call condinit_mod_wave(hp,hm)
!        call condinit_mod_random_3(hp,hm)

      IF(istart.EQ.0) THEN  
         
         oldtime = 0.
         time = 0. 
         ioutt = 1
         ioutf = 1

         call energie_mod

!*************RIPRISTINA******
         call out_field

!**********TEST CORRENTE_RES**********
         call corrente_res

!         call out_field         
!*****************************
         if( mpime == root ) then
            open(unit = 19, POSITION='append', status='unknown',
     &           file='Energie.dat')
            open(unit = 29, POSITION='append', status='unknown',
     &           file='int_rhs.dat')
            
         end if
         
         call out_inv
         
         if( mpime == root ) then
            close (unit=19)
            close (unit=29)
         endif
   
        nsortie=nsortie+1
        
        ENDIF



      IF(istart.EQ.1) THEN  

         time = oldtime
!         call out_field
           uu=hm

!******** TAU NE 0 ************
!           call B_func
!
!           uu_temp = uu
!           uu = d2_uu
!           call poisson_mod
!           uu  =uu_temp
!******************************
!******** TAU eq  0 ************

           call poisson_mod
!******************************
           call helm_mod(hp)

           psi_temp = psi
           hpy_temp = hpy
           psiy_temp = psiy
!****************************************
           psi = psi_temp
           hpy = hpy_temp
           psiy = psiy_temp

 
      ENDIF

!------------------------------------------------------
! boucle sur le temps:
! metodo di Eulero per i primi 2 passi
!------------------------------------------------------

        call nvtxStartRange('first_timestep',1)
        time = dt + oldtime
        oldtime = time

        write(*,*) 'time=', time

        uu =hm 
        cur_old=cur
        if (de.eq.0) then 
          psi = hp
          call corrente_res 
        else
!          call corrente_res 
          cur=(hp-psi)/de2    !!!!!DA RIPRISTINARE SUBITO!!!!!!!!!!!!!!
        endif

!        call rhs_xy_lin(hp,hp_pascal,rhshpxy,rhshmxy,rhspascalxy,
!     &           rhstparexy,rhstperpexy)

!        WRITE(*,*) 'BEFORE RHS_XY',mpime
        call rhs_xy(hp,rhshpxy,rhshmxy)
!        WRITE(*,*) 'AFTER RHS_XY',mpime
!        WRITE(*,*) 'BEFORE RHS_Z',mpime
        call rhs_z(rhshpz,rhshmz)
!        WRITE(*,*) 'AFTER RHS_Z',mpime


        flxp1 = rhshpxy +rhshpz
        flxm1 = rhshmxy +rhshmz
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        call helm_mod(flxp1)
!        dpsi_dt = psi

!        flxpascal1 = dpsi_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call out_fieldrhs(flxp1,flxm1)
!        call out_field_rhs(rhshpxy,rhshmxy)

        hp = hp  + flxp1 * dt
        hm = hm  + flxm1 * dt


!!!!!!!TEST KH!!!!!!!!!!!!!!!!!!!
         uu = hm
         call poisson_mod
!         call out_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call filterx_t(hp)
        call filterx_t(hm)

        call filtery_t(hp)
        call filtery_t(hm) 


!        write(*,*) 'BEFORE FILTERZ HACK'

        if (nz .ge.4) then

        call filterz_t(hp)
        call filterz_t(hm) 

        endif

!        write(*,*) 'AFTER FILTERZ HACK'


      call nvtxEndRange()

!--------------------------------------------------

        call nvtxStartRange('second_timestep',2)
        time = dt + oldtime

        oldtime = time

        uu = hm

!******** TAU NE 0 ************
!          call B_func

!           uu_temp = uu
!           uu = d2_uu
!           call poisson_mod
!           uu  =uu_temp
!******************************
!******** TAU eq 0************

           call poisson_mod
!******************************


        cur_old=cur
        if (de.eq.0) then
           psi = hp
        call corrente_res
        else
           call helm_mod(hp)
!           call corrente_res 
           cur = (hp -psi)/de2     !!!!DA RIPRISTINARE  SUBITO!!!!!
           psi_temp = psi
           hpy_temp = hpy
           psiy_temp = psiy
!****************************************

           psi = psi_temp
           hpy = hpy_temp
           psiy = psiy_temp

        endif

!        call heat_flux_mod   UNSTABLE

!        call rhs_xy_lin(hp,hp_pascal,rhshpxy,rhshmxy,rhspascalxy,
!     &           rhstparexy,rhstperpexy)
        call rhs_xy(hp,rhshpxy,rhshmxy)
        call rhs_z(rhshpz,rhshmz)

        flxp2 = rhshpxy +rhshpz
        flxm2 = rhshmxy +rhshmz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        flxpascal2 = rhspascalxy +rhspascalz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        call helm_mod(flxp2)
!        dpsi_dt = psi

!        flxpascal2 = dpsi_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call out_field_rhs(flxp2,flxm2)
!        call out_field_rhs(rhshpxy,rhshmxy)

        hp = hp  + flxp2 * dt
        hm = hm  + flxm2 * dt

!!!!!!!TEST KH!!!!!!!!!!!!!!!!!!!
         uu = hm
         call poisson_mod
!         call out_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call filterx_t(hp)
        call filterx_t(hm)

        call filtery_t(hp)
        call filtery_t(hm) 

        if (nz .ge.4) then

        call filterz_t(hp)
        call filterz_t(hm) 

        endif
           
      call nvtxEndRange()
!------------------------------------------------------
!      inizio Adams-Bashforth
!-----------------------------------------------------

        ir = 1
        call nvtxStartRange('main_loop',3)
        do  it = 1, nstep 

        call nvtxStartRange('timestep',4)
           time = it * dt + oldtime
     
!          call heat_flux_mod        

           uu = hm

!******** TAU NE 0 ************
!          call B_func

!           uu_temp = uu
!           uu = d2_uu
!           call poisson_mod
!           uu  =uu_temp
!******************************
!******** TAU eq 0 ************
           call poisson_mod
!******************************

           cur_old=cur
           if (de.eq.0) then
              psi = hp
              
              call corrente_res

           else
              call helm_mod(hp)
!              call corrente_res 
              cur = (hp -psi)/de2    !!!!!!DA RIPRISTINARE SUBITO
              psi_temp = psi
              hpy_temp = hpy
              psiy_temp = psiy
!*****************************************
              psi = psi_temp
              hpy = hpy_temp
              psiy = psiy_temp

           endif

!           call heat_flux_mod  UNSTABLE       

!           call rhs_xy_lin(hp,hp_pascal,rhshpxy,rhshmxy,rhspascalxy,
!     &         rhstparexy,rhstperpexy)

           call rhs_xy(hp,rhshpxy,rhshmxy)
           call rhs_z(rhshpz,rhshmz)

!*************TEST DPSIDT****************           
           flxp_av = rhshpxy +rhshpz
!           call helm_mod(flxp_av)
!           dpsi_dt = psi
!****************************************

               if(ir.EQ.1) then
                  flxp3 = rhshpxy +rhshpz
                  flxm3 = rhshmxy +rhshmz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  flxpascal3 = rhspascalxy +rhspascalz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  call helm_mod(flxp3)    
!                  dpsi_dt = psi
!                  flxpascal3 = dpsi_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               elseif(ir.EQ.2) then
                  flxp1 = rhshpxy +rhshpz
                  flxm1 = rhshmxy +rhshmz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  flxpascal1 = rhspascalxy +rhspascalz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  call helm_mod(flxp1)    
!                  dpsi_dt = psi
!                  flxpascal1 = dpsi_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               elseif(ir.EQ.3) then
                  flxp2 = rhshpxy +rhshpz
                  flxm2 = rhshmxy +rhshmz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  flxpascal2 = rhspascalxy +rhspascalz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  call helm_mod(flxp2)    
!                  dpsi_dt = psi
!                  flxpascal2 = dpsi_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               endif
 
               if(ir.EQ.1) then
                 hp = hp+ab1*flxp1+ab2*flxp2+ab3*flxp3
                 hm = hm+ab1*flxm1+ab2*flxm2+ab3*flxm3 

 
!***********TEST ANNA*******************
                 flxp_av = (ab1*flxp1+ab2*flxp2+ab3*flxp3)/dt   
!                 call helm_mod(flxp_av)
!                 dpsi_dt = psi 
!***************************************
                 elseif(ir.EQ.2) then
                 hp = hp+ab1*flxp2+ab2*flxp3+ab3*flxp1
                 hm = hm+ab1*flxm2+ab2*flxm3+ab3*flxm1

!***********TEST ANNA*******************
                 flxp_av = (ab1*flxp2+ab2*flxp3+ab3*flxp1)/dt
!                 call helm_mod(flxp_av)
!                 dpsi_dt = psi
!***************************************
                 elseif(ir.EQ.3) then
                 hp = hp+ab1*flxp3+ab2*flxp1+ab3*flxp2
                 hm = hm+ab1*flxm3+ab2*flxm1+ab3*flxm2

!***********TEST ANNA*******************
                 flxp_av = (ab1*flxp3+ab2*flxp1+ab3*flxp2)/dt
!                 call helm_mod(flxp_av)
!                 dpsi_dt = psi
!***************************************
                 ir = 0
               endif
 
               ir = ir + 1 


               call filterx_t(hp)
               call filterx_t(hm)


               call filtery_t(hp)
               call filtery_t(hm) 


               if (nz .ge.4) then

               call filterz_t(hp)
               call filterz_t(hm)

       
               endif
    
 
!               call filterx_t(qpare)        
!               call filtery_t(qpare)
!               call filterz_t(qpare)

         if( mpime == root ) then
         if(mod(it,10).eq.0) then
         open(unit = 70, status='unknown',
     &      file='tempo.dat')  
            write(70,*) it
         end if
         end if

         if( mpime == root ) then
            close (unit=70)
         end if
            
!------------------------------------------------------------
!*******************************

        if(ioutt.ge.noutt) then
           nsortie=nsortie+1
           uu = hm 


           if( mpime == root ) then
              open(unit = 19, POSITION='append', status='unknown',
     &             file='Energie.dat')
              open(unit = 29, POSITION='append', status='unknown',
     &             file='int_rhs.dat')

           end if

           call energie_mod

           call out_inv
           if( mpime == root ) then
              close (unit=19)
              close (unit=29)

           end if
           ioutt = 0
        endif

        if(ioutf.ge.noutf) then
           uu = hm
           psi_temp = psi
           hpy_temp = hpy
           psiy_temp = psiy 
           call helm_mod(flxp_av)
           dpsi_dt = psi 
           psi = psi_temp
           hpy = hpy_temp
           psiy = psiy_temp
           e_para = (e_para+dpsi_dt)/abs_b

           call filterx_t(e_para)
           call filtery_t(e_para)

           if (nz .ge.4) then
             call filterz_t(e_para)
           endif   

           call out_field

!           call out_field_rhs(rhshpxy+rhshpz,rhshmxy+rhshmz)
!           call out_field_rhs(rhshpxy,rhshmxy)

           ioutf = 0
        endif

        ioutt = ioutt + 1
        ioutf = ioutf + 1
      
      call nvtxEndRange()
      enddo
      call nvtxEndRange()
!********************************   
!    PARTE NUOVA (per il riavvio)

        oldtime = time

        call gpm_write(hp, hm, 11)

!*************test riavvio************************
!************************************

!	close (unit = 12)
!	close (unit = 14)

!101     format(1x, 1e11.4)
!102     format(1x, 2(1x, 1e12.5))
104     format(1x, 2(1x, 1e22.9))
103     format(1x, 5(1x, 1e22.9))
105     format(1x, 4(1x, 1e12.5))
109     format(1x, 1(1x, 1e22.9))
110     format(1x, 6(1x, 1e22.9))
129     format(1x, 3(1x, 1e22.9))
!106     format(1x, 6(1x, 1e12.5))
!107     format(1x, 7(1x, 1e12.5))
!108     format(1x, 6(1x, 1e12.5))

! .. deallocate memory
        deallocate(psi, phi, cur, uu)
        deallocate(flxm1, flxm2, flxm3)
        deallocate(flxp1, flxp2, flxp3, flxp_av)
        deallocate(hp, hm)
        deallocate(dpsi_dt)
        deallocate(psi_temp)
        deallocate(psiy_temp)
        deallocate(hpy_temp)
        deallocate(rhshpxy)
        deallocate(rhshmxy)
        deallocate(rhshpz)
        deallocate(rhshmz)
        deallocate(cur_old)
        deallocate(fact_aa,fact_bb,fact_cc,fact_ww,ipv_my)
        deallocate(fact_aa_1,fact_bb_1,fact_cc_1)
        deallocate(fact_ww_1,ipv_my_1)
        deallocate(e_para,abs_b)
!$acc exit data delete(aux_alfa_2_G, aux_beta_2_G, aux_gamma_1_G)
!        deallocate(d2_uu)
!        deallocate(uu_temp)

        call hangup()
        call nvtxEndRange()
        stop
        end 
