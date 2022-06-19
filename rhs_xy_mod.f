      subroutine rhs_xy(hp1,hp1_pascal,rhshpxy1,rhshmxy1,rhspascalxy1,
     &            rhstparexy1,rhstperpexy1)
      
      include 'par.inc'

      dimension d1(nx),d2(nx),dd1(nx),dd2(nx),ddd1(nx),dd3(ny)
      dimension ddd2(nx),ddd3(nx)
      dimension ddd4(nx),ddd5(nx)
      dimension ddd6(nx),ddd7(nx)
      dimension f1(nx),f2(nx),ff1(nx),ff2(nx),fff1(nx),ff3(ny)
      dimension fff2(nx),fff3(nx)
      dimension fff4(nx),fff5(nx)
      dimension fff6(nx),fff7(nx)
      dimension aux(nx,nyl,nzl),aux_t(nxl,ny,nzl)
      
      real(8), dimension (nx,nyl,nzl) ::  
     &     rhshpxy1,rhshmxy1,rhspascalxy1,hp1,
     &     rhstparexy1,rhstperpexy1,
     &     hp1_pascal

      real(8), dimension (nx,nyl,nzl) ::  
     &     phix,psix,hpx,hmx,curx,cury,
     &     tparex,tperpex,
     &     tparey,tperpey,
     &     pascalx,pascaly,
     &     qparex,qperpex,
     &     qparey,qperpey,
     &     hp_pascalx,hp_pascaly


c      dimension zmagn(nx,nyl)



c     derivo (rispetto a x)

      do iz = 1,nzl
         do iy = 1,nyl
            do ix = 1,nx
               f1(ix) = phi(ix,iy,iz)
               f2(ix) = psi(ix,iy,iz)
               ff1(ix) = cur(ix,iy,iz)
               ff2(ix) = hp1(ix,iy,iz)
               fff1(ix) = uu(ix,iy,iz)
               fff2(ix) = pascal(ix,iy,iz)
               fff3(ix) = hp1_pascal(ix,iy,iz)
               fff4(ix) = tpare(ix,iy,iz)
               fff5(ix) = tperpe(ix,iy,iz)
               fff6(ix) = qpare(ix,iy,iz)
               fff7(ix) = qperpe(ix,iy,iz)
               
            enddo
            CALL der1x_free(f1,d1)
            CALL der1x_free(f2,d2)
            CALL der1x_free(ff1,dd1)
            CALL der1x_free(ff2,dd2)
            CALL der1x_free(fff1,ddd1)
            CALL der1x_free(fff2,ddd2)
            CALL der1x_free(fff3,ddd3)
            CALL der1x_free(fff4,ddd4)
            CALL der1x_free(fff5,ddd5)
            CALL der1x_free(fff6,ddd6)
            CALL der1x_free(fff7,ddd7)
            do ix = 1,nx
!***************** Bikley jet **************************
!               phix(ix,iy,iz) = d1(ix)+ phieq/(dcosh(x(ix)))**2.0d0
!*******************************************************
!***************** Vortex sheet **************************
               phix(ix,iy,iz) = d1(ix)+ phieq*tanh(x(ix)/eq_l)
!*******************************************************
!***************Harris pinch****************************
               psix(ix,iy,iz) = d2(ix) - psoeq * dtanh(x(ix)/eq_l)
     &                 + asym*yl/zl
               curx(ix,iy,iz) = dd1(ix) - psoeq  
     &          * 2.0d0 * dtanh(x(ix)/eq_l)
     &          *(1.0d0 - (dtanh(x(ix)/eq_l))**2.0d0)/(eq_l**2.)
               hpx(ix,iy,iz) = dd2(ix) - psoeq * dtanh(x(ix)/eq_l)
     &              + asym*yl/zl
     &          - de2 * psoeq * 2.0d0 * dtanh(x(ix)/eq_l) * 
     &          (1.0d0 - dtanh(x(ix)/eq_l)*dtanh(x(ix)/eq_l))
     &         /(eq_l**2.)
!********************************************************
!************* Bikley jet ******************************* 
!               hmx(ix,iy,iz) = ddd1(ix) - 2.0d0*phieq
!     &              *(1.0d0-2.0d0*(dsinh(x(ix)))**2.0d0)
!     &              /(dcosh(x(ix)))**4.0d0          
!********************************************************
!************* vortex sheet ******************************* 
               hmx(ix,iy,iz) = ddd1(ix) - phieq
     &          * 2.0d0 * dtanh(x(ix)/eq_l)
     &          *(1.0d0 - (dtanh(x(ix)/eq_l))**2.0d0)/(eq_l**2.) 
!********************************************************
               pascalx(ix,iy,iz) = ddd2(ix)  
               hp_pascalx(ix,iy,iz) = ddd3(ix)  
               tparex(ix,iy,iz) = ddd4(ix)  
               tperpex(ix,iy,iz) = ddd5(ix)  
               qparex(ix,iy,iz) = ddd6(ix)  
               qperpex(ix,iy,iz) = ddd7(ix)  

            enddo
         enddo
      enddo


!---------derivo cur rispetto a y--------------

        aux = cur

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        cury = aux

!--------------------------------------------------
!---------derivo pascal rispetto a y--------------

        aux = pascal

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        pascaly = aux

!--------------------------------------------------
!---------derivo hp_pascal rispetto a y--------------

        aux = hp1_pascal

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        hp_pascaly = aux

!--------------------------------------------------
!---------derivo T_e_|| rispetto a y--------------

        aux = tpare

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        tparey = aux

!-------------------------------------------------
!---------derivo T_e_perp rispetto a y--------------
        
        aux = tperpe
     
        call trasponi_yx(aux,aux_t,1)
        
        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        tperpey = aux

!-------------------------------------------------
!---------derivo Q_e_|| rispetto a y--------------

        aux = qpare

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        qparey = aux

!-------------------------------------------------
!---------derivo T_e_perp rispetto a y--------------

        aux = qperpe

        call trasponi_yx(aux,aux_t,1)

        do iz = 1,nzl
           do ix = 1,nxl
              do iy = 1,ny
                 ff3(iy) = aux_t(ix,iy,iz)
              enddo
              CALL der1y(ff3,dd3)
              do iy = 1,ny
                 aux_t(ix,iy,iz) = dd3(iy)
              enddo
           enddo
        enddo

        call trasponi_yx(aux,aux_t,-1)

        qperpey = aux

!-------------------------------------------------

        rhshpxy1 = - (phix * hpy - phiy * hpx)
!     &             - (0.3d0*0.3d0) * (psix * hmy - psiy * hmx)
     &             - (0.0) * (psix * hmy - psiy * hmx)
     &             - 0.0d0 * (psix * tparey - psiy * tparex)
     &             - eta * cur_old

        rhshmxy1 = - (phix * hmy - phiy * hmx)
     &             - 1.0*2./beta_e * (psix * cury - psiy * curx)

!        rhspascalxy1 = - (psix * pascaly - psiy * pascalx)

        rhstparexy1 = - (phix * tparey - phiy * tparex)
     &             - 4./beta_e * (psix * cury - psiy * curx)
     &             + 0.0*(psix * qparey - psiy * qparex)

!*****************TEST*****************************
!        rhstparexy1 = 2.*rhshmxy1
!**************************************************


        rhstperpexy1 = - (phix * tperpey - phiy * tperpex)
     &             + 0.0*(psix * qperpey - psiy * qperpex)

!******************** OK - Borgogno version ********
        rhspascalxy1 = - (phix * hp_pascaly - phiy * hp_pascalx)
     &             - (pascalx * hmy - pascaly * hmx)
!****************************************************
!     &             - eta * cur_old


        e_para = phix * psiy - phiy * psix
        abs_b = sqrt(psix**2. + psiy**2. + 1.0d0)

!      do iz = 1,nzl

!        rhshpxy1(1:nx,1:nyl,iz) =  - (phix(1:nx,1:nyl) * hpy(1:nx,1:nyl)
!     &          - phiy(1:nx,1:nyl) * hpx(1:nx,1:nyl))
!     &       - rhos2 * (psix(1:nx,1:nyl) * hmy(1:nx,1:nyl)
!     &          - psiy(1:nx,1:nyl) * hmx(1:nx,1:nyl)) 
!     &          - eta * cur_old(1:nx,1:nyl,iz)
!        rhshmxy1(1:nx,1:nyl,iz) = - (phix(1:nx,1:nyl) * hmy(1:nx,1:nyl)
!     &       - phiy(1:nx,1:nyl) * hmx(1:nx,1:nyl)) 
!     &       - (psix(1:nx,1:nyl) * cury(1:nx,1:nyl)  
!     &       - psiy(1:nx,1:nyl) * curx(1:nx,1:nyl)) 

!      enddo



c ********** Parte utile per l'energia magnetica

c         do iy = 1,nyl
c            do ix = 1, nx
c               zmagn(ix,iy) = (psix(ix,iy)**2. + psiy(ix,iy)**2.)
c     &                         * dxyz(ix)
c            enddo
c         enddo
c         Ene_z(iz) = sum(zmagn)

c      Emag = sum(Ene_z)

c      CALL PARALLEL_SUM_REAL(Emag,1)


      end 
