      subroutine rhs_xy(hp1,rhshpxy1,rhshmxy1)
      use nvtx
      use gpu_func 
      include 'par.inc'

      dimension d1(nx,nyl),d2(nx,nyl),dd1(nx,nyl),dd2(nx,nyl),
     & ddd1(nx,nyl),dd3(ny)
      dimension f1(nx,nyl),f2(nx,nyl),ff1(nx,nyl),ff2(nx,nyl),
     & fff1(nx,nyl),ff3(ny)
      dimension aux(nx,nyl,nzl),aux_t(nxl,ny,nzl)
      
      real(8), dimension (nx,nyl,nzl) ::  
     &     rhshpxy1,rhshmxy1,hp1

      real(8), dimension (nx,nyl,nzl) ::  
     &     phix,psix,hpx,hmx,curx,cury

      integer(4) ierr



c      dimension zmagn(nx,nyl)



c     derivo (rispetto a x)
       call nvtxStartRange('rhs_xy',13)
!$acc enter data copyin(phi,psi,cur,hp1,uu,x) 
!$acc& create(f1(1:nx,1:nyl),f2(1:nx,1:nyl),ff1(1:nx,1:nyl),
!$acc& ff2(1:nx,1:nyl),fff1(1:nx,1:nyl))
!$acc& create(d1(1:nx,1:nyl),d2(1:nx,1:nyl),dd1(1:nx,1:nyl),
!$acc& dd2(1:nx,1:nyl),ddd1(1:nx,1:nyl))
!$acc& create(phix(1:nx,1:nyl,1:nzl),curx(1:nx,1:nyl,1:nzl))
!$acc& create(psix(1:nx,1:nyl,1:nzl))
!$acc& create(hpx(1:nx,1:nyl,1:nzl),hmx(1:nx,1:nyl,1:nzl))       
      do iz = 1,nzl
!$acc parallel loop present(f1,f2,ff1,ff2,fff1,phi,psi,cur,hp1,uu) 
!$acc& collapse(2)      
         do iy = 1,nyl
            do ix = 1,nx
               f1(ix,iy) = phi(ix,iy,iz)
               f2(ix,iy) = psi(ix,iy,iz)
               ff1(ix,iy) = cur(ix,iy,iz)
               ff2(ix,iy) = hp1(ix,iy,iz)
               fff1(ix,iy) = uu(ix,iy,iz)
               
            enddo
         enddo

            CALL der1x_y(f1,d1,use_gpu)
            CALL der1x_y(f2,d2,use_gpu)
            CALL der1x_y(ff1,dd1,use_gpu)
            CALL der1x_y(ff2,dd2,use_gpu)
            CALL der1x_y(fff1,ddd1,use_gpu)
       call nvtxStartRange('after der1x_y',13)
         ierr=cudaDeviceSynchronize() 
!$acc parallel loop collapse(2)
!$acc& present(phix,d1,x,psix,d2,curx,dd1,hpx,dd2,hmx,ddd1)         
        do iy = 1,nyl
            do ix = 1,nx
!***************** Bikley jet **************************
               phix(ix,iy,iz) = 1.0*d1(ix,iy)+ phieq/(dcosh(x(ix)))**2.0d0
!*******************************************************
!***************** Vortex sheet **************************
!               phix(ix,iy,iz) = d1(ix)+ phieq*tanh(x(ix)/eq_l)
!*******************************************************
!***************Harris pinch****************************
               psix(ix,iy,iz) = 1.0*d2(ix,iy) - psoeq
     &                 * dtanh(x(ix)/eq_l)
     &                 + asym*yl/zl
               curx(ix,iy,iz) = 1.0*dd1(ix,iy) - psoeq  
     &          * 2.0d0 * dtanh(x(ix)/eq_l)
     &          *(1.0d0 - (dtanh(x(ix)/eq_l))**2.0d0)/(eq_l**2.)
               hpx(ix,iy,iz) = 1.0*dd2(ix,iy) - psoeq 
     &              * dtanh(x(ix)/eq_l)
     &              + asym*yl/zl
     &          - de2 * psoeq * 2.0d0 * dtanh(x(ix)/eq_l) * 
     &          (1.0d0 - dtanh(x(ix)/eq_l)*dtanh(x(ix)/eq_l))
     &         /(eq_l**2.)
!********************************************************
!************* Bikley jet ******************************* 
               hmx(ix,iy,iz) = 1.0*ddd1(ix,iy) - 2.0d0*phieq
     &              *(1.0d0-2.0d0*(dsinh(x(ix)))**2.0d0)
     &              /(dcosh(x(ix)))**4.0d0          
!********************************************************
!************* vortex sheet ******************************* 
!               hmx(ix,iy,iz) = ddd1(ix) - phieq
!     &          * 2.0d0 * dtanh(x(ix)/eq_l)
!     &          *(1.0d0 - (dtanh(x(ix)/eq_l))**2.0d0)/(eq_l**2.) 
!********************************************************

            enddo
         enddo
         call nvtxEndRange()
      enddo
!$acc exit data delete(phi,psi,cur,hp1,uu,x)
!$acc& delete(f1,f2,ff1,ff2,fff1)
!$acc& delete(d1,d2,dd1,dd2,ddd1)
!$acc& copyout(phix(1:nx,1:nyl,1:nzl),curx(1:nx,1:nyl,1:nzl))
!$acc& copyout(psix(1:nx,1:nyl,1:nzl))
!$acc& copyout(hpx(1:nx,1:nyl,1:nzl),hmx(1:nx,1:nyl,1:nzl))       

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

        rhshpxy1 = - (phix * hpy - phiy * hpx)
!     &             - (0.3d0*0.3d0) * (psix * hmy - psiy * hmx)
     &             - (rhos2) * (psix * hmy - psiy * hmx)
     &             - eta * cur_old

        rhshmxy1 = - (phix * hmy - phiy * hmx)
     &             - 1.0*2./beta_e * (psix * cury - psiy * curx)

!        rhspascalxy1 = - (psix * pascaly - psiy * pascalx)


!*****************TEST*****************************
!        rhstparexy1 = 2.*rhshmxy1
!**************************************************





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

        call nvtxEndRange()
      end 
