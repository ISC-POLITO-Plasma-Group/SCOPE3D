      subroutine energie_mod
      use nvtx
      include 'par.inc'

      dimension zmagn(nx,nyl)
      dimension fx(nx), fy(ny), f1x(nx), f1y(ny)
      dimension Ene_z(nzl)
      dimension aux(nx,nyl,nzl), aux_t(nxl,ny,nzl)
      dimension E_aux(nx,nyl,nzl)
      real*8 deta,aux_1

      call nvtxStartRange('energie_mod',7)
c*******energia magnetica
      
      do iz = 1,nzl

         do iy = 1, nyl
            do ix = 1, nx
               fx(ix) = psi(ix,iy,iz) 
            enddo
            CALL der1x_free(fx, f1x)
            do ix = 1, nx
               f1x(ix) = f1x(ix) - psoeq * dtanh(x(ix)/eq_l)/eq_l  
     &              + asym*yl/zl
               E_aux(ix,iy,iz) = f1x(ix) * f1x(ix)
            enddo
         enddo
  
      enddo

      aux = psi 
      call trasponi_yx(aux,aux_t,1)

      do iz = 1,nzl
 
         do ix = 1, nxl
            do iy = 1, ny
               fy(iy) = aux_t(ix,iy,iz)
            enddo
            CALL der1y(fy, f1y)
            do iy = 1, ny
               aux_t(ix,iy,iz) =  f1y(iy) 
            enddo
         enddo

      enddo

      call trasponi_yx(aux,aux_t,-1)

      zmagn = 0.0

      do iz = 1,nzl

         do iy = 1, nyl
            do ix = 1, nx
               zmagn(ix,iy) = (E_aux(ix,iy,iz) + 
     &          aux(ix,iy,iz)*aux(ix,iy,iz)) * dxyz(ix)
            enddo
         enddo
         
         Ene_z(iz) = sum(zmagn)

      enddo
      
      Emag = sum(Ene_z) 

      CALL PARALLEL_SUM_REAL(Emag,1)

      Emag = Emag/beta_e  

c*******energia cinetica elettronica 

      if (de.eq.0) then 
         deta= eta
      else
         deta=de
      endif
      do iz = 1, nzl
         do iy = 1, nyl
            do ix = 1, nx
        aux_1 = cur(ix,iy,iz) + psoeq*(1.0d0 - (tanh(x(ix)/eq_l))**2.)
     &    /(eq_l**2.)
          zmagn(ix,iy) = deta * deta * aux_1*aux_1*dxyz(ix)
            enddo
         enddo
         
         Ene_z(iz) = sum(zmagn)

      enddo

      Eke  = sum(Ene_z) 
      CALL PARALLEL_SUM_REAL(Eke,1)

      Eke  = Eke/beta_e 


c*******energia cinetica ionica

      do iz = 1,nzl
         do iy = 1, nyl
            do ix = 1, nx
               zmagn(ix,iy) = - phi(ix,iy,iz) * uu(ix,iy,iz) *
     &              dxyz(ix)
            enddo
         enddo
         
         Ene_z(iz) = sum(zmagn)
         
        enddo

        Ekp  = sum(Ene_z) 
        CALL PARALLEL_SUM_REAL(Ekp,1)

        Ekp  =Ekp/2.d0

C*******energia potenziale elettronica

        do iz =1, nzl
           do iy = 1, nyl
              do ix = 1, nx
                 zmagn(ix,iy) = uu(ix,iy,iz) * uu(ix,iy,iz)
     &                * dxyz(ix)
              enddo
           enddo
           Ene_z(iz) = sum(zmagn)
           
        enddo
        
        Epe  = sum(Ene_z) 

        CALL PARALLEL_SUM_REAL(Epe,1)

        Epe = Epe/2.0d0

c*******energia totale
        
        Etot = Emag + Eke + Ekp + Epe  
        call nvtxEndRange()
        return
        end

