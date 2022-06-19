      subroutine helm_mod(hp)
      use nvtx
c  risolve Helmholtz: nabla^2 psi + elambda psi = g2 (psi output)

      include 'par.inc' 
        
      dimension q1(ny) 
      dimension F0_my(nxl,ny),f1(nx,nyl)
      dimension aux_aa(nx-1),aux_cc(nx-1)
      dimension aux_bb(nx),aux_f1(nx)
      dimension aux_ww(nx-2),f1_my(nx,nyl)
      dimension q2(ny), F0_my_red(nx,nyl)
      dimension hp(nx,nyl,nzl)
      dimension aux(nx,nyl,nzl), aux_t(nxl,ny,nzl)
      dimension aux_t2(nxl,ny,nzl),aux_t3(nxl,ny,nzl)
      dimension ff4(ny),dd4(ny),f4(ny),d4(ny)

      call nvtxStartRange('helm_mod',11)
      aux = hp

      call trasponi_yx(aux,aux_t,1)


!---------derivo F_t (rispetto a y)---------------

      aux_t2 = aux_t

      do iz = 1,nzl
         do ix = 1,nxl
            do iy = 1, ny
               ff4(iy) = aux_t2(ix,iy,iz)
            enddo
            CALL der1y(ff4,dd4)
            do iy = 1, ny
               aux_t2(ix,iy,iz) = dd4(iy)
            enddo
         enddo
      enddo

!--------------------------------------------------


      do iz = 1,nzl

         do ix = 1,nxl
            do iy = 1,ny
               q1(iy) = -aux_t(ix,iy,iz)/de2
            enddo

            call drfftf(ny, q1, wsavey)

            do iy = 1,ny
               F0_my(ix,iy) = q1(iy)/ny
            enddo
         enddo

         aux_t(1:nxl,1:ny,iz) = F0_my(1:nxl,1:ny)

      enddo

      call trasponi_yx(aux,aux_t,-1)

      do iz = 1, nzl

         F0_my_red(1:nx,1:nyl) = aux(1:nx,1:nyl,iz)

         do my = 1,nyl
            do ix = 2, nx-1
        f1(ix,my) = alfa(ix) * F0_my_red(ix-1,my) + F0_my_red(ix,my) 
     &           + beta(ix) * F0_my_red(ix+1,my) 
            enddo
         enddo

!     Condizioni al contorno - caso k = 0
         
         f1(1,1)  = 0.0 
         f1(nx,1) = 0.0

!     Condizioni al contorno - caso k noneq 0

         do my = 2,nyl
            f1(1,my) =  0.0
            f1(nx,my) = 0.0
         enddo
         
         
!     Soluzione del sistema lineare A x = y
         
         do my = 1,nyl
            do ix = 1,nx-1  
               aux_aa(ix) = fact_aa_1(ix,my)
               aux_bb(ix) = fact_bb_1(ix,my)
               aux_cc(ix) = fact_cc_1(ix,my)
               ipv(ix) = ipv_my_1(ix,my)
               aux_f1(ix) = f1(ix,my)
            enddo
            aux_bb(nx) = fact_bb_1(nx,my)
            ipv(nx) = ipv_my_1(nx,my)
            aux_f1(nx) = f1(nx,my)
            
            do ix = 1,nx-2
               aux_ww(ix) = fact_ww_1(ix,my)
            enddo
            
            
            CALL  DGTTRS(TRANS,nx,NRHS,aux_aa,aux_bb,aux_cc,
     &           aux_ww,ipv,aux_f1,LDB,INFO)
            
            if (info > 0 .or. info < 0) then
               write(*,*) 'Problemi soluzione Helmoltz, info:', info
               stop
            endif
                        
            do ix = 1,nx
               f1_my(ix,my) = aux_f1(ix) 
            enddo
         
         enddo

         aux(1:nx,1:nyl,iz) = f1_my(1:nx,1:nyl)
   
      enddo

      call trasponi_yx(aux,aux_t,1)
         
      do iz = 1,nzl 

         do ix = 1,nxl
            do my = 1,ny
               q2(my) =  aux_t(ix,my,iz) 
            enddo
            call drfftb(ny, q2, wsavey)
            do iy = 1,ny
               aux_t(ix,iy,iz) = q2(iy)
            enddo
         enddo
         
      enddo


!---------derivo psi_t (rispetto a y)--------------

      aux_t3 = aux_t

      do iz = 1,nzl
         do ix = 1,nxl
            do iy = 1, ny
               f4(iy) = aux_t3(ix,iy,iz)
            enddo
            CALL der1y(f4,d4)
            do iy = 1, ny
               aux_t3(ix,iy,iz) = d4(iy)
            enddo
         enddo
      enddo

!--------------------------------------------------


      call trasponi_yx(aux,aux_t,-1) 

      psi = aux

 
!-----------------------------------------------------

         call trasponi_yx(aux,aux_t2,-1)

         hpy = aux

         call trasponi_yx(aux,aux_t3,-1)

         psiy = aux

!-----------------------------------------------------

      call nvtxEndRange() 
      return        
      end  
