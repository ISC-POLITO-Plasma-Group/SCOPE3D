      subroutine filter_yz_mod(g)


      include 'par.inc'   

      real*8 :: g(nx,nyl,nzl)
      dimension q3(nz),q4(nz)
      dimension q1(ny),q2(ny)
      dimension F0_my(nxl,ny)
      dimension F0_mz(nyl,nz)
      dimension aux(nx,nyl,nzl),aux_t(nxl,ny,nzl)
      dimension aux_t2(nxl2,nyl,nz)
      dimension out(nyl,nz)

      aux = g

      call trasponi_yx(aux,aux_t,1)

       write(*,*) 'AFTER TRASPONI_YX 1'

      do iz = 1, nzl
         do ix = 1,nxl
            do iy = 1,ny
               q1(iy) = aux_t(ix,iy,iz)
            enddo
            call drfftf(ny, q1, wsavey)
            do iy = 1,ny
               F0_my(ix,iy) = q1(iy)/ny
            enddo
         enddo
         aux_t(1:nxl,1:ny,iz) = F0_my(1:nxl,1:ny)
      enddo
      
      call trasponi_yx(aux,aux_t,-1)
       write(*,*) 'AFTER TRASPONI_YX -1'

      call trasponi_zx(aux,aux_t2,1)
       write(*,*) 'AFTER TRASPONI_ZX 1'
      
      do ix = 1, nxl2
         do iy = 1,nyl
            do iz = 1,nz
               q3(iz) = aux_t2(ix,iy,iz)
            enddo
            call drfftf(nz, q3, wsavez)
            do iz = 1,nz
               F0_mz(iy,iz) = q3(iz)/nz
            enddo
         enddo
         aux_t2(ix,1:nyl,1:nz) = F0_mz(1:nyl,1:nz)
      enddo      
      
      
      do ix = 1, nxl2
         
!     m = 1
         
         do iy = 2,nyl-1,2

!!!!!!BISOGNA COSTRUIRE INDICE GLOBALE y e fare attenzione agli indici!!!!!!
            
            iyg = iy + iylg - 1 

!     m = 1, n .ne. 1
            
            out(iy,1)=aux_t2(ix,iy,1)*work4(ky(iyg))
            out(iy+1,1)=aux_t2(ix,iy+1,1)*work4(ky(iyg))

         enddo

!     m = 1, n = 1
         
         if (mpime == root) then
            out(1,1)=aux_t2(ix,1,1)
         else
            iyg = 1 + iylg - 1 
            out(1,1)=aux_t2(ix,1,1)*work4(ky(iyg))
         endif

!     m = 1, n = nyl

         iyg_nyl = nyl + iylg - 1
         
          if (iyg_nyl .eq. ny) then
            out(nyl,1)=0.0
         else
            out(nyl,1)=aux_t2(ix,nyl,1)*work4(ky(iyg_nyl))
         endif

!     m .ne. 1
         
         do iy = 2, nyl-1, 2
            iyg = iy + iylg - 1
            do iz = 2, nz-1, 2
 
               out(iy,iz)     = aux_t2(ix,iy,iz)*work5(kz(iz))*
     &         work4(ky(iyg)) 
               out(iy+1,iz)   = aux_t2(ix,iy+1,iz)*work5(kz(iz))*
     &         work4(ky(iyg))  
               out(iy,iz+1)   = aux_t2(ix,iy,iz+1)*work5(kz(iz))*
     &         work4(ky(iyg))  
               out(iy+1,iz+1) = aux_t2(ix,iy+1,iz+1)*work5(kz(iz))*
     &         work4(ky(iyg))             
            enddo
         enddo

!     m .ne. 1, n = 1
         if (mpime == root) then
            do iz = 2, nz-1, 2
               out(1,iz)     = aux_t2(ix,1,iz)*work5(kz(iz))
               out(1,iz+1)   = aux_t2(ix,1,iz+1)*work5(kz(iz))   
            enddo
         else
            iyg = 1 + iylg - 1
            do iz = 2, nz-1, 2
               out(1,iz)     = aux_t2(ix,1,iz)*work5(kz(iz))*
     &         work4(ky(iyg))
               out(1,iz+1)   = aux_t2(ix,1,iz+1)*work5(kz(iz))*
     &         work4(ky(iyg))   
            enddo
         endif

!     m .ne. 1, n = nyl

         iyg = nyl + iylg - 1

         if (iyg .eq. ny) then

            do iz=1,nz 
               out(nyl,iz)=0.0
            enddo
            
         else 
         
            do iz = 2, nz-1, 2          
               out(nyl,iz)     = aux_t2(ix,nyl,iz)*work5(kz(iz))*
     &         work4(ky(iyg))
               out(nyl,iz+1)   = aux_t2(ix,nyl,iz+1)*work5(kz(iz))*
     &         work4(ky(iyg)) 
            enddo
            
         endif

         do iy = 1,nyl 
            out(iy,nz)=0.0
         enddo

         do iy = 1, nyl
            do iz = 1, nz
               q4(iz) = out(iy,iz)
            enddo
            call drfftb(nz, q4, wsavez)
            do iz = 1, nz
               aux_t2(ix,iy,iz) = q4(iz)
            enddo
         enddo
         
      enddo
         
      call trasponi_zx(aux,aux_t2,-1)
       write(*,*) 'AFTER TRASPONI_ZX -1'


      call trasponi_yx(aux,aux_t,1)
      
      do iz=1,nzl
         do ix=1,nxl
            do iy=1,ny
               q2(iy)= aux_t(ix,iy,iz)
            enddo
            call drfftb(ny, q2, wsavey)            
            do iy = 1, ny
               aux_t(ix,iy,iz) = q2(iy)
            enddo               
         enddo
      enddo
      
      call trasponi_yx(aux,aux_t,-1)

      g = aux
      
!--------------------------------------------------
      
      return
      end
