      subroutine corrente_res
      use nvtx
c     This routine gives the laplacian of the magnetic flux function
      
      include 'par.inc'   


      dimension q1(ny),F00(nx),q3(ny),q5(ny)
      dimension aux_aa(nx1),aux_cc(nx)
      dimension aux_bb(nx1),aux_f1(nx)
      dimension aux_ww(nx-2)
      dimension q2(nx)
      dimension d2psidy2(nx,nyl,nzl),d2psidx2(nx,nyl,nzl)
      dimension aux(nx,nyl,nzl), aux_t(nxl,ny,nzl)
      dimension aux2(nx,nyl,nzl), aux_t2(nxl,ny,nzl)

      call nvtxStartRange('corrente_res',9)
      aux = psi

      call trasponi_yx(aux,aux_t,1)

      do iz=1,nzl
         do ix = 1,nxl
            do iy = 1,ny
               q1(iy) = aux_t(ix,iy,iz)
            enddo
            call der2y(q1, q3)
!            call der1y(q1, q5)
            do iy = 1,ny
               aux_t(ix,iy,iz) = q3(iy)
!               aux_t2(ix,iy,iz) = q5(iy)
            enddo
         enddo
      enddo

      call trasponi_yx(aux,aux_t,-1)    

      d2psidy2 = aux

      aux = psi

      call trasponi_yx(aux,aux_t,1)

      do iz=1,nzl
         do ix = 1,nxl
            do iy = 1,ny
               q1(iy) = aux_t(ix,iy,iz)
            enddo
            call der1y(q1, q3)
            do iy = 1,ny
               aux_t(ix,iy,iz) = q3(iy)
            enddo
         enddo
      enddo

      call trasponi_yx(aux,aux_t,-1)

!      call trasponi_yx(aux2,aux_t2,-1)    

      psiy = aux
      hpy = psiy

      do iz=1,nzl
         do iy = 1,nyl
            do ix = 1 ,nx
               F00(ix) = psi(ix,iy,iz)
            enddo
            do ix = 2, nx1
               q2(ix) = aa_2(ix) * F00(ix-1) + bb_2(ix) * F00(ix)
     &              + cc_2(ix) * F00(ix+1)
            enddo

            q2(1)  = aasecder1 * F00(1) + bbsecder1 * F00(2)
     &             + ccsecder1 * F00(3) + ddsecder1 * F00(4)
            q2(nx) = aasecdernx * F00(nx) + bbsecdernx * F00(nx-1)
     &             + ccsecdernx * F00(nx-2) + ddsecdernx * F00(nx-3)
     
            aux_aa = fact_alfa_2
            aux_bb = fact_beta_2
            aux_cc = fact_gamma_2
            aux_f1 = q2
            aux_ww = fact_ww_2 
     
            ipv = ipv2

            CALL  DGTTRS(TRANS,nx,NRHS,aux_aa,aux_cc,aux_bb, 
     &         aux_ww,ipv,aux_f1,LDB,INFO)
            if (info > 0 .or. info < 0) then
               write(*,*) 'Problemi soluzione, info:', info
               stop
            endif
            do ix = 1,nx
               d2psidx2(ix,iy,iz) = aux_f1(ix)
            enddo
         enddo
      enddo

      cur = - d2psidx2 - d2psidy2

      call nvtxEndRange()
      return
      end
