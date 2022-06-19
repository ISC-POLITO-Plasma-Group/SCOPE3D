      subroutine B_func

c     This routine gives the laplacian of the magnetic flux function
      
      include 'par.inc'   


      dimension q1(ny),F00(nx),q3(ny)
      dimension aux_aa(nx1),aux_cc(nx)
      dimension aux_bb(nx1),aux_f1(nx)
      dimension aux_ww(nx-2)
      dimension q2(nx)
      dimension d2psidy2(nx,nyl,nzl),d2psidx2(nx,nyl,nzl)
      dimension aux(nx,nyl,nzl), aux_t(nxl,ny,nzl)


      aux = uu

      call trasponi_yx(aux,aux_t,1)

      do iz=1,nzl
         do ix = 1,nxl
            do iy = 1,ny
               q1(iy) = aux_t(ix,iy,iz)
            enddo
            call der2y(q1, q3)
            do iy = 1,ny
               aux_t(ix,iy,iz) = q3(iy)
            enddo
         enddo
      enddo

      call trasponi_yx(aux,aux_t,-1)    

      d2psidy2 = aux

      do iz=1,nzl
         do iy = 1,nyl
            do ix = 1 ,nx
               F00(ix) = uu(ix,iy,iz)
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

      d2_uu =  uu - tau*(d2psidx2 + d2psidy2)
    
       write(*,*) 'TAU, BETA', tau, beta_e,mpime

 
      return
      end
