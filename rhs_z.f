      subroutine rhs_z(rhshpz,rhshmz)
      use nvtx	
      include 'par.inc'                               
      dimension q4(nx,nyl,nzl),q4_t(nxl2,nyl,nz)
      dimension q5(nz),q6(nz)
      real*8 :: rhshpz(nx,nyl,nzl)
      real*8 :: rhshmz(nx,nyl,nzl)
      real*8 zi

!      q4 = -(phi - rhos2 * uu)  Bz = 1
!      q4 = phi - rhos2 * uu 

      call nvtxStartRange('rhs_z.f',13)
      if (nz .ge. 4) then

      q4 = phi  
       write(*,*),'before trasponi_zx rhs_z 1'

      call trasponi_zx(q4,q4_t,1)
       write(*,*),'after trasponi_zx rhs_z 1'

      do j=1,nyl
         do i = 1, nxl2
            do k=1,nz
               q5(k) = q4_t(i,j,k)
            enddo

!***********************************************
      call nvtxStartRange('der1z',22)
      call drfftf(nz, q5, wsavez)

      do k = 2, nz-1, 2
         zi    = gradz * float(k) / float(2 * nz)
         q6(k)   = - q5(k+1) * zi 
         q6(k+1) =   q5(k)  * zi 
      enddo

      q6(1)  = 0.0
      q6(nz) = 0.0

      call drfftb(nz, q6, wsavez)
      call nvtxEndRange()
!***********************************************

            do k = 1,nz
               q4_t(i,j,k) = q6(k)
            enddo
         enddo
      enddo

      call trasponi_zx(q4,q4_t,-1)

      e_para = e_para-q4

      rhshpz = q4

      q4 = -0.0*uu

      call trasponi_zx(q4,q4_t,1)

      do j=1,nyl
         do i = 1, nxl2
            do k=1,nz
               q5(k) = q4_t(i,j,k)
            enddo

!***********************************************
      
      call nvtxStartRange('der1z',22)
      call drfftf(nz, q5, wsavez)

      do k = 2, nz-1, 2
         zi    = gradz * float(k) / float(2 * nz)
         q6(k)   = - q5(k+1) * zi
         q6(k+1) =   q5(k)  * zi
      enddo

      q6(1)  = 0.0
      q6(nz) = 0.0

      call drfftb(nz, q6, wsavez)
      call nvtxEndRange()
!***********************************************

            do k = 1,nz
               q4_t(i,j,k) = q6(k)
            enddo
         enddo
      enddo

      call trasponi_zx(q4,q4_t,-1)

      rhshpz = rhshpz + q4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      q4 = cur  Bz = 1
      q4 = -2./beta_e*cur 

      call trasponi_zx(q4,q4_t,1)

      do j=1,nyl
         do i = 1, nxl2
            do k=1,nz
               q5(k) = q4_t(i,j,k)
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call nvtxStartRange('der1z',22)            
      call drfftf(nz, q5, wsavez)

      do k = 2, nz-1, 2
         zi    = gradz * float(k) / float(2 * nz)
         q6(k)   = - q5(k+1) * zi 
         q6(k+1) =   q5(k)  * zi 
      enddo

      q6(1)  = 0.0
      q6(nz) = 0.0

      call drfftb(nz, q6, wsavez)
      call nvtxEndRange()
!***********************************************

            do k = 1,nz
               q4_t(i,j,k) = q6(k)
            enddo
         enddo
      enddo

      call trasponi_zx(q4,q4_t,-1)

      rhshmz = q4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

     else

      rhshpz = 0.0
      rhshmz = 0.0

      endif
      call nvtxEndRange()
      return
      end


