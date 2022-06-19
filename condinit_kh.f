      subroutine condinit_kh(hp1,hm1,hp1_pascal)
      use nvtx 
c  Condizione iniziale caso random

      include 'par.inc'


      real*8 :: hp1(nx,nyl,nzl),hm1(nx,nyl,nzl)
      real*8 :: hp1_pascal(nx,nyl,nzl)
      dimension dx2_pascal(nx)
      dimension aux_der1x(nx),F0(nx)
      double precision xs

      real*8 :: ampl,ph_diff1,ph_diff2
      real(4) :: dran0
      integer*4 :: idum

      call nvtxStartRange('condinit_kh',6)
c  Initial state

        idum = 164235

        sk1 = xl * pi / yl

	IF (istart.EQ.0) THEN

        cur = 0.0d0

        ampl = 0.0
        ph_diff1 = 0.0
        ph_diff2 = 0.0  

!        do i_m = 1,8
!!!!!!!!! SINGLE MODE  !!!!!!!!!!!!!
        do i_m = 1,16
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        do i_n = 2,2
        do i_n = 1,1
        
         if( mpime == root ) then

               ampl = 1.0

               CALL rann(dran0,idum)
               ph_diff1 = dran0*pi
               CALL rann(dran0,idum)
               ph_diff2 = dran0*pi

               write(*,*) 'ph_diff ROOT',ph_diff1,ph_diff2,i_m,i_n,mpime

         endif

         call bcast_real(ampl,1,root)

         call bcast_real(ph_diff1,1,root)
         call bcast_real(ph_diff2,1,root)

           do iy = 1, nyl

              iyg = iy + iylg - 1 ! costruisco l'indice globale

              yy = y(iyg)

              do iz = 1, nzl
                 izg = iz + izlg - 1 ! costruisco l'indice globale

                 zz = z(izg)    

                 do ix = 1, nx
                    xx = x(ix)

                       sky = xl * i_m * pi / yl
                       skz = xl * (i_n-1) * pi / zl

                       uu(ix,iy,iz) = uu(ix,iy,iz) + 
     &              pso*cos(sky*yy+skz*zz+ph_diff1)*
     &                (1.0d0 + tanh(x(ix))) / 
     &                (cosh(x(ix)))**2.0d0

                 enddo
              enddo
           enddo
      
           enddo
           enddo
 
           write(*,*) 'inizio poisson'
           
           call poisson_mod
           
           write(*,*) 'fine poisson'
          

           psi = 0.0d0
           cur = 0.0d0 
           hp1 = psi + de2 * cur
           hm1 = uu

 
        ELSEIF(istart.EQ.1) THEN
           
           call gpm_read(hp1, hm1, hp1_pascal, 10)
   
!           hp1 = 0.01 * hp1
!           hm1 = 0.01 * hm1
!           phi = 0.01 * phi
!           psi = 0.01 * psi

           if(irestart_zero.eq.1) then 
              oldtime =0.0d0
	      call out_field
           endif
           
           
        ENDIF 
      call nvtxEndRange()  
        return
        end
     
!************************************************************************
!*     subroutine RANN
!*     IN : idum (integer flag for generating random numbers)
!*     OUT: dran0 (dimension random number)
!************************************************************************
      subroutine rann(dran0,idum)
      implicit real*4 (d)
      integer*4 idum,ia,im,iq,ir,mask,m

      parameter(ia=16807,im=2147483647,dam=1./im)
      parameter(iq=127773,ir=2836,mask=123459876)

      idum=ieor(idum,mask)
      m=idum/iq
      idum=ia*(idum-m*iq)-ir*m
      if(idum.lt.0.) idum=idum+im
      dran0=dam*idum
      idum=ieor(idum,mask)
      
      return
      end
 
