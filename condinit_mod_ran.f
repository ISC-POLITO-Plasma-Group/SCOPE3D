      subroutine condinit_mod_ran(hp1,hm1,hp1_pascal)

c  Condizione iniziale caso random

      include 'par.inc'


      real*8 :: hp1(nx,nyl,nzl),hm1(nx,nyl,nzl)
      dimension aux_der1x(nx),F0(nx)
      double precision xs

      real*8 :: ampl,ph_diff1,ph_diff2
      real(4) :: dran0
      integer*4 :: idum

c  Initial state

        idum = 164235

        sk1 = 2.0 * pi / yl

	IF (istart.EQ.0) THEN

        cur = 0.0d0

        ampl = 0.0
        ph_diff1 = 0.0
        ph_diff2 = 0.0  

        do i_m = 1,1

        do i_n = 1,1
        
         if( mpime == root ) then

               ampl = 1.0

               CALL rann(dran0,idum)
               ph_diff1 = dran0*pi
               CALL rann(dran0,idum)
               ph_diff2 = dran0*pi

              write(*,*) 'ph_diff ROOT',ph_diff1,phi_diff2,i_m,i_n,mpime

         endif

         call bcast_real(ampl,1,root)

         call bcast_real(ph_diff1,1,root)
         call bcast_real(ph_diff2,1,root)

!         write(*,*) 'ph_diff',ph_diff1,ph_diff2,i_m,i_n,mpime


           if( mpime == root ) then

c**************IF PSI_EQ NE 1/cosh^2*****************              
             hel = ((float(i_n-1)/float(i_m))*(yl/zl)) / psoeq * eq_l
!!!!!!!!!!!!!!!!!Bz<0 (PegSchep MODEL) !!!!!!!!!!!!!!!!!!!!
             xs = -dlog((1.0d0-hel)/(1.0d0+hel))/2.*eq_l 
!!!!!!!!!!!!!!!!!Bz>0!!!!!!!!!!!!!!!!!!!!
!             xs = -xs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
             write(*,*) 'XS',xs,hel,mpime
c*****************************************************

c**************IF ASYM NE 0**************************
             if (asym.ne.0) then
                hel = ((float(1)/float(1))*(yl/zl)) / psoeq
                xs = -dlog((1.0d0-hel)/(1.0d0+hel))/2.
                write(*,*) 'XS',xs,mpime
             endif
c****************************************************

           endif

           call bcast_real(xs,1,root)

!           write(*,*) 'XS',xs,mpime


           if (de.eq.0) then
              delta= (eta/sk1)**(1.0d0/3.0d0)
              delta2=delta*delta
           else
              delta=de
              delta2=de2
           endif

           do iy = 1, nyl

              iyg = iy + iylg - 1 ! costruisco l'indice globale

              yy = y(iyg)

              do iz = 1, nzl
                 izg = iz + izlg - 1 ! costruisco l'indice globale

                 zz = z(izg)    

                 do ix = 1, nx
                    xx = x(ix)

                    eexx_1 = 2.0 * abs(xx-xs) / delta

                    if (eexx_1 .LT. 1.0d0)  then

                       sky = xl * i_m * pi / yl
                       skz = xl * (i_n-1) * pi / zl

                       cur(ix,iy,iz) = cur(ix,iy,iz) + 
     &               pso*sqrt(2.0d0/pi/delta2)* 
     &               cos(sky*yy+skz*zz+ph_diff1)*
     &                      (exp(-((x(ix)-xs)**2.0d0)/
     &                      (4.0*delta2)))


                     else
                       
                          cur(ix,iy,iz) = cur(ix,iy,iz)
                       
                     endif

                 enddo
              enddo
           enddo
      
           enddo
           enddo
           
           uu = cur
           call poisson_mod
           
           do iz = 1,nzl
              do iy = 1, nyl
                 do ix = 1, nx
                    psi(ix,iy,iz) = - phi(ix,iy,iz)
                    psiy(ix,iy,iz) = - phiy(ix,iy,iz)                    


                 enddo
              enddo
           enddo
           
           phi=0.0d0
           uu=0.0d0  
           hp1 = psi + de2 * cur
           hm1 = uu


        ELSEIF(istart.EQ.1) THEN
           
           call gpm_read(hp1, hm1,10)
   
!           hp1 = 0.01 * hp1
!           hm1 = 0.01 * hm1
!           phi = 0.01 * phi
!           psi = 0.01 * psi

           if(irestart_zero.eq.1) then 
              oldtime =0.0d0
	      call out_field
           endif
           
           
        ENDIF 
        
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
 
