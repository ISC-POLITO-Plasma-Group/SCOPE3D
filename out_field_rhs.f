       subroutine out_field_rhs(flxp1,flxm1)

       include 'par.inc'

       real*8 :: flxp1(nx,nyl,nzl),flxm1(nx,nyl,nzl)

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Flux_F.dat',form='unformatted')
       end if

       call outwrite(12,flxp1)

       if( mpime == root ) then
	  close (unit=12)
       end if

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Flux_U.dat',form='unformatted')
       end if

       call outwrite(12,flxm1)

       if( mpime == root ) then
	  close (unit=12)
       end if

        
       return
       end

