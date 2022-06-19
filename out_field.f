       subroutine out_field
       use nvtx
       include 'par.inc'
      
       call nvtxStartRange('out_field',8)
       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Psi.dat',form='unformatted')
       end if

       call outwrite(12,psi)

       if( mpime == root ) then
	  close (unit=12)
       end if

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Phi.dat',form='unformatted')
       end if

       call outwrite(12,phi)

        write(*,*) 'SCRITTO PHI',myrank

       if( mpime == root ) then
	  close (unit=12)
       end if

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='U.dat',form='unformatted')
       end if

       call outwrite(12,uu)

       if( mpime == root ) then
	  close (unit=12)
       end if

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='J.dat',form='unformatted')
       end if

       call outwrite(12,cur)

       if( mpime == root ) then
	  close (unit=12)
       end if

       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Dpsidt.dat',form='unformatted')
       end if

       call outwrite(12,dpsi_dt)

       if( mpime == root ) then
	  close (unit=12)
       end if
       
       if( mpime == root ) then
          open(unit = 12,POSITION='append',status='unknown',
     &      file='Epara.dat',form='unformatted')
       end if

       call outwrite(12,e_para)

       if( mpime == root ) then
          close (unit=12)
       end if


       call nvtxEndRange()
       return
       end

