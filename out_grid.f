       subroutine out_grid

       include 'par.inc'

       if( mpime == root ) then
         write(18,110) x
       end if

110    format(1x, 1e12.5)

       return
       end

