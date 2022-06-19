       subroutine out_inv

       include 'par.inc'


       if( mpime == root ) then
         write(19,110) time, Emag, Eke, Ekp, Epe, Etpare, Etperpe, Etot
         write(29,129) time, dint_drhsf_dz,dint_drhsu_dz 
       end if

110    format(1x, 8(1x, 1e22.9))
129     format(1x, 3(1x, 1e22.9))

       return
       end

