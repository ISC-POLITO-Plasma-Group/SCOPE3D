          subroutine filterx_t(g_filt)
          use nvtx
          include 'par.inc'

          dimension g_filt(nx,nyl,nzl), filtrox(nx), TW(nx)

c  filtro 
          call nvtxStartRange('filterx_t',14)
          do k = 1,nzl          

          do j = 1, nyl

            do i = 1, nx
              TW(i) = g_filt(i,j,k)
            enddo
 
! **    x Filtering

!     Soluzione del sistema lineare tridiagonale fattorizzato LU
!     Differenze finite 3x5

            do i = 3, nx-2
             filtrox(i) = a_filt * TW(i) + b_filt * (TW(i+1) 
     &              + TW(i-1)) + c_filt * (TW(i+2) + TW(i-2))
            enddo

!   Boundary conditions

        filtrox(1) = a11_filt * TW(1) + a12_filt * TW(2)
     &    + a13_filt * TW(3) + a14_filt * TW(4) + a15_filt * TW(5)

        filtrox(2) = a21_filt * TW(1) + a22_filt * TW(2)
     &    + a23_filt * TW(3) + a24_filt * TW(4) + a25_filt * TW(5)

        filtrox(nx-1) = a21_filt * TW(nx)   + a22_filt * TW(nx-1) 
     &                 + a23_filt * TW(nx-2) + a24_filt * TW(nx-3)
     &                 + a25_filt * TW(nx-4)

        filtrox(nx)   = a11_filt * TW(nx) + a12_filt * TW(nx-1) 
     &      + a13_filt * TW(nx-2) + a14_filt * TW(nx-3) 
     &                             + a15_filt * TW(nx-4)

            CALL  DGTTRS(TRANS, nx, NRHS, filt_m_t, filt_d_t
     &       , filt_p_t, filt_2_t,jpv_t, filtrox, NDB_t, INFO_t)

          if (info > 0 .or. info < 0) then
           write(*,*) 'Problemi soluzione, info:', info_t
           stop
          endif

! **    Filtered function

            do i = 1, nx
              g_filt(i,j,k) = filtrox(i)
            enddo
       
         enddo

         enddo
        call nvtxEndRange()
        return
        end 
