        implicit real*8(a-h,o-z)

        parameter(nx=1024, ny=96, nz=96, nx1 = nx - 1)

        parameter(nprocx=4096)

        CHARACTER*1 TRANS

        integer nzl ! number of z planes per processor
        integer nzg ! number of global z planes = nz
        integer izlg ! global index of the local first plane
        integer nyl, nyg, iylg
        integer nzp ! number of z local planes per processor as a function of the processor
        integer nyp ! number of y local planes per processor as a function of the processor
        integer izlgp ! global index of the local first plane as a function of the processor
        integer iylgp
        integer nxl,nxl2
        common/dimensions/nzl, nyl, nxl, nxl2, nzg, nyg, izlg, iylg
     &   , nzp(nprocx), nyp(nprocx), izlgp(nprocx), iylgp(nprocx)


        common/boite1/istart,nstep,ioutt,ioutf,mm,inc,jump,irestart_zero  
     &     ,kbc_p,kbc_h  
        common/boite2/xl,yl,zl,dy,dz,pi,grady,gradz,
     &   c0yz,c0xz,c0xy,dyz

        common/boite3/dt,time,oldtime,omega,omegax,alfa1,
     &   ggp,ggm

        common/boite4/de,rhos,tau,eta,beta_e,ee,cc,sk1,sk,zmxt,zmyt
        common/boite5/gpmin, gpmax,gp_inv,g_2inv,new_Cas3D
        common/boite6/gm1, gm2, gm3, gm4, gmmin, gmmax,gm_inv
        common/boite7/Emag, Eke, Ekp, Epe, Etpare, Etperpe,Etot
	common/boite8/t0,gridalfa,gridbeta,gridgamma
	common/boite9/de2,rhos2,pso,psoeq,phieq,eq_l,asym	
	common/boite10/front_l,front_r,delta_ps	

        common/LUint1/nrhs,ndb,ndb_t,ldb,info,info_t
        common/LUint3/jpv(nx1), jpv_t(nx)

        common/solvpoiss2/alfa(nx1), beta(nx1)

        common/firstderG/aa_1_G(nx),bb_1_G(nx),cc_1_G(nx),
     &      d_1_G(nx),alfa_1_G(nx-1),beta_1_G(nx-1),
     &      aux_alfa_1_G(nx-1),aux_beta_1_G(nx-1),
     &      aux_gamma_1_G(nx),ww_1_G(nx-2),ipv_1_G(nx),
     &      aux_alfa_2_G(nx),aux_beta_2_G(nx)
     
        common/secderivx1/alfa_2(nx1),beta_2(nx1)
        common/secderivx2/fact_alfa_2(nx1),ipv2(nx)
        common/secderivx5/fact_ww_2(nx-2)
        common/secderivx3/fact_beta_2(nx1),fact_gamma_2(nx)
        common/secderivx4/aa_2(nx1),bb_2(nx1),cc_2(nx1)
        common/secderivx5/aasecder1,bbsecder1,ccsecder1,ddsecder1
        common/secderivx6/aasecdernx,bbsecdernx,ccsecdernx,ddsecdernx

        common/filtro1/a_filt, b_filt, c_filt, a11_filt, a12_filt 
     &        , a13_filt, a14_filt, a15_filt, a21_filt, a22_filt 
     &        , a23_filt, a24_filt, a25_filt
        common/filtro2/filt_d(nx1), filt_p(nx1-1), filt_m(nx1-1)
     &             , filt_2(nx1-2), filt_2_t(nx-2) 
     &           ,  filt_d_t(nx), filt_p_t(nx-1), filt_m_t(nx-1)

        common/Fields1/x(nx),y(ny),z(nz),yk(ny/2+1), zk(nz/2+1)
        common/wave_numb/ky(ny),kz(nz)

        common/Fields3/work4(ny/2),work5(nz/2)
        common/Fields5/wsavey(2*ny+15),wsavez(2*nz+15)

        common/Fields7/dx(nx),dxy(nx),dxz(nx),dxyz(nx)
	common/caratteri/TRANS
        

        integer mpime, nproc, root, group,nprocy,nprocz
        integer nprow, npcol, myrow, mycol, iprow, ipcol
        common/parallel/mpime, nproc, root, group, nprow, npcol
     &       ,nprocy,nprocz
     &       , myrow, mycol, iprow(nprocx), ipcol(nprocx)


	real*8, pointer :: phi(:,:,:), psi(:,:,:)
        real*8, pointer :: cur(:,:,:), uu(:,:,:)
        real*8, pointer :: cur_old(:,:,:), dpsi_dt(:,:,:)
        real*8, pointer :: hpy(:,:,:), hmy(:,:,:)
        real*8, pointer :: psiy(:,:,:), phiy(:,:,:)
        real*8, pointer :: e_para(:,:,:),abs_b(:,:,:)
        real*8, pointer :: d2_uu(:,:,:)

        common/Fields/phi, psi, cur, uu, cur_old, dpsi_dt,
     &   hpy, hmy, psiy, phiy, e_para, abs_b,d2_uu 


        real*8, pointer :: fact_aa(:,:),fact_bb(:,:)
        real*8, pointer :: fact_cc(:,:),fact_ww(:,:)
        common/solvpoiss1/fact_aa,fact_bb,fact_cc,fact_ww

        real*8, pointer :: fact_aa_1(:,:),fact_bb_1(:,:)
        real*8, pointer :: fact_cc_1(:,:),fact_ww_1(:,:)
        common/solvhelm1/fact_aa_1,fact_bb_1,fact_cc_1,fact_ww_1


        real*8, pointer :: ipv_my(:,:),ipv_my_1(:,:)
        common/LUint2/ipv(nx),ipv_my,ipv_my_1


!        real*8, allocatable :: phi(:,:,:), psi(:,:,:)
!        real*8, allocatable :: cur(:,:,:), uu(:,:,:)
!        real*8, allocatable :: cur_old(:,:,:)
