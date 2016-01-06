 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_dft_grid
    use m_definitions
    use m_atoms
    use m_shells
    use m_basis
    use m_rysquad,only : stiel


    type grid_type
        integer                  ::  itype
        integer                  ::  ngrid_type
        real(dp),allocatable    ::  xt(:,:) !size = (3,ngird_type)
        real(dp),allocatable    ::  wt(:)   !size = (ngird_type)
    end type grid_type

    !For atomic type
    type grid_set
        integer                        ::  ngt ! ng=ntype
        type(grid_type),allocatable   ::  gts(:) 
    end type grid_set

    !For all atoms
    type dft_grid
        integer                        ::  ngt    !ng=natom
        type(grid_type),allocatable   ::  gts(:) 
    end type dft_grid

    contains
!====================================================================
    !¡°log-squared¡± quadrature 
    ! see Gill, P. M. W.; Chien, S. -H. J Comput Chem 2003, 24, 732.
    subroutine log2_stiel(n,x,w)
       implicit none
       integer,intent(in) :: n
       real(dp),intent(out) :: x(1:n),w(1:n)
       !==
       call stiel(0.0_dp,1.0_dp,4.3_dp,func,x,w)
       contains
       function func(x,del)
            real(dp),intent(in) :: x,del
            real(dp) :: func
            real(dp) :: eps = 1.0d-15
            !Attention : here eps = 1.0d-15
            if(abs(x) < eps) then
                func =  1.1929d3
            else
                func = log(x)*log(x)
            end if
        end function func
    end subroutine log2_stiel
!====================================================================
!
!see Chien S H, Gill P M W. J. Comput. Chem., 2006, 27(6):730-739.
!SG-0: A small standard grid for DFT quadrature on large systems 
    subroutine setup_grid_sg_0(mol,grids)
        implicit none
        type(geom),intent(in)           ::  mol
        type(grid_set),intent(inout)    ::  grids
        !==
        integer     :: itype
        grids%ngt = mol%ntype
        allocate( grids%gts( grids%ngt ) )
        do itype = 1,grids%ngt
            call set_sg_0( NINT(mol%ztype(itype)),grids%gts(itype))
        end do
    end subroutine setup_grid_sg_0
    !
    !=====
    subroutine set_sg_0(zatom,gridt)
        implicit none
        integer,intent(in)                ::    zatom
        type(grid_type),intent(out)        ::    gridt  
        !===
        real(dp),save               ::  x_log2_23(23),w_log2_23(23)
        real(dp),save               ::  x_log2_26(26),w_log2_26(26)
        logical,save               ::  init = .true.
        !
        real(dp)                     ::  x_log2_23_tmp(23),w_log2_23_tmp(23)
        real(dp)                     ::  x_log2_26_tmp(26),w_log2_26_tmp(26)
        integer,pointer            ::  ld_indx(:,:)
        integer                     ::  ii,jj,ngrid_type,ngrid_rad
        integer                     ::  il,ir,ij,coun
        real(dp)                     ::  rsf
        !==
        if(init) then
            call log2_stiel(23,x_log2_23,w_log2_23)
            w_log2_23(:) = w_log2_23(:) * 4.0_dp * pi
            call log2_stiel(26,x_log2_26,w_log2_26)
            w_log2_26(:) = w_log2_26(:) * 4.0_dp * pi
            init = .false.
        end if
        ld_indx => lebedev_partition_sg_0(zatom)
        rsf = radsf(zatom)
        ngrid_type = sum( ld_indx(:,1) * ld_indx(:,2) )
        gridt%ngrid_type = ngrid_type
        allocate( gridt%xt(3,ngrid_type), gridt%wt(ngrid_type) )
        ngrid_rad = sum(ld_indx(:,2))
        select case(ngrid_rad)
            case(23)
                    !ri = - R*ln(xi)
                    !wi = -(ai/xi)*R^3
                    !FIXME{ x_log2_23(:) /= x_log2_23(1:23)}
                    x_log2_23_tmp(:) = -rsf * log(x_log2_23(:))
                    w_log2_23_tmp(:) = rsf**3 * (w_log2_23(:) / x_log2_23(:))
                    coun = 23
                    il = 1
                    ir = 0
                    do ii = 1,size(ld_indx,1)
                        do jj = 1,ld_indx(ii,2)
                            ir = ir + ld_indx(ii,1)
                            call lebedev_x(gridt%xt(:,il:ir),gridt%wt(il:ir),ld_indx(ii,1))
                            gridt%xt(:,il:ir) = gridt%xt(:,il:ir) * x_log2_23_tmp(coun)
                            gridt%wt(il:ir) = gridt%wt(il:ir) * w_log2_23_tmp(coun)
                            il = ir + 1
                            coun = coun - 1
                        end do
                    end do
            case(26)
                    !ri = - R*ln(xi)
                    !wi = -(ai/xi)*R^3
                    x_log2_26_tmp(:) = -rsf * log(x_log2_26(:))
                    w_log2_26_tmp(:) = rsf**3 * (w_log2_26(:) / x_log2_26(:))
                    coun = 26
                    il = 1
                    ir = 0
                    do ii = 1,size(ld_indx,1)
                        do jj = 1,ld_indx(ii,2)
                            ir = ir + ld_indx(ii,1)
                            call lebedev_x(gridt%xt(:,il:ir),gridt%wt(il:ir),ld_indx(ii,1))
                            gridt%xt(:,il:ir) = gridt%xt(:,il:ir) * x_log2_26_tmp(coun)
                            gridt%wt(il:ir) = gridt%wt(il:ir) * w_log2_26_tmp(coun)
                            il = ir + 1
                            coun = coun - 1
                        end do
                    end do
            case default
                    write(stdout,*) 'Ngrid_rad: ',ngrid_rad
                    stop 'SG0 : Radial Points Error '
        end select
        !================================
    end subroutine set_sg_0
    !
    !
    !The x-point Lebedev grid
    subroutine lebedev_x(roots,weights,npoint)
        implicit none
        real(dp),intent(out)    ::  roots(:,:),weights(:)
        integer,intent(in)     ::  npoint
        !===
        integer     ::  ndum
        ndum = assert_eq(size(roots,2),size(weights),npoint,'arg lebedev_x')
        select case(npoint)
            case(6)
                call ld0006(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(14)
                call ld0014(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(26)
                call ld0026(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(38)
                call ld0038(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(50)
                call ld0050(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(74)
                call ld0074(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(86)
                call ld0086(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(110)
                call ld0110(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(146)
                call ld0146(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(170)
                call ld0170(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case(230)
                call ld0230(roots(1,:),roots(2,:),roots(3,:),weights(:),ndum)
            case default
                 write(stdout,*) 'grid points: ',npoint
                 stop 'Lebedev grid is not available'
        end select
    end subroutine lebedev_x
    !
    !Lebedev Partition = x^y
    !The notation xy indicates that the x-point Lebedev grid is used at y successive radial points
    function lebedev_partition_sg_0(zatom)
        implicit none
        integer,intent(in)      ::  zatom
        !===
        integer,pointer         :: lebedev_partition_sg_0(:,:)
        !FIXME,here only few element
        integer,target,save     ::  lp_H_xy(11,2)
        integer,target,save     ::  lp_Li_xy(11,2)
        integer,target,save     ::  lp_Be_xy(12,2)
        integer,target,save     ::  lp_B_xy(7,2)
        integer,target,save     ::  lp_C_xy(13,2)
        integer,target,save     ::  lp_N_xy(10,2)
        integer,target,save     ::  lp_O_xy(11,2)
        integer,target,save     ::  lp_F_xy(10,2)
        integer,target,save     ::  lp_Na_xy(8,2)
        integer,target,save     ::  lp_Mg_xy(13,2)
        integer,target,save     ::  lp_Al_xy(15,2)
        integer,target,save     ::  lp_Si_xy(11,2)
        integer,target,save     ::  lp_P_xy(11,2)
        integer,target,save     ::  lp_S_xy(12,2)
        integer,target,save     ::  lp_Cl_xy(12,2)
        !
        logical,save             ::  init = .true.
        if(init) then
            lp_H_xy(:,1) = (/6, 14, 26, 38, 74, 110, 146, 86, 50, 38, 14/)
            lp_H_xy(:,2) = (/6, 3,  1,  1,  1,  1,   6,   1,  1,  1,  1/)

            lp_Li_xy(:,1) = (/6, 14, 26, 38, 74, 110, 146, 86, 50, 38, 14/)
            lp_Li_xy(:,2) = (/6, 3,  1,  1,  1,  1,   6,   1,  1,  1,  1/)

            lp_Be_xy(:,1) = (/6, 14, 26, 38, 74, 86, 110, 146, 50, 38, 14, 6/)
            lp_Be_xy(:,2) = (/4, 2,  1,  2,  1,  1,  2,   5,   1,  1,  1,  2/)

            lp_B_xy(:,1) = (/6, 26, 38, 86, 146, 38, 6/)
            lp_B_xy(:,2) = (/4, 4,  3,  3,  6,   1,  2/)

            lp_C_xy(:,1) = (/6, 14, 26, 38, 50, 86, 110, 146, 170, 146, 86, 38, 14/)
            lp_C_xy(:,2) = (/6, 2,  1,  2,  2,  1,  1,   1,   2,   2,   1,  1,  1/)

            lp_N_xy(:,1) = (/6, 14, 26, 38, 74, 110, 170, 146, 86, 50/)
            lp_N_xy(:,2) = (/6, 3,  1,  2,  2,  1,   2,   3,   1,  2/)

            lp_O_xy(:,1) = (/6, 14, 26, 38, 50, 86, 110, 86, 50, 38, 6/)
            lp_O_xy(:,2) = (/5, 1,  2,  1,  4,  1,  5,   1,  1,  1,  1/)

            lp_F_xy(:,1) = (/6, 38, 50, 74, 110, 146, 110, 86, 50, 6/)
            lp_F_xy(:,2) = (/4, 2,  4,  2,  2,   2,   2,   3,  1,  1/)

            lp_Na_xy(:,1) = (/6, 14, 26, 38, 50, 110, 74, 6/)
            lp_Na_xy(:,2) = (/6, 2,  3,  1,  2,  8,   2,  2/)

            lp_Mg_xy(:,1) = (/6, 14, 26, 38, 50, 74, 110, 146, 110, 86, 38, 14, 6/)
            lp_Mg_xy(:,2) = (/5, 2,  2,  2,  2,  1,  2,   4,   1,   1,  2,  1,  1/)

            lp_Al_xy(:,1) = (/6, 14, 26, 38, 50, 74, 86, 146, 170, 110, 86, 74, 26, 14, 6/)
            lp_Al_xy(:,2) = (/6, 2,  1,  2,  2,  1,  1,  2,   2,   2,   1,  1,  1,  1,  1/)

            lp_Si_xy(:,1) = (/6, 14, 38, 50, 74, 110, 146, 170, 86, 50, 6/)
            lp_Si_xy(:,2) = (/5, 4,  4,  3,  1,  2,   1,   3,   1,  1,  1/)

            lp_P_xy(:,1) = (/6, 14, 38, 50, 74, 110, 146, 170, 86, 50, 6/)
            lp_P_xy(:,2) = (/5, 4,  4,  3,  1,  2,   1,   3,   1,  1,  1/)

            lp_S_xy(:,1) = (/6, 14, 26, 38, 50, 74, 110, 170, 146, 110, 50, 6/)
            lp_S_xy(:,2) = (/4, 1,  8,  2,  1,  2,  1,   3,   1,   1,   1,  1/)

            lp_Cl_xy(:,1) = (/6, 14, 26, 38, 50, 74, 110, 170, 146, 110, 50, 6/)
            lp_Cl_xy(:,2) = (/4, 7,  2,  2,  1,  1,  2,   3,   1,   1,   1,  1/)
    
            init = .false.
        end if
        select case(zatom)
            case(1)
                lebedev_partition_sg_0 => lp_H_xy
            case(3)
                lebedev_partition_sg_0 => lp_Li_xy
            case(4)
                lebedev_partition_sg_0 => lp_Be_xy
            case(5)
                lebedev_partition_sg_0 => lp_B_xy
            case(6)
                lebedev_partition_sg_0 => lp_C_xy
            case(7)
                lebedev_partition_sg_0 => lp_N_xy
            case(8)
                lebedev_partition_sg_0 => lp_O_xy
            case(9)
                lebedev_partition_sg_0 => lp_F_xy
            case(11)
                lebedev_partition_sg_0 => lp_Na_xy
            case(12)
                lebedev_partition_sg_0 => lp_Mg_xy
            case(13)
                lebedev_partition_sg_0 => lp_Al_xy
            case(14)
                lebedev_partition_sg_0 => lp_Si_xy
            case(15)
                lebedev_partition_sg_0 => lp_P_xy
            case(16)
                lebedev_partition_sg_0 => lp_S_xy
            case(17)
                lebedev_partition_sg_0 => lp_Cl_xy
            case default
                 write(stdout,*) 'Element : ',element_name(real(zatom,dp))
                 stop 'Lebedev partition is not available'
        end select
    end function lebedev_partition_sg_0
    !
    !Radial Scale Factors
    function radsf(zatom)
        implicit none
        integer,intent(in)  ::  zatom
        !===
        real(dp)        :: radsf
        select case(zatom)
            case(1)
               radsf = 1.30_dp
            case(3)
                radsf = 1.95_dp
            case(4)
                radsf = 2.20_dp
            case(5)
                radsf = 1.45_dp
            case(6)
                radsf = 1.20_dp
            case(7)
                radsf = 1.10_dp
            case(8)
                radsf = 1.10_dp
            case(9)
                radsf = 1.20_dp
            case(11)
                radsf = 2.30_dp
            case(12)
                radsf = 2.20_dp
            case(13)
                radsf = 2.10_dp
            case(14)
                radsf = 1.30_dp
            case(15)
                radsf = 1.30_dp
            case(16)
                radsf = 1.10_dp
            case(17)
                radsf = 1.45_dp
            case default
                 write(stdout,*) 'Element : ',element_name(real(zatom,dp))
                 stop 'Radial scale factor is not available'
        end select
    end function radsf
!====================================================================
    subroutine destroy_grid_set(grids)
        implicit none
        type(grid_set),intent(inout)    :: grids
        !===
        integer         ::  i
        if(allocated(grids%gts)) then
            do i=1,size(grids%gts)
                call destroy_grid_type(grids%gts(i))
            end do
        end if
    end subroutine destroy_grid_set
    !
    subroutine destroy_dft_grid(dgrid)
        implicit none
        type(dft_grid),intent(inout)    :: dgrid
        !===
        integer         ::  i
        if(allocated(dgrid%gts)) then
            do i=1,size(dgrid%gts)
                call destroy_grid_type(dgrid%gts(i))
            end do
        end if
    end subroutine destroy_dft_grid
    !
    subroutine destroy_grid_type(gt)
        implicit none
        type(grid_type),intent(inout) :: gt
        !==
        if(allocated(gt%xt)) deallocate(gt%xt)
        if(allocated(gt%wt)) deallocate(gt%wt)
    end subroutine destroy_grid_type
!====================================================================
!Partitionning scheme of Axel Becke, J. Chem. Phys. 88, 2547 (1988).
    subroutine setup_partition(mol,dgrid,print_)
        implicit none
        type(geom),intent(in)         :: mol
        type(dft_grid),intent(inout)     :: dgrid
        logical,optional,intent(in)  :: print_
        !==
        real(dp),parameter      :: threshold = 1.0e-10
        integer                 :: i,itype,iatom,coun
        integer                 :: natom,ngrid_type
        type(grid_set)          :: grids_type
        real(dp)                 :: s_becke(mol%natom,mol%natom)
        real(dp)                 :: tmp(mol%natom)
        logical                 :: non_diagtri(mol%natom,mol%natom)
        real(dp),allocatable    :: r_grid(:,:),w_grid(:)
        !===
        call setup_grid_sg_0(mol,grids_type)

        !
        natom = mol%natom
        non_diagtri(:,:) = .true.
        forall(i=1:natom) non_diagtri(i,i) = .false.
        if(present(print_)) then
            write(stdout,'(/)')
            write(stdout,*),"=========== Grid Data per Atom (Initiazation Parameter) ==========="
            write(stdout,'(2x,a8,5x,a8)'),"Element","Ngrid"
            coun = 0
            do iatom = 1,natom
                 itype = mol%iatom_to_type(iatom)
                 write(stdout,'(5x,a2,8x,i8)'),element_name(real(mol%zatom(iatom),dp)),grids_type%gts(itype)%ngrid_type
                 coun = coun + grids_type%gts(itype)%ngrid_type
            end do
            write(stdout,'(/,a30,5x,i8)'),"Total Number of Grid Points:",coun
        end if
        dgrid%ngt = natom
        allocate( dgrid%gts( natom ) )
        do iatom = 1,natom
            itype = mol%iatom_to_type(iatom)
            ngrid_type = grids_type%gts(itype)%ngrid_type
            allocate(  r_grid(3,ngrid_type),               &
                      & w_grid(ngrid_type) ,                &
                      & dgrid%gts(iatom)%xt(3,ngrid_type),  &
                      & dgrid%gts(iatom)%wt(ngrid_type)     )
            r_grid(:,:) = grids_type%gts(itype)%xt(:,:)
            forall(i=1:ngrid_type) r_grid(:,i) = r_grid(:,i) + mol%x(:,iatom)
            w_grid(:)  = grids_type%gts(itype)%wt(:)
            do i = 1,ngrid_type
                call cutoff_function(r_grid(:,i),w_grid(i),iatom)
            end do
            dgrid%gts(iatom)%ngrid_type = ngrid_type
            dgrid%gts(iatom)%xt(:,:) = r_grid(:,:)
            dgrid%gts(iatom)%wt(:)   = w_grid(:)
            !!!!
            !do i = 1,ngrid_type
            !    write(stdout,'(2f16.10)'),w_grid(i),grids_type%gts(itype)%wt(i)
            !end do
            deallocate(r_grid,w_grid)
        end do
        call destroy_grid_set( grids_type )
        contains
        !==================
        subroutine cutoff_function(xt,wt,iatom_)
            implicit none
            real(dp),intent(in)              ::  xt(:)
            real(dp),intent(inout)           ::  wt
            integer,intent(in)              ::  iatom_
            !character(len = 50),intent(in)  ::  scheme
            !==
            integer       ::  ia,ja,it
            !upper section
            forall(ia=1:natom,ja=1:natom,(ia<ja .and. ia/=ja)) 
                s_becke(ia,ja) = (  norm2(xt(:)-mol%x(:,ia)) - norm2(xt(:)-mol%x(:,ja)) ) / mol%dist(ia,ja)
            end forall
            !odd symmetry
            forall(ia=1:natom,ja=1:natom,(ia<ja .and. ia/=ja)) 
                s_becke(ja,ia) = -s_becke(ia,ja)
            end forall
            where(non_diagtri) 
                s_becke = 0.5_dp * (1.0_dp - smooth_step(smooth_step(smooth_step(s_becke))))
            end where
            tmp(:) = product(s_becke,dim=2,mask = non_diagtri)
            wt = wt * tmp(iatom_) / sum(tmp)
            !write(stdout,'(f16.10)'),tmp(iatom_) / sum(tmp)
        end subroutine cutoff_function
        !
        pure function smooth_step(mu)
            implicit none
            real(dp),intent(in)     ::  mu(:,:)
            !==
            real(dp) :: smooth_step(size(mu,1),size(mu,2))
            smooth_step = 0.5_dp * mu * (3.0_dp - mu * mu)
        end  function smooth_step
    end subroutine setup_partition
!====================================================================

end module m_dft_grid
