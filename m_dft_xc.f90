 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_dft_xc
    use m_definitions
    use m_atoms
    use m_shells
    use m_basis
    use m_dft_grid


    type gvalue_shell
        integer                            ::  ngs
        real(dp),allocatable              ::  gaussian_r(:) !(nv)
        real(dp),allocatable              ::  gaussian_grad(:,:) !(3,nv)
        real(dp),allocatable              ::  gaussian_lapl(:,:) !(3,nv)
    end type gvalue_shell

    type gvalue_basis
        integer                                      :: ngb
        type(gvalue_shell),allocatable              ::  gshells(:)
    end type gvalue_basis
    contains
!====================================================================
    subroutine calc_density_value(mol,shellset,basis,p_mat,xx,den_r)
        implicit none
        type(geom),intent(in)	::  mol
        type(shell_set),intent(in)	::  shellset
        type(basis_set),intent(in) 	::  basis
        real(dp),intent(in)		::  p_mat(:,:)
        real(dp),intent(in)		::  xx(3)
        real(dp),intent(out)        ::  den_r
        !==
        real(dp)                     ::  bfp_mat(size(p_mat,1),size(p_mat,2))
        integer			::  i,j,nc
        integer			::  i1_l,i1_r,i2_l,i2_r,nbf1,nbf2
        integer			::  ibasis1,ibasis2,ishell1,ishell2
        do i=1,basis%nbas
            do j=1,i
                ibasis1 = basis%bf(i)%ibas
                ibasis2 = basis%bf(j)%ibas
                ishell1 = basis%bf(i)%ishell
                ishell2 = basis%bf(j)%ishell
                nc = basis%basis_pair_data( ibasis1,ibasis2 )%nc12

                i1_l = basis%bf(i)%ibf_left
                i1_r = basis%bf(i)%ibf_right
                i2_l = basis%bf(j)%ibf_left
                i2_r = basis%bf(j)%ibf_right

                call calc_bfpair( basis%bf(i),&
				     basis%bf(j),&
                                   nc,&
				     shellset%shell_pair_data( ishell1,ishell2 ),&
				     basis%basis_pair_data( ibasis1,ibasis2 ),&
                                   mol,&
                                   xx,&
				     bfp_mat( i1_l : i1_r, i2_l : i2_r )	)
                if( i/=j ) then
                    bfp_mat( i2_l : i2_r, i1_l : i1_r ) = transpose( bfp_mat( i1_l : i1_r, i2_l : i2_r ) )
                end if
	    end do
	end do
        den_r = sum( p_mat * bfp_mat )
    end subroutine calc_density_value
    !===========================
    subroutine calc_bfpair(bf1,bf2,nc,shp,bfp,mol,xx,bf_r)
        implicit none
        type(basis_function),intent(in) ::  bf1,bf2
        type(basis_pair),intent(in)	    ::  bfp
        type(shell_pair),intent(in)	    ::  shp
        type(geom),intent(in)	    ::  mol
        integer,intent(in)             ::  nc
        real(dp),intent(in)		    ::  xx(:)
        real(dp),intent(out)            ::  bf_r(:,:)                               
        !==
        real(dp)         ::      zeta(nc),coeff12(nc,bfp%isp12)
        real(dp)         ::      fact(bfp%isp12)
        real(dp)         ::      vxa(3),vxb(3),vxp2(nc)
        real(dp)         ::      vxa_power(3,0:bf1%am),vxb_power(3,0:bf2%am)
        integer         ::      isp1,isp2,ibf1,ibf2
        integer         ::      amm1,amm2,i,j,ij
        integer         ::      imx,imy,jmx,jmy
        integer         ::      bf1_nx,bf1_ny,bf1_nz
        integer         ::      bf2_nx,bf2_ny,bf2_nz
        !

        vxa(:) = xx(:) - mol%x(:,bf1%iatom)
        vxb(:) = xx(:) - mol%x(:,bf2%iatom)
        forall(i=0:bf1%am) vxa_power(:,i) = vxa(:)**i
        forall(j=0:bf2%am) vxb_power(:,j) = vxb(:)**j
        forall(i=1:nc) vxp2(i) = sum( (xx(:) - bfp%p(i,:))**2 )
        !
        zeta(:)   =  shp%zeta(:)
        zeta(:)   =  exp( -zeta(:) * vxp2(:))
        coeff12(:,:) =  bfp%c12(:,:)
        forall(i=1:bfp%isp12)  fact(i) = sum( zeta(:) * coeff12(:,i) )
        !
        isp1 = bf1%isp
        isp2 = bf2%isp
        ibf1 = 1
        do i=1,isp1
            amm1 = bf1%am
            if(isp1 == 2 .and. i == 1) amm1 = 0
                do imx=0,amm1
                    bf1_nx = amm1 - imx
                    do imy=0,imx
                        bf1_ny = imx - imy
                        bf1_nz = imy

                        ibf2 = 1
                        do j=1,isp2
                            ij = i+(j-1)*isp1
                            amm2 = bf2%am
                            if(isp2 == 2 .and. j == 1) amm2 = 0
                            do jmx=0,amm2
                                bf2_nx = amm2 - jmx
                                do jmy=0,jmx
                                    bf2_ny = jmx - jmy
                                    bf2_nz = jmy

                                    bf_r(ibf1,ibf2) = fact(ij) * gaussian_r( bf1_nx,bf1_ny,bf1_nz,&
                                                                           & bf2_nx,bf2_ny,bf2_nz )
            
                                end do
                            end do
                        end do
                        ibf1 = ibf1 + 1
		    end do
	        end do
        end do
        contains
        !============================
        pure function gaussian_r(n1_x,n1_y,n1_z,n2_x,n2_y,n2_z)
            implicit none
            integer,intent(in)  ::  n1_x,n1_y,n1_z,n2_x,n2_y,n2_z
            !===
            real(dp)    ::  gaussian_r
            gaussian_r = 1.0_dp
            !Bra
            if(n1_x /= 0) gaussian_r = gaussian_r * vxa_power(1,n1_x)
            if(n1_y /= 0) gaussian_r = gaussian_r * vxa_power(2,n1_y)
            if(n1_z /= 0) gaussian_r = gaussian_r * vxa_power(3,n1_z)
            !Ket
            if(n2_x /= 0) gaussian_r = gaussian_r * vxb_power(1,n2_x)
            if(n2_y /= 0) gaussian_r = gaussian_r * vxb_power(2,n2_y)
            if(n2_z /= 0) gaussian_r = gaussian_r * vxb_power(3,n2_z)
            !
        end function
    end subroutine calc_bfpair
end module m_dft_xc
