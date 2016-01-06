 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
program main
        use m_definitions
        use m_atoms
        use m_dft_grid
        use m_dft_xc
        implicit none

        character(len=100) 			::  geom_path,filename
        type(geom) 				::  mol
        type(grid_type)                     ::  gridt
        type(dft_grid)                      ::  grids
        !==
        integer                             ::  iatom,igrid
        real(dp)				::  start,endt
        real(dp)                             ::  su,rho
        logical                             ::  print_

        filename = "out_grid.txt"
        open(stdout,file = filename,status = "old")

        geom_path = "O2.txt"
        call readatoms(geom_path,mol,print_)

        call setup_partition(mol,grids,print_)
        su = 0.0_dp
        do iatom=1,mol%natom
            do igrid=1,grids%gts(iatom)%ngrid_type
                su = su + grids%gts(iatom)%wt(igrid) * gaussian_r(grids%gts(iatom)%xt(:,igrid))
            end do
        end do
        write(stdout,'(a,2x,f20.10)'),"su:",su
        write(stdout,'(a,2x,f20.10)'),"exact:",2*pi**1.5 

        call destroy_geom(mol)
        call destroy_dft_grid(grids)

        contains
        pure function gaussian_r(xx)
            implicit none
            real(dp),intent(in)    ::  xx(3)
            !==
            real(dp) :: gaussian_r
            gaussian_r = exp(-sum((xx-mol%x(:,1))**2)) +exp(-sum((xx-mol%x(:,2))**2))
        end function
        !
end program
