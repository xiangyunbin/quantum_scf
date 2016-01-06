 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*

program main
        use m_definitions
        use m_atoms
        use m_shells
        use m_basis
        use m_eri
        use m_hamilton
        use m_dft_grid
        use m_dft_xc
        implicit none

        character(len=100) 			::  basis_path,basis_name
        character(len=100) 			::  geom_path,filename
        type(geom) 				::  mol
        type(shell_set) 			::  shells
        type(basis_set) 			::  basis
        type(basis_quartet_set)		::  basisqs
        type(grid_type)                     ::  gridt
        type(dft_grid)                      ::  grids
        !==
        integer,parameter                 ::  iter_max = 25
        real(dp),parameter                 ::  threshold = 1.0d-8
        integer				::  nbas,ne,ndocc
        integer                             :: iter
        real(dp)				::  time0,time1,time2,time3,time4
        real(dp),allocatable               ::  S(:,:),K(:,:),V(:,:)
        real(dp),allocatable               ::  H(:,:),F(:,:)
        real(dp),allocatable               ::  D(:,:),S_minus_sqrt(:,:)
        real(dp)                             ::  ehf,enuc,ehf_old,err
        real(dp)                             ::  su,rho
        logical                             ::  print_

        call cpu_time(time0)

        filename = "out.txt"
        open(stdout,file = filename,status = "old")

        geom_path = "ch4.txt"
        call readatoms(geom_path,mol,print_)

        basis_path = "/home/xiang/gfort/scf/basis"
        basis_name = "STO-3G"
        call read_shell_set(basis_path,basis_name,mol,shells,print_)

        call init_basis_set(mol,shells,basis,print_)

        nbas = basis%nbf
        allocate( S( nbas, nbas),&
		  & K( nbas, nbas),&
	  	  & V( nbas, nbas),&
                 & H(nbas,nbas),&
                 & F(nbas,nbas),&
                 & D(nbas,nbas), &
                 & S_minus_sqrt(nbas,nbas)&
                 )

        call do_overlap_kinetic_matrix(mol,shells,basis,S,K,print_)
        call do_nucleus_matrix(mol,shells,basis,V,print_)  
        call cpu_time(time1)

        call precalc_eri(shells,basis,basisqs,print_)
        call cpu_time(time2)

        call nucleus_nucleus_energy(mol,enuc)
        call seutp_S_minus_sqrt(S,S_minus_sqrt,print_)
        ne = sum(mol%zatom(:))
        ndocc = ne/2
        H = K + V
        call generalized_eigen_solver(H,S_minus_sqrt,D,ndocc)
        ehf = sum( D * (H + H) )
        write(stdout,'(10x,a,17x,a)'),"Energy", "Error"

        do iter=1,iter_max
            ehf_old = ehf
            call setup_fock_matirx(D,basis,basisqs,F)
            F = F + H
            call generalized_eigen_solver(F,S_minus_sqrt,D,ndocc)
            ehf = sum( D * (H + F) ) 
            err = (ehf - ehf_old)/ehf
            write(stdout,'(2f12.6)'),ehf + enuc,err
            if( abs(err) < threshold ) exit
            if( iter == iter_max-1 ) stop 'SCF : Not Convergence! '
        end do
        call cpu_time(time3)

        write(stdout,'(/)')
        write(stdout,'(a,2x,f12.6)'),"1e Hamiltonian",(time1-time0)
        write(stdout,'(a,2x,f12.6)'),"Pre ERI :",(time2-time1)
        write(stdout,'(a,2x,f12.6)'),"SCF :",(time3-time2)
        write(stdout,'(a,2x,f12.6)'),"Total Time :",(time3-time0)


        call destroy_geom(mol)
        call destroy_shell_set(shells)
        call destroy_basis_set(basis)
        call destroy_basis_quartet_set(basisqs)
        call destroy_dft_grid(grids)
        deallocate(S,K,V,H,F,D,S_minus_sqrt)
end program
