 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*

module m_shells
	use m_definitions
	use m_tools
	use m_atoms
	implicit none

	type shell
	    integer 		 	:: am
	    logical			:: bool_sp
	    integer  			:: itype
	    integer		        :: ncontr
	    real(dp),allocatable	:: alpha(:)
	    real(dp),allocatable    :: coeff(:,:)
	end type shell

	type shell_pair
	    real(dp),allocatable		:: o2a(:)				!1/(2A)
	    real(dp),allocatable		:: o2b(:)				!1/(2B)
	    real(dp),allocatable		:: zeta(:)				!(A+B)
	    real(dp),allocatable	        :: oozeta(:)			        !1/(A+B)
	    real(dp),allocatable	        :: oo2zeta(:)			        !1/2(A+B)
	    real(dp),allocatable	        :: rho(:)				!AB/(A+B)
	    real(dp),allocatable	        :: boab(:)				!B/(A+B)
	    real(dp),allocatable	        :: aoab(:)				!A/(A+B)
	end type shell_pair

	type shell_set
            character(len=100)                ::      basis_set_name
            integer                           ::      nshell
            integer,allocatable   		   ::     nbf_shell(:)
            integer,allocatable		   ::     nbf_atype(:)
            integer,allocatable		   ::     ishell_atype(:)
            integer,allocatable		   ::     nshell_atype(:)
            type(shell),allocatable          ::     shells(:)
            type(shell_pair),allocatable     ::     shell_pair_data(:,:)
	end type shell_set

	contains
!=========================================================================
	subroutine read_shell_set(basis_path,basis_name,mol,basis,print_)
		implicit none
		character(len=100),intent(in)	:: basis_path
		character(len=100),intent(in)	:: basis_name
		type(geom),intent(in)	        :: mol
		type(shell_set),intent(out)	        :: basis
               logical,optional,intent(in)       :: print_
		!===
		character(len=100)            	:: basis_filename
		real(dp),allocatable          	:: alpha(:),coeff1(:),coeff2(:)
		logical                       	:: file_exists
		integer                       	:: basisfile
		integer				:: itype,ibf_file,ig,ishell,ishell_l
		integer				:: am_tmp,nbf_file,ng
		integer				:: ii,jj
		if(TRIM(basis_name) == 'none') return
               basis%basis_set_name = basis_name
		basis%nshell = 0
		do itype = 1,mol%ntype
			basis_filename=ADJUSTL(TRIM(basis_path)//'/'&
						&//TRIM(ADJUSTL(element_name(mol%ztype(itype))))//'_'&
						&//TRIM(basis_name))
			inquire(file=TRIM(basis_filename),exist=file_exists)
			if(.not.file_exists) then
				write(stdout,'(a,a)') ' Looking for file ',TRIM(basis_filename)
		 		stop'Basis set file not found'
			end if
			open(newunit=basisfile,file=TRIM(basis_filename),status='old')
			read(basisfile,*) nbf_file
			if(nbf_file<1) stop'ERROR in basis set file'
			basis%nshell = basis%nshell + nbf_file
			close(basisfile)
		end do
		allocate(       basis%shells( basis%nshell ),    &
				& basis%nbf_shell( basis%nshell ), &
				& basis%nbf_atype( mol%ntype ),    &
				& basis%ishell_atype( mol%ntype ), &
				& basis%nshell_atype( mol%ntype ) )
		ishell = 1
		do itype = 1,mol%ntype
			basis_filename=ADJUSTL(TRIM(basis_path)//'/'&
						&//TRIM(ADJUSTL(element_name(mol%ztype(itype))))//'_'&
						&//TRIM(basis_name))
			open(newunit=basisfile,file=TRIM(basis_filename),status='old')
			read(basisfile,*) nbf_file
			ishell_l = ishell
			basis%ishell_atype(itype) = ishell_l
			do ibf_file = 1,nbf_file
				read(basisfile,*) ng,am_tmp
				if(ng<1) stop'ERROR in basis set file'
				allocate(alpha(ng),coeff1(ng),coeff2(ng))
                              !
                              !Inverted order read
				if(am_tmp < 10) then
					do ig=ng,1,-1
						read(basisfile,*) alpha(ig),coeff1(ig)
					end do
				else
					do ig=ng,1,-1
						read(basisfile,*) alpha(ig),coeff1(ig),coeff2(ig)
					end do
				end if
				basis%nbf_shell(ishell) = number_basis_function_am('CART',am_tmp)
				call init_shell(basis%shells(ishell),ng,am_tmp,itype,alpha,coeff1,coeff2)
				ishell = ishell + 1
				deallocate(alpha,coeff1,coeff2)
			end do
			basis%nshell_atype(itype) = ishell-ishell_l
			basis%nbf_atype(itype) = sum( basis%nbf_shell(ishell_l : ishell-1) )
			close(basisfile)
		end do
		allocate(basis%shell_pair_data(basis%nshell,basis%nshell))
		do ii=1,basis%nshell
			do jj=1,basis%nshell
				call do_shellpair( basis%shells(ii) , basis%shells(jj) , basis%shell_pair_data(ii,jj) )
			end do
		end do
               if(present(print_))  then
                    write(stdout,'(/)')
                    write(stdout,'(/,a)'), '=========== Shell List ============'
                    call print_shell_set(basis)
               end if
	end subroutine 
!=========================================================================
	subroutine init_shell(sh,ncontr,am,itype,alpha,coeff1,coeff2)
		implicit none
		type(shell),intent(out)	:: sh
		integer,intent(in)		:: ncontr,am,itype
		real(dp),intent(in)		:: alpha(:)
		real(dp),intent(in)		:: coeff1(:),coeff2(:)
		if(am < 0) stop 'Data error in init_shell'
		allocate(sh%alpha(ncontr))
		if(am == 10) then
			allocate(sh%coeff(ncontr,2))
			sh%bool_sp = .true.
		else 
			allocate(sh%coeff(ncontr,1))
			sh%bool_sp = .false.
		end if
		sh%ncontr = ncontr
		sh%am = am
		sh%itype = itype
		sh%alpha(:) = alpha(:)
		if(am /= 10) then
			call renorm(coeff1,sh%coeff(:,1),am)
		else
			call renorm(coeff1,sh%coeff(:,1),0)
			call renorm(coeff2,sh%coeff(:,2),1)
			sh%am = 1
		end if
		contains
		!===
               !Normalization
		subroutine renorm(coeff,norm_coeff,am)
			use m_tools,only : double_factorial
			implicit none
			real(dp),intent(in)		:: coeff(:)
			integer,intent(in)		:: am
			real(dp),intent(out)	        :: norm_coeff(:)
			real(dp)			:: tmp(ncontr)
			real(dp),parameter      :: sqrt_Pi_cubed = 5.568327996831707_dp
			integer		        :: i
			if( any(alpha(:) < 0.0) ) stop 'ERROR for alpha coefficent '
			tmp(:) = 2.0_dp * alpha(:)
			tmp(:) = tmp(:)**(am+1) * sqrt(tmp(:))
			tmp(:) = sqrt(2.0_dp ** am * tmp(:) / (sqrt_Pi_cubed * double_factorial(2*am-1)))
			norm_coeff(:) = tmp(:) * coeff(:)
		end subroutine
	end subroutine 
!=========================================================================
	subroutine do_shellpair(sh1,sh2,shellp)
		implicit none
		type(shell),intent(in)		:: sh1,sh2
		type(shell_pair),intent(out)	        :: shellp
		integer				:: i,j,nc1,nc2,ntol
		real(dp)				:: tmp( sh1%ncontr,sh2%ncontr )
		nc1 = sh1%ncontr
		nc2 = sh2%ncontr
		ntol = nc1*nc2
		!FIXME
		allocate(       shellp%o2a(ntol),		&
		        &       shellp%o2b(ntol),		&
		        &       shellp%zeta(ntol),		&
		        &       shellp%oozeta(ntol),	        &
			&	shellp%oo2zeta(ntol),	        &
			&	shellp%rho(ntol),		&
			&	shellp%aoab(ntol),		&
			&	shellp%boab(ntol)		)
		!FIXME
		forall(i=1:nc1,j=1:nc2) 
			shellp%o2a(i+(j-1)*nc1) = sh1%alpha(i)					 !A
			shellp%o2b(i+(j-1)*nc1) = sh2%alpha(j)					 !B
			shellp%zeta(i+(j-1)*nc1) = sh1%alpha(i) + sh2%alpha(j) 			 !(A+B)
			shellp%rho(i+(j-1)*nc1) = sh1%alpha(i) * sh2%alpha(j)		        !(A*B)
		end forall
		shellp%oozeta = 1.0_dp / shellp%zeta							!1/(A+B)
		shellp%rho = shellp%rho * shellp%oozeta					        !A*B/(A+B)
		shellp%oo2zeta = 0.5_dp * shellp%oozeta						!1/2(A+B)
		shellp%aoab = shellp%o2a *  shellp%oozeta						!A/(A+B)
		shellp%boab = shellp%o2b *  shellp%oozeta					    	!B/(A+B)
		shellp%o2a = 0.5_dp / shellp%o2a							!1/2A
		shellp%o2b = 0.5_dp / shellp%o2b							!1/2B
	end subroutine
!=========================================================================
	subroutine destroy_shell_set(shellset)
		implicit none
		type(shell_set),intent(inout) :: shellset
		integer		        :: i,j
		if(allocated(shellset%nbf_shell))  deallocate(shellset%nbf_shell)
		if(allocated(shellset%nbf_atype))  deallocate(shellset%nbf_atype)
		if(allocated(shellset%ishell_atype))  deallocate(shellset%ishell_atype)
		if(allocated(shellset%nshell_atype))  deallocate(shellset%nshell_atype)
		if(allocated(shellset%shells)) then
			do i=1,size(shellset%shells)
				call destroy_shell(shellset%shells(i))
			end do
			deallocate(shellset%shells)
		end if
		if(allocated(shellset%shell_pair_data)) then
			do i=1,size(shellset%shell_pair_data,1)
				do  j=1,size(shellset%shell_pair_data,2)
					call destroy_shell_pair(shellset%shell_pair_data(i,j))
				end do
			end do
			deallocate(shellset%shell_pair_data)
		end if
	end subroutine
	!=========
	subroutine destroy_shell(sh)
		implicit none
		type(shell),intent(inout) :: sh
		if(allocated(sh%alpha))	 	deallocate(sh%alpha)
		if(allocated(sh%coeff)) 	deallocate(sh%coeff)
		end subroutine
		!=========
		subroutine destroy_shell_pair(shellp)
		type(shell_pair),intent(inout) :: shellp
		if(allocated(shellp%o2a)) then
			deallocate(	shellp%o2a,  	&
					  &	shellp%o2b,   	&
					  &	shellp%zeta,   	&
					  &	shellp%oozeta,  &
					  &	shellp%oo2zeta, &
					  &	shellp%rho,   	&
					  &	shellp%boab,   	&
					  &	shellp%aoab   	)
		end if
	end subroutine
!=========================================================================
	function number_basis_function_am(gaussian_type,am)
		implicit none
	 	character(len=4),intent(in) :: gaussian_type
	 	integer,intent(in)          :: am
	 	integer                     :: number_basis_function_am
		!=====

	 	select case(gaussian_type)
	 	case('CART')
		   	select case(am)
		   	case(0)
				number_basis_function_am = 1
		   	case(1)
				number_basis_function_am = 3
			case(2)
				number_basis_function_am = 6
			case(3)
			 	number_basis_function_am = 10
			case(4)
			 	number_basis_function_am = 15
			case(5)
			 	number_basis_function_am = 21
			case(6)
			 	number_basis_function_am = 28
			case(7)
			 	number_basis_function_am = 36
			case(10) ! stands for SP orbitals
				number_basis_function_am = 4 
			case default
				write(stdout,*) 'am=',am
			 	stop'number_basis_function_am: not implemented'
			end select
	 	case('PURE')
			if(am/=10) then
				number_basis_function_am = 2 * am + 1
			else ! stands for SP orbitals
				number_basis_function_am = 4 
		   	endif
		end select
	end function
!=========================================================================
	subroutine print_shell_set(basis)
		implicit none
		type(shell_set),intent(in) :: basis
		integer					   :: i,ig,nam
		write(stdout,'(/,a,x,i3)') "Number of shells :",basis%nshell
               write(stdout,'(/,a,x,a)') "Basis set :",basis%basis_set_name
		write(stdout,'(/,a,x)') "The information of shells below :"
		do i=1,basis%nshell
			write(stdout,'(a,i3,2x,a,i3,2x,a,i3,2x,a,l)') "ishell=",i,"mom=",basis%shells(i)%am,&
									 & "ishell_type=",basis%shells(i)%itype,&
									 & "bool_sp=",basis%shells(i)%bool_sp
			write(stdout,'(2x,a,5x,a)') "alpha= ","coeff= "
			do ig=1,basis%shells(i)%ncontr
				if( size(basis%shells(i)%coeff(:,:),2) == 1 ) then
					write(stdout,'(x,2f12.6)'),basis%shells(i)%alpha(ig),&
									& basis%shells(i)%coeff(ig,1)  
				else
					write(stdout,'(x,3f12.6)'),basis%shells(i)%alpha(ig),  &
									& basis%shells(i)%coeff(ig,1),&
								        & basis%shells(i)%coeff(ig,2)  
				end if
			end do
			write(stdout,'(/)')
		end do
	end subroutine
!=========================================================================
end module m_shells
