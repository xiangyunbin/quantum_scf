 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_atoms
	use m_definitions

	integer,parameter,private		:: nelement_max = 54
	character(len=2),parameter,private :: element_list(nelement_max) =                       &
  	(/' H',                                                                                'He', &  !  2
         'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 10
         'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 18
         ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
         'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe'  /) ! 54

        type geom
        integer :: natom
        integer :: ntype							
        real(dp),allocatable  :: zatom(:)			
        real(dp),allocatable  :: ztype(:)			
        real(dp),allocatable  :: x(:,:)				
        real(dp),allocatable  :: dist2(:,:) 
        real(dp),allocatable  :: dist(:,:) 			
        integer,allocatable  :: iatom_to_type(:)	       
        integer,allocatable  :: itype_number(:)
        end type geom

contains
!=========================================================================

	function element_number(element_name)
		implicit none
		character(len=2),intent(in)	:: element_name
		integer						:: element_number
		element_number = 1
		do while( ADJUSTL(element_name) /= ADJUSTL(element_list(element_number)) )
			if( element_number == nelement_max ) then
				write(stdout,'(a,a)')    ' Input symbol ',element_name
		 		write(stdout,'(a,i3,a)') ' Element symbol is not one of first ',nelement_max,' elements'
		 		stop'element symbol not understood'
			end if
			element_number = element_number + 1
		end do
	end function
!=========================================================================

	function element_name(zatom)
		implicit none
		real(dp),intent(in)	:: zatom
		character(len=2) 	:: element_name

		if( NINT(zatom) > nelement_max ) then
			write(stdout,'(a,i3,a)') 'Element symbol is not one of first ',nelement_max,' elements'
			stop 'element symbol not understood'
		end if
		element_name = element_list(NINT(zatom))
	end function
!=========================================================================

	subroutine readatoms(filename,mol,print_)
		implicit none
		character(len = 50),intent(in) 	:: filename
		type(geom),intent(out) 		:: mol
               logical,optional,intent(in)       :: print_
               !
		integer 				:: i,ios,N
		character(len = 2) 			:: symbol
		open(7,file = filename,status = "old",iostat = ios)
		if(ios /= 0) then
			write(stdout,'(a,a)') ' Looking for file ',TRIM(filename)
		 	stop'Molecular geometry file not found'
		end if
		read(7,*) N
		mol%natom = N
		allocate(mol%zatom(N),mol%x(3,N))
		do i = 1,N
			read(7,*) symbol,mol%x(1:3,i)
			mol%x(:,i) = mol%x(:,i) / bohr_A
			mol%zatom(i) = element_number(symbol)
		end do
		call do_distmatrix(mol)
		call do_screen(mol)
		close(7)
               if(present(print_)) then
                    write(stdout,'(/,a)'), '=========== Atom List ============'
                    call print_atoms(mol)
               end if
		contains
			subroutine do_distmatrix(mol)
				implicit none
				type(geom),intent(inout) 	:: mol
				integer 			  		:: i,j
				allocate( mol%dist2(mol%natom,mol%natom) )
                              allocate( mol%dist(mol%natom,mol%natom) )
                              mol%dist2(:,:) = 0.0_dp
                              mol%dist(:,:) = 0.0_dp
				do i=1,mol%natom
					mol%dist2(i,i) = 0.0_dp
					do j=i+1,mol%natom
                                            mol%dist2(i,j) = SUM( (mol%x(:,i) - mol%x(:,j))**2 )
                                            mol%dist2(j,i) = mol%dist2(i,j)
                                            mol%dist(i,j) = sqrt(mol%dist2(i,j))
                                            mol%dist(j,i) = mol%dist(i,j)
					end do
				end do
			end subroutine
			!==============
			subroutine do_screen(mol)
				implicit none
				type(geom),intent(inout) 	:: mol
				!==
				integer    					:: tmp(nelement_max)
				integer 					:: i,iatom
				integer					:: ind_atom(nelement_max)
				logical					:: mark(nelement_max)
				tmp(:) = 0
				ind_atom(:) = (/(i,i=1,nelement_max)/)
				do i=1,mol%natom
					iatom = NINT(mol%zatom(i))
					tmp(iatom) = tmp(iatom) + 1  
				end do
				mark = (tmp > 0)
				mol%ntype = count(mark)
				allocate(mol%ztype(mol%ntype))
				allocate(mol%itype_number(mol%ntype))
				allocate(mol%iatom_to_type(mol%natom))

				mol%ztype(:) = pack(ind_atom,mark)
				mol%itype_number(:) = pack(tmp(:),mark)
				do i=1,mol%ntype
					where( NINT(mol%zatom(:)) == NINT(mol%ztype(i)) ) mol%iatom_to_type = i
				end do
			end subroutine 
	end subroutine readatoms

!=========================================================================
	subroutine destroy_geom(mol)
		implicit none
		type(geom),intent(inout) :: mol
		!===
		if(allocated(mol%zatom)) 			deallocate(mol%zatom)
		if(allocated(mol%ztype))			deallocate(mol%ztype)
		if(allocated(mol%x)) 				deallocate(mol%x)
		if(allocated(mol%dist2)) 			deallocate(mol%dist2)
               if(allocated(mol%dist)) 			deallocate(mol%dist)
		if(allocated(mol%iatom_to_type))	        deallocate(mol%iatom_to_type)
		if(allocated(mol%itype_number))		deallocate(mol%itype_number)
	end subroutine destroy_geom
!=========================================================================

	subroutine nucleus_nucleus_energy(mol,energy)
		implicit none
		type(geom),intent(in) 	:: mol
		real(dp),intent(out) 	:: energy
		integer              	:: iatom,jatom
		energy = 0.0_dp
		do iatom=1,mol%natom
			do jatom=iatom+1,mol%natom
				energy = energy + mol%zatom(iatom) * mol%zatom(jatom) / SQRT( mol%dist2(iatom,jatom) )
			end do
		end do
	end subroutine nucleus_nucleus_energy

!=========================================================================

	subroutine print_atoms(mol) 
		implicit none
		type(geom),intent(in) :: mol
		integer 			  :: i
		write(stdout,'(a,x,i3)') "NAtoms :",mol%natom
		write(stdout,'(a,x,i3)') "NAtomic Type :",mol%ntype
		write(stdout,'(/,a,2x,a)') "Element","Number"
		do i=1,mol%ntype
			write(stdout,'(x,a3,6x,i3)') element_name(mol%ztype(i)),mol%itype_number(i)
		end do
		write(stdout,'(/,a)') "Input orientation:" 
		write(stdout,'(a,2x,a,4x,a)') "Element","TypeIndex","Coordination"
               write(stdout,'(a)') "-------------------------------------------"
		do i=1,mol%natom
			write(stdout,'(a3,8x,i3,4x,3f12.6)') element_name(mol%zatom(i)),mol%iatom_to_type(i),mol%x(:,i)
		end do
	end subroutine 
!=========================================================================
        function element_covalent_radius(zatom)
            implicit none
            integer,intent(in)  ::  zatom
            real(dp)              ::  element_covalent_radius
            !===
            ! Data from Cambridge Structural Database
            !
            ! Values are first given in picometer
            ! They will be converted in bohr just after
            select case(zatom)
            case(1)
                element_covalent_radius =  31.
            case( 2)
                element_covalent_radius =  28.
            case( 3)
                element_covalent_radius = 128.
            case( 4)
                element_covalent_radius =  96.
            case( 5)
                element_covalent_radius =  84.
            case( 6)
                element_covalent_radius =  73.
            case( 7)
                element_covalent_radius =  71.
            case( 8)
                element_covalent_radius =  66.
            case( 9)
                element_covalent_radius =  57.
            case(10) ! Ne.
                element_covalent_radius =  58.
            case(11)
                element_covalent_radius = 166.
            case(12)
                element_covalent_radius = 141.
            case(13)
                element_covalent_radius = 121.
            case(14)
                element_covalent_radius = 111.
            case(15)
                element_covalent_radius = 107.
            case(16)
                element_covalent_radius = 105.
            case(17)
                element_covalent_radius = 102.
            case(18) ! Ar.
                element_covalent_radius = 106.
            case(19)
                element_covalent_radius = 203.
            case(20)
                element_covalent_radius = 176.
            case(21)
                element_covalent_radius = 170.
            case(22)
                element_covalent_radius = 160.
            case(23)
                element_covalent_radius = 153.
            case(24)
                element_covalent_radius = 139.
            case(25)
                element_covalent_radius = 145.
            case(26)
                element_covalent_radius = 145.
            case(27)
                element_covalent_radius = 140.
            case(28)
                element_covalent_radius = 124.
            case(29)
                element_covalent_radius = 132.
            case(30)
                element_covalent_radius = 122.
            case(31)
                element_covalent_radius = 120.
            case(32)
                element_covalent_radius = 119.
            case(34)
                element_covalent_radius = 120.
            case(35)
                element_covalent_radius = 120.
            case(36) ! Kr.
                element_covalent_radius = 116.
            case default
                stop 'radius not available'
            end select
            ! pm to bohr conversion
            element_covalent_radius = element_covalent_radius / bohr_A * 0.01_dp
        end function element_covalent_radius
!=========================================================================
end module m_atoms
