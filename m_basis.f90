 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_basis
	use m_definitions
	use m_tools
	use m_atoms
	use m_shells

	type basis_function
		integer 	        :: am												
               logical	        :: bool_sp 										
               integer	        :: isp												
		integer	        :: ibas												
		integer	        :: iatom											
		integer	        :: ishell											
		integer	        :: nbf												
		integer 	        :: ibf_left,ibf_right								
	end type basis_function
	!==
	type basis_pair
		integer		:: nc12									
		integer		:: ibas,jbas							
        	integer		:: isp1,isp2,isp12					
            	integer		:: am1,am2,amp							
               real(dp)		:: AB2,AB(3)									
		real(dp),allocatable	:: p(:,:),pa(:,:),pb(:,:)		!P,and P-A,P-B
		real(dp),allocatable	:: c12(:,:)           			!C_1*C_2*exp(-k12*|A-B|^2)
	end type basis_pair
	!==
	type basis_set
		integer					:: nbf					
		integer					:: nbas					
		type(basis_function),allocatable	        :: bf(:)							
               type(basis_pair),allocatable	        :: basis_pair_data(:,:)		
	end type basis_set

	contains
!=========================================================================
	subroutine init_basis_set(mol,shellset,basis,print_)
		implicit none
		type(geom),intent(in)	    :: mol
		type(shell_set),intent(in)	    :: shellset
		type(basis_set),intent(out)     :: basis
               logical,optional,intent(in)    :: print_
		!==
		integer			:: iatom,itype,ishell_l,ibf,ibas,itmp
		integer			:: nbf_type,nshell_atype,nbf,nbas
		integer			:: i,j
		nbf = 0
		nbas = 0
		do iatom=1,mol%natom
			itype = mol%iatom_to_type(iatom)
			nbas = nbas + shellset%nshell_atype(itype)
			nbf = nbf + shellset%nbf_atype(itype)
		end do
		basis%nbf = nbf
		basis%nbas = nbas
		allocate( basis%bf( nbas ) )
		allocate( basis%basis_pair_data( nbas,nbas ) )
		ibas = 0
		ibf = 1
		do iatom=1,mol%natom
			itype = mol%iatom_to_type(iatom)
			ishell_l = shellset%ishell_atype(itype)
			nshell_atype = shellset%nshell_atype(itype)
			do i=0,nshell_atype-1
				ibas = ibas + 1
				itmp =  ishell_l + i 	
				basis%bf( ibas )%ibas 		= 	ibas
				basis%bf( ibas )%am 		= 	shellset%shells( itmp )%am
				basis%bf( ibas )%bool_sp 	= 	shellset%shells( itmp )%bool_sp
				basis%bf( ibas )%iatom		= 	iatom
				basis%bf( ibas )%ishell 	= 	itmp
				basis%bf( ibas )%ibf_left	= 	ibf
				basis%bf( ibas )%nbf            =       shellset%nbf_shell( itmp )
				ibf = ibf + basis%bf(ibas)%nbf
				basis%bf(ibas)%ibf_right 	= 	ibf - 1
				if( basis%bf( ibas )%bool_sp .eqv. .true. ) then
					basis%bf( ibas )%isp  =  2
				else
					basis%bf( ibas )%isp  =  1
				end if
			end do
		end do
		!FIXME : the symmetrical problem
		do i=1,nbas
			do j=1,i
			    call do_basispair( basis%bf(i) , basis%bf(j) , basis%basis_pair_data(i,j) )
			end do
		end do
               if(present(print_)) then
                    write(stdout,'(/)')
                    write(stdout,'(/,a)'), '=========== Basis List ============'
                    call print_basis(basis)
               end if
		contains
		subroutine do_basispair(bf1,bf2,basisp)
			implicit none
			type(basis_function),intent(in) 	  ::  bf1,bf2
			type(basis_pair),intent(out)          ::  basisp
			integer				  ::  i,j,k,l,ij
			integer				  ::  nc1,nc2,nc12,ishell1,ishell2
			integer				  ::  isp1,isp2,isp12
			real(dp),allocatable			  ::  coeff1(:,:),coeff2(:,:),rho(:)
			real(dp)				  ::  AB2,A(3),B(3)
			!===
			ishell1 = bf1%ishell
			ishell2 = bf2%ishell
			nc1 = shellset%shells( ishell1 )%ncontr
			nc2 = shellset%shells( ishell2 )%ncontr
			nc12 = nc1 * nc2
			A(:) = mol%x(:,bf1%iatom)
			B(:) = mol%x(:,bf2%iatom)
			isp1 = bf1%isp
			isp2 = bf2%isp
			isp12 = isp1 * isp2
			basisp%AB2 =  mol%dist2(bf1%iatom,bf2%iatom)
                      basisp%AB(:) =  A(:)-B(:)
			basisp%isp1 = isp1
			basisp%isp2 = isp2
			basisp%isp12 = isp12
			basisp%am1 = bf1%am
			basisp%am2 = bf2%am
			basisp%amp = bf1%am + bf2%am
			allocate(  coeff1(nc1,isp1), &
				&  coeff2(nc2,isp2), &
			        &  rho(nc12)         )
			coeff1(:,:) = shellset%shells( ishell1 )%coeff(:,:)
			coeff2(:,:) = shellset%shells( ishell2 )%coeff(:,:)
			rho(:) = shellset%shell_pair_data( ishell1,ishell2 )%rho(:)
			allocate( basisp%p( nc12,3 ), 	   &
			        & basisp%pa( nc12,3 ),     &
				& basisp%pb( nc12,3 ),     &
				& basisp%c12( nc12,isp12 ) )
			basisp%nc12 = nc12
			basisp%ibas = bf1%ibas
			basisp%jbas = bf2%ibas
			do i=1,3
				basisp%p(:,i) = shellset%shell_pair_data( ishell1,ishell2 )%aoab(:) * A(i) &
							 &+ shellset%shell_pair_data( ishell1,ishell2 )%boab(:) * B(i)
				basisp%pa(:,i) = basisp%p(:,i) - A(i)
				basisp%pb(:,i) = basisp%p(:,i) - B(i)
			end do
			do i=1,isp1
				do j=1,isp2
					ij = i+(j-1)*isp1
					forall(k=1:nc1,l=1:nc2) 
						basisp%c12(k+(l-1)*nc1,ij) = coeff1(k,i) * coeff2(l,j)
					end forall
					basisp%c12(:,ij) = basisp%c12(:,ij) * exp( -rho(:)*basisp%AB2)
				end do
			end do
			deallocate(coeff1,coeff2,rho)
		end subroutine
	end subroutine 

!=========================================================================
!
!Obara and Saika, JCP  87 3963 (1986)
	subroutine do_overlap_kinetic_matrix(mol,shellset,basis,overlap_mat,kinitic_mat,print_)
		implicit none
		type(geom),intent(in)	:: 	mol
		type(shell_set),intent(in)	:: 	shellset
		type(basis_set),intent(in) 	:: 	basis
		real(dp),intent(out)		::      overlap_mat(:,:),kinitic_mat(:,:)
               logical,optional,intent(in) ::     print_
		!==
		integer			::      i,j,nc
		integer			::	i1_l,i1_r,i2_l,i2_r,nbf1,nbf2
		integer			:: 	ibasis1,ibasis2,ishell1,ishell2
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
				

				call overlap_kinetic_recurrence( basis%bf(i),&
								 basis%bf(j),&
								 nc,&
								 shellset%shell_pair_data( ishell1,ishell2 ),&
								 basis%basis_pair_data( ibasis1,ibasis2 ),&
								 overlap_mat( i1_l : i1_r, i2_l : i2_r ),&
								 kinitic_mat( i1_l : i1_r, i2_l : i2_r ) &
								)
				if( i /= j ) then
					overlap_mat( i2_l : i2_r, i1_l : i1_r ) = transpose( overlap_mat( i1_l : i1_r, i2_l : i2_r ) )
					kinitic_mat( i2_l : i2_r, i1_l : i1_r ) = transpose( kinitic_mat( i1_l : i1_r, i2_l : i2_r ) )
				end if
			end do
		end do
               if(present(print_)) then
                    write(stdout,'(/,a)'), '=========== Overlap Matrix ============'
                    write(stdout,'(/)')
                    do i=1,basis%nbf
                        !FIXME
                        write(stdout,'(i3,4x,100f12.6)'),i,overlap_mat(i,:)
                    end do
                    write(stdout,'(/)')
                    write(stdout,'(/,a)'), '=========== Kinitic Matrix ============'
                    write(stdout,'(/)')
                    do i=1,basis%nbf
                        write(stdout,'(i3,4x,100f12.6)'),i,kinitic_mat(i,:)
                    end do
               end if
	end subroutine do_overlap_kinetic_matrix
!=========================================================================
!
        subroutine overlap_kinetic_recurrence(bf1,bf2,nc,shp,bfp,s_ab,k_ab)
		implicit none
		type(basis_function),intent(in)		:: bf1,bf2
		type(basis_pair),intent(in)			:: bfp
		type(shell_pair),intent(in)			:: shp
		integer,intent(in)				:: nc
		real(dp),intent(out)				:: s_ab(1:bf1%nbf,1:bf2%nbf),k_ab(1:bf1%nbf,1:bf2%nbf)			
		!===
		real(dp)             				:: oozeta_ab(nc),boab(nc),aoab(nc)
		real(dp)                                    :: ksi_ab(nc),fact(nc),coeff12(nc,bfp%isp12)
 		real(dp)             				:: p(nc,3),pa(nc,3),pb(nc,3),ab2
 		real(dp)             				:: s_tmp_x(nc,0:bf1%am,0:bf2%am)
 		real(dp)            				:: s_tmp_y(nc,0:bf1%am,0:bf2%am)
		real(dp)             				:: s_tmp_z(nc,0:bf1%am,0:bf2%am)
		real(dp)             				:: k_tmp_x(nc,0:bf1%am,0:bf2%am)
		real(dp)             				:: k_tmp_y(nc,0:bf1%am,0:bf2%am)
		real(dp)             				:: k_tmp_z(nc,0:bf1%am,0:bf2%am)
		integer              				:: i,j,ij
		integer				        :: ix,iy,iz
		integer              				:: ixa,iya,iza
		integer              				:: ixb,iyb,izb
		integer              				:: ixap,iyap,izap
		integer              				:: ixbp,iybp,izbp
		integer					:: am1,am2,amm1,amm2
		integer				 	:: bf1_nx,bf1_ny,bf1_nz
		integer				 	:: bf2_nx,bf2_ny,bf2_nz
		integer				 	:: ibf1,ibf2,isp1,isp2
		integer					:: imx,imy,jmx,jmy
			

		am1 			= 	bf1%am
		am2 			= 	bf2%am
		oozeta_ab(:) 	        = 	shp%oozeta(:)
		ksi_ab(:) 		= 	shp%rho(:)
		aoab(:) 		= 	shp%aoab(:)
		boab(:) 		= 	shp%boab(:)
		coeff12(:,:) 	        = 	bfp%c12(:,:)
		ab2 			= 	bfp%AB2
		p(:,:)  		= 	bfp%p(:,:)
		pa(:,:) 		= 	bfp%pa(:,:)
		pb(:,:) 		= 	bfp%pb(:,:)
		fact(:) 		= 	shp%oo2zeta(:)
		!
		! direction X
		!
		s_tmp_x(:,0,0) =  ( pi*oozeta_ab(:) )**1.5_dp 
		k_tmp_x(:,0,0) = ksi_ab(:) * ( 3.0_dp - 2.0_dp * ksi_ab(:) * ab2 ) * s_tmp_x(:,0,0)

		call do_recurrence_x()
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

					    !
					    ! direction Y
					    !
					    s_tmp_y(:,0,0) = s_tmp_x(:,bf1_nx,bf2_nx)
					    k_tmp_y(:,0,0) = k_tmp_x(:,bf1_nx,bf2_nx)
					    call do_recurrence_y(bf1_ny,bf2_ny)

					    !
					    ! direction Z
	    				    !
					    s_tmp_z(:,0,0) = s_tmp_y(:,bf1_ny,bf2_ny)
					    k_tmp_z(:,0,0) = k_tmp_y(:,bf1_ny,bf2_ny)
					    call do_recurrence_z(bf1_nz,bf2_nz)
					    !
					    s_ab(ibf1,ibf2) = sum( s_tmp_z(:,bf1_nz,bf2_nz) * coeff12(:,ij) )
					    k_ab(ibf1,ibf2) = sum( k_tmp_z(:,bf1_nz,bf2_nz) * coeff12(:,ij) )

					    ibf2 = ibf2 + 1
					end do
				    end do
			        end do
			        ibf1 = ibf1 + 1
		        end do
	            end do
		end do
		contains
		!==============================================
		subroutine do_recurrence_x()
			implicit none
			do ix=1,am1+am2
			        do ixa=0,MIN(am1,ix)
					ixb=ix-ixa
					if(ixb>am2) cycle

					if(ixa>0) then
						ixap=ixa-1

						s_tmp_x(:,ixap+1,ixb) = pa(:,1) * s_tmp_x(:,ixap,ixb)
						if(ixap>0)  s_tmp_x(:,ixap+1,ixb) = s_tmp_x(:,ixap+1,ixb) + fact(:) * ixap * s_tmp_x(:,ixap-1,ixb)
						if(ixb>0)   s_tmp_x(:,ixap+1,ixb) = s_tmp_x(:,ixap+1,ixb) + fact(:) * ixb  * s_tmp_x(:,ixap,ixb-1)

						k_tmp_x(:,ixap+1,ixb) = pa(:,1) * k_tmp_x(:,ixap,ixb) + 2.0_dp * ksi_ab(:) * s_tmp_x(:,ixap+1,ixb)
						if(ixap>0)  k_tmp_x(:,ixap+1,ixb) = k_tmp_x(:,ixap+1,ixb) + fact(:) * ixap * k_tmp_x(:,ixap-1,ixb) &
									  & -boab(:) * ixap * s_tmp_x(:,ixap-1,ixb)
						if(ixb>0)   k_tmp_x(:,ixap+1,ixb) = k_tmp_x(:,ixap+1,ixb) + fact(:) * ixb  * k_tmp_x(:,ixap,ixb-1)

					else
					   	ixbp=ixb-1

					   	s_tmp_x(:,ixa,ixbp+1) = pb(:,1) * s_tmp_x(:,ixa,ixbp)
					   	if(ixbp>0) s_tmp_x(:,ixa,ixbp+1) = s_tmp_x(:,ixa,ixbp+1) + fact(:) * ixbp * s_tmp_x(:,ixa,ixbp-1)
					   	if(ixa>0)  s_tmp_x(:,ixa,ixbp+1) = s_tmp_x(:,ixa,ixbp+1) + fact(:) * ixa  * s_tmp_x(:,ixa-1,ixbp)

					   	k_tmp_x(:,ixa,ixbp+1) = pb(:,1) * k_tmp_x(:,ixa,ixbp) +  2.0_dp * ksi_ab(:) * s_tmp_x(:,ixa,ixbp+1)
					   	if(ixbp>0) k_tmp_x(:,ixa,ixbp+1) = k_tmp_x(:,ixa,ixbp+1) + fact(:) * ixbp * k_tmp_x(:,ixa,ixbp-1) &
								     & -aoab(:) * ixbp * s_tmp_x(:,ixa,ixbp-1)
					   	if(ixa>0)  k_tmp_x(:,ixa,ixbp+1) = k_tmp_x(:,ixa,ixbp+1) + fact(:) * ixa  * k_tmp_x(:,ixa-1,ixbp)
					 endif
				 enddo
			enddo
		end subroutine do_recurrence_x
		!=============================================
		subroutine do_recurrence_y(am1y,am2y)
			implicit none
			integer,intent(in)		::   am1y,am2y
			!===
			do iy=1,am1y+am2y

				do iya=0,MIN(am1y,iy)
					iyb=iy-iya
					 if(iyb>am2y) cycle

					 if(iya>0) then
						 iyap=iya-1

						 s_tmp_y(:,iyap+1,iyb) = pa(:,2) * s_tmp_y(:,iyap,iyb)
						 if(iyap>0)  s_tmp_y(:,iyap+1,iyb) = s_tmp_y(:,iyap+1,iyb) + fact(:) * iyap * s_tmp_y(:,iyap-1,iyb)
						 if(iyb>0)   s_tmp_y(:,iyap+1,iyb) = s_tmp_y(:,iyap+1,iyb) + fact(:) * iyb  * s_tmp_y(:,iyap,iyb-1)

						 k_tmp_y(:,iyap+1,iyb) = pa(:,2) * k_tmp_y(:,iyap,iyb) + 2.0_dp * ksi_ab(:) * s_tmp_y(:,iyap+1,iyb)
						 if(iyap>0)  k_tmp_y(:,iyap+1,iyb) = k_tmp_y(:,iyap+1,iyb) + fact(:) * iyap * k_tmp_y(:,iyap-1,iyb) &
									  & -boab(:) * iyap * s_tmp_y(:,iyap-1,iyb)
						 if(iyb>0)   k_tmp_y(:,iyap+1,iyb) = k_tmp_y(:,iyap+1,iyb) + fact(:) * iyb  * k_tmp_y(:,iyap,iyb-1)

					 else
						 iybp=iyb-1

					   	 s_tmp_y(:,iya,iybp+1) = pb(:,2) * s_tmp_y(:,iya,iybp)
					   	 if(iybp>0) s_tmp_y(:,iya,iybp+1) = s_tmp_y(:,iya,iybp+1) + fact(:) * iybp * s_tmp_y(:,iya,iybp-1)
					   	 if(iya>0)  s_tmp_y(:,iya,iybp+1) = s_tmp_y(:,iya,iybp+1) + fact(:) * iya  * s_tmp_y(:,iya-1,iybp)

					   	 k_tmp_y(:,iya,iybp+1) = pb(:,2) * k_tmp_y(:,iya,iybp) +  2.0_dp * ksi_ab(:) * s_tmp_y(:,iya,iybp+1)
					   	 if(iybp>0) k_tmp_y(:,iya,iybp+1) = k_tmp_y(:,iya,iybp+1) + fact(:) * iybp * k_tmp_y(:,iya,iybp-1) &
								    &  -aoab(:) * iybp * s_tmp_y(:,iya,iybp-1)
					   	 if(iya>0)  k_tmp_y(:,iya,iybp+1) = k_tmp_y(:,iya,iybp+1) + fact(:) * iya  * k_tmp_y(:,iya-1,iybp)
					endif
				enddo
			enddo
		end subroutine do_recurrence_y
		!===========================================
		subroutine do_recurrence_z(am1z,am2z)
			implicit none
			integer,intent(in)		::   am1z,am2z
			!===
			do iz=1,am1z+am2z
				do iza=0,MIN(am1z,iz)
					izb=iz-iza
					 if(izb>am2z) cycle

					 if(iza>0) then
					   	izap=iza-1

					   	s_tmp_z(:,izap+1,izb) = pa(:,3) * s_tmp_z(:,izap,izb)
					   	if(izap>0)  s_tmp_z(:,izap+1,izb) = s_tmp_z(:,izap+1,izb) + fact(:) * izap * s_tmp_z(:,izap-1,izb)
					   	if(izb>0)   s_tmp_z(:,izap+1,izb) = s_tmp_z(:,izap+1,izb) + fact(:) * izb  * s_tmp_z(:,izap,izb-1)

					   	k_tmp_z(:,izap+1,izb) = pa(:,3) * k_tmp_z(:,izap,izb) + 2.0_dp * ksi_ab(:) * s_tmp_z(:,izap+1,izb)
					   	if(izap>0)  k_tmp_z(:,izap+1,izb) = k_tmp_z(:,izap+1,izb) + fact(:) * izap * k_tmp_z(:,izap-1,izb) &
								&  -boab(:) * izap * s_tmp_z(:,izap-1,izb)
					   	if(izb>0)   k_tmp_z(:,izap+1,izb) = k_tmp_z(:,izap+1,izb) + fact(:) * izb  * k_tmp_z(:,izap,izb-1)

					 else
					   	izbp=izb-1

					   	s_tmp_z(:,iza,izbp+1) = pb(:,3) * s_tmp_z(:,iza,izbp)
					   	if(izbp>0) s_tmp_z(:,iza,izbp+1) = s_tmp_z(:,iza,izbp+1) + fact(:) * izbp * s_tmp_z(:,iza,izbp-1)
					   	if(iza>0)  s_tmp_z(:,iza,izbp+1) = s_tmp_z(:,iza,izbp+1) + fact(:) * iza  * s_tmp_z(:,iza-1,izbp)

					   	k_tmp_z(:,iza,izbp+1) = pb(:,3) * k_tmp_z(:,iza,izbp) +  2.0_dp * ksi_ab(:) * s_tmp_z(:,iza,izbp+1)
					   	if(izbp>0) k_tmp_z(:,iza,izbp+1) = k_tmp_z(:,iza,izbp+1) + fact(:) * izbp * k_tmp_z(:,iza,izbp-1) &
						        & -aoab(:) * izbp * s_tmp_z(:,iza,izbp-1)
					   	if(iza>0)  k_tmp_z(:,iza,izbp+1) = k_tmp_z(:,iza,izbp+1) + fact(:) * iza  * k_tmp_z(:,iza-1,izbp)
					endif
				enddo
			enddo
		end subroutine do_recurrence_z
	end subroutine overlap_kinetic_recurrence
!=========================================================================
!
	subroutine do_nucleus_matrix(mol,shellset,basis,nucleus_mat,print_)
		implicit none
		type(geom),intent(in)	    :: 	mol
		type(shell_set),intent(in)	    :: 	shellset
		type(basis_set),intent(in) 	    :: 	basis
		real(dp),intent(out)		    ::  nucleus_mat(:,:)
               logical,optional,intent(in)    :: print_
		!==
		integer			    ::  i,j,nc
		integer			    ::	i1_l,i1_r,i2_l,i2_r,nbf1,nbf2
		integer			    :: 	ibasis1,ibasis2,ishell1,ishell2
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
				

				call nucleus_recurrence( basis%bf(i),&
							    basis%bf(j),&
							    nc,&
						           shellset%shell_pair_data( ishell1,ishell2 ),&
							    basis%basis_pair_data( ibasis1,ibasis2 ),&
							    mol,&
							    nucleus_mat( i1_l : i1_r, i2_l : i2_r )& 
                                                        )

				if( i /= j ) then
					nucleus_mat( i2_l : i2_r, i1_l : i1_r ) = transpose( nucleus_mat( i1_l : i1_r, i2_l : i2_r ) )
				end if
			end do
		end do
                if(present(print_)) then
                    write(stdout,'(/,a)'), '=========== Nucleus Matrix ============'
                    write(stdout,'(/)')
                    do i=1,basis%nbf
                        write(stdout,'(i3,4x,50f12.6)'),i,nucleus_mat(i,:)
                    end do
                    write(stdout,'(/)')
                end if
	end subroutine do_nucleus_matrix
!=========================================================================
!see Obara and Saika, JCP  87 3963 (1986)
	subroutine nucleus_recurrence(bf1,bf2,nc,shp,bfp,mol,v_ab)
		use m_tools,only : boys_function
		implicit none
		type(basis_function),intent(in)		:: bf1,bf2
		type(basis_pair),intent(in)			:: bfp
		type(shell_pair),intent(in)			:: shp
		type(geom),intent(in)		        :: mol
		integer,intent(in)				:: nc
		real(dp),intent(out)				:: v_ab(1:bf1%nbf,1:bf2%nbf)		
		!===
		real(dp)             				:: zeta_ab(nc),ksi_ab(nc),ab2,fact(nc)
 		real(dp)             				:: p(nc,3),pa(nc,3),pb(nc,3),pc(nc,3)
 		real(dp)           				:: v_tmp_x_m(nc,0:bf1%am,0:bf2%am)
 		real(dp)             				:: v_tmp_y_m(nc,0:bf1%am,0:bf2%am)
 		real(dp)             				:: v_tmp_z_m(nc,0:bf1%am,0:bf2%am)
 		real(dp)             				:: v_tmp_x_mp1(nc,0:bf1%am,0:bf2%am)
 		real(dp)             				:: v_tmp_y_mp1(nc,0:bf1%am,0:bf2%am)
 		real(dp)             				:: v_tmp_z_mp1(nc,0:bf1%am,0:bf2%am)
		real(dp)             				:: bigu(nc),coeff12(nc,bfp%isp12)
		real(dp)					:: fmu(nc,0:bf1%am+bf2%am),vtmp
 		integer              				:: ix,iy,iz
 		integer              				:: ixa,iya,iza
 		integer              				:: ixb,iyb,izb
 		integer              				:: ixap,iyap,izap
 		integer             				:: ixbp,iybp,izbp
 		integer              				:: ixam,iyam,izam
 		integer              				:: ixbm,iybm,izbm
		integer				 	:: mm,am1,am2,amm1,amm2
		integer				 	:: bf1_nx,bf1_ny,bf1_nz
		integer				 	:: bf2_nx,bf2_ny,bf2_nz
		integer			     		:: ibf1,ibf2,isp1,isp2
		integer				 	:: imx,imy,jmx,jmy		
		integer					:: i,j,ij,iatom
		!===
		am1 			= 	bf1%am
		am2 			= 	bf2%am
		zeta_ab(:) 	        = 	shp%zeta(:)
		ksi_ab(:) 		= 	shp%rho(:)
		coeff12(:,:) 	        = 	bfp%c12(:,:)
		ab2 			= 	bfp%AB2
		p(:,:)  		= 	bfp%p(:,:)
		pa(:,:) 		= 	bfp%pa(:,:)
		pb(:,:) 		= 	bfp%pb(:,:)
		fact(:) 		= 	shp%oo2zeta(:)


		forall(i=1:bfp%isp12) coeff12(:,i) = 2.0_dp * (pi / zeta_ab(:)) * coeff12(:,i)
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


					    v_ab(ibf1,ibf2) = 0.0_dp
					    !
                                         !
					    do iatom = 1,mol%natom
					        call do_single_atom_nucleus(iatom,vtmp)
						 v_ab(ibf1,ibf2) = v_ab(ibf1,ibf2)  + vtmp
					    end do

					    ibf2 = ibf2 + 1
					end do
				    end do
				end do
			        ibf1 = ibf1 + 1
		            end do
		        end do
		end do
		contains
!=========================
		subroutine do_single_atom_nucleus(ia,vv)
			implicit none
			integer,intent(in)		::   ia
			real(dp),intent(out)	       ::   vv
			integer			::   inc
			!===
			forall(inc=1:nc) pc(inc,:)  =   p(inc,:) - mol%x(:,ia)
			bigu(:)			    =   zeta_ab(:) * ( pc(:,1)**2 + pc(:,2)**2 + pc(:,3)**2)
			call boys_function(fmu,nc,am1+am2,bigu)
		
			v_tmp_x_mp1(:,:,:) =  0.0_dp
 			v_tmp_y_mp1(:,:,:) =  0.0_dp
 			v_tmp_z_mp1(:,:,:) =  0.0_dp

			do mm = am1 + am2,0,-1
				v_tmp_x_m(:,0,0) = fmu(:,mm) 
				ixam=0
	   			ixbm=0
				call do_recurrence_x(bf1_nx,bf2_nx)

	   			! direction Y
	   			v_tmp_y_m(:,0,0) = v_tmp_x_m(:,ixam,ixbm)
				iyam=0
	   			iybm=0
				call do_recurrence_y(bf1_ny,bf2_ny)
	   			! direction Z
				v_tmp_z_m(:,0,0) = v_tmp_y_m(:,iyam,iybm)
				izam=0
				izbm=0
				call do_recurrence_z(bf1_nz,bf2_nz)

	  			v_tmp_x_mp1(:,0:ixam,0:ixbm) =  v_tmp_x_m(:,0:ixam,0:ixbm)
	   			v_tmp_y_mp1(:,0:iyam,0:iybm) =  v_tmp_y_m(:,0:iyam,0:iybm)
	   			v_tmp_z_mp1(:,0:izam,0:izbm) =  v_tmp_z_m(:,0:izam,0:izbm)
			end do 
			vv =  - real(mol%zatom(ia),dp) * sum( v_tmp_z_m(:,izam,izbm) * coeff12(:,ij) )

		end subroutine do_single_atom_nucleus
!=============================================================
		subroutine do_recurrence_x(am1x,am2x)
			implicit none
			integer,intent(in)		::  am1x,am2x
			 do ix=1,am1 + am2 - mm
		 		do ixa=0,MIN(am1x,ix)
			   		ixb=ix-ixa
			   		if(ixb>am2x) cycle
			   		ixam=MAX(ixam,ixa)
			   		ixbm=MAX(ixbm,ixb)

			   		if(ixa>0) then
				 		ixap=ixa-1
				 		v_tmp_x_m(:,ixap+1,ixb) = pa(:,1) * v_tmp_x_m(:,ixap,ixb) - pc(:,1) * v_tmp_x_mp1(:,ixap,ixb)
				 		if(ixap>0) v_tmp_x_m(:,ixap+1,ixb) = v_tmp_x_m(:,ixap+1,ixb) + &
										fact(:) * ixap * ( v_tmp_x_m(:,ixap-1,ixb) -  v_tmp_x_mp1(:,ixap-1,ixb) )
				 		if(ixb>0)  v_tmp_x_m(:,ixap+1,ixb) = v_tmp_x_m(:,ixap+1,ixb) + &
										fact(:) * ixb  * ( v_tmp_x_m(:,ixap,ixb-1) -  v_tmp_x_mp1(:,ixap,ixb-1) )
			   		else
				 		ixbp=ixb-1
				 		v_tmp_x_m(:,ixa,ixbp+1) = pb(:,1) * v_tmp_x_m(:,ixa,ixbp) - pc(:,1) * v_tmp_x_mp1(:,ixa,ixbp)
				 		if(ixbp>0) v_tmp_x_m(:,ixa,ixbp+1) = v_tmp_x_m(:,ixa,ixbp+1) + &
										fact(:) * ixbp * ( v_tmp_x_m(:,ixa,ixbp-1) -  v_tmp_x_mp1(:,ixa,ixbp-1) )
				 		if(ixa>0)  v_tmp_x_m(:,ixa,ixbp+1) = v_tmp_x_m(:,ixa,ixbp+1) + &
										fact(:) * ixa  * ( v_tmp_x_m(:,ixa-1,ixbp) -  v_tmp_x_mp1(:,ixa-1,ixbp) )
			   		endif
		 		enddo
   			enddo
		end subroutine do_recurrence_x
!=========================
		subroutine do_recurrence_y(am1y,am2y)
			implicit none
			integer,intent(in)		::  am1y,am2y
			do iy=1,am1 + am2 - mm
     			    do iya=0,MIN(am1y,iy)
       				iyb=iy-iya
       				if(iyb>am2y) cycle
       				iyam=MAX(iyam,iya)
				    iybm=MAX(iybm,iyb)

				    if(iya>0) then
						iyap=iya-1
					 	v_tmp_y_m(:,iyap+1,iyb) = pa(:,2) * v_tmp_y_m(:,iyap,iyb) - pc(:,2) * v_tmp_y_mp1(:,iyap,iyb)
					 	if(iyap>0) v_tmp_y_m(:,iyap+1,iyb) = v_tmp_y_m(:,iyap+1,iyb) + &
											fact(:) * iyap * ( v_tmp_y_m(:,iyap-1,iyb) -  v_tmp_y_mp1(:,iyap-1,iyb) )
					 	if(iyb>0)  v_tmp_y_m(:,iyap+1,iyb) = v_tmp_y_m(:,iyap+1,iyb) + &
											fact(:) * iyb  * ( v_tmp_y_m(:,iyap,iyb-1) -  v_tmp_y_mp1(:,iyap,iyb-1) )
				   else
					 	iybp=iyb-1
					 	v_tmp_y_m(:,iya,iybp+1) = pb(:,2) * v_tmp_y_m(:,iya,iybp) - pc(:,2) * v_tmp_y_mp1(:,iya,iybp)
					 	if(iybp>0) v_tmp_y_m(:,iya,iybp+1) = v_tmp_y_m(:,iya,iybp+1) + &
										fact(:) * iybp * ( v_tmp_y_m(:,iya,iybp-1) -  v_tmp_y_mp1(:,iya,iybp-1) )
					 	if(iya>0)  v_tmp_y_m(:,iya,iybp+1) = v_tmp_y_m(:,iya,iybp+1) + &
										fact(:) * iya  * ( v_tmp_y_m(:,iya-1,iybp) -  v_tmp_y_mp1(:,iya-1,iybp) )
				  endif
		   		enddo
		        enddo
		end subroutine do_recurrence_y
!=========================
		subroutine do_recurrence_z(am1z,am2z)
			implicit none
				integer,intent(in)		::  am1z,am2z
				do iz=1,am1 + am2 - mm
					do iza=0,MIN(am1z,iz)
						izb=iz-iza
					   	if(izb>am2z) cycle
					   	izam=MAX(izam,iza)
					   	izbm=MAX(izbm,izb)

					   	if(iza>0) then
							izap=iza-1
							v_tmp_z_m(:,izap+1,izb) = pa(:,3) * v_tmp_z_m(:,izap,izb) - pc(:,3) * v_tmp_z_mp1(:,izap,izb)
						 	if(izap>0) v_tmp_z_m(:,izap+1,izb) = v_tmp_z_m(:,izap+1,izb) + &
											fact(:) * izap * ( v_tmp_z_m(:,izap-1,izb) -  v_tmp_z_mp1(:,izap-1,izb) )
						 	if(izb>0)  v_tmp_z_m(:,izap+1,izb) = v_tmp_z_m(:,izap+1,izb) + &
											fact(:) * izb  * ( v_tmp_z_m(:,izap,izb-1) -  v_tmp_z_mp1(:,izap,izb-1) )
					 	else
							izbp=izb-1
						 	v_tmp_z_m(:,iza,izbp+1) = pb(:,3) * v_tmp_z_m(:,iza,izbp) - pc(:,3) * v_tmp_z_mp1(:,iza,izbp)
						 	if(izbp>0) v_tmp_z_m(:,iza,izbp+1) = v_tmp_z_m(:,iza,izbp+1) + &
											fact(:) * izbp * ( v_tmp_z_m(:,iza,izbp-1) -  v_tmp_z_mp1(:,iza,izbp-1) )
						 	if(iza>0)  v_tmp_z_m(:,iza,izbp+1) = v_tmp_z_m(:,iza,izbp+1) + &
											fact(:) * iza  * ( v_tmp_z_m(:,iza-1,izbp) -  v_tmp_z_mp1(:,iza-1,izbp) )
						endif
					enddo
   				enddo
		end subroutine do_recurrence_z
	end subroutine
!=========================================================================
	subroutine destroy_basis_set(basis)
		implicit none
		type(basis_set),intent(inout)   :: basis
		integer			        :: i,j
		if(allocated(basis%bf)) deallocate(basis%bf)
		if(allocated(basis%basis_pair_data)) then
			do i=1,size(basis%basis_pair_data,1)
				do j=1,i
					deallocate( basis%basis_pair_data(i,j)%p,  &
					          & basis%basis_pair_data(i,j)%pa, &
				        	  & basis%basis_pair_data(i,j)%pb, &
			        		  & basis%basis_pair_data(i,j)%c12 )
				end do
			end do
			deallocate(basis%basis_pair_data)
		end if
	end subroutine  destroy_basis_set
!=========================================================================
	subroutine print_basis(basis)
		implicit none
		type(basis_set),intent(inout) :: basis
		integer						  :: ib	
		write(stdout,'(/,a,x,i3)') "Total number of basis functions :",basis%nbf
		write(stdout,'(/,a,x)') "The information of basis functions below :"
		do ib=1,basis%nbas
			write(stdout,'(/,i3,3x,a,i3,3x,a,i3,3x,a,i3,5x,a,i3)') ib,'iatom=',basis%bf(ib)%iatom,&
									& 'ishell=',basis%bf(ib)%ishell,&
                                                                   & 'ibas=',basis%bf(ib)%ibas,'am=',basis%bf(ib)%am
			write(stdout,'(6x,a,l,3x,a,i3,3x,a,i3)') 'bool_sp=',basis%bf(ib)%bool_sp,'ibf_left=',&
									& basis%bf(ib)%ibf_left,'ibf_right=',basis%bf(ib)%ibf_right
		end do 
	end subroutine

end module
