 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_eri
        use m_definitions
        use m_tools
        use m_atoms
        use m_shells
        use m_basis
        use m_rysroot


        type eri_quartet
            integer                     ::  nbf1,nbf2,nbf3,nbf4
            real(dp),pointer            ::  eri_data(:) 
        end type

        type basis_quartet_set
            integer			         ::	 nsize_quartet
            integer                          ::     nbas
            type(eri_quartet),allocatable   ::     eri_buffer(:)		
        end type basis_quartet_set

	contains
!====================================================================
        subroutine precalc_eri(shells,basis,basisqs,print_)
                implicit none
                type(shell_set),intent(in)				::  shells
                type(basis_set),intent(in)				::  basis
                type(basis_quartet_set),intent(out)		::  basisqs
                logical,optional,intent(in)                     ::  print_
                !==
                integer         ::      s1,s2,s3,s4,s4_max
                integer         ::      nbas,nbas_pair,nsize
                integer         ::      nbf1,nbf2,nbf3,nbf4,nbf_tot
                integer         ::      iquartet
                nbas           = 	basis%nbas
                nbas_pair      =	nbas * (nbas + 1)/2
                basisqs%nbas   =     nbas
                basisqs%nsize_quartet   =   nbas_pair * (nbas_pair + 1)/2
                allocate( basisqs%eri_buffer( basisqs%nsize_quartet ) )
                iquartet = 1
                if(present(print_))  write(stdout,'(/,a)'), '=========== ERI DATA ============'
                do s1=1,nbas
                    nbf1      =     basis%bf(s1)%nbf
                    do s2=1,s1
                        nbf2      =     basis%bf(s2)%nbf
                        do s3=1,s1
                            nbf3      =     basis%bf(s3)%nbf
                            if( s1 == s3) then
                                s4_max = s2					
                            else
                                s4_max = s3 
                            end if
                            do s4=1,s4_max
                                nbf4      =     basis%bf(s4)%nbf
                                basisqs%eri_buffer(iquartet)%nbf1 = nbf1
                                basisqs%eri_buffer(iquartet)%nbf2 = nbf2
                                basisqs%eri_buffer(iquartet)%nbf3 = nbf3
                                basisqs%eri_buffer(iquartet)%nbf4 = nbf4
                            
                                nsize = nbf1 * nbf2 * nbf3 * nbf4
                                allocate( basisqs%eri_buffer(iquartet)%eri_data(nsize) )

                                call do_single_eri( basis%bf(s1),&
                                                 &  basis%bf(s2),&
                                                 &  basis%bf(s3),&
                                                 &  basis%bf(s4),&
                                                 &  basisqs%eri_buffer(iquartet)%eri_data)

                                iquartet = iquartet + 1
                            end do
		        end do
                    end do
                end do
                contains
                subroutine do_single_eri(bf1,bf2,bf3,bf4,buffer)
                    implicit none
                    type(basis_function),intent(in)     ::  bf1,bf2,bf3,bf4
                    real(dp)                             ::   buffer(:) 
                    !===
                    integer         ::      ns,am_total,order
                    integer         ::      ibas1,ibas2,ibas3,ibas4
                    integer         ::      ish1,ish2,ish3,ish4
                    integer         ::      am1,am2,am3,am4
                    integer         ::      np,ip,ipp

                    ibas1 	 =     basis%bf(s1)%ibas
                    ibas2 	 =     basis%bf(s2)%ibas
                    ibas3 	 =     basis%bf(s3)%ibas
                    ibas4 	 =     basis%bf(s4)%ibas

                    ish1  	 =     basis%bf(s1)%ishell
                    ish2  	 =     basis%bf(s2)%ishell
                    ish3  	 =     basis%bf(s3)%ishell
                    ish4  	 =     basis%bf(s4)%ishell

                    am1 = bf1%am
                    am2 = bf2%am
                    am3 = bf3%am
                    am4 = bf4%am
                    am_total = am1 + am2 + am3 + am4
                    if(am_total == 0) then
                        call eri_ssss( basis%basis_pair_data(ibas1,ibas2),&
                                    &  basis%basis_pair_data(ibas3,ibas4), &
                                    &  shells%shell_pair_data(ish1,ish2), &
                                    &  shells%shell_pair_data(ish3,ish4),&
                                    &  buffer)
                        if(present(print_)) write(stdout,'(i3,4x,i3,4x,i3,4x,i3,4x,f12.6)'),ibas1,ibas2,ibas3,ibas4,buffer(:)
                    else
                        order = mod(am_total,2)
                        order = (am_total - order)/2 + 1
                        call eri_rys(  basis%basis_pair_data(ibas1,ibas2),&
                                    &  basis%basis_pair_data(ibas3,ibas4), &
                                    &  shells%shell_pair_data(ish1,ish2), &
                                    &  shells%shell_pair_data(ish3,ish4),&
                                    &  order,&
                                    &  buffer)
                        if(present(print_)) then
                            np = nsize/8 + 1
                            ipp = 1
                            do ip = 1,np
                                if(np == 1) then
                                    write(stdout,'(i3,4x,i3,4x,i3,4x,i3,4x,8f12.6)'),&
                                                    ibas1,ibas2,ibas3,ibas4,buffer(1:nsize)
                                else
                                    if(ip == 1) then
                                           write(stdout,'(i3,4x,i3,4x,i3,4x,i3,4x,8f12.6)'),&
                                                    ibas1,ibas2,ibas3,ibas4,buffer(ipp : ipp + 7)
                                            ipp = ipp + 8
                                    else if (ip == np) then
                                        write(stdout,'(28x,8f12.6)'),&
                                                        buffer(ipp : nsize)
                                    else
                                        write(stdout,'(28x,8f12.6)'),&
                                                        buffer(ipp : ipp +7)
                                        ipp = ipp + 8
                                    end if
                                end if

                            end do

                        end if
                    end if
                end subroutine do_single_eri
	end subroutine precalc_eri
!============================================================
        subroutine eri_ssss(bfp1,bfp2,shp1,shp2,buffer)
            implicit none
            type(basis_pair),intent(in)     ::      bfp1,bfp2
            type(shell_pair),intent(in)     ::      shp1,shp2
            real(dp),intent(out)             ::      buffer(:)
            !==
            real(dp),parameter     ::      const_eri = 2.0_dp * pi**2.5_dp
            real(dp)                ::      rho(bfp1%nc12 * bfp2%nc12)
            real(dp)                ::      tt(bfp1%nc12 * bfp2%nc12)
            real(dp)                ::      vboys(bfp1%nc12 * bfp2%nc12,1) 
            real(dp)                ::      pa(bfp1%nc12 * bfp2%nc12,3)   
            real(dp)                ::      qc(bfp1%nc12 * bfp2%nc12,3)
            real(dp)                ::      pq(bfp1%nc12 * bfp2%nc12,3)
            integer                ::      i,j,ij,nc12,nc34

            nc12 = bfp1%nc12
            nc34 = bfp2%nc12
            forall(i=1:nc12,j=1:nc34)
                  rho(i+(j-1)*nc12)  =  shp1%zeta(i) * shp2%zeta(j)/( shp1%zeta(i) + shp2%zeta(j) )
                  pq(i+(j-1)*nc12,:) =  bfp1%p(i,:) - bfp2%p(j,:)
            end forall
            tt(:) = rho(:) * ( pq(:,1)**2 + pq(:,2)**2 + pq(:,3)**2)
            call boys_function( vboys,nc12 * nc34,0,tt(:))
            buffer(1) = 0.0_dp
            do i=1,nc12
                do j=1,nc34
                    ij = i+(j-1)*nc12
                    buffer(1) = buffer(1) + bfp1%c12(i,1) * bfp2%c12(j,1)  * vboys(ij,1) / &
                                       ( shp1%zeta(i) * shp2%zeta(j) * sqrt( shp1%zeta(i) + shp2%zeta(j) ))
                end do
            end do
            buffer(1) = const_eri * buffer(1)
        end subroutine eri_ssss
!====================================================================
!Rys quadrature
        subroutine eri_rys(bfp1,bfp2,shp1,shp2,norder,buffer)
            implicit none
            type(basis_pair),intent(in)     ::      bfp1,bfp2
            type(shell_pair),intent(in)     ::      shp1,shp2
            integer,intent(in)              ::      norder
            real(dp),intent(out)             ::      buffer(:)
            !===
            real(dp),parameter              ::  threshold = 1.0e-10
            real(dp) 				 ::  rho(bfp1%nc12 * bfp2%nc12)	
            real(dp) 				 ::  tt(bfp1%nc12 * bfp2%nc12)
            real(dp)				 ::  pa(bfp1%nc12 * bfp2%nc12,3)
            real(dp)				 ::  qc(bfp1%nc12 * bfp2%nc12,3)
            real(dp)				 ::  pq(bfp1%nc12 * bfp2%nc12,3)
            real(dp)                         ::  roots(norder),weights(norder)
            real(dp)                         ::  Ir(0:bfp1%amp,0:bfp2%amp,3)  !Rys Matrix,3 direction
            real(dp)                         ::  ab_x0(3,0:bfp1%am2),cd_x0(3,0:bfp2%am2)
            real(dp)                         ::  tmp(size(buffer))
            logical                         ::  is_nonzero(size(buffer))
            integer                         ::  nc12,nc34,am1,am2,am3,am4
            integer                         ::  ip,jp,i_nc,j_nc,inc
            integer                         ::  iorder
            
            nc12 = bfp1%nc12
            nc34 = bfp2%nc12
            am1 = bfp1%am1
            am2 = bfp1%am2
            am3 = bfp2%am1
            am4 = bfp2%am2
            forall(ip=0:am2) ab_x0(:,ip) = bfp1%AB**ip
            forall(ip=0:am4) cd_x0(:,ip) = bfp2%AB**ip
            forall(ip=1:nc12,jp=1:nc34)
                rho(ip+(jp-1)*nc12)  =  shp1%zeta(ip) * shp2%zeta(jp)/( shp1%zeta(ip) + shp2%zeta(jp) )
                pa(ip+(jp-1)*nc12,:) =  bfp1%pa(ip,:)
                qc(ip+(jp-1)*nc12,:) =  bfp2%pa(jp,:)
                pq(ip+(jp-1)*nc12,:) =  bfp1%p(ip,:) - bfp2%p(jp,:)
            end forall
            tt(:) = rho(:) * ( pq(:,1)**2 + pq(:,2)**2 + pq(:,3)**2)
            
            !
            is_nonzero(:) = .true.
            buffer(:) = 0.0_dp
            do i_nc=1,nc12
                do j_nc=1,nc34
                    inc = i_nc+(j_nc-1)*nc12
                    call rysquad(norder,tt(inc),roots,weights)
                    do iorder = 1,norder
                        call rys_recurrence(roots(iorder),weights(iorder),tmp)
                        where(is_nonzero) buffer = buffer + tmp
                    end do
                    if(inc == 1) is_nonzero = (abs(tmp) > threshold)
                end do
            end do
            contains
            !===================================================
            subroutine rys_recurrence(root,weight,buf)  
                implicit none
                real(dp),intent(in)     ::  root,weight
                real(dp),intent(out)    ::  buf(:)
                !===
                integer                 ::  ii,jj,kk,ll
                integer                 ::  ij_bra,kl_ket
                integer                 ::  isp1,isp2,isp3,isp4
                integer                 ::  amm1,amm2,amm3,amm4
                integer                 ::  imx,imy,imz
                integer                 ::  jmx,jmy,jmz
                integer                 ::  kmx,kmy,kmz
                integer                 ::  lmx,lmy,lmz
                integer                 ::  bf1_nx,bf1_ny,bf1_nz
                integer                 ::  bf2_nx,bf2_ny,bf2_nz
                integer                 ::  bf3_nx,bf3_ny,bf3_nz
                integer                 ::  bf4_nx,bf4_ny,bf4_nz
                integer                 ::  icount


                call build_Ir(root)
                icount = 1
                ![ab|..] path
                isp1 = bfp1%isp1
                isp2 = bfp1%isp2
                do ii=1,isp1
                  !
                  ![ab|..] -> a
                  amm1 = am1
                  if(isp1 == 2 .and. ii == 1) amm1 = 0
                  do imx=0,amm1
                    bf1_nx = amm1 - imx
                    do imy=0,imx
                      bf1_ny = imx - imy
		       bf1_nz = imy
                      do jj=1,isp2
                        ij_bra = ii+(jj-1)*isp1
                        !
                        ![ab|..] -> b
			  amm2 = am2
                        if(isp2 == 2 .and. jj == 1) amm2 = 0
                        do jmx=0,amm2
                          bf2_nx = amm2 - jmx
                          do jmy=0,jmx
                            bf2_ny = jmx - jmy
                            bf2_nz = jmy
                            !
                            ![ab|cd] path
                            isp3 = bfp2%isp1
                            isp4 = bfp2%isp2
                            do kk=1,isp3
                              !
                              ![ab|cd] -> c
                              amm3 = am3
                              if(isp3 == 2 .and. kk == 1) amm3 = 0
                              do kmx=0,amm3
                                bf3_nx = amm3 - kmx
                                do kmy=0,kmx
                                  bf3_ny = kmx - kmy
		                   bf3_nz = kmy
                                  do ll=1,isp4
                                    kl_ket = kk+(ll-1)*isp3
                                    !
                                    ![ab|cd] ->d
			              amm4 = am4
                                    if(isp4 == 2 .and. ll == 1) amm4 = 0
                                    do lmx=0,amm4
                                      bf4_nx = amm4 - lmx
                                      do lmy=0,lmx
                                        bf4_ny = lmx - lmy
                                        bf4_nz = lmy

                                        !
                                        if( is_nonzero(icount)  .eqv. .false.) then
                                            icount = icount + 1
                                            cycle
                                        end if
                                        !
                                        buf(icount) = &
                                             & rys_I(1,bf1_nx,bf2_nx,bf3_nx,bf4_nx) * &
					        & rys_I(2,bf1_ny,bf2_ny,bf3_ny,bf4_ny) * &
					        & rys_I(3,bf1_nz,bf2_nz,bf3_nz,bf4_nz) * &
                                             & weight
                                        !inc : global variable
                                        buf(icount) = 2.0_dp * bfp1%c12(i_nc,ij_bra) * bfp2%c12(j_nc,kl_ket) * &
                                                    & sqrt(rho(inc)/pi) * buf(icount)
                                        icount = icount + 1
                                        !
                                    end do
                                  end do
                                end do
                              end do
                            end do
                           end do
                           ![|cd]

                          end do
                        end do
                      end do
                    end do
                  end do
                end do
                ![ab|cd]
            end subroutine rys_recurrence
            !===================================================
            subroutine  build_Ir(rt)
                real(dp),intent(in) :: rt
                !==
                real(dp) :: rootab,zeta_ab
                real(dp) :: xab(3),xxa(3),xxb(3),bxba(3),axab(3)
                real(dp) :: zeta_12,zeta_34
                real(dp) :: B00,B10,B10p,C00(3),C00p(3)
                integer :: am_ij,am_kl,n,m
                !===
                am_ij       =   am1 + am2
                am_kl       =   am3 + am4
                ! inc,i_nc,j_nc : global variable
                ! Rys coefficients 
                zeta_12     =   shp1%zeta(i_nc)
                zeta_34     =   shp2%zeta(j_nc)
                zeta_ab     =   1.0_dp / (zeta_12 + zeta_34)
                xab(:)      =   pq(inc,:)
                xxa(:)      =   pa(inc,:)
                xxb(:)      =   qc(inc,:)
                bxba(:)     =   -zeta_34 * xab(:)
                axab(:)     =   zeta_12 * xab(:)
                rootab      =   rt * zeta_ab
                B00         =   0.5_dp*rootab
                B10         =   (0.5_dp - B00 * zeta_34)/zeta_12
                B10p        =   (0.5_dp - B00 * zeta_12)/zeta_34
                C00(:)      =   xxa(:) + bxba(:) * rootab
                C00P(:)     =   xxb(:) + axab(:) * rootab

		
                Ir(0,0,:) = pi / sqrt(zeta_12 * zeta_34)
                if( am_ij .gt. 0 ) Ir(1,0,:) = C00(:) * Ir(0,0,:)
                if( am_kl .gt. 0 ) Ir(0,1,:) = C00p(:) * Ir(0,0,:)
                if(am_ij .gt. 1) then
                    do n=1,am_ij-1
                        Ir(n+1,0,:) = n * B10 * Ir(n-1,0,:) + C00(:) * Ir(n,0,:)
                    end do
		 end if
                if(am_kl .gt. 1) then
                    do m=1,am_kl-1
                        Ir(0,m+1,:) = m * B10p * Ir(0,m-1,:) + C00p(:) * Ir(0,m,:)
                        end do
                    end if
                if( am_ij .ne. 0 .and. am_kl .ne. 0 ) then
                    do n=0,am_ij-1
                        Ir(n+1,1,:) = (n+1) * B00 * Ir(n,0,:) + C00p(:) * Ir(n+1,0,:)
                        do m=1,am_kl-1
                            Ir(n+1,m+1,:) = m * B10p * Ir(n+1,m-1,:) + &
			        & (n+1) * B00 * Ir(n,m,:) + C00p(:) * Ir(n+1,m,:)
                        end do
                    end do
                end if
            end subroutine build_Ir
            !===============================================
            function rys_I(id,am_i,am_j,am_k,am_l)
                implicit none
                integer,intent(in)  ::  id,am_i,am_j,am_k,am_l
                !===
                real(dp)         :: rys_I,ixjxm
                integer         ::  n,m,ind

                rys_I = 0.0_dp
                do m=0,am_l
                    ind = am_k + m
                    ixjxm = 0.0_dp
                    do n=0,am_j
                        ixjxm = ixjxm + bico(am_j,n) * &
				& Ir(am_i + n,am_k + m,id) * ab_x0(id,am_j - n)
                    end do
                    rys_I = rys_I + bico(am_l,m) * &
		        & ixjxm * cd_x0(id,am_l - m)
                end do
               
            end function rys_I

       end subroutine eri_rys
!====================================================================
        subroutine destroy_basis_quartet_set(basisqs)
            implicit none
            type(basis_quartet_set),intent(inout)      ::  basisqs
            !===
            integer             ::  i
            if(allocated(basisqs%eri_buffer)) then
                do i=1,size(basisqs%eri_buffer)
                    deallocate(basisqs%eri_buffer(i)%eri_data)
                end do
                deallocate(basisqs%eri_buffer)
            end if
        end subroutine destroy_basis_quartet_set
!====================================================================
end module m_eri


