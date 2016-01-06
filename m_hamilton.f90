 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_hamilton
    use m_definitions
    use m_atoms
    use m_basis
    use m_shells
    use m_tools,only : tred2,tqli_mat
    use m_eri

    contains
!=============================================================================
    subroutine setup_fock_matirx(dmat,basis,basisqs,hmat)
        implicit none
        real(dp),intent(in)                   ::  dmat(:,:)
        type(basis_set),intent(in)            ::  basis
        type(basis_quartet_set),intent(in)    ::  basisqs
        real(dp),intent(out)                  ::   hmat(:,:)
        !==
        real(dp),parameter  ::  threshold = 1.0e-10
        real(dp),pointer    ::  buffer(:)
        real(dp)             ::  G(size(dmat,1),size(dmat,2))
        integer             ::  n,nbas
        integer             ::  iquartet
        integer             ::  s1,s2,s3,s4,s4_max
        integer             ::  bf1_first,bf2_first,bf3_first,bf4_first
        integer             ::  n1,n2,n3,n4
        integer             ::  bf1,bf2,bf3,bf4
        integer             ::  f1,f2,f3,f4,f1234
        real(dp)             ::  s12_deg,s34_deg,s1234_deg,s12_34_deg
        real(dp)             ::  value_scal_by_deg

        nbas = basisqs%nbas
        iquartet = 1
        G(:,:) = 0.0_dp
        do s1 = 1,nbas
            bf1_first   =   basis%bf(s1)%ibf_left
            n1          =   basis%bf(s1)%nbf
            do s2 = 1,s1
                bf2_first   =   basis%bf(s2)%ibf_left
                n2          =   basis%bf(s2)%nbf
                do s3 = 1,s1
                    bf3_first   =   basis%bf(s3)%ibf_left
                    n3          =   basis%bf(s3)%nbf
                    if(s1 == s3) then
                        s4_max = s2
                    else
                        s4_max = s3
                    end if
                    do s4 = 1,s4_max
                        bf4_first   =   basis%bf(s4)%ibf_left
                        n4          =   basis%bf(s4)%nbf

                        s12_deg = 2.0_dp
                        if(s1 == s2) s12_deg = 1.0_dp
                        !
                        s34_deg = 2.0_dp
                        if(s3 == s4) s34_deg = 1.0_dp
                        !
                        if(s1 == s3) then
                            if(s2 == s4) then
                                s12_34_deg = 1.0_dp
                            else
                                s12_34_deg = 2.0_dp
                            end if
                        else
                            s12_34_deg = 2.0_dp
                        end if
                        !
                        s1234_deg = s12_deg * s34_deg * s12_34_deg
                        buffer => basisqs%eri_buffer(iquartet)%eri_data(:)
                        f1234 = 1
                        do f1 = 0,n1-1
                            bf1 = f1 + bf1_first
                            do f2 = 0,n2-1
                                bf2 = f2 + bf2_first
                                do f3 = 0,n3-1
                                    bf3 = f3 + bf3_first
                                    do f4 = 0,n4-1
                                        bf4 = f4 + bf4_first
                                        if(abs(buffer(f1234)) < threshold) then
                                            f1234 = f1234 + 1
                                            cycle
                                        end if
                                        
                                        value_scal_by_deg = s1234_deg * buffer(f1234)
                                        G( bf1,bf2 ) = G( bf1,bf2 ) + dmat( bf3,bf4 ) * value_scal_by_deg
					   G( bf3,bf4 ) = G( bf3,bf4 ) + dmat( bf1,bf2 ) * value_scal_by_deg
					   G( bf1,bf3 ) = G( bf1,bf3 ) - 0.25_dp * dmat( bf2,bf4 ) * value_scal_by_deg
                                        G( bf2,bf4 ) = G( bf2,bf4 ) - 0.25_dp * dmat( bf1,bf3 ) * value_scal_by_deg
					   G( bf1,bf4 ) = G( bf1,bf4 ) - 0.25_dp * dmat( bf2,bf3 ) * value_scal_by_deg
					   G( bf2,bf3 ) = G( bf2,bf3 ) - 0.25_dp * dmat( bf1,bf4 ) * value_scal_by_deg
                                        f1234 = f1234 + 1
                                    end do
                                end do
                            end do
                        end do !(ab|cd) loop
                        iquartet = iquartet + 1
                    end do

                end do
            end do
        end do !basis loop
        hmat = 0.5_dp * ( G + transpose(G))
        !
    end subroutine setup_fock_matirx
!=============================================================================
    subroutine seutp_S_minus_sqrt(S,S_minus_sqrt,print_)
       implicit none
       real(dp),intent(in)      ::      S(:,:)
       real(dp),intent(inout)   ::      S_minus_sqrt(:,:)
       !===
       real(dp)                 ::      tmp(size(S,1),size(S,2))
       real(dp)                 ::      d(size(S,1)),e(size(S,1))
       integer                 ::      i,j,n
       logical,optional       ::      print_
       n = size(S,1)
       S_minus_sqrt = S
       call tred2(S_minus_sqrt,d,e)
       call tqli_mat(d,e,S_minus_sqrt)
       !S^(-1/2) = E * d^(-1/2) * transpose(E)
       forall(i=1:n) tmp(i,:) =  S_minus_sqrt(:,i)/sqrt(d(i))
       S_minus_sqrt = matmul( S_minus_sqrt , tmp)
       if(present(print_)) then
           write(stdout,'(/,a)'), '=========== S^(-1/2) ============'
           write(stdout,'(/)')
           do i=1,n
                write(stdout,'(i3,4x,100f12.6)'),i, S_minus_sqrt(i,:)
           end do
           write(stdout,'(/)')
        end if
    end subroutine 
!=============================================================================
!Solve generalized eigenvalue problem  : HC = SCE
!{S^(-1/2) * H * S^(-1/2)} * S^(1/2) * C = S^(1/2) * C * E
!=> H' = S^(-1/2) * H * S^(-1/2) , X = S^(1/2) * C
!=> H'X = XE
!=> C =  S^(-1/2) * X
    subroutine generalized_eigen_solver(H,S_minus_sqrt,D,ndocc,print_)
        implicit none
        real(dp),intent(in)     ::      H(:,:),S_minus_sqrt(:,:)
        real(dp),intent(inout)  ::      D(:,:)
        integer,intent(in)     ::      ndocc
        logical,optional       ::      print_
        !===
        real(dp)            ::  Hp(size(H,1),size(H,2))
        real(dp)            ::  dd(size(H,1)),ee(size(H,1))
        integer            ::  i
        Hp = matmul(S_minus_sqrt,matmul(H,S_minus_sqrt))
        call tred2(Hp,dd,ee)
        call tqli_mat(dd,ee,Hp)
        call eigsort(dd,Hp)
        D = matmul(S_minus_sqrt,Hp)
        !D = C * C'
        D = matmul(D(:,1:ndocc),transpose(D(:,1:ndocc)))
        if(present(print_)) then
            write(stdout,'(/,a)'), '=========== Density Matrix ============'
            do i=1,size(dd)
                 !FIXME
                 write(stdout,'(i3,4x,100f12.6)'),i,D(i,:)
           end do
           write(stdout,'(/)')
           write(stdout,'(a)'), '=== Energies ==='
           write(stdout,'(2x,a,7x,a)'), '#','[Ha]'
           do i=1,size(dd)
                write(stdout,'(i3,4x,f12.6)'),i,dd(i)
                if(i == ndocc) write(stdout,'(a)'),'---------------------------'
           end do
           write(stdout,'(/)')
        end if
    end subroutine
!=============================================================================
end module m_hamilton

	
