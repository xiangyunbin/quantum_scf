 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
 module m_definitions
    use,intrinsic :: iso_fortran_env, only: OUTPUT_UNIT,ERROR_UNIT

    integer,parameter  :: dp=kind(0.d0)
    integer,parameter  :: sp=kind(0.0)
    integer,parameter  :: stdout = OUTPUT_UNIT
    integer,parameter  :: stderr = ERROR_UNIT

    !
    ! Physical constants
    ! Values from NIST CODATA 2010
    real(dp),parameter    :: Ha_eV        = 27.21138505_dp
    real(dp),parameter    :: bohr_A       = 0.52917721092_dp
    real(dp),parameter    :: c_speedlight = 137.035999074_dp
    !
    ! Mathematical constants
    real(dp),parameter    :: pi         =  3.14159265358979323_dp
    real(dp),parameter    :: pi2        =  pi**2
    real(dp),parameter    :: sqrtpi     =  sqrt(pi)
    real(dp),parameter    :: twopi      =  2.0_dp * pi
    real(dp),parameter    :: sqrt2      =  sqrt(2.0_dp)
    complex(dp),parameter :: im        =  (0.0_dp,1.0_dp)

    !
    integer,parameter :: npar_arth = 16,npara2_arth = 8,npar_poly = 8
    !
    !
    interface assert_eq
        module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
    end interface assert_eq
    !
    !
    interface swap
        module procedure swap_d,swap_v
    end interface swap
    !
    !
    interface arth
        module procedure  arth_r,arth_i
    end interface arth
    !
    !
    interface imaxloc
        module procedure imaxloc_r,imaxloc_i
    end interface imaxloc
    !
    !
    interface iminloc
        module procedure iminloc_r,iminloc_i
    end interface iminloc
    !
    !
    interface outerdiff
        module procedure outerdiff_d,outerdiff_i
    end interface outerdiff

    !
    interface poly
        module procedure poly_r,poly_rv
    end interface

    !
    interface reallocate
	module procedure reallocate_rv,reallocate_rm,&
				reallocate_iv,reallocate_im
    end interface
    contains
!=========================================================================
    !swaps corresponding elements
    !single vaule
    subroutine swap_d(a,b)
        real(dp),intent(inout) :: a,b
        real(dp) :: dum
        dum = a
        a = b
        b = dum
    end subroutine
    !vector value
    subroutine swap_v(a,b)
        real(dp),intent(inout) :: a(:),b(:)
        real(dp) :: dum(size(a))
        dum = a
        a = b
        b = dum
    end subroutine
    !
!=========================================================================
!
!routines for skew operations on matrices
!
    !sets the diagonal elements of a matrix
    subroutine put_diag(mat,diag)
        real(dp),intent(in) :: diag(:)
        real(dp),intent(inout) :: mat(:,:)
        integer :: i
        n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'put_diag')
        do i=1,n
            mat(i,i)=diag(i)
        end do
    end subroutine
    !
    ! gets diagonal of matrix
    function get_diag(mat)
        real(dp),intent(inout) :: mat(:,:)
        real(dp) :: get_diag(min(size(mat,1),size(mat,2)))
        integer :: i
        do i=1,min(size(mat,1),size(mat,2))
            get_diag(i) = mat(i,i)
        end do
    end function
    !
    !multiplies vector into diagonal of a matrix
    subroutine diagmult(mat,diag)
        real(dp),intent(in) :: diag(:)
        real(dp),intent(inout) :: mat(:,:)
        integer :: i
        n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'put_diag')
        do i=1,n
            mat(i,i)=diag(i)*mat(i,i)
        end do
    end subroutine
    !
    !adds vector to diagnonal of a matirx
    subroutine diagadd(mat,diag)
        real(dp),intent(in) :: diag(:)
        real(dp),intent(inout) :: mat(:,:)
        integer :: i
        n = assert_eq(size(diag),min(size(mat,1),size(mat,2)),'put_diag')
        do i=1,n
            mat(i,i)=diag(i) + mat(i,i)
        end do
    end subroutine

    !returns an upper triangular mask
    function upper_triangle(j,k,extra)
        integer,intent(in) :: j,k
        integer,optional,intent(in) :: extra
        logical,dimension(j,k) :: upper_triangle
        !
        integer :: n,jj,kk
        n=0
        if(present(extra)) n=extra
        do jj=1,j
            do kk=1,k
                upper_triangle(jj,kk)=(jj-kk < n)
            end do
        end do
    end function upper_triangle
    !
    !returns a lower triangular mask
    function lower_triangle(j,k,extra)
        integer,intent(in) :: j,k
        integer,optional,intent(in) :: extra
        logical,dimension(j,k) :: lower_triangle
        !
        integer :: n,jj,kk
        n=0
        if(present(extra)) n=extra
        do jj=1,j
            do kk=1,k
                lower_triangle(jj,kk)=(kk-jj < n)
            end do
        end do
    end function lower_triangle

    !returns a unit matrix
    subroutine unit_matrix(mat)
        real(dp),intent(inout) :: mat(:,:)
        integer :: i,n
        n = min(size(mat,1),size(mat,2))
        mat(:,:) = 0.0_dp
        forall(i=1:n) mat(i,i) = 1.0_dp
    end subroutine

    !returns a diagonal matrix for any constant value
    subroutine diag_matrix(mat,val)
        real(dp),intent(out) :: mat(:,:)
        real(dp),intent(in) :: val
        integer :: i,n
        n = min(size(mat,1),size(mat,2))
        mat(:,:) = 0.0_dp
        forall(i=1:n) mat(i,i) = val
    end subroutine
!=========================================================================
    !
    function assert_eq2(n1,n2,string)
        character(len = *),intent(in) :: string
        integer,intent(in) :: n1,n2
        integer :: assert_eq2
        if(n1 == n2) then
            assert_eq2 = n1
        else
            write(stdout,*) 'error : an assert_eq failed with this tag',string
            stop 'program terminated by assert_eq2'
        end if
    end function

    function assert_eq3(n1,n2,n3,string)
        character(len = *),intent(in) :: string
        integer,intent(in) :: n1,n2,n3
        integer :: assert_eq3
        if(n1 == n2 .and. n2 == n3) then
            assert_eq3 = n1
        else
            write(stdout,*) 'error : an assert_eq failed with this tag',string
            stop 'program terminated by assert_eq3'
        end if
    end function

    function assert_eq4(n1,n2,n3,n4,string)
        character(len = *),intent(in) :: string
        integer,intent(in) :: n1,n2,n3,n4
        integer :: assert_eq4
        if(n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq4 = n1
        else
            write(stdout,*) 'error : an assert_eq failed with this tag',string
            stop 'program terminated by assert_eq4'
        end if
    end function

   function assert_eqn(nn,string)
        character(len = *),intent(in) :: string
        integer,dimension(:),intent(in) :: nn
        integer :: assert_eqn
        if(all(nn(2:) == nn(1))) then
            assert_eqn = nn(1)
        else
            write(stdout,*) 'error : an assert_eq failed with this tag',string
            stop 'program terminated by assert_eqn'
        end if
    end function
!=========================================================================
    function ifirstloc(mask)
        logical,dimension(:),intent(in) :: mask
        integer :: ifirstloc
        integer,dimension(1) :: loc
        loc = maxloc(merge(1,0,mask))
        ifirstloc = loc(1)
        if(.not.mask(ifirstloc)) ifirstloc = size(mask)+1
    end function ifirstloc
    !
    !
    function imaxloc_r(arr)
        real(dp),dimension(:),intent(in) :: arr
        integer :: imaxloc_r
        integer,dimension(1) :: imax
        imax = maxloc(arr(:))
        imaxloc_r = imax(1)
    end function imaxloc_r
    !
    !
    function imaxloc_i(arr)
        integer,dimension(:),intent(in) :: arr
        integer :: imaxloc_i
        integer,dimension(1) :: imax
        imax = maxloc(arr(:))
        imaxloc_i = imax(1)
    end function imaxloc_i
    !
    !
    function iminloc_r(arr)
        real(dp),dimension(:),intent(in) :: arr
        integer :: iminloc_r
        integer,dimension(1) :: imin
        imin = minloc(arr(:))
        iminloc_r = imin(1)
    end function iminloc_r
    !
    function iminloc_i(arr)
        integer,dimension(:),intent(in) :: arr
        integer :: iminloc_i
        integer,dimension(1) :: imin
        imin = minloc(arr(:))
        iminloc_i = imin(1)
    end function iminloc_i

!=========================================================================
    function outerand(a,b)
        logical,dimension(:),intent(in) :: a,b
        logical,dimension(size(a),size(b)) :: outerand
        outerand = spread(a,dim = 2,ncopies = size(b)) .and.&
                    spread(b,dim = 1,ncopies = size(a))
    end function outerand
    !
    !
    function outerprod(a,b)
        real(dp),dimension(:),intent(in) :: a,b
        real(dp),dimension(size(a),size(b)) :: outerprod
        outerprod = spread(a,dim = 2,ncopies = size(b)) * &
                    spread(b,dim = 1,ncopies = size(a))
    end function outerprod
    !
    !
    function outerdiv(a,b)
        real(dp),dimension(:),intent(in) :: a,b
        real(dp),dimension(size(a),size(b)) :: outerdiv
        outerdiv = spread(a,dim = 2,ncopies = size(b)) / &
                    spread(b,dim = 1,ncopies = size(a))
    end function outerdiv
    !
    !
    function outersum(a,b)
        real(dp),dimension(:),intent(in) :: a,b
        real(dp),dimension(size(a),size(b)) :: outersum
        outersum = spread(a,dim = 2,ncopies = size(b)) + &
                    spread(b,dim = 1,ncopies = size(a))
    end function outersum
!=========================================================================
    function outerdiff_d(a,b)
        real(dp),dimension(:),intent(in) :: a,b
        real(dp),dimension(size(a),size(b)) :: outerdiff_d
        outerdiff_d = spread(a,dim = 2,ncopies = size(b)) - &
                    spread(b,dim = 1,ncopies = size(a))
    end function outerdiff_d
    
    function outerdiff_i(a,b)
        integer,dimension(:),intent(in) :: a,b
        integer,dimension(size(a),size(b)) :: outerdiff_i
        outerdiff_i = spread(a,dim = 2,ncopies = size(b)) - &
                    spread(b,dim = 1,ncopies = size(a))
    end function outerdiff_i
!=========================================================================
    !
    !see matlab function : linspace
    function linspace(a,b,n) 
        real(dp),intent(in) :: a,b
        integer,intent(in) :: n
        !===
        real(dp) :: dx
        real(dp) :: linspace(n)
        integer :: i
        if(n == 1) then
            linspace(1) = a
            return
        else
            dx = (b-a)/(n-1)
            linspace(1) = a
            do i = 1,n-1
                linspace(i+1) = linspace(i) + dx 
            end do
            return
        end if
    end function linspace
    !
    !
    function arth_r(first,increment,n)
        implicit none
        real(dp),intent(in) :: first,increment
        integer,intent(in) :: n
        real(dp) :: arth_r(n)
        integer :: k,k2
        real(dp) :: temp
        if(n > 0) arth_r(1) = first
        if(n <= npar_arth) then
            do k = 2,n
                arth_r(k) = arth_r(k-1) + increment
            end do
        else
            do k=2,npara2_arth
	        arth_r(k)=arth_r(k-1)+increment
	    end do
	    temp=increment*npara2_arth
	    k=npara2_arth
	    do
	        if (k >= n) exit
		 k2=k+k
		 arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
		 temp=temp+temp
		 k=k2
	     end do
        end if
    end function arth_r
    !
    !
    function arth_i(first,increment,n)
        implicit none
        integer,intent(in) :: first,increment,n
        integer :: arth_i(n)
        integer :: k,k2
        real(dp) :: temp
        if(n > 0) arth_i(1) = first
        if(n <= npar_arth) then
            do k = 2,n
                arth_i(k) = arth_i(k-1) + increment
            end do
        else
            do k=2,npara2_arth
	        arth_i(k)=arth_i(k-1)+increment
	    end do
	    temp=increment*npara2_arth
	    k=npara2_arth
	    do
	        if (k >= n) exit
		 k2=k+k
		 arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
		 temp=temp+temp
		 k=k2
	     end do
        end if
    end function arth_i
    !
    !
    !
    function poly_r(x,coeffs)
        real(dp),intent(in) :: x
        real(dp),intent(in) :: coeffs(:)
        real(dp) :: poly_r
        real(dp) :: pow
        real(dp),allocatable :: vec(:)
        integer :: i,n,nn
        !===
        n = size(coeffs)
        if(n<=0) then
            poly_r = 0.0_dp
        else if(n < npar_poly) then
            poly_r = coeffs(n)
            do i=n-1,1,-1
                poly_r=x*poly_r + coeffs(i)
            end do
        else
            allocate(vec(n+1))
            pow = x
            vec(1:n) = coeffs
            do 
                vec(n+1) = 0.0_dp
                nn=ishft(n+1,-1)
                vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                if (nn == 1) exit
	         pow=pow*pow
		 n=nn
            end do
            poly_r=vec(1)
	    deallocate(vec)
        end if
    end function poly_r
    !
    !
     function poly_rv(x,coeffs)
        real(dp),intent(in) :: x(:)
        real(dp),intent(in) :: coeffs(:)
        real(dp) :: poly_rv(size(x))
        integer :: i,n,m
        !===
        m = size(coeffs)
        n = size(x)
        if(m<=0) then
            poly_rv = 0.0_dp
        else if(m < n .or. m < npar_poly) then
            poly_rv = coeffs(m)
            do i=m-1,1,-1
                poly_rv=x*poly_rv + coeffs(i)
            end do
        else
            do i=1,n
                poly_rv(i) = poly_r(x(i),coeffs)
            end do
        end if
    end function poly_rv
!=========================================================================
	function reallocate_rv(p,n)
		real(dp),dimension(:),pointer :: p,reallocate_rv
		integer,intent(in) :: n
		integer :: nold,ierr
		allocate(reallocate_rv(n),stat=ierr)
		if(ierr /= 0) stop 'reallocate_rv: problem in attempt to allocate memory'
		if(.not.associated(p)) return
		nold = size(p)
		reallocate_rv(1:min(nold,n)) = p(1:min(nold,n))
		deallocate(p)
	end function reallocate_rv

	function reallocate_iv(p,n)
		integer,dimension(:),pointer :: p,reallocate_iv
		integer,intent(in) :: n
		integer :: nold,ierr
		allocate(reallocate_iv(n),stat=ierr)
		if(ierr /= 0) stop 'reallocate_iv: problem in attempt to allocate memory'
		if(.not.associated(p)) return
		nold = size(p)
		reallocate_iv(1:min(nold,n)) = p(1:min(nold,n))
		deallocate(p)
	end function reallocate_iv

	function reallocate_rm(p,n,m)
		real(dp),dimension(:,:),pointer :: p,reallocate_rm
		integer,intent(in) :: n,m
		integer :: nold,mold,ierr
		allocate(reallocate_rm(n,m),stat=ierr)
		if(ierr /= 0) stop 'reallocate_rm: problem in attempt to allocate memory'
		if(.not.associated(p)) return
		nold = size(p,1)
		mold = size(p,2)
		reallocate_rm(1:min(nold,n),1:min(mold,m)) = &
			p(1:min(nold,n),1:min(mold,m))
		deallocate(p)
	end function reallocate_rm

	function reallocate_im(p,n,m)
		integer,dimension(:,:),pointer :: p,reallocate_im
		integer,intent(in) :: n,m
		integer :: nold,mold,ierr
		allocate(reallocate_im(n,m),stat=ierr)
		if(ierr /= 0) stop 'reallocate_rm: problem in attempt to allocate memory'
		if(.not.associated(p)) return
		nold = size(p,1)
		mold = size(p,2)
		reallocate_im(1:min(nold,n),1:min(mold,m)) = &
			p(1:min(nold,n),1:min(mold,m))
		deallocate(p)
	end function reallocate_im
end module m_definitions
