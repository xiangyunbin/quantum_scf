 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*

module m_rysquad
    use m_definitions
    
    real(dp),allocatable,private :: a(:),b(:)
    integer,private :: j_index

    contains
!====================================================================
    !rys quadrature for large T( 50 <= T)
    !Nmax < 20
     subroutine rysq_asym(t,x,w)
        implicit none
        real(dp),intent(in) :: t
        real(dp),intent(out) :: x(:),w(:)
        !===
        integer :: n,i
        real(dp) :: ot,ot2,o2t,o2t2,et,tmp,amu,expmt
        real(dp) :: alpha(2*size(x) - 1),beta(2*size(x)-1)
        real(dp) :: anu(2*size(x))
        real(dp) :: a(size(x)),b(size(x))
        n = size(x)
        ot = 1.0_dp / t
        ot2 = ot / t
        o2t = ot / 2.0_dp
        o2t2 = o2t / t
        alpha(1) = o2t
        do i = 2,2*n-1
            alpha(i) = (4.0_dp * i - 3.0_dp)*o2t
            beta(i) = (i-1.0_dp)*(2.0_dp*i-3.0_dp)*o2t2
        end do
        amu = boys_val(0,t) 
        anu(1) = 1.0_dp
        do i = 2,2*n-1
            if(i == 2) then
                anu(2) = 1.0_dp - 1.5_dp * ot
            else
                anu(i) = (1.0_dp - (2.0_dp * i - 2.5_dp) * ot)*anu(i-1) &
                        -(i-2.0_dp)*(i-1.5_dp)*ot2*anu(i-2)
            end if
        end do
        et = exp(-t)
        expmt = -et * o2t
        do i = 2*n,2,-1
            anu(i) = expmt * anu(i-1)
        end do
        anu(1) = amu
        call orthog(anu,alpha,beta,a,b)
        call gausscof(a,b,amu,x,w)
     end subroutine
!====================================================================
    !rys quadrature for small T(0<= T <= 5)
    !Miller's algorithm
    subroutine rysq_small(t,x,w)
        implicit none
        real(dp),intent(in) :: t
        real(dp),intent(out) :: x(:),w(:)
        !===
        integer :: n,i,k
        real(dp) :: ri,si,miu0,i2,i4,scal
        real(dp) :: ot,ot2,o2t,o2t2
        real(dp) :: alpha(2*size(x) - 1),beta(2*size(x) - 1)
        real(dp),allocatable :: anu(:)
        real(dp) :: a(size(x)),b(size(x))
        n = size(x)
        if(t < 0.01_dp) then
            call gauss_shift_jacobi(x,w)
            return
        else 
            ot = 1.0_dp / t
            ot2 = ot / t
            o2t = ot / 2.0_dp
            o2t2 = o2t / t
            alpha(1) = o2t
            alpha(1) = 1.0_dp/3.0_dp
            do i = 2,2*n-1
                i2 = 2.0_dp * i
                alpha(i) = ( 2.0_dp*(i-1.0_dp)*(i-0.5_dp) - 0.25_dp)/((i2-0.5_dp)*(i2 - 2.5_dp))
                beta(i) = ((i - 1.0_dp)**2 * (i - 1.5_dp)**2 )&
                    /((i2 - 2.5_dp)**2 * (i2 - 3.5_dp) * (i2 - 1.5_dp ))
            end do
            if(n < 5) then
                k = 2*n + 10
            else
                k = 2*n
            end if
            allocate(anu(k))
            anu(k) = 1.0_dp
            anu(k - 1) = 0.0_dp
            miu0 = boys_val(0,t) 
            do i = k-1,2,-1
                i2 = 2.0_dp * (i-1.0_dp)
                i4 = 4.0_dp * (i-1.0_dp)
                ri = (i2 + 1.0_dp) * o2t + (i2 + 1.0_dp)/&
                        ((i4 + 3.0_dp) * (i4 - 1.0_dp))
                si = ((i4 + 1.0_dp)*(i4 - 1.0_dp)**2 *(i4 - 3.0_dp))/&
                        (i2 * (i2 + 1.0_dp)*(i2 - 1.0_dp)**2)
                anu(i-1) = - (ri*si) * anu(i) + si*anu(i+1)
            end do
            scal = miu0 / anu(1)
            anu(1:2*n) = anu(1:2*n) * scal
            call orthog(anu(1:2*n),alpha,beta,a,b)
            call gausscof(a,b,miu0,x,w)
            deallocate(anu)
        end if
    end subroutine
!====================================================================
    !rys quadrature for intermediate T(5< T < 60)
    !LU-decomposition algorithm
    subroutine rysq_mid(t,x,w)
        implicit none
        real(dp),intent(in) :: t
        real(dp),intent(out) :: x(:),w(:)
        !===
        integer :: n,i,k
        real(dp) :: ri,r1,miu0,i2,i2p
        real(dp) :: o2t,u1
        real(dp) :: alpha(2*size(x) - 1),beta(2*size(x) - 1)
        real(dp) :: a(size(x)),b(size(x)) 
        real(dp),allocatable :: osn(:),xn(:),oogam(:),anu(:)
        n = size(x)
        alpha(1) = 1.0_dp/3.0_dp
        do i = 2,2*n-1
            i2 = i - 0.5_dp
            i2p = 2.0_dp * i - 1.5_dp
            alpha(i) = (2.0_dp * (i-1.0_dp) * i2 - 0.25_dp)/&
                            ((i2p + 1.0_dp)*(i2p - 1.0_dp)) 
            beta(i) = ((i-1.0_dp)*(i2 - 0.5_dp)*(i - 1.5_dp)*(i2 - 1.0_dp))/&
                            ((i2p - 1.0)**2 * (i2p - 2.0_dp)*i2p)
        end do
        if(n < 5) then
             k = 2*n + 20
        elseif(n < 10) then
             k = 2*n + 15
        else
             k = 2*n + 10 
        end if
        allocate(osn(k))
        allocate(xn(k))
        allocate(oogam(k))
        allocate(anu(k))
        o2t = 0.5_dp / t
        r1 = 3.0_dp * o2t + 1.0_dp/7.0_dp
        osn(1) = 7.5_dp
        miu0 = boys_val(0,t)
        oogam(1) = -1.0_dp / (r1 * osn(1))
        xn(1) =  - miu0
        do i = 2,k
            i2 = 2.0_dp*i + 1.0_dp
            i2p = 2.0_dp*i + 0.5_dp
            ri = i2 * o2t + 0.25_dp * i2 /((i2p + 1.0_dp)*(i2p - 1.0_dp))
            osn(i) = 4.0_dp * i2p * (i2p - 1.0_dp)**2 *(i2p - 2.0_dp)/&
                    (i * i2 * (i2 - 2.0_dp) * (i - 0.5_dp))
            oogam(i) = - ri * osn(i) + oogam(i-1)*osn(i-1)
            oogam(i) = 1.0_dp / oogam(i)
            !LU
            xn(i) = oogam(i-1)*xn(i-1)
        end do
        !LU
        anu(k) = -xn(k-1)*oogam(k-1)
        do i = k-1,2,-1
            anu(i) = -(anu(i+1) * osn(i-1) + xn(i-1))*oogam(i-1)
        end do
        anu(1) = miu0
        call orthog(anu(1:2*n),alpha,beta,a,b)
        call gausscof(a,b,miu0,x,w)
        deallocate(osn,xn,oogam,anu)
    end subroutine
!====================================================================
    function boys_val(n,t)
        implicit none
        integer,intent(in)   :: n
        real(dp),intent(in)  :: t
        !=====
        real(dp)  :: boys_val
        integer,parameter  :: maxfac=100
        real(dp),parameter :: eps=1.0d-17
        integer :: i,m,k
        integer :: m2
        real(dp) :: t2,num,sum,term1,et,tt
        real(dp),parameter :: kk = 0.5_dp * SQRT( pi ) 
        real(dp),save :: df(2*maxfac)=0.0_dp
        !=====
        if(n>maxfac) stop ' boys function Fm(t) for a too high m value'
        if( ABS(df(1))<1.d-10 ) then
            df(1:3) = 1.0_dp
            do i=4,2*maxfac
                df(i) = (i-2) * df(i-2)
            enddo
        endif
        if( t > 20.0 ) then ! For big t's do upward recursion 
            t2 = 2 * t
            et = exp(-t)
            tt  = sqrt(t)
            boys_val = kk *erf(tt) / tt
            do m=0,n-1
                boys_val = ( (2*m+1) * boys_val - et ) / t2
            enddo
        else
        !   For smaller t's compute F with highest n using
        !   asymptotic series (see I. Shavitt in
        !   Methods in Computational Physics, ed. B. Alder eta l,
        !   vol 2, 1963, page 8)
            et = exp(-t)
            t2 = 2.0_dp * t
            m2 = 2 * n
            num = df(m2 + 1)
            sum = 1.0_dp /(m2 + 1)
            do i=1,maxfac-1
                num = num * t2
                term1 = num / df(m2 + 2*i + 3)
                sum = sum + term1
                if(ABS(term1) < eps) exit
            enddo
            boys_val = sum * et 
        endif
    end function
!====================================================================
    subroutine rys_jacobi(t,x,w)
        implicit none
        real(dp),intent(in) :: t
        real(dp),intent(out) :: x(:),w(:)
        !===
        integer :: n
        n = size(x)
        if(t < 5.0_dp) then
            call rysq_small(t,x,w)

        elseif(t < 50.0_dp) then
            call rysq_mid(t,x,w)

        else
            call rysq_asym(t,x,w)
        end if
    end subroutine
!====================================================================
    !QL algorithm
    subroutine tqli(d,e,z)
        implicit none
        real(dp),intent(inout) :: d(:),e(:)
        real(dp),optional,intent(inout) :: z(:)
        !===
        integer :: i,iter,l,m,n,ndum
        real(dp),parameter :: eps = 1.0e-8
        real(dp) :: b,c,dd,f,g,p,r,s
        real(dp) :: ff
        n = size(d)
        if(present(z)) ndum = n
        e(:) = eoshift(e(:),1)
        do l = 1,n
            iter = 0
            iterate : do
                do m = l,n-1
                    dd = abs(d(m)) + abs(d(m+1))
                    if(abs(e(m)) < eps * dd) exit
                end do
                if(m == l) exit iterate
                if(iter == 30) stop 'too many iterations in tqli'
                iter = iter + 1
                g =(d(l+1) - d(l))/(2.0_dp * e(l))
                r = pythag(g,1.0_dp)
                g = d(m) - d(l) + e(l)/(g+sign(r,g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp
                do i = m-1,l,-1
                    f = s*e(i)
                    b = c*e(i)
                    r = pythag(f,g)
                    e(i+1) = r
                    if( r == 0.0_dp ) then
                        d(i+1) = d(i+1) - p
                        e(m) = 0.0_dp
                        cycle iterate
                    end if
                    s = f/r
                    c = g/r
                    g = d(i+1) - p
                    r = (d(i) - g)*s + 2.0_dp*c*b
                    p = s*r
                    d(i+1) = g + p
                    g = c*r - b
                    if(present(z)) then
                        ff = z(i+1)
                        z(i+1) = s*z(i) + c*ff
                        z(i) = c*z(i) - s*ff
                    end if
                end do
                d(l) = d(l) - p
                e(l) = g
                e(m) = 0.0_dp
            end do iterate
        end do
        contains
        function pythag(a,b)
            implicit none
            real(dp) :: a,b
            real(dp) :: pythag
            real(dp) :: absa,absb
            absa = abs(a)
            absb = abs(b)
            if(absa > absb) then
                pythag = absa * sqrt(1.0_dp + (absb/absa)**2)
            else
                if(absb == 0.0_dp) then
                    pythag = 0.0_dp
                else
                    pythag = absb * sqrt(1.0_dp + (absa/absb)**2)
                end if
            end if
         end function
    end subroutine
!====================================================================
    subroutine eigsort_v(d,v)
        implicit none
        real(dp),intent(inout) :: d(:),v(:)
        integer :: i,j,n
        integer :: imin(1)
        n = size(d)
        do i = 1,n-1
            imin = minloc(d(i:n))
            j = imin(1) + i - 1
            if( j /= i) then
                call swap(d(i),d(j))
                call swap(v(i),v(j))
            end if
        end do
    end subroutine
!===============================================
    !Gaussian quadrature from the jacobi matrix
    subroutine gausscof(a,b,amu0,x,w)
        implicit none
        real(dp),intent(inout) :: a(:),b(:)
        real(dp),intent(in) :: amu0
        real(dp),intent(out) :: x(:),w(:)

        real(dp) :: z(size(a))
        integer :: n
        n = size(a)
        b(2:n) = sqrt(b(2:n))
        z(:) = 0.0_dp
        z(1) = 1.0_dp
        call tqli(a,b,z)
        call eigsort_v(a,z)
        x = a
        w = amu0 * z(:)**2
    end subroutine gausscof

!====================================================================
    !wheeler's algorithm
    !Wheeler, J.C. 1974, ¡°Modified Moments and Gaussian Quadratures,¡± Rocky Mountain Journal of Mathematics, vol. 4, pp. 287¨C296
    subroutine orthog(anu,alpha,beta,a,b)
        implicit none
        real(dp),intent(in) :: anu(:),alpha(:),beta(:)
        real(dp),intent(out) :: a(:),b(:)
        integer :: k,n,ndum
        real(dp),dimension(2*size(a)+1,2*size(a)+1) :: sig
        n = size(a)
        ndum = 2*n
        sig(1,3:2*n) = 0.0_dp
        sig(2,2:2*n + 1) = anu(1:2*n)
        a(1) = alpha(1) + anu(2)/anu(1)
        b(1)=0.0_dp
        do k=3,n+1 
            sig(k,k:2*n-k+3)=sig(k-1,k+1:2*n-k+4)+(alpha(k-1:2*n-k+2) &
                -a(k-2))*sig(k-1,k:2*n-k+3)-b(k-2)*sig(k-2,k:2*n-k+3) &
                +beta(k-1:2*n-k+2)*sig(k-1,k-1:2*n-k+2)
            a(k-1)=alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)
            b(k-1)=sig(k,k)/sig(k-1,k-1)
        end do
    end subroutine orthog

!====================================================================
    !
    !Shift Gauss-Legendre quadrature
    !interval [0,1]    weight unity
    subroutine gauss_shift_legendre(x,w) 
        implicit none
        real(dp),intent(inout) :: x(:),w(:)
        integer :: i,n
        real(dp) :: a(size(x)),b(size(x))
        real(dp) :: amu0,tmp
        n = size(x)
        a(:) = 0.5_dp
        do i = 2,n
            tmp = 1.0_dp/(i-1)**2
            b(i) = 1.0_dp/(4.0_dp *(4.0_dp - tmp))
        end do
        amu0 = 1.0_dp
        call gausscof(a,b,amu0,x,w)
    end subroutine gauss_shift_legendre

!====================================================================
    !
    !Shift Gauss-Jacobi quadrature(p=1/2,q=1/2)
    !interval [0,1]    weight 1/sqrt(x)
    subroutine gauss_shift_jacobi(x,w) 
        implicit none
        real(dp),intent(inout) :: x(:),w(:)
        integer :: i,n
        real(dp) :: a(size(x)),b(size(x))
        real(dp) :: amu0,tmp
        n = size(x)
        a(1) = 1.0_dp/3.0_dp
        do i = 2,n
            a(i) = ( 2.0_dp*(i-1.0_dp)*(i-0.5_dp) - 0.25_dp)/((2.0_dp * i-0.5_dp)*(2.0_dp *i - 2.5_dp))
            b(i) = ((i - 1.0_dp)**2 * (i - 1.5_dp)**2 )&
                    /((2.0_dp * i - 2.5_dp)**2 * (2.0_dp * i - 3.5_dp) * (2.0_dp * i - 1.5_dp ))
        end do
        amu0 = 1.0_dp
        call gausscof(a,b,amu0,x,w)
    end subroutine gauss_shift_jacobi

!====================================================================
    !
    !Gauss-LOG quadrature
    subroutine gauss_log(x,w) 
        implicit none
        real(dp),intent(inout) :: x(:),w(:)
        integer :: i,n
        real(dp) :: a(size(x)),b(size(x))
        real(dp) :: anu(2*size(x)),alpha(2*size(x) - 1),beta(2*size(x) - 1)
        real(dp) :: amu0,tmp
        n = size(x)
        alpha(:) = 0.5_dp
        beta(1) = 0.0_dp
        do i = 2,2*n-1
            tmp = 1.0_dp/(i-1)**2
            beta(i) = 0.25_dp/(4.0_dp - tmp)
        end do
        anu(1) = 1.0_dp
        do i = 2,2*n-1
            if(i == 2) then
                anu(2) = -0.25_dp
            else
                tmp = (i - 1.0_dp) * (i - 2.0_dp) / (i * (2.0_dp*i - 3.0_dp))
                anu(i) = -0.5_dp * anu(i-1) * tmp
            endif
        end do
        call orthog(anu,alpha,beta,a,b)
        amu0 = 1.0_dp
        call gausscof(a,b,amu0,x,w)
    end subroutine gauss_log


!====================================================================
    !
    !Integration : DE rule
    subroutine derule(func,a,b,hmax,s,n)
        implicit none
        real(dp),intent(in) :: a,b,hmax
        real(dp),intent(inout) :: s
        integer,intent(in) :: n
        interface
            function func(x,del)
                use m_definitions
                real(dp),intent(in) :: x,del
                real(dp) :: func
            end function
        end interface
        real(dp) :: del,fact,q,sum,t,twoh
        integer :: it,j
        if( n == 1 ) then
            s = hmax * (b - a) * func(0.5_dp*(b+a),0.5_dp*(b-a))/2.0_dp
        else
            it = 2.0_dp**(n-2)
            twoh = hmax/it
            t = 0.5_dp*twoh
            sum = 0.0_dp
            do j = 1,it
                q = exp(-2.0_dp*sinh(t))
                del = (b-a)*q/(1.0_dp + q)
                fact = q/(1.0_dp+q)**2 * cosh(t)
                sum = sum + fact * (func(a+del,del) + func(b-del,del))
                t = t + twoh
            end do
            s = 0.5*s + (b-a) * twoh * sum
        end if
    end subroutine
!====================================================================
    subroutine stiel(x1,x2,hmax,func,x,w)
        implicit none
        real(dp),intent(in) :: x1,x2,hmax
        real(dp),intent(inout) :: x(:),w(:)
        interface 
            function func(x,del)
                use m_definitions
                real(dp),intent(in) :: x,del
                real(dp) :: func
            end function
        end interface
        !
        real(dp) :: amu0,c,oldc
        integer :: i,n
        n = size(x)
        if(allocated(a)) deallocate(a,b)
        allocate(a(n),b(n))
        oldc = 1.0_dp 
        do i = 1,n
            j_index = i
            c = quad(pp)
            b(i) = c/oldc
            a(i) = quad(ppx)/c
            oldc = c
        end do
        amu0 = b(1)
        call gausscof(a,b,amu0,x,w)
        !===
        contains
        !
        function pp(x,del)
            real(dp),intent(in):: x,del
            real(dp) :: pp
            pp = p(x)
            pp = pp * func(x,del) * pp
        end function
        !
        function ppx(x,del)
            real(dp),intent(in):: x,del
            real(dp) :: ppx
            ppx = p(x)
            ppx = ppx * func(x,del) * ppx * x
        end function
        !===
        function quad(fun)
            implicit none
            interface
                function fun(x,del)
                    use m_definitions
                    real(dp),intent(in):: x,del
                    real(dp) :: fun
                end function
            end interface
            real(dp) :: quad
            real(dp),parameter :: eps = 1.0e-9_dp
            integer,parameter :: nmax = 15
            integer :: k
            real(dp) :: olds,sums
            olds = 0.0_dp
            sums = 0.0_dp
            do k = 1,nmax
                call derule(fun,x1,x2,hmax,sums,k)
                if( k > 3) then
                    if(abs(sums - olds) <= eps * abs(olds)) then
                        quad = sums
                        return
                    end if
                end if
                olds = sums
            end do
            stop 'no convergence in quad'
        end function
    end subroutine
!====================================================================
    function p(x)
        real(dp),intent(in) :: x
        real(dp) :: p
        real(dp) :: pval,pj,pjm1
        integer :: i
        if(j_index == 1) then
            p = 1.0_dp
            return
        else
            pjm1 = 0.0_dp
            pj = 1.0_dp
            do i = 1,j_index-1
                pval = (x - a(i))*pj - b(i)*pjm1
                pjm1 = pj
                pj = pval
            end do
            p = pval
            return
        end if
    end function
!====================================================================
    subroutine rys_stiel(n,t,x,w)
       implicit none
       integer,intent(in) :: n
       real(dp),intent(in) :: t
       real(dp),intent(out) :: x(1:n),w(1:n)
       !==
       call stiel(0.0_dp,1.0_dp,4.3_dp,func,x,w)
       contains
       function func(x,del)
            real(dp),intent(in) :: x,del
            real(dp) :: func
            real(dp) :: eps = 1.0d-15
            if(abs(x) < eps) then
                func =  exp(-t*x) /(2.0_dp*sqrt(eps))
            else
                func = exp(-t*x) /(2.0_dp*sqrt(x))
            end if
        end function
    end subroutine
!====================================================================
end module m_rysquad
