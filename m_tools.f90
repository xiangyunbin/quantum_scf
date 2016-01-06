 !*
 !* Copyright (C) 2016 XiangYunBin <1084066694@qq.com>
 !*
module m_tools
	use m_definitions

	contains
!=========================================================================
	subroutine boys_function(fnt,nc,n,t)
		implicit none
 		integer,intent(in)   :: nc,n
		real(dp),intent(in)  :: t(1:nc)
		real(dp),intent(out) :: fnt(1:nc,0:n)
		!=====
		integer,parameter  :: maxfac=100
		real(dp),parameter :: eps=1.0d-17
		integer :: i,m,k
		integer :: m2
		real(dp) :: t2(nc),num(nc),summ(nc)
		real(dp) :: term1(nc),term2(nc),et(nc),tt(nc)
		logical :: isbig(nc),unfinished(nc)
		real(dp),parameter :: kk = 0.5_dp * SQRT( pi ) 
		real(dp),save :: df(2*maxfac)=0.0_dp
		!=====

		if(n>maxfac) stop' boys function Fm(t) for a too high m value'

		if( ABS(df(1))<1.d-10 ) then
			!   write(stdout,*) 'initialize df'
			df(1:3) = 1.0_dp
			do i=4,2*maxfac
		 		df(i) = (i-2) * df(i-2)
	   		enddo
 		endif

		isbig(:) = ( t(:) > 20.0 )
		unfinished(:) = .not. isbig(:)
		t2(:) = 2.0_dp * t(:)
   		et(:) = exp(-t(:))

		! For big t's do upward recursion 
 		where( isbig )  
   			tt(:)  = sqrt(t(:))
   			fnt(:,0) = kk *erf(tt(:)) / tt(:)
		end where
		do m=0,n-1
			where(isbig )
     			fnt(:,m+1) = ( (2*m+1) * fnt(:,m) - et(:) ) / t2(:)
			end where
   		enddo
		!   For smaller t's compute F with highest n using
		!   asymptotic series (see I. Shavitt in
		!   Methods in Computational Physics, ed. B. Alder eta l,
		!   vol 2, 1963, page 8)
		m2 = 2 * n
		num(:) = df(m2+1)
		summ(:) = 1.0_dp / ( m2 + 1 )
		do i=1,maxfac-1
			where(unfinished)
				num(:) = num(:) * t2(:)
				term1(:) = num(:) / df(m2 + 2*i + 3)
				summ(:) = summ(:) + term1(:)
				unfinished = ( term1(:) > eps )
			end where
			if( .not. any( unfinished ) ) exit
		enddo
		where(.not. isbig)
		    fnt(:,n) = summ(:) * et(:) 
		end where
		!
		! And then do downward recursion 
		do m=n-1,0,-1
			where(.not. isbig)
				fnt(:,m)= ( t2(:) * fnt(:,m+1) + et(:) ) / ( 2 * m + 1 )
			end where
		enddo

		

        end subroutine boys_function

!=========================================================================

    function double_factorial(intin)
         implicit none
         integer,intent(in) :: intin
         real(dp) :: double_factorial
        !=====
         ! just hard coded for some small integers

         select case (intin)
         case(-1) 
           double_factorial = 1.0_dp
         case( 0) 
           double_factorial = 1.0_dp
         case( 1) 
           double_factorial = 1.0_dp
         case( 2) 
           double_factorial = 2.0_dp
         case( 3) 
           double_factorial = 3.0_dp
         case( 4) 
           double_factorial = 8.0_dp
         case( 5) 
           double_factorial = 15.0_dp
         case( 6) 
           double_factorial = 48.0_dp
         case( 7) 
           double_factorial = 105.0_dp
         case( 8) 
           double_factorial = 384.0_dp
         case( 9) 
           double_factorial = 945.0_dp
         case(10) 
           double_factorial = 3840.0_dp
         case(11) 
           double_factorial = 10395.0_dp
         case(12) 
           double_factorial = 46080.0_dp
         case(13) 
           double_factorial = 135135.0_dp
         case(14) 
           double_factorial = 645120.0_dp
         case(15)
           double_factorial = 2027025.0_dp
         case(16) 
           double_factorial = 10321920.0_dp
         case(17) 
           double_factorial = 34459425.0_dp
         case(18) 
           double_factorial = 185794560.0_dp
         case(19) 
           double_factorial = 654729075.0_dp
         case(20) 
           double_factorial = 3715891200.0_dp
         case(21) 
           double_factorial = 13749310575.0_dp
         case(22)
           double_factorial = 81749606400.0_dp
         case(23) 
           double_factorial = 316234143225.0_dp
         case(25) 
           double_factorial = 7905853580625.0_dp
         case(27) 
           double_factorial = 213458046676875.0_dp
         case(29) 
           double_factorial = 6190283353629375.0_dp
         case(31) 
           double_factorial = 191898783962510625.0_dp
         case default
           write(stdout,*) 'integer =',intin
           !stop'double factorial not coded for this integer value'
           write(stdout,*) 'double factorial not coded for this integer value'
           double_factorial = 1
         end select

    end function double_factorial
!=========================================================================
    function factrl(intin)
         implicit none
         integer,intent(in) :: intin
         real(dp) :: factrl
         real(dp),save :: a(0:30)
         logical,save:: init = .true.
         integer :: i
        !==
         if(init) then
	        init = .false.
	        a(0) = 1.0_dp
	        do i = 1,30
		        a(i)=i*a(i-1)
	        end do
	        init = .false.
         end if
         if(intin < 0 .or. intin > 30) stop 'factrl out of range'
         factrl = a(intin)
         return
    end function factrl
!=========================================================================
    function bico(n,k)
         integer,intent(in) :: n,k
         real(dp) :: bico
        !FIXME
         integer,parameter:: para_max = 15
         real(dp),save :: bico_mat(0:para_max,0:para_max)
         integer :: i,j
         logical,save:: init = .true.
        !===
         if(init) then
	         do i=0,para_max
		        do j=0,i
			        bico_mat(i,j) = nint(factrl(i) / (factrl(j) * factrl(i-j) ) )
		        end do
	         end do
  	         init = .false.	
         end if
         !if((n > 15) .or. (n < 0) .or. (k > n) .or. (k < 0)) stop 'bad args in bico'
         bico = bico_mat(n,k)
         return
    end function bico
!=========================================================================
    !None useful
    subroutine cross_product(u1,u2,u3)
         implicit none
         real(dp),intent(in)  :: u1(3),u2(3)
         real(dp),intent(out) :: u3(3)
        !=====

         u3(1) = u1(2) * u2(3) - u1(3) * u2(2)
         u3(2) = u1(3) * u2(1) - u1(1) * u2(3)
         u3(3) = u1(1) * u2(2) - u1(2) * u2(1)
    end subroutine cross_product

!=========================================================================
    function capitalize(str)
         implicit none
         character(*), intent(in) :: str
         character(LEN(str))      :: capitalize
        !=====
         character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
         character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
         integer :: ic, ii
        !=====

         capitalize = str
         do ii=1,LEN_TRIM(str)
           ic = INDEX(low,str(ii:ii))
           if (ic > 0) capitalize(ii:ii) = cap(ic:ic)
         end do
    end function capitalize

!=========================================================================
    subroutine tred2(a,d,e,novectors)
            implicit none
            real(dp),intent(inout)          ::    a(:,:)
            real(dp),intent(out)            ::    d(:),e(:)
            logical,optional,intent(in)   ::    novectors
            !===
            integer                         ::    i,j,l,n
            real(dp)                         ::    f,g,h,hh,scale
            real(dp)                         ::    gg(size(a,1)) 
            logical                         ::    yesvec
            n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
            if (present(novectors)) then
                yesvec=.not. novectors
            else
                yesvec=.true.
            end if
            do i=n,2,-1
		l=i-1
		h=0.0
		if (l > 1) then
			scale=sum(abs(a(i,1:l)))
			if (scale == 0.0) then
				e(i)=a(i,l)
			else
				a(i,1:l)=a(i,1:l)/scale
				h=sum(a(i,1:l)**2)
				f=a(i,l)
				g=-sign(sqrt(h),f)
				e(i)=scale*g
				h=h-f*g
				a(i,l)=f-g
				if (yesvec) a(1:l,i)=a(i,1:l)/h
				do j=1,l
				    e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
					+dot_product(a(j+1:l,j),a(i,j+1:l)))/h
				end do
				f=dot_product(e(1:l),a(i,1:l))
				hh=f/(h+h)
				e(1:l)=e(1:l)-hh*a(i,1:l)
				do j=1,l
					a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
				end do
			end if
		else
			e(i)=a(i,l)
		end if
		d(i)=h
	end do
	if (yesvec) d(1)=0.0
	e(1)=0.0
	do i=1,n
		if (yesvec) then
			l=i-1
			if (d(i) /= 0.0) then
				gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
				a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
			end if
			d(i)=a(i,i)
			a(i,i)=1.0
			a(i,1:l)=0.0
			a(1:l,i)=0.0
		else
			d(i)=a(i,i)
		end if
	end do
    end subroutine tred2
!=========================================================
    subroutine tqli_mat(d,e,z)
            implicit none
            real(dp),intent(inout)              ::    d(:),e(:)
            real(dp),optional,intent(inout)    ::    z(:,:)
            !===
            integer         ::      i,iter,l,m,n,ndum
            real(dp)         ::      b,c,dd,f,g,p,r,s
            real(dp)         ::      ff(size(e))
	    n=assert_eq(size(d),size(e),'tqli: n')
	    if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
	    e(:)=eoshift(e(:),1)
	    do l=1,n
		iter=0
		iterate: do
			do m=l,n-1
				dd=abs(d(m))+abs(d(m+1))
				if (abs(e(m))+dd == dd) exit
			end do
			if (m == l) exit iterate
			if (iter == 30) stop 'too many iterations in tqli'
			iter=iter+1
			g=(d(l+1)-d(l))/(2.0_sp*e(l))
			r=pythag(g,1.0_dp)
			g=d(m)-d(l)+e(l)/(g+sign(r,g))
			s=1.0
			c=1.0
			p=0.0
			do i=m-1,l,-1
				f=s*e(i)
				b=c*e(i)
				r=pythag(f,g)
				e(i+1)=r
				if (r == 0.0) then
					d(i+1)=d(i+1)-p
					e(m)=0.0
					cycle iterate
				end if
				s=f/r
				c=g/r
				g=d(i+1)-p
				r=(d(i)-g)*s+2.0_sp*c*b
				p=s*r
				d(i+1)=g+p
				g=c*r-b
				if (present(z)) then
					ff(1:n)=z(1:n,i+1)
					z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
					z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
				end if
			end do
			d(l)=d(l)-p
			e(l)=g
			e(m)=0.0
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
        end function pythag
   end subroutine tqli_mat
!=========================================================
    subroutine eigsort(d,v)
        implicit none
        real(dp),intent(inout) :: d(:),v(:,:)
        integer :: i,j,n
        n = size(d)
        do i = 1,n-1
            j = iminloc(d(i:n)) + i - 1
            if( j /= i) then
                call swap(d(i),d(j))
                call swap(v(:,i),v(:,j))
            end if
        end do
    end subroutine
!=========================================================

end module m_tools
