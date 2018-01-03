program landau_lifchitz_2d_main
use landau_lifchitz_2d_mod

real(DP),parameter::pi=acos(-1.)
real(DP)::alpha,theta,dt,e,r,dummy,a(3,3),a1l(3,3),a1fixed(3,3),a2l(3,3),a2fixed(3,3),a3l1(3,3),a3fixed1(3,3),a3l2(3,3),a3fixed2(3,3),a3l3(3,3),a3fixed3(3,3),fl(3,3),ffixed(3,3)
real(DP),allocatable::m0(:,:),m(:,:),c(:),f(:),v(:,:),b1(:,:),b2(:,:),x(:,:),mlength(:),e_matrix(:,:)
real(DP),allocatable::a1(:,:),a2(:,:),a3(:,:),a123(:,:)
integer::i,j,n,l,info,q,out1_2d,oute_2d,num_elem,in1=8,i1,i2,i3,t
integer,allocatable::ipiv(:),elem(:,:)

print *,"alpha value ( suggested value: 0.1, 0.5, 0.8, 1 )"
read *,alpha

print *,"theta value (between 0.5 to 1)"
read *,theta

a=0
a(1,2)=1
a(1,3)=1
a(2,1)=-1
a(2,3)=1
a(3,1)=-1
a(3,2)=-2

open(newunit=out1_2d, file='plot_lle_data_2d',action='write')
open(newunit=oute_2d, file='plot_energy_data_2d',action='write')

!-------------------------*get data from circle.msh-------------------------

call getdata_x("circle.msh",n,x,num_elem,elem)

!----------------------------------------------------------------

allocate(m0(3,1:n),m(3,1:n),b1(3,1:n),b2(3,1:n))
allocate(a1(1:2*n,1:2*n),a2(1:2*n,1:2*n),a3(1:2*n,1:2*n),a123(1:2*n,1:2*n))
allocate(f(1:2*n),v(3,1:n),c(1:2*n),e_matrix(1:n,1:n),mlength(1:n),ipiv(2*n))

!----------------------------------------------1.compute m0-----------------------------------------------------------------
 	do q=1,n
 	r=sqrt(x(1,q)**2+x(2,q)**2)
 	
	 	if (r == 0) then
		m0(1,q)=-1.*x(2,q)*pi/2
		m0(2,q)=x(1,q)*pi/2
		m0(3,q)=cos(pi*r/2)
		else
		m0(1,q)=-1.*x(2,q)*sin(pi*r/2)/r
		m0(2,q)=x(1,q)*sin(pi*r/2)/r
		m0(3,q)=cos(pi*r/2)
		end if
		
	end do

m=m0
write(out1_2d,"(f9.5,f9.5,f9.5)") (m(1,i),m(2,i),m(3,i),i=1,n)
write(out1_2d,*)"------------------------------"

dt=1./(n-150)

!loop for how many times we gonna update m
do t=1,n

!----------------------------------------------2.compute b1,b2(check)------------------------------------------------------------------
do i=1,n
z=0
z=m(2,i)*m(3,i)
a(1,1)=z
a(2,2)=z
a(3,3)=z
b1(:,i)=matmul(a,m(:,i))
b2(:,i)=cross_prod(m(:,i),b1(:,i))
end do
!print *,b1
!print *,b2

!-------------------------------compute f----------------------------------------------
f=0

do i=1,num_elem
i1=elem(1,i)
i2=elem(2,i)
i3=elem(3,i)

ffixed(1,1)=(x(2,i3)-x(2,i2))**2+(x(1,i3)-x(1,i2))**2
ffixed(2,2)=(x(2,i1)-x(2,i3))**2+(x(1,i1)-x(1,i3))**2
ffixed(3,3)=(x(2,i1)-x(2,i2))**2+(x(1,i1)-x(1,i2))**2
ffixed(1,2)=(x(2,i2)-x(2,i3))*(x(2,i3)-x(2,i1))+(x(1,i3)-x(1,i2))*(x(1,i1)-x(1,i3))
ffixed(2,1)=ffixed(1,2)
ffixed(1,3)=(x(2,i2)-x(2,i3))*(x(2,i1)-x(2,i2))+(x(1,i3)-x(1,i2))*(x(1,i2)-x(1,i1))
ffixed(3,1)=ffixed(1,3)
ffixed(2,3)=(x(2,i3)-x(2,i1))*(x(2,i1)-x(2,i2))+(x(1,i3)-x(1,i1))*(x(1,i1)-x(1,i2))
ffixed(3,2)=ffixed(2,3)

fl=ffixed/(2*jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3)))
!print *,fl
f(i1)=f(i1)+dot_prod(m(:,i1),b1(:,i1))*fl(1,1)+dot_prod(m(:,i2),b1(:,i1))*fl(1,2)+dot_prod(m(:,i3),b1(:,i1))*fl(1,3)
f(i2)=f(i2)+dot_prod(m(:,i1),b1(:,i2))*fl(1,2)+dot_prod(m(:,i2),b1(:,i2))*fl(2,2)+dot_prod(m(:,i3),b1(:,i2))*fl(2,3)
f(i3)=f(i3)+dot_prod(m(:,i1),b1(:,i3))*fl(1,3)+dot_prod(m(:,i2),b1(:,i3))*fl(3,2)+dot_prod(m(:,i3),b1(:,i3))*fl(3,3)

f(i1+n)=f(i1+n)+dot_prod(m(:,i1),b2(:,i1))*fl(1,1)+dot_prod(m(:,i2),b2(:,i1))*fl(1,2)+dot_prod(m(:,i3),b2(:,i1))*fl(1,3)
f(i2+n)=f(i2+n)+dot_prod(m(:,i1),b2(:,i2))*fl(1,2)+dot_prod(m(:,i2),b2(:,i2))*fl(2,2)+dot_prod(m(:,i3),b2(:,i2))*fl(2,3)
f(i3+n)=f(i3+n)+dot_prod(m(:,i1),b2(:,i3))*fl(1,3)+dot_prod(m(:,i2),b2(:,i3))*fl(3,2)+dot_prod(m(:,i3),b2(:,i3))*fl(3,3)

end do
!print *,f
f=-(1.+(alpha**2))*f

!-------------------------------compute a1------------------------------------------------

a1=0
a1fixed=(1./24)*reshape((/2,1,1,1,2,1,1,1,2/),(/3,3/))
do i=1,num_elem
i1=elem(1,i)
i2=elem(2,i)
i3=elem(3,i)
!print *,jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3))
a1l=jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3))*a1fixed

a1(i1,i1)=a1(i1,i1)+dot_prod(b1(:,i1),b1(:,i1))*a1l(1,1)
a1(i1,i2)=a1(i1,i2)+dot_prod(b1(:,i2),b1(:,i1))*a1l(1,2)
a1(i1,i3)=a1(i1,i3)+dot_prod(b1(:,i3),b1(:,i1))*a1l(1,3)
a1(i2,i1)=a1(i2,i1)+dot_prod(b1(:,i1),b1(:,i2))*a1l(2,1)
a1(i2,i2)=a1(i2,i2)+dot_prod(b1(:,i2),b1(:,i2))*a1l(2,2)
a1(i2,i3)=a1(i2,i3)+dot_prod(b1(:,i3),b1(:,i2))*a1l(2,3)
a1(i3,i1)=a1(i3,i1)+dot_prod(b1(:,i1),b1(:,i3))*a1l(3,1)
a1(i3,i2)=a1(i3,i2)+dot_prod(b1(:,i2),b1(:,i3))*a1l(3,2)
a1(i3,i3)=a1(i3,i3)+dot_prod(b1(:,i3),b1(:,i3))*a1l(3,3)

a1(i1+n,i1)=a1(i1+n,i1)+dot_prod(b1(:,i1),b2(:,i1))*a1l(1,1)
a1(i1+n,i2)=a1(i1+n,i2)+dot_prod(b1(:,i2),b2(:,i1))*a1l(1,2)
a1(i1+n,i3)=a1(i1+n,i3)+dot_prod(b1(:,i3),b2(:,i1))*a1l(1,3)
a1(i2+n,i1)=a1(i2+n,i1)+dot_prod(b1(:,i1),b2(:,i2))*a1l(2,1)
a1(i2+n,i2)=a1(i2+n,i2)+dot_prod(b1(:,i2),b2(:,i2))*a1l(2,2)
a1(i2+n,i3)=a1(i2+n,i3)+dot_prod(b1(:,i3),b2(:,i2))*a1l(2,3)
a1(i3+n,i1)=a1(i3+n,i1)+dot_prod(b1(:,i1),b2(:,i3))*a1l(3,1)
a1(i3+n,i2)=a1(i3+n,i2)+dot_prod(b1(:,i2),b2(:,i3))*a1l(3,2)
a1(i3+n,i3)=a1(i3+n,i3)+dot_prod(b1(:,i3),b2(:,i3))*a1l(3,3)

a1(i1,i1+n)=a1(i1,i1+n)+dot_prod(b2(:,i1),b1(:,i1))*a1l(1,1)
a1(i1,i2+n)=a1(i1,i2+n)+dot_prod(b2(:,i2),b1(:,i1))*a1l(1,2)
a1(i1,i3+n)=a1(i1,i3+n)+dot_prod(b2(:,i3),b1(:,i1))*a1l(1,3)
a1(i2,i1+n)=a1(i2,i1+n)+dot_prod(b2(:,i1),b1(:,i2))*a1l(2,1)
a1(i2,i2+n)=a1(i2,i2+n)+dot_prod(b2(:,i2),b1(:,i2))*a1l(2,2)
a1(i2,i3+n)=a1(i2,i3+n)+dot_prod(b2(:,i3),b1(:,i2))*a1l(2,3)
a1(i3,i1+n)=a1(i3,i1+n)+dot_prod(b2(:,i1),b1(:,i3))*a1l(3,1)
a1(i3,i2+n)=a1(i3,i2+n)+dot_prod(b2(:,i2),b1(:,i3))*a1l(3,2)
a1(i3,i3+n)=a1(i3,i3+n)+dot_prod(b2(:,i3),b1(:,i3))*a1l(3,3)

a1(i1+n,i1+n)=a1(i1+n,i1+n)+dot_prod(b2(:,i1),b2(:,i1))*a1l(1,1)
a1(i1+n,i2+n)=a1(i1+n,i2+n)+dot_prod(b2(:,i2),b2(:,i1))*a1l(1,2)
a1(i1+n,i3+n)=a1(i1+n,i3+n)+dot_prod(b2(:,i3),b2(:,i1))*a1l(1,3)
a1(i2+n,i1+n)=a1(i2+n,i1+n)+dot_prod(b2(:,i1),b2(:,i2))*a1l(2,1)
a1(i2+n,i2+n)=a1(i2+n,i2+n)+dot_prod(b2(:,i2),b2(:,i2))*a1l(2,2)
a1(i2+n,i3+n)=a1(i2+n,i3+n)+dot_prod(b2(:,i3),b2(:,i2))*a1l(2,3)
a1(i3+n,i1+n)=a1(i3+n,i1+n)+dot_prod(b2(:,i1),b2(:,i3))*a1l(3,1)
a1(i3+n,i2+n)=a1(i3+n,i2+n)+dot_prod(b2(:,i2),b2(:,i3))*a1l(3,2)
a1(i3+n,i3+n)=a1(i3+n,i3+n)+dot_prod(b2(:,i3),b2(:,i3))*a1l(3,3)

end do


a1=alpha*a1
!print *,a1

!-------------------------------compute a2----------------------------------------------

a2=0
e_matrix=0

do i=1,num_elem
i1=elem(1,i)
i2=elem(2,i)
i3=elem(3,i)

a2fixed(1,1)=(x(2,i3)-x(2,i2))**2+(x(1,i3)-x(1,i2))**2
a2fixed(2,2)=(x(2,i1)-x(2,i3))**2+(x(1,i1)-x(1,i3))**2
a2fixed(3,3)=(x(2,i1)-x(2,i2))**2+(x(1,i1)-x(1,i2))**2
a2fixed(1,2)=(x(2,i2)-x(2,i3))*(x(2,i3)-x(2,i1))+(x(1,i3)-x(1,i2))*(x(1,i1)-x(1,i3))
a2fixed(2,1)=a2fixed(1,2)
a2fixed(1,3)=(x(2,i2)-x(2,i3))*(x(2,i1)-x(2,i2))+(x(1,i3)-x(1,i2))*(x(1,i2)-x(1,i1))
a2fixed(3,1)=a2fixed(1,3)
a2fixed(2,3)=(x(2,i3)-x(2,i1))*(x(2,i1)-x(2,i2))+(x(1,i3)-x(1,i1))*(x(1,i1)-x(1,i2))
a2fixed(3,2)=a2fixed(2,3)

a2l=a2fixed/(2*jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3)))

e_matrix(i1,i1)=e_matrix(i1,i1)+dot_prod(m(:,i1),m(:,i1))*a2l(1,1)
e_matrix(i1,i2)=e_matrix(i1,i2)+dot_prod(m(:,i2),m(:,i1))*a2l(1,2)
e_matrix(i1,i3)=e_matrix(i1,i3)+dot_prod(m(:,i3),m(:,i1))*a2l(1,3)
e_matrix(i2,i1)=e_matrix(i2,i1)+dot_prod(m(:,i1),m(:,i2))*a2l(2,1)
e_matrix(i2,i2)=e_matrix(i2,i2)+dot_prod(m(:,i2),m(:,i2))*a2l(2,2)
e_matrix(i2,i3)=e_matrix(i2,i3)+dot_prod(m(:,i3),m(:,i2))*a2l(2,3)
e_matrix(i3,i1)=e_matrix(i3,i1)+dot_prod(m(:,i1),m(:,i3))*a2l(3,1)
e_matrix(i3,i2)=e_matrix(i3,i2)+dot_prod(m(:,i2),m(:,i3))*a2l(3,2)
e_matrix(i3,i3)=e_matrix(i3,i3)+dot_prod(m(:,i3),m(:,i3))*a2l(3,3)

a2(i1,i1)=a2(i1,i1)+dot_prod(b1(:,i1),b1(:,i1))*a2l(1,1)
a2(i1,i2)=a2(i1,i2)+dot_prod(b1(:,i2),b1(:,i1))*a2l(1,2)
a2(i1,i3)=a2(i1,i3)+dot_prod(b1(:,i3),b1(:,i1))*a2l(1,3)
a2(i2,i1)=a2(i2,i1)+dot_prod(b1(:,i1),b1(:,i2))*a2l(2,1)
a2(i2,i2)=a2(i2,i2)+dot_prod(b1(:,i2),b1(:,i2))*a2l(2,2)
a2(i2,i3)=a2(i2,i3)+dot_prod(b1(:,i3),b1(:,i2))*a2l(2,3)
a2(i3,i1)=a2(i3,i1)+dot_prod(b1(:,i1),b1(:,i3))*a2l(3,1)
a2(i3,i2)=a2(i3,i2)+dot_prod(b1(:,i2),b1(:,i3))*a2l(3,2)
a2(i3,i3)=a2(i3,i3)+dot_prod(b1(:,i3),b1(:,i3))*a2l(3,3)

a2(i1,i1+n)=a2(i1,i1+n)+dot_prod(b2(:,i1),b1(:,i1))*a2l(1,1)
a2(i1,i2+n)=a2(i1,i2+n)+dot_prod(b2(:,i2),b1(:,i1))*a2l(1,2)
a2(i1,i3+n)=a2(i1,i3+n)+dot_prod(b2(:,i3),b1(:,i1))*a2l(1,3)
a2(i2,i1+n)=a2(i2,i1+n)+dot_prod(b2(:,i1),b1(:,i2))*a2l(2,1)
a2(i2,i2+n)=a2(i2,i2+n)+dot_prod(b2(:,i2),b1(:,i2))*a2l(2,2)
a2(i2,i3+n)=a2(i2,i3+n)+dot_prod(b2(:,i3),b1(:,i2))*a2l(2,3)
a2(i3,i1+n)=a2(i3,i1+n)+dot_prod(b2(:,i1),b1(:,i3))*a2l(3,1)
a2(i3,i2+n)=a2(i3,i2+n)+dot_prod(b2(:,i2),b1(:,i3))*a2l(3,2)
a2(i3,i3+n)=a2(i3,i3+n)+dot_prod(b2(:,i3),b1(:,i3))*a2l(3,3)

a2(i1+n,i1)=a2(i1+n,i1)+dot_prod(b1(:,i1),b2(:,i1))*a2l(1,1)
a2(i1+n,i2)=a2(i1+n,i2)+dot_prod(b1(:,i2),b2(:,i1))*a2l(1,2)
a2(i1+n,i3)=a2(i1+n,i3)+dot_prod(b1(:,i3),b2(:,i1))*a2l(1,3)
a2(i2+n,i1)=a2(i2+n,i1)+dot_prod(b1(:,i1),b2(:,i2))*a2l(2,1)
a2(i2+n,i2)=a2(i2+n,i2)+dot_prod(b1(:,i2),b2(:,i2))*a2l(2,2)
a2(i2+n,i3)=a2(i2+n,i3)+dot_prod(b1(:,i3),b2(:,i2))*a2l(2,3)
a2(i3+n,i1)=a2(i3+n,i1)+dot_prod(b1(:,i1),b2(:,i3))*a2l(3,1)
a2(i3+n,i2)=a2(i3+n,i2)+dot_prod(b1(:,i2),b2(:,i3))*a2l(3,2)
a2(i3+n,i3)=a2(i3+n,i3)+dot_prod(b1(:,i3),b2(:,i3))*a2l(3,3)

a2(i1+n,i1+n)=a2(i1+n,i1+n)+dot_prod(b2(:,i1),b2(:,i1))*a2l(1,1)
a2(i1+n,i2+n)=a2(i1+n,i2+n)+dot_prod(b2(:,i2),b2(:,i1))*a2l(1,2)
a2(i1+n,i3+n)=a2(i1+n,i3+n)+dot_prod(b2(:,i3),b2(:,i1))*a2l(1,3)
a2(i2+n,i1+n)=a2(i2+n,i1+n)+dot_prod(b2(:,i1),b2(:,i2))*a2l(2,1)
a2(i2+n,i2+n)=a2(i2+n,i2+n)+dot_prod(b2(:,i2),b2(:,i2))*a2l(2,2)
a2(i2+n,i3+n)=a2(i2+n,i3+n)+dot_prod(b2(:,i3),b2(:,i2))*a2l(2,3)
a2(i3+n,i1+n)=a2(i3+n,i1+n)+dot_prod(b2(:,i1),b2(:,i3))*a2l(3,1)
a2(i3+n,i2+n)=a2(i3+n,i2+n)+dot_prod(b2(:,i2),b2(:,i3))*a2l(3,2)
a2(i3+n,i3+n)=a2(i3+n,i3+n)+dot_prod(b2(:,i3),b2(:,i3))*a2l(3,3)
end do

a2=theta*dt*(1.+(alpha**2))*a2
!print *,a2
!-------------------------------compute a3----------------------------------------------

a3=0
a3fixed1=(1./600)*reshape((/30,10,10,10,10,5,10,5,10/),(/3,3/))
a3fixed2=(1./600)*reshape((/10,10,5,10,30,10,5,10,10/),(/3,3/))
a3fixed3=(1./600)*reshape((/10,5,10,5,10,10,10,10,30/),(/3,3/))

do i=1,num_elem
i1=elem(1,i)
i2=elem(2,i)
i3=elem(3,i)

a3l1=jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3))*a3fixed1
a3l2=jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3))*a3fixed2
a3l3=jacobian(x(1:2,i1),x(1:2,i2),x(1:2,i3))*a3fixed3

a3(i1,i1)=a3(i1,i1)+dc_prod(b1(:,i1),b1(:,i1),m(:,i1))*a3l1(1,1)+dc_prod(b1(:,i1),b1(:,i1),m(:,i2))*a3l2(1,1)+dc_prod(b1(:,i1),b1(:,i1),m(:,i3))*a3l3(1,1)
a3(i1,i2)=a3(i1,i2)+dc_prod(b1(:,i2),b1(:,i1),m(:,i1))*a3l1(1,2)+dc_prod(b1(:,i2),b1(:,i1),m(:,i2))*a3l2(1,2)+dc_prod(b1(:,i2),b1(:,i1),m(:,i3))*a3l3(1,2)
a3(i1,i3)=a3(i1,i3)+dc_prod(b1(:,i3),b1(:,i1),m(:,i1))*a3l1(1,3)+dc_prod(b1(:,i3),b1(:,i1),m(:,i2))*a3l2(1,3)+dc_prod(b1(:,i3),b1(:,i1),m(:,i3))*a3l3(1,3)
a3(i2,i1)=a3(i2,i1)+dc_prod(b1(:,i1),b1(:,i2),m(:,i1))*a3l1(2,1)+dc_prod(b1(:,i1),b1(:,i2),m(:,i2))*a3l2(2,1)+dc_prod(b1(:,i1),b1(:,i2),m(:,i3))*a3l3(2,1)
a3(i2,i2)=a3(i2,i2)+dc_prod(b1(:,i2),b1(:,i2),m(:,i1))*a3l1(2,2)+dc_prod(b1(:,i2),b1(:,i2),m(:,i2))*a3l2(2,2)+dc_prod(b1(:,i2),b1(:,i2),m(:,i3))*a3l3(2,2)
a3(i2,i3)=a3(i2,i3)+dc_prod(b1(:,i3),b1(:,i2),m(:,i1))*a3l1(2,3)+dc_prod(b1(:,i3),b1(:,i2),m(:,i2))*a3l2(2,3)+dc_prod(b1(:,i3),b1(:,i2),m(:,i3))*a3l3(2,3)
a3(i3,i1)=a3(i3,i1)+dc_prod(b1(:,i1),b1(:,i3),m(:,i1))*a3l1(3,1)+dc_prod(b1(:,i1),b1(:,i3),m(:,i2))*a3l2(3,1)+dc_prod(b1(:,i1),b1(:,i3),m(:,i3))*a3l3(3,1)
a3(i3,i2)=a3(i3,i2)+dc_prod(b1(:,i2),b1(:,i3),m(:,i1))*a3l1(3,2)+dc_prod(b1(:,i2),b1(:,i3),m(:,i2))*a3l2(3,2)+dc_prod(b1(:,i2),b1(:,i3),m(:,i3))*a3l3(3,2)
a3(i3,i3)=a3(i3,i3)+dc_prod(b1(:,i3),b1(:,i3),m(:,i1))*a3l1(3,3)+dc_prod(b1(:,i3),b1(:,i3),m(:,i2))*a3l2(3,3)+dc_prod(b1(:,i3),b1(:,i3),m(:,i3))*a3l3(3,3)

!-------------------------

a3(i1,i1+n)=a3(i1,i1+n)+dc_prod(b2(:,i1),b1(:,i1),m(:,i1))*a3l1(1,1)+dc_prod(b2(:,i1),b1(:,i1),m(:,i2))*a3l2(1,1)+dc_prod(b2(:,i1),b1(:,i1),m(:,i3))*a3l3(1,1)
a3(i1,i2+n)=a3(i1,i2+n)+dc_prod(b2(:,i2),b1(:,i1),m(:,i1))*a3l1(1,2)+dc_prod(b2(:,i2),b1(:,i1),m(:,i2))*a3l2(1,2)+dc_prod(b2(:,i2),b1(:,i1),m(:,i3))*a3l3(1,2)
a3(i1,i3+n)=a3(i1,i3+n)+dc_prod(b2(:,i3),b1(:,i1),m(:,i1))*a3l1(1,3)+dc_prod(b2(:,i3),b1(:,i1),m(:,i2))*a3l2(1,3)+dc_prod(b2(:,i3),b1(:,i1),m(:,i3))*a3l3(1,3)
a3(i2,i1+n)=a3(i2,i1+n)+dc_prod(b2(:,i1),b1(:,i2),m(:,i1))*a3l1(2,1)+dc_prod(b2(:,i1),b1(:,i2),m(:,i2))*a3l2(2,1)+dc_prod(b2(:,i1),b1(:,i2),m(:,i3))*a3l3(2,1)
a3(i2,i2+n)=a3(i2,i2+n)+dc_prod(b2(:,i2),b1(:,i2),m(:,i1))*a3l1(2,2)+dc_prod(b2(:,i2),b1(:,i2),m(:,i2))*a3l2(2,2)+dc_prod(b2(:,i2),b1(:,i2),m(:,i3))*a3l3(2,2)
a3(i2,i3+n)=a3(i2,i3+n)+dc_prod(b2(:,i3),b1(:,i2),m(:,i1))*a3l1(2,3)+dc_prod(b2(:,i3),b1(:,i2),m(:,i2))*a3l2(2,3)+dc_prod(b2(:,i3),b1(:,i2),m(:,i3))*a3l3(2,3)
a3(i3,i1+n)=a3(i3,i1+n)+dc_prod(b2(:,i1),b1(:,i3),m(:,i1))*a3l1(3,1)+dc_prod(b2(:,i1),b1(:,i3),m(:,i2))*a3l2(3,1)+dc_prod(b2(:,i1),b1(:,i3),m(:,i3))*a3l3(3,1)
a3(i3,i2+n)=a3(i3,i2+n)+dc_prod(b2(:,i2),b1(:,i3),m(:,i1))*a3l1(3,2)+dc_prod(b2(:,i2),b1(:,i3),m(:,i2))*a3l2(3,2)+dc_prod(b2(:,i2),b1(:,i3),m(:,i3))*a3l3(3,2)
a3(i3,i3+n)=a3(i3,i3+n)+dc_prod(b2(:,i3),b1(:,i3),m(:,i1))*a3l1(3,3)+dc_prod(b2(:,i3),b1(:,i3),m(:,i2))*a3l2(3,3)+dc_prod(b2(:,i3),b1(:,i3),m(:,i3))*a3l3(3,3)

!---------------------------

a3(i1+n,i1)=a3(i1+n,i1)+dc_prod(b1(:,i1),b2(:,i1),m(:,i1))*a3l1(1,1)+dc_prod(b1(:,i1),b2(:,i1),m(:,i2))*a3l2(1,1)+dc_prod(b1(:,i1),b2(:,i1),m(:,i3))*a3l3(1,1)
a3(i1+n,i2)=a3(i1+n,i2)+dc_prod(b1(:,i2),b2(:,i1),m(:,i1))*a3l1(1,2)+dc_prod(b1(:,i2),b2(:,i1),m(:,i2))*a3l2(1,2)+dc_prod(b1(:,i2),b2(:,i1),m(:,i3))*a3l3(1,2)
a3(i1+n,i3)=a3(i1+n,i3)+dc_prod(b1(:,i3),b2(:,i1),m(:,i1))*a3l1(1,3)+dc_prod(b1(:,i3),b2(:,i1),m(:,i2))*a3l2(1,3)+dc_prod(b1(:,i3),b2(:,i1),m(:,i3))*a3l3(1,3)
a3(i2+n,i1)=a3(i2+n,i1)+dc_prod(b1(:,i1),b2(:,i2),m(:,i1))*a3l1(2,1)+dc_prod(b1(:,i1),b2(:,i2),m(:,i2))*a3l2(2,1)+dc_prod(b1(:,i1),b2(:,i2),m(:,i3))*a3l3(2,1)
a3(i2+n,i2)=a3(i2+n,i2)+dc_prod(b1(:,i2),b2(:,i2),m(:,i1))*a3l1(2,2)+dc_prod(b1(:,i2),b2(:,i2),m(:,i2))*a3l2(2,2)+dc_prod(b1(:,i2),b2(:,i2),m(:,i3))*a3l3(2,2)
a3(i2+n,i3)=a3(i2+n,i3)+dc_prod(b1(:,i3),b2(:,i2),m(:,i1))*a3l1(2,3)+dc_prod(b1(:,i3),b2(:,i2),m(:,i2))*a3l2(2,3)+dc_prod(b1(:,i3),b2(:,i2),m(:,i3))*a3l3(2,3)
a3(i3+n,i1)=a3(i3+n,i1)+dc_prod(b1(:,i1),b2(:,i3),m(:,i1))*a3l1(3,1)+dc_prod(b1(:,i1),b2(:,i3),m(:,i2))*a3l2(3,1)+dc_prod(b1(:,i1),b2(:,i3),m(:,i3))*a3l3(3,1)
a3(i3+n,i2)=a3(i3+n,i2)+dc_prod(b1(:,i2),b2(:,i3),m(:,i1))*a3l1(3,2)+dc_prod(b1(:,i2),b2(:,i3),m(:,i2))*a3l2(3,2)+dc_prod(b1(:,i2),b2(:,i3),m(:,i3))*a3l3(3,2)
a3(i3+n,i3)=a3(i3+n,i3)+dc_prod(b1(:,i3),b2(:,i3),m(:,i1))*a3l1(3,3)+dc_prod(b1(:,i3),b2(:,i3),m(:,i2))*a3l2(3,3)+dc_prod(b1(:,i3),b2(:,i3),m(:,i3))*a3l3(3,3)

!--------------------

a3(i1+n,i1+n)=a3(i1+n,i1+n)+dc_prod(b2(:,i1),b2(:,i1),m(:,i1))*a3l1(1,1)+dc_prod(b2(:,i1),b2(:,i1),m(:,i2))*a3l2(1,1)+dc_prod(b2(:,i1),b2(:,i1),m(:,i3))*a3l3(1,1)
a3(i1+n,i2+n)=a3(i1+n,i2+n)+dc_prod(b2(:,i2),b2(:,i1),m(:,i1))*a3l1(1,2)+dc_prod(b2(:,i2),b2(:,i1),m(:,i2))*a3l2(1,2)+dc_prod(b2(:,i2),b2(:,i1),m(:,i3))*a3l3(1,2)
a3(i1+n,i3+n)=a3(i1+n,i3+n)+dc_prod(b2(:,i3),b2(:,i1),m(:,i1))*a3l1(1,3)+dc_prod(b2(:,i3),b2(:,i1),m(:,i2))*a3l2(1,3)+dc_prod(b2(:,i3),b2(:,i1),m(:,i3))*a3l3(1,3)
a3(i2+n,i1+n)=a3(i2+n,i1+n)+dc_prod(b2(:,i1),b2(:,i2),m(:,i1))*a3l1(2,1)+dc_prod(b2(:,i1),b2(:,i2),m(:,i2))*a3l2(2,1)+dc_prod(b2(:,i1),b2(:,i2),m(:,i3))*a3l3(2,1)
a3(i2+n,i2+n)=a3(i2+n,i2+n)+dc_prod(b2(:,i2),b2(:,i2),m(:,i1))*a3l1(2,2)+dc_prod(b2(:,i2),b2(:,i2),m(:,i2))*a3l2(2,2)+dc_prod(b2(:,i2),b2(:,i2),m(:,i3))*a3l3(2,2)
a3(i2+n,i3+n)=a3(i2+n,i3+n)+dc_prod(b2(:,i3),b2(:,i2),m(:,i1))*a3l1(2,3)+dc_prod(b2(:,i3),b2(:,i2),m(:,i2))*a3l2(2,3)+dc_prod(b2(:,i3),b2(:,i2),m(:,i3))*a3l3(2,3)
a3(i3+n,i1+n)=a3(i3+n,i1+n)+dc_prod(b2(:,i1),b2(:,i3),m(:,i1))*a3l1(3,1)+dc_prod(b2(:,i1),b2(:,i3),m(:,i2))*a3l2(3,1)+dc_prod(b2(:,i1),b2(:,i3),m(:,i3))*a3l3(3,1)
a3(i3+n,i2+n)=a3(i3+n,i2+n)+dc_prod(b2(:,i2),b2(:,i3),m(:,i1))*a3l1(3,2)+dc_prod(b2(:,i2),b2(:,i3),m(:,i2))*a3l2(3,2)+dc_prod(b2(:,i2),b2(:,i3),m(:,i3))*a3l3(3,2)
a3(i3+n,i3+n)=a3(i3+n,i3+n)+dc_prod(b2(:,i3),b2(:,i3),m(:,i1))*a3l1(3,3)+dc_prod(b2(:,i3),b2(:,i3),m(:,i2))*a3l2(3,3)+dc_prod(b2(:,i3),b2(:,i3),m(:,i3))*a3l3(3,3)

end do
a3=-1.*a3
!print *,a3

!-------------------------compute energy-----------------------------------

!e=(1./2)*sum(sum(e_matrix,1),1)
e=sum(sum(e_matrix,1),1)
write(oute_2d,*)e


!----------------------------------------------5.compute a123 and solve a123*c=f----------------------------------------------------------------------------

c=0
a123=0
a123=a1+a2+a3
c=f
call dgesv(2*n,1,a123,2*n,ipiv,c,2*n,info)

!print *,"info:"
!print *,info


do i=1,n
!----------------------------------------------compute v-----------------------------------------------------------------------------------------------
v(:,i)=c(i)*b1(:,i)+c(n+i)*b2(:,i)
!----------------------------------------------compute mlength( length of m )-------------------------------------------------------------------
mlength(i)=sqrt(dot_prod(m(:,i)+dt*v(:,i),m(:,i)+dt*v(:,i)))
!----------------------------------------------update m----------------------------------------------------------------------------------------
m(:,i)=(m(:,i)+dt*v(:,i))/mlength(i)
end do


write(out1_2d,"(f9.5,f9.5,f9.5)") (m(1,i),m(2,i),m(3,i),i=1,n)
write(out1_2d,*)"------------------------------"




end do

deallocate(m0,m,b1,b2,x,elem)
deallocate(a1,a2,a3,a123)
deallocate(f,v,c,mlength,ipiv)

close(unit=out1_2d)
close(unit=oute_2d)

end
