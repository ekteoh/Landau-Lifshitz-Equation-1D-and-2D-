module landau_lifchitz_2d_mod
implicit none
   integer,parameter::DP=selected_real_kind(15)

contains

real(DP) function dot_prod(a,b)
    real(DP),intent(in)::a(3),b(3)
    real(DP)::s
    integer::i
    
    s=0
	    do i=1,size(a)
	    s=s+a(i)*b(i)
	    end do
    dot_prod=s
    end function
    
function cross_prod(a,b)
    real(DP),intent(in)::a(3),b(3)
    real(DP)::cross_prod(3)
    
	   cross_prod(1)=a(2)*b(3)-a(3)*b(2)
	   cross_prod(2)=a(3)*b(1)-a(1)*b(3)
	   cross_prod(3)=a(1)*b(2)-a(2)*b(1)

    end function
    
real(DP) function dc_prod(a,b,c)
    real(DP),intent(in)::a(3),b(3),c(3)
    real(DP)::s(3)
	s=cross_prod(b,c)
    dc_prod=dot_prod(a,s)

    end function
      
    
real(DP) function jacobian(a,b,c)
    real(DP),intent(in)::a(2),b(2),c(2)
	jacobian=(c(1)-b(1))*(a(2)-c(2))-(a(1)-c(1))*(c(2)-b(2))

    end function
    

subroutine getdata_x(input_file,n,x,num_elem,elem)
	character(len=11),intent(in)::input_file
	real(DP),allocatable,intent(out)::x(:,:)
	real(DP)::dummy
	integer,intent(out)::num_elem,n
	integer,allocatable,intent(out)::elem(:,:)
	integer::l,in1=99
	
	open(unit=in1,file=input_file)

	do l=1,8
	read(in1,*)
	end do
	read(in1,*)n
	
	allocate(x(2,n))
	
	do l=1,n
	read(in1,*) dummy,x(1,l),x(2,l)
	end do
	
	do l=1,2
	read(in1,*)
	end do
	
	read(in1,*)num_elem
	allocate(elem(3,num_elem))
	
	do l=1,num_elem
	read(in1,*) dummy,dummy,dummy,dummy,dummy,elem(1,l),elem(2,l),elem(3,l)
	end do
    
    close(unit=in1)

end




end
