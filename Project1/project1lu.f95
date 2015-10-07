program project1_lu
use F90library
use constants
implicit none
real(4)::time=0.d0,values(2)=0.d0 !these variables are used to measure how long does it take to carry out all the computations in the program
real(8),allocatable::f(:),matrix(:,:)!these variables can be dynamically dimensioned
real(8)::d=0.d0,emax=-1.d3,eps=0.d0,h=0.d0,temp=2.d0,uclosed=0.d0!these variables are used to compare the numerical solution with analytical one
integer(4)::i=0,j=0,n=0
integer(4),allocatable::indx(:)

write(*,*) 'Enter the number of points:'
read(*,*) n ; h=1.d0/dfloat(n+1)
write(*,*) 'The number of points is:',n,', press any key to continue.'!reading the number of grid points
read(*,*)
allocate(f(n),indx(n),matrix(n,n)) ; indx=0
do i=1,n
	f(i)=h**2*1.d2*dexp(-1.d1*dfloat(i)/dfloat(n+1))
	if((i>1).and.(i<n))then
		matrix(i,i)=2.d0 ; matrix(i,i-1)=-1.d0 ; matrix(i,i+1)=-1.d0
	end if
end do
matrix(1,1)=2.d0 ; matrix(1,2)=-1.d0 ; matrix(n,n)=2.d0 ; matrix(n,n-1)=-1.d0!allocating and initialazing the matrix and the source term

call lu_decompose(matrix,n,indx,d)
call lu_linear_equation(matrix,n,indx,f)!these two subroutines will solve the equation system using the LU method
call etime(values,time)!finding out how long did it take to run the algorithm, below showing miscellaneuos results
open(10,file='output_lu.dat')
write(10,*) '#          x          v(x)          u(x)          log(x)          epsilon'
do i=n-1,1,-1
	uclosed=1.d0-(1.d0-dexp(-1.d1))*dfloat(i)/dfloat(n+1)-dexp(-1.d1*dfloat(i)/dfloat(n+1)) ; eps=dlog10(dabs((f(i)-uclosed)/uclose&
	&d))
	write(10,*) dfloat(i)/dfloat(n+1),f(i),uclosed,dlog10(dfloat(i)/dfloat(n+1)),eps!data exit and recursively finding the maximum error
	emax=dmax1(eps,emax)
end do
close(10)
write(*,*) 'The time elpased during runtime was:',time,'and the maximum relative error is:',emax
deallocate(f,indx,matrix)!deallocation

end program project1_lu
