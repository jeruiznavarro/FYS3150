program project1_regular
implicit none
real(4)::time=0.d0,values(2)=0.d0 !these variables are used to measure how long does it take to carry out all the computations in the program
real(8),allocatable::aux(:),f(:),u(:)!these variables can be dynamically dimensioned
real(8)::emax=-1.d3,eps=0.d0,h=0.d0,temp=2.d0,uclosed=0.d0!these variables are used to compare the numerical solution with analytical one
integer(8)::i=0,j=0,n=0

write(*,*) 'Enter the number of points:'
read(*,*) n ; h=1.d0/dfloat(n+1)
write(*,*) 'The number of points is:',n,', press any key to continue.'!reading the number of grid points
read(*,*)
allocate(aux(n),f(n),u(n)) ; aux=0.d0 ; u=0.d0!allocating and initialazing variables with dynamic memory
do i=1,n
	f(i)=h**2*1.d2*dexp(-1.d1*dfloat(i)/dfloat(n+1))!initializing the source term
end do

u(1)=f(1)/temp ; j=j+1
do i=2,n
	aux(i)=-1.d0/temp ; j=j+1
	temp=2.d0+aux(i) ; j=j+1
	u(i)=(f(i)+u(i-1))/temp ; j=j+2
end do!forward substitution part of the algorithm

!el contador j del numero de floats a lo mejor es asi y no +1 por linea porque en cada linea puede haber mas de un float, ademas, a
!lo mejor hay que contar las operaciones al calcular uclosed, eps y emax 

open(10,file='output_reg.dat')
write(10,*) '#          x          v(x)          u(x)          log(x)          epsilon'
do i=n-1,1,-1
	u(i)=u(i)-u(i+1)*aux(i+1) ; j=j+2!backwards substition part of the algorithm
	uclosed=1.d0-(1.d0-dexp(-1.d1))*dfloat(i)/dfloat(n+1)-dexp(-1.d1*dfloat(i)/dfloat(n+1)) ; eps=dlog10(dabs((u(i)-uclosed)/uclose&
	&d))
	write(10,*) dfloat(i)/dfloat(n+1),u(i),uclosed,dlog10(dfloat(i)/dfloat(n+1)),eps!data exit and recursively finding the maximum error
	emax=dmax1(eps,emax)
end do
close(10)
call etime(values,time)!finding out how long did it take to run the algorithm, below showing miscellaneuos results
write(*,*) 'The number of floats is:',j,'the time elpased during runtime was:',time,'and the maximum relative error is:',emax
deallocate(aux,f,u)!deallocation

end program project1_regular
