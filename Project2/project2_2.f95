program project2_2
implicit none
real(4)::time=0.d0,values(2)=0.d0 !these variables are used to measure how long does it take to carry out all the computations in the program
real(8),allocatable::evecs(:,:),matrix(:,:),lambda(:,:)!these variables can be dynamically dimensioned
real(8)::amax=0.d0,c=0.d0,eik=0.d0,eil=0.d0,esum=0.d0,h=0.d0,mik=0.d0,mil=0.d0,mkk=0.d0,mll=0.d0,omega=0.d0,rhomax=0.d0,s=0.d0,t=0.&
&d0,tau=0.d0,x=1.d9!these variables are used mostly as temporal storage for values in the algorithm
integer(8)::coulomb=1,i=0,j=0,k=0,l=0,m=0,n=0,nn=0

write(*,*) 'Enter the width of the well, the number of points, the maximum value of rho and whether you want coulomb interaction or&
& not (1 and 0 respectively).'
read(*,*) omega,n,rhomax,coulomb!reading initializing parameters
h=rhomax/dfloat(n)!discretization step
allocate(evecs(n-1,n-1),matrix(n-1,n-1),lambda(n-1,2))
evecs=0.d0 ; matrix=0.d0 ; lambda=1.d9
matrix(1,1)=2.d0+omega**2*h**4+dble(coulomb)*h ; matrix(1,2)=-1.d0 ; matrix(n-1,n-1)=2.d0+omega**2*h**4*(dfloat(n-1))**2+dble(coulo&
&mb)*h/dfloat(n-1) ; matrix(n-1,n-2)=-1.d0
evecs(1,1)=1.d0 ; evecs(n-1,n-1)=1.d0
do i=2,n-2
	matrix(i,i)=2.d0+omega**2*h**4*(dfloat(i))**2+dble(coulomb)*h/dfloat(i) ; matrix(i,i-1)=-1.d0 ; matrix(i,i+1)=-1.d0
	evecs(i,i)=1.d0
end do!allocating and initialazing both matrices, the one with the eigenvectors and the one containing the information of the potential

do
	do i=1,n-1
		do j=i+1,n-1
			if(dabs(matrix(i,j))>amax)then
				amax=dabs(matrix(i,j)) ; k=i ; l=j!finding the non-diagonal element with the biggest value
			end if
		end do
	end do
	if(amax**2<=1.d-9)exit!once the desired tolerance is reached the algorithm is finished
	if(matrix(k,l)/=0.d0)then
		amax=0.d0 ; tau=(matrix(l,l)-matrix(k,k))/(2.d0*matrix(k,l))
		if(tau>=0.d0)then
			t=1.d0/(tau+dsqrt(1.d0+tau**2))
		else
			t=-1.d0/(-tau+dsqrt(1.d0+tau**2))
		end if
		c=1.d0/dsqrt(1.d0+t**2) ; s=c*t!computing the values of the sine, cosine and tangent of the angle of the rotation
	else
		c=1.d0 ; s=0.d0
	end if
	mkk=matrix(k,k) ; mll=matrix(l,l)
	matrix(k,k)=c**2*mkk-2.d0*c*s*matrix(k,l)+s**2*mll ; matrix(l,l)=s**2*mkk+2.d0*c*s*matrix(k,l)+c**2*mll
	matrix(k,l)=0.d0 ; matrix(l,k)=0.d0!performing the rotation on isolated elements
	do i=1,n-1
		if((i/=k).and.(i/=l))then
			mik=matrix(i,k) ; mil=matrix(i,l)
			matrix(i,k)=c*mik-s*mil ; matrix(k,i)=matrix(i,k) ; matrix(i,l)=s*mik+c*mil ; matrix(l,i)=matrix(i,l)!performing the rotation 				in the appropiate rows and columns
		end if
		eik=evecs(k,i) ; eil=evecs(l,i)
		evecs(k,i)=c*eik-s*eil ; evecs(l,i)=c*eil+s*eik!the eigenvectors need to be rotated as well
	end do
	m=m+1
end do
call etime(values,time)!finding out how long did it take to run the algorithm

do i=1,n-1
	esum=esum+evecs(1,i)**2!computing the norm of the eigenvectors
	if(matrix(i,i)/h**2<x)then
		x=matrix(i,i)/h**2 ; nn=i!finding the smallest eigenvalue
	end if
end do
lambda(1,1)=x ; lambda(1,2)=dfloat(nn) ; x=1.d9
do i=1,n-1
	evecs(1,i)=evecs(1,i)/dsqrt(esum)!the eigenvectors are normalized
end do
esum=0.d0
do i=2,n-1
	do j=1,n-1
		esum=esum+evecs(i,j)**2!computing the norm of the eigenvectors
		if((matrix(j,j)/h**2<x).and.(matrix(j,j)/h**2>lambda(i-1,1)))then
			x=matrix(j,j)/h**2 ; nn=j!finding the smallest eigenvalue
		end if
	end do
	lambda(i,1)=x ; lambda(i,2)=dfloat(nn) ; x=1.d9
	do j=1,n-1
		evecs(i,j)=evecs(i,j)/dsqrt(esum)!the eigenvectors are normalized
	end do
	esum=0.d0
end do

open(10,file='output2.dat')!results are here
do i=1,n-1
	write(10,*) h*dfloat(i),evecs(idnint(lambda(1,2)),i)**2
	write(10,*) h*dfloat(i),evecs(idnint(lambda(2,2)),i)**2
end do
write(10,*) '#',lambda(1,1),lambda(2,1)
write(10,*) '#	The time elpased during runtime was:',time,'and the number of rotations performed was:',m
write(*,*) 'The time elpased during runtime was:',time,'and the number of rotations performed was:',m
write(*,*) lambda(1,1),lambda(2,1)
call flush(10)
close(10)
deallocate(evecs,matrix,lambda)!deallocation

end program project2_2
