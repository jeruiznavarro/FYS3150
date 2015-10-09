!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Running this program will test the algorithm for a simple 3x3 matrix with knonw results, it will check that the eigenvectors are always orthogonal, that the eigenvalues are correct after the iterations are over and that the biggest element of the matrix in any given iteration of the algorithm is always smaller than the biggest eigenvalue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program project2_test
implicit none
logical::success=.true.!if true a message will be printed at the end saying that all test were successful
real(8)evecs(3,3),matrix(3,3),lambda(3,2)!these variables can be dynamically dimensioned
real(8)::amax=0.d0,c=0.d0,eik=0.d0,eil=0.d0,esum=0.d0,mik=0.d0,mil=0.d0,mkk=0.d0,mll=0.d0,s=0.d0,t=0.d0,tau=0.d0,x=1.d9,ortog=0.d0!these variables are used mostly as temporal storage for values in the algorithm
integer(8)::i=0,j=0,k=0,l=0,m=0,nn=0

do j=1,3
	do i=1,3
		matrix(i,j)=i*j!the eigenvalues for this matrix are 0,0 and 14
		if(i==j)then
			evecs(i,j)=1.d0
		else
			evecs(i,j)=0.d0
		end if
	end do
end do!allocating and initialazing both matrices, the one with the eigenvectors and the one containing the information of the potential

do
	do j=1,3
		do i=j+1,3
			if(dabs(matrix(i,j))>amax)then
				amax=dabs(matrix(i,j)) ; k=i ; l=j!finding the non-diagonal element with the biggest value
				if(amax>=1.401d1)then
					write(*,*) 'The biggest element of the matrix is too big, it can not be bigger than the biggest eigenvalue!'
					success=.false.
				end if
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
	do i=1,3
		if((i/=k).and.(i/=l))then
			mik=matrix(i,k) ; mil=matrix(i,l)
			matrix(i,k)=c*mik-s*mil ; matrix(k,i)=matrix(i,k) ; matrix(i,l)=s*mik+c*mil ; matrix(l,i)=matrix(i,l)!performing the rotation 				in the appropiate rows and columns
		end if
		eik=evecs(k,i) ; eil=evecs(l,i)
		evecs(k,i)=c*eik-s*eil ; evecs(l,i)=c*eil+s*eik!the eigenvectors need to be rotated as well
	end do
	ortog=evecs(1,1)*evecs(2,1)+evecs(1,2)*evecs(2,2)+evecs(1,3)*evecs(2,3)
	if(dabs(ortog)>=1.d-4)then
		write(*,*) 'Eigenvectors 1 and 2 are not orthogonal!'
		success=.false.
	end if
	ortog=evecs(1,1)*evecs(3,1)+evecs(1,2)*evecs(3,2)+evecs(1,3)*evecs(3,3)
	if(dabs(ortog)>=1.d-4)then
		write(*,*) 'Eigenvectors 1 and 3 are not orthogonal!'
		success=.false.
	end if
	ortog=evecs(3,1)*evecs(2,1)+evecs(3,2)*evecs(2,2)+evecs(3,3)*evecs(2,3)
	if(dabs(ortog)>=1.d-4)then
		write(*,*) 'Eigenvectors 3 and 2 are not orthogonal!'
		success=.false.
	end if
	m=m+1
end do
do i=1,3
	do j=1,3
		write(*,*) evecs(i,j)
	end do
end do


do i=1,3
	esum=esum+evecs(1,i)**2!computing the norm of the eigenvectors
	if(matrix(i,i)<x)then
		x=matrix(i,i) ; nn=i!finding the smallest eigenvalue
	end if
end do
lambda(1,1)=x ; lambda(1,2)=dfloat(nn) ; x=1.d9
if(dabs(lambda(1,1))>=1.d-2)then
	write(*,*) 'First eigenvalue is wrong! It is:',lambda(1,1)!eigenvalues for this matrix are 14 and 0
	success=.false.
end if
do i=1,3
	evecs(1,i)=evecs(1,i)/dsqrt(esum)!the eigenvectors are normalized
end do
esum=0.d0
do i=2,3
	do j=1,3
		esum=esum+evecs(i,j)**2!computing the norm of the eigenvectors
		if((matrix(j,j)<x).and.(matrix(j,j)>lambda(i-1,1)))then
			x=matrix(j,j) ; nn=j!finding the smallest eigenvalue
		end if
	end do
	lambda(i,1)=x ; lambda(i,2)=dfloat(nn) ; x=1.d9
	if((dabs(lambda(i,1))>=1.d-2).and.(i==2))then
		write(*,*) 'Second eigenvalue is wrong! It is:',lambda(i,1)!eigenvalues for this matrix are 14 and 0
		success=.false.
	end if
	if((dabs(lambda(i,1)-1.4d1)>=1.d-2).and.(i==3))then
		write(*,*) 'Third eigenvalue is wrong! It is:',lambda(i,1)!eigenvalues for this matrix are 14 and 0
		success=.false.
	end if
	do j=1,3
		evecs(i,j)=evecs(i,j)/dsqrt(esum)!the eigenvectors are normalized
	end do
	esum=0.d0
end do
if(success)write(*,*) 'No errors during runtime!'
end program project2_test