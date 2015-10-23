program test_legendre
	use constants
	use F90library
	implicit none
	character(100)::argument=' '
	real(8),parameter::pi=dacos(-1.d0)
	real(8)::answer=6.4d1
	real(8)::legendre_integral=0.d0
	real(8)::lower_integration_limit=-pi*5.d-1
	real(8)::upper_integration_limit=pi*5.d-1
	real(8),allocatable::weights(:)
	real(8)::x_1=0.d0
	real(8)::x_2=0.d0
	real(8)::x_3=0.d0
	real(8)::y_1=0.d0
	real(8)::y_2=0.d0
	real(8)::y_3=0.d0
	real(8),allocatable::zeros(:)
	integer(8)::i=0
	integer(8)::j=0
	integer(8)::k=0
	integer(8)::l=0
	integer(8)::m=0
	integer(8)::n=0
	integer(8)::number_of_iterations=0
	integer(4)::number_of_iterations_4=0
	call get_command_argument(1,argument)
	read(argument,*) number_of_iterations
	number_of_iterations_4=number_of_iterations
	allocate(weights(number_of_iterations),zeros(number_of_iterations))
	do i=1,number_of_iterations
		weights(i)=0.d0
		zeros(i)=0.d0
	end do
	call gauleg(lower_integration_limit,upper_integration_limit,zeros,weights,number_of_iterations_4)
	do i=1,number_of_iterations
		do j=1,number_of_iterations
			do k=1,number_of_iterations
				do l=1,number_of_iterations
					do m=1,number_of_iterations
						do n=1,number_of_iterations
							x_1=zeros(i)
							x_2=zeros(j)
							x_3=zeros(k)
							y_1=zeros(l)
							y_2=zeros(m)
							y_3=zeros(n)
							legendre_integral=legendre_integral+weights(i)*weights(j)*weights(k)*weights(l)*weights(m)*weights(n)*dcos(x_1)*dcos(x_2)*dcos(x_3)*dcos(y_1)*dcos(y_2)*dcos(y_3)
						end do
					end do
				end do
			end do
		end do
	end do
	write(*,*) 'The value for the laguerre integral is:',legendre_integral,' which is',legendre_integral*1.d2/answer,'% of the closed form value.'
	deallocate(weights,zeros)
end program test_legendre


REAL(8) FUNCTION pythag(a,b)
	REAL(8)  :: a,b
	REAL(8)  :: absa,absb
	absa=ABS(a)
	absb=ABS(b)
	IF(absa > absb) THEN
		pythag=absa*sqrt(1.+(absb/absa)**2)
	ELSE
		IF(absb == 0.) THEN
			pythag=0.
		ELSE
			pythag=absb*sqrt(1.+(absa/absb)**2)
		ENDIF
	ENDIF
END FUNCTION pythag