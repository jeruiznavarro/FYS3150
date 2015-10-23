program test_laguerre
	implicit none
	character(100)::argument=' '
	real(8)::answer=1.5625d-2
	real(8)::laguerre_integral=0.d0
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
	call get_command_argument(1,argument)
	read(argument,*) number_of_iterations
	allocate(weights(number_of_iterations),zeros(number_of_iterations))
	do i=1,number_of_iterations
		weights(i)=0.d0
		zeros(i)=0.d0
	end do
	call gauss_laguerre(zeros,weights,number_of_iterations,2.d0)
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
							laguerre_integral=laguerre_integral+weights(i)*weights(j)*weights(k)*weights(l)*weights(m)*weights(n)*dcos(x_1)*dcos(x_2)*dcos(x_3)*dcos(y_1)*dcos(y_2)*dcos(y_3)
						end do
					end do
				end do
			end do
		end do
	end do
	write(*,*) 'The value for the laguerre integral is:',laguerre_integral,' which is',laguerre_integral*1.d2/answer,'% of the closed form value.'
	deallocate(weights,zeros)
end program test_laguerre