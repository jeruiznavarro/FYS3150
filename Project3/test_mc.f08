program test_mc
	implicit none
	character(100)::argument=' '
	character(8)::datechar='--------'
    character(10)::timechar='----------'
	real(8),parameter::answer=2.d0/3.d0
	real(8)::dran_u
	real(8)::lower_integration_limit=-1.d0
	real(8)::mc_integral=0.d0
	real(8)::mc_sample=0.d0
	real(8)::mc_squared_sample_sum=0.d0
	real(8)::standard_deviation=0.d0
	real(8)::time_1=0.d0
	real(8)::time_2=0.d0
	real(8)::upper_integration_limit=1.d0
	real(8)::x=0.d0
	integer(8)::i=0
	integer(8)::number_of_mc_iterations=0
	integer(8)::seed(8)=0
	call date_and_time(date=datechar,time=timechar,values=seed)
    call dran_ini(1000*seed(8)+3*seed(7)*seed(6)/10)
	call get_command_argument(1,argument)
	read(argument,*) number_of_mc_iterations
	call cpu_time(time_1)
	do i=1,number_of_mc_iterations
		x=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		mc_sample=x**2
		mc_integral=mc_integral+mc_sample
		mc_squared_sample_sum=+mc_squared_sample_sum+mc_sample**2
	end do
	call cpu_time(time_2)
	mc_integral=mc_integral/dble(number_of_mc_iterations)
	standard_deviation=dsqrt(mc_squared_sample_sum/dble(number_of_mc_iterations)-mc_integral**2)/dble(number_of_mc_iterations)
	write(*,*) 'The value for the mc integral is:',mc_integral*(upper_integration_limit-lower_integration_limit),' which is',mc_integral*(upper_integration_limit-lower_integration_limit)*1.d2/answer,' % of the closed form value, and with a standard deviation of',standard_deviation,' . The time spent during integration was',time_2-time_1
end program test_mc