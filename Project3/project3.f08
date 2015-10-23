program project3
	use constants
	use F90library
	implicit none
	logical::file_exists=.true.
	character(100)::argument=' '
	character(8)::datechar='--------'
    character(10)::timechar='----------'
	real(8),parameter::pi=dacos(-1.d0)
	real(8),parameter::answer=5.d0*pi**2/2.56d2
	real(8)::calculate_module
	real(8)::dran_u
	real(8)::legendre_integral=0.d0
	real(8)::laguerre_integral=0.d0
	real(8)::lower_integration_limit=-2.d0
	real(8)::mc_integral=0.d0
	real(8)::mc_sample=0.d0
	real(8)::mc_squared_sample_sum=0.d0
	real(8)::mcis_integral=0.d0
	real(8)::phi_1=0.d0
	real(8)::phi_2=0.d0
	real(8)::proximity_threshold=0.d0
	real(8)::r_1=0.d0
	real(8)::r_2=0.d0
	real(8)::r_12=0.d0
	real(8)::random_try=0.d0
	real(8)::standard_deviation=0.d0
	real(8)::standard_deviation_is=0.d0
	real(8)::theta_1=0.d0
	real(8)::theta_2=0.d0
	real(8)::time_1=0.d0
	real(8)::time_2=0.d0
	real(8)::upper_integration_limit=2.d0
	real(8)::wavefunction
	real(8),allocatable::weights(:)
	real(8),allocatable::weights_for_phi(:)
	real(8),allocatable::weights_for_theta(:)
	real(8)::x_1=0.d0
	real(8)::x_2=0.d0
	real(8)::x_3=0.d0
	real(8)::y_1=0.d0
	real(8)::y_2=0.d0
	real(8)::y_3=0.d0
	real(8),allocatable::zeros(:)
	real(8),allocatable::zeros_for_phi(:)
	real(8),allocatable::zeros_for_theta(:)
	integer(8)::i=0
	integer(8)::j=0
	integer(8)::k=0
	integer(8)::l=0
	integer(8)::m=0
	integer(8)::n=0
	integer(8)::number_of_iterations=0
	integer(4)::number_of_iterations_4=0
	integer(8)::number_of_mc_iterations=0
	integer(8)::seed(8)=0
	wavefunction(r_1,r_2,r_12)=dexp(-4.d0*(r_1+r_2))/r_12!the hydrogenic inspired wavefunction is here
	calculate_module(x_1,x_2,x_3)=dsqrt(x_1**2+x_2**2+x_3**2)
	inquire(file='output.dat',exist=file_exists)
	if(file_exists)then
		open(10,file='output.dat',status='old',position='append',action='write')
	else
		open(10,file='output.dat',status='new',action='write')
	end if
	call date_and_time(date=datechar,time=timechar,values=seed)
    call dran_ini(1000*seed(8)+3*seed(7)*seed(6)/10)!starting the random number generator with a seed taken from the computer clock
	call get_command_argument(1,argument)
	read(argument,*) number_of_iterations
	number_of_iterations_4=number_of_iterations
	call get_command_argument(2,argument)
	read(argument,*) number_of_mc_iterations
	call get_command_argument(3,argument)
	read(argument,*) proximity_threshold!reading the three arguments for the program
	write(10,*) 'Run from ',datechar,' at ',timechar,' with',number_of_iterations,' zeros,',number_of_mc_iterations,' mc cycles and',proximity_threshold,'proximity threshold.'
	allocate(weights(number_of_iterations),zeros(number_of_iterations),weights_for_phi(number_of_iterations),zeros_for_phi(number_of_iterations),weights_for_theta(number_of_iterations),zeros_for_theta(number_of_iterations))!allocating arrays
	do i=1,number_of_iterations
		weights(i)=0.d0
		zeros(i)=0.d0
		weights_for_phi(i)=0.d0
		zeros_for_phi(i)=0.d0
		weights_for_theta(i)=0.d0
		zeros_for_theta(i)=0.d0
	end do

	!#####	GAUSS-LEGENDRE INTEGRATION STARTS HERE	#####

	call cpu_time(time_1)
	call gauleg(lower_integration_limit,upper_integration_limit,zeros,weights,number_of_iterations_4)!getting the legendre weights and integration points (zeros of the polynomials)
	do i=1,number_of_iterations
		do j=1,number_of_iterations
			do k=1,number_of_iterations
				do l=1,number_of_iterations
					do m=1,number_of_iterations
						do n=1,number_of_iterations
							r_1=calculate_module(zeros(i),zeros(j),zeros(k))
							r_2=calculate_module(zeros(l),zeros(m),zeros(n))
							r_12=calculate_module(zeros(i)-zeros(l),zeros(j)-zeros(m),zeros(k)-zeros(n))
							if(r_12<=proximity_threshold)cycle!if the electrons are too close the integral will blow up, so we skip dangerous integration points
							legendre_integral=legendre_integral+weights(i)*weights(j)*weights(k)*weights(l)*weights(m)*weights(n)*wavefunction(r_1,r_2,r_12)!integrating
						end do
					end do
				end do
			end do
		end do
	end do
	call cpu_time(time_2)
	write(*,*) 'The value for the legendre integral is:',legendre_integral,' which is',legendre_integral*1.d2/answer,' % of the closed form value. The time spent during integration was',time_2-time_1,' seconds.'
	write(10,*) 'The value for the legendre integral is:',legendre_integral,' which is',legendre_integral*1.d2/answer,' % of the closed form value. The time spent during integration was',time_2-time_1,' seconds.'
	do i=1,number_of_iterations
		weights(i)=0.d0
		zeros(i)=0.d0
	end do

	!#####	GAUSS-LAGUERRE INTEGRATION STARTS HERE	#####

	call cpu_time(time_1)
	call gauss_laguerre(zeros,weights,number_of_iterations,2.d0)!getting the laguerre weights and integration points (zeros of the polynomials) for the radial variable
	call gauleg(0.d0,2.d0*pi,zeros_for_phi,weights_for_phi,number_of_iterations_4)!getting the legendre weights and integration points (zeros of the polynomials) for the azimuthal variable
	call gauleg(0.d0,pi,zeros_for_theta,weights_for_theta,number_of_iterations_4)!getting the legendre weights and integration points (zeros of the polynomials) for the polar variable
	do i=1,number_of_iterations
		do j=1,number_of_iterations
			do k=1,number_of_iterations
				do l=1,number_of_iterations
					do m=1,number_of_iterations
						do n=1,number_of_iterations
							r_1=4.d0*zeros(i)
							r_2=4.d0*zeros(l)!the "4" factors are here to account for the change of variable when using the laguerre polynomials (exp(-x) instead of our function exp(-4x))
							phi_1=zeros_for_phi(j)
							phi_2=zeros_for_phi(m)
							theta_1=zeros_for_theta(k)
							theta_2=zeros_for_theta(n)if((r_1==r_2).and.((phi_1==phi_2).and.(theta_1==theta_2)))cycle!sometimes if the value for r_12 should be 0 there can be NaN values if the radicand is practically 0 but negative
							r_12=dsqrt(r_1**2+r_2**2-2.d0*r_1*r_2*(dcos(theta_1)*dcos(theta_2)+dsin(theta_1)*dsin(theta_2)*dcos(phi_1-phi_2)))
							if((r_12<=proximity_threshold).or.(isnan(r_12)))cycle!if the electrons are too close the integral will blow up, so we skip dangerous integration points
							laguerre_integral=laguerre_integral+weights(i)*weights_for_phi(j)*weights_for_theta(k)*weights(l)*weights_for_phi(m)*weights_for_theta(n)*dsin(theta_1)*dsin(theta_2)/(r_12*2.56d2)!integrating, the exponential term of the wavefunction and the r_1**2 and r_2**2 factors are already absorved by the laguerre weights
						end do
					end do
				end do
			end do
		end do
	end do
	call cpu_time(time_2)
	write(*,*) 'The value for the laguerre integral is:',laguerre_integral,' which is',laguerre_integral*1.d2/answer,' % of the closed form value. The time spent during integration was',time_2-time_1,' seconds.'
	write(10,*) 'The value for the laguerre integral is:',laguerre_integral,' which is',laguerre_integral*1.d2/answer,' % of the closed form value. The time spent during integration was',time_2-time_1,' seconds.'
	deallocate(weights,zeros,weights_for_phi,zeros_for_phi,weights_for_theta,zeros_for_theta)

	!#####	MONTE CARLO INTEGRATION STARTS HERE	#####

	call cpu_time(time_1)
	do i=1,number_of_mc_iterations
		x_1=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		x_2=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		x_3=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		y_1=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		y_2=2.d0*upper_integration_limit*dran_u()+lower_integration_limit
		y_3=2.d0*upper_integration_limit*dran_u()+lower_integration_limit!uniform random selection of the variables for each iteration of the integration
		r_1=calculate_module(x_1,x_2,x_3)
		r_2=calculate_module(y_1,y_2,y_3)
		r_12=calculate_module(x_1-y_1,x_2-y_2,x_3-y_3)
		if(r_12<=proximity_threshold)cycle!if the electrons are too close the integral will blow up, so we skip dangerous integration points
		mc_sample=wavefunction(r_1,r_2,r_12)!sampling the wavefunction
		mc_integral=mc_integral+mc_sample
		mc_squared_sample_sum=+mc_squared_sample_sum+mc_sample**2
	end do
	call cpu_time(time_2)
	mc_integral=mc_integral/dble(number_of_mc_iterations)
	standard_deviation=dsqrt((mc_squared_sample_sum/dble(number_of_mc_iterations)-mc_integral**2)/dble(number_of_mc_iterations))
	write(*,*) 'The value for the mc integral is:',mc_integral*(upper_integration_limit-lower_integration_limit)**6,' which is',mc_integral*(upper_integration_limit-lower_integration_limit)**6*1.d2/answer,' % of the closed form value, and with a standard deviation of',standard_deviation*(upper_integration_limit-lower_integration_limit)**6,' . The time spent during integration was',time_2-time_1,' seconds.'
	write(10,*) 'The value for the mc integral is:',mc_integral*(upper_integration_limit-lower_integration_limit)**6,' which is',mc_integral*(upper_integration_limit-lower_integration_limit)**6*1.d2/answer,' % of the closed form value, and with a standard deviation of',standard_deviation*(upper_integration_limit-lower_integration_limit)**6,' . The time spent during integration was',time_2-time_1,' seconds.'

	!#####	MONTE CARLO INTEGRATION WITH IMPORTANCE SAMPLING STARTS HERE	#####

	call cpu_time(time_1)
	do i=1,number_of_mc_iterations
10		random_try=dran_u()
		if(random_try==1.d0)go to 10
		r_1=-2.5d-1*dlog(1.d0-random_try)
20		random_try=dran_u()
		if(random_try==1.d0)go to 20
		r_2=-2.5d-1*dlog(1.d0-random_try)!exponential random selection of the radial variables for each iteration of the integration, to avoid log(0), there are two loops to cycle over the generation of numbers
		phi_1=2.d0*pi*dran_u()
		phi_2=2.d0*pi*dran_u()
		theta_1=pi*dran_u()
		theta_2=pi*dran_u()!uniform random selection of the angular variables for each iteration of the integration
		r_12=dsqrt(r_1**2+r_2**2-2.d0*r_1*r_2*(dcos(theta_1)*dcos(theta_2)+dsin(theta_1)*dsin(theta_2)*dcos(phi_1-phi_2)))
		if(r_12<=proximity_threshold)cycle!if the electrons are too close the integral will blow up, so we skip dangerous integration points
		mc_sample=r_1**2*r_2**2*dsin(theta_1)*dsin(theta_2)/r_12!sampling the wavefunction
		mcis_integral=mcis_integral+mc_sample
		mc_squared_sample_sum=+mc_squared_sample_sum+mc_sample**2
	end do
	call cpu_time(time_2)
	mcis_integral=mcis_integral/dble(number_of_mc_iterations)
	standard_deviation_is=dsqrt((mc_squared_sample_sum/dble(number_of_mc_iterations)-mcis_integral**2)/dble(number_of_mc_iterations))
	write(*,*) 'The value for the mc integral with importance sampling is:',2.5d-1*pi**4*mcis_integral,' which is',2.5d-1*pi**4*mcis_integral*1.d2/answer,'% of the closed form value, and with a standard deviation of',2.5d-1*pi**4*standard_deviation_is,' . The time spent during integration was',time_2-time_1,' seconds.'
	write(10,*) 'The value for the mc integral with importance sampling is:',2.5d-1*pi**4*mcis_integral,' which is',2.5d-1*pi**4*mcis_integral*1.d2/answer,'% of the closed form value, and with a standard deviation of',2.5d-1*pi**4*standard_deviation_is,' . The time spent during integration was',time_2-time_1,' seconds.'
	write(10,"(/)")
	close(10)
end program project3


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