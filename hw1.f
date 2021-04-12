	parameter (mx=5000)
	real H,D,k_s,nu,rho,L,g
	real f, Re
	real inc
	data eps/1.0E-3/ L/50./
	real Tf,dt
	real V(mx),T(mx)
	character string*5,filename*9,len*3
	integer index

c..Arrange the output file name
      write(string,'(f5.3)') L/1000
      read(string,'(2x,a3)') len
      filename = 'L-'//len//'.dat'
	open(1,file=filename)

c..Final time and time step
	Tf = 60
	dt = 0.02
	index = Tf/dt+1
c..Flow parameters and dimensions
	H   = 8 
	D   = 0.3
	k_s = 1.0E-4
	nu  = 1.0E-3
	rho = 1000
	g   = 9.81
c..Initial Velocity and Time
	V(1) = 0.
	T(1) = 0.
c..Loop to calculate steady-state velocities
	do i=1,index
		Re  = (rho*V(i)*D)/nu
		f   = 0.25/((log10((k_s/(3.7*D))+(5.74/(Re**0.9))))**2)
		inc = H-((1+f*(L/D))*((V(i)**2)/(2*g)))
		V(i+1) = V(i) + dt*inc*(g/L)
		T(i+1) = T(i) + dt
		write(1,'(2f15.10)') V(i),T(i) ! get output data
	enddo	
c..Output the steadt-state vel. and corresponding time
	print*, "Steady vel","    time"
	print*, V(index),T(index)

	close(1)
c..Find the convergence time and speed
	k=0
	eps1=1
	do while(eps1.gt.eps)
		k=k+1
		eps1 = (V(index)-V(k))/V(index)
	enddo
	print*, k, T(k), V(k)

	stop
	end