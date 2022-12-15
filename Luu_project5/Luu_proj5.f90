!	This is Project # 5
!
!
!	Purpose: The purpose of this program is to use the Runge-Kutta algorithm
!			 to deal with orbital mechanics. 
!
! TO USER: PLEASE READ THE READ ME FILE !!!!
!
!
!	Programmer: 		Date: 			Version:		Description:		
! - - - - - - - -	   - - - - - -	   - - - - - -		- - - - - - - 
!	Ken Luu			   11/26/2017		1.0				Created from the depths of Mordor

! CODING BEGINS ! 

	Module Params 
		
		Implicit None
		real :: M,m1,m2
		real, parameter :: G = 1.0! Newton's gravitational constant 
	
	End Module Params 
	
	Program Narsil 
 	
 ! Global Variables :
 Use Params 
 ! - - - - - - - - - !
 	
 	
! Local Variables : 
! - - - - - - - - - - - !	
	
	Implicit None
	
	integer :: n 		! PLEASE READ THIS! this is the number that symbolizes how many ODEs we have, USER can change if needed 	
	real, allocatable :: xvec(:)			! coordinates and velocities 	
	real, allocatable :: F(:)				! derivatives
	real t									! time 
	real dt									! time step 
	integer i,j 							! this is for the do loop that is dependent on the Number of steps 
	real :: E 								! the total energy of the system 
	real :: fhold							!function holding the derivatives
	real :: var 							! variables 
	real :: x1,x2,y1,y2,vx1,vy1,vx2,vy2		! coordinates and velocity 
	real :: r1,r2,r12						! radii functions 
	!real :: A,B,C,D,W 						! holding values for the declared function called "Fvel"
	integer :: Nsteps
	real :: L 								! angular momentum 
! - - - - - - - - - - - - !	
 	
 	print*," "
 	print*," "
 	Print*,"WELCOME TO PROJECT #5"
 	Print*,"If USER is confused, please look at the attached: Readme.txt " 
 	Print*,"Written by Ken Luu" 
 	print*," "
 	print*," "
 	
 	Print*,"Enter the stationary mass M:"					! stationary mass
 	Read*, M
 	
 	Print*,"Enter the mass m1 (mass of the planet):"		! mass of the planet	
 	Read*, m1
 	
 	Print*,"Enter the mass m2 (mass of the second planet/moon):" 	! mass of the planet or the moon
 	Read*, m2 
 	
 	Print*,"Enter the initial position(x1):" 				! initial position for the first mass (x-direction)
 	Read*, x1
 	
 	Print*,"Enter the initial position(y1):"				! initial position for the first mass (y-direction)
 	Read*, y1
 	
 	Print*,"Enter the initial position(x2):" 				! initial position for the second mass (planet or moon)
 	Read*, x2
 	
 	Print*,"Enter the initial position(y2):"				! initial position for the second mass (planet or moon)
 	Read*, y2
 	
 	Print*,"Enter the initial velocity(vx1):" 				! initial velocity for the first mass (x direction)
 	Read*, vx1
 	
 	Print*,"Enter the initial velocity(vy1):"				! initial velocity for the first mass (y direction)	
 	Read*, vy1
 	
 	Print*,"Enter the initial velocity(vx2):"				! initial velocity for the second mass (x direction)
 	Read*, vx2
 	
 	Print*,"Enter the initial velocity(vy2):"				! initial velocity for the second mass (y direction)
 	Read*, vy2
 	
 	Print*,"Please enter the dt value(time step): "			! time step 
 	Read*, dt
 	
 	Print*,"Please enter the number of steps:"				! the number of steps (lattice points)
 	Read*, Nsteps
 
 	Print*,"Enter the number of dependent variables  you want to solve for:" !in our case it is 8
 	read*, n
 
  allocate( xvec(1:n) )										 		! fills in the variables 
  allocate( F(1:n) )												! derivatives 
 
 	open(unit=10,file="onemass.dat",status="unknown")
 	open(unit=11,file="twomass.dat",status="unknown")
 	open(unit=12,file="energycons.dat",status="unknown")
 	open(unit=13,file="angmoment.dat",status="unknown")
 
 t = 0.0 														! initialization of the time step

! - - - - - - - - - - - - - - - - - - - - - - !
 Call Variables(n,x1,x2,y1,y2,vx1,vy1,vx2,vy2,xvec) 
 
! fvelocity = fvel(A,B,C,D,W)
 
  do i = 1,Nsteps ! this loop is with respect to time 
		
	do j = 1, n	 ! number of dependent variables 
		
		r1 = sqrt( (xvec(1))**2 + (xvec(3))**2 ) 
		r2 = sqrt( (xvec(2))**2 + (xvec(4))**2 )
		r12 = sqrt( (xvec(1)-xvec(2))**2 + (xvec(3)-xvec(4))**2 )
 
 		call Eval(n,r1,r2,r12,F,xvec) 
 			
 		var = xvec(j)
 			
		call Runge(j,var,n,dt,t,r1,r2,r12,F,xvec) 

  	 enddo

 call Energy(r1,r2,r12,n,xvec,E)
 call Aangmoment(xvec,n,r1,r2,L)

 write(10,*) xvec(1), xvec(3)
 write(11,*) xvec(2), xvec(4) 
 write(12,*) t, E
 write(13,*) t, L
		
	t = t + dt 									! time changing over time 

 enddo		 
 
 close(10)
 close(11)
 close(12)	
 close(13)
 
 end program Narsil 	
 	
Subroutine Variables(n,x1,x2,y1,y2,vx1,vy1,vx2,vy2,xvec)

! This subroutine is used to make the arbitrary vector "xvec" fill up the values that are 
! being declared 	
	
	USE Params
	implicit none 
	integer :: n
	real :: xvec(1:n) 
	real :: x1,y1,x2,y2,vx1,vy1,vx2,vy2
	
	xvec(1) = x1  							! x position of the planet
	xvec(2) = x2							! x position of the second mass
	xvec(3) = y1							! y position of the first mass
	xvec(4) = y2 							! y position of the second mass
	xvec(5) = vx1 							! velocity in the x direction of the first mass
	xvec(6) = vy1							! velocity in the y direction of the first mass
	xvec(7) = vx2							! velocity in the x direction of the second mass 
    xvec(8) = vy2							! velocity in the y direction of the second mass 
	
	return
	end subroutine Variables 

Subroutine EVAL(n,r1,r2,r12,F,xvec)

! SUBROUTINE EVAL: used for evaluating the derivatives
! each F vector fills up a derivative of the given ODEs which we can then pass on
! to solve in the Runge-Kutta 

! Global Variables: 
	USE Params 
! - - - - - - - - - - !
	
	Implicit None
! Passed Variables : 
	integer:: n
	real :: xvec(1:n)						! coordinate and velocity 
	real :: F(1:n) 						! derivatives 
	real :: r1,r2,r12					! radii functions
! - - - - - - - - - - ! 
! - - - - - - - - - 
! Functions: 
! - - - - - - - - - - - ! 

	F(1) =  xvec(5) 
	F(2) = 	xvec(7) 
	F(3) =  xvec(6)
	F(4) =  xvec(8)  
	F(5) = -(G*M/(r1**3))*xvec(1) - ((G*m2/(r12**3))*(xvec(1)-xvec(2)))
	F(6) = -(G*M/(r1**3))*xvec(3) - ((G*m2/(r12**3))*(xvec(3)-xvec(4)))
	F(7) = -(G*M/(r2**3))*xvec(2) - ((G*m1/(r12**3))*(xvec(2)-xvec(1)))
	F(8) = -(G*M/(r2**3))*xvec(4) - ((G*m1/(r12**3))*(xvec(4)-xvec(3)))
	
	return
	end subroutine EVAL

	Subroutine RUNGE(j,var,n,dt,t,r1,r2,r12,F,xvec)

 !Uses the 4th-order Runge-Kutta algorithm to evaluate the given ODEs ( 8 of them)
 !hence the variables given for the vectors, the numbers can be changed depending 
 !on the number of ODEs 

! Global Variables: 	
	USE Params 
! - - - - - - - - - - - !
	
	Implicit None
! Passed Variables: 
	integer :: n
	real :: xvec(1:n)									! the independent variables, coordinates of position and velocity
	real :: t,dt										! time, the independent variable 
	real :: Var 										! the variables: x1,y1,x2,y2,vx1,vy1,vx2,vy2
	real :: F(1:n)				    					! derivatives 
	real :: r1,r2,r12									! radii functions 

! - - - - - - - - - - -	!
	
! Local Variables: 
 integer :: j, l 									! dependent variable index	
 real :: k1(1:n)									! increments 
 real :: k2(1:n)									! increments 
 real :: k3(1:n)									! increments 
 real :: k4(1:n) 									! increments 
 real :: fhold 										! function holding the derivatives 
 	
	fhold(var,t) = F(j) 							! setting a "holding" function equalling the derivatives 
 	
! - - - - - - - - - !  
 
 do l = 1,n

 k1(l) = F(l)*dt
 k2(l) = fhold(xvec(l)+0.5*k1(l),t+0.5*dt)*dt
 k3(l) = fhold(xvec(l)+0.5*k2(l),t+0.5*dt)*dt
 k4(l) = fhold(xvec(l)+k3(l),t+dt)*dt

 enddo


 xvec(j) = xvec(j) + (1./6.)*( k1(j) + 2.0*k2(j) +2.*k3(j) + k4(j) ) ! our differential equations
	
return

end subroutine RUNGE 
		
! testing the conservation of energy 

Subroutine Energy(r1,r2,r12,n,xvec,E) 

! E is the energy that we are testing for 
! this subroutine makes sure that the energy is conserved 

USE Params

implicit none 
real :: E 									! energy
integer :: n 								! dependent variables
real :: xvec(1:n)							! holding the variables
real :: r1,r2,r12							! radii functions


E = 0.5*m1*(xvec(5)**2+xvec(6)**2)+0.5*m2*(xvec(7)**2 + xvec(8)**2)-G*(m1*M/r1) - G*(m2*M/r2)-G*(m1*m2/r12)  !total energy

return
end subroutine Energy

Subroutine Aangmoment(xvec,n,r1,r2,L)
USE Params

implicit none 
	! this is for the angular momentum 
	! We know that angular momentum: L = rxp (cross product) 
	! v is just the magnitude of the velocity 
	! m1 is the first mass, m2 the is second mass 
	! r1 is the radius with respect to the first mass
	! r2 is the radius with respect to the second mass						

real :: L 											! angular momentum 
real :: r1,r2										! radii 
integer :: n 										! index for the variables 
real :: xvec(1:n) 									! our arbitrary vectors holding our variables 

L = m1*(xvec(1)*xvec(6) - xvec(3)*xvec(5)) + m2*(xvec(2)*xvec(8) - xvec(4)*xvec(7)) ! angular momentum 

return
end subroutine Aangmoment 

