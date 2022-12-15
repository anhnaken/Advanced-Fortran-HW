!THIS IS MY ADVANCED PROJECT 
!
!
!
! Purpose: The purpose of this program is calculate neutron flux from a reactor via a triple integral 
!			using Boole's method and the Monte Carlo Method. 
!
!	Programmer: 		Date:			Description:
!	- - - - - -			- - - - -		- - - - - - -
!	Ken Luu				09/28/17		Created from the depths of Mordor
!
!	Quantity:			Units:			Description:
!	Input:
!	- - - - - 			- - - - -		- - - - - - - 
!	N					unitless		dimension of our loops
!	Msamp				unitless		dimension of the Monte Carlo method 
!	W					meters			width of the box	
!	Height				meters			height of the box
!	D					meters 			depth of the box
!	x0					meters			initial point chosen starting from a corner in the x-direction
!	y0					meters			initial point chosen starting from a corner in the y-direction
!	
!
!	Output:
!	- - - - -
!   QuadFlux			Neutrons/sm^3	this is the flux via Boole's method
!	ApproxFlux			Neutrons/sm^3   this is the approximated flux from the equation
!	Stndev				unit-less		this is the standard deviation of the Monte Carlo
!	MCQuad				Neutron/sm^3	this is the flux via the Monte Carlo method
!

!* * * * * * * * * * * * * * CODING BEGINS * * * * * * * * * * * * * * * * * * * *!


!THIS is for the HOLLOW sphere, if you look at the subroutine that does the triple integral and 
!the monte carlo you can see the if statement that deals with the hollow sphere 


Program Quadrature
	
	Implicit None
	Real:: f_xyz, x, y, z, dx, dy, dz, SumBoole, a, QuadBoole, array1, array2, array3
	Real:: W, Height, D, x0, y0 !these are dimensions and starting point relating to the reactor
	Real:: QuadFlux !the integral we want 
	Integer:: N, Nz, Ny, Nx !these are declared for our do loops and lattice points (N)
	Integer, Parameter:: Np = 500 !this is the dimension of our array
	Real:: array(0:Np) ! this is the general array for our Boole's Subroutine
	Real:: arrayx(0:Np)
	Real:: arrayy(0:Np)
	Real:: arrayz(0:Np)
	Real:: ApproxFlux, fxyz, h, pi, QuadZYX, MCQuad, sum, stndev, sum2, standard_deviation, R
	
	
	! Monte Carlo variable declarations 
	
	Integer:: seed,msamp 
	
	
	
	
	Write(*,*)"COMPUTATIONAL ADVANCED PROJECT #2"
	Write(*,*)"This program was written by Ken Luu"
	
	
! ***** NOTE ABOUT USER INPUTS *****
! - - - - - - - - - - - - - - - - - !
! N: this is our lattice points do our do loops, I fixed this to be 100 !
! Height: I fixed this value at 1.	!
! W: I fixed this value at 2.	!
! D: I fixed this value at 3.	!
! y0: I fixed this value at 1.	!
! x0: this was my varying value starting from 0 !
! seed: any "random" numbers that I wish to choose
! msamp: I fixed this value at 200 because I did not want my lattice points for the Monte Carlo !
!		 to be the same lattice points as my triple integral!
!
	
	
 
	
	
	
	Write(*,*)"Enter the N(Lattice points) value:" 
	Read(*,*) N
	Write(*,*) N
	
	
	Write(*,*)"Enter the H(height) value:"
	Read(*,*) Height
	
	
	Write(*,*)"Enter the W(width) value:"
	Read(*,*) W
	
	
	Write(*,*)"Enter the D(depth) value:"
	Read(*,*) D
	
	
	Write(*,*)"Enter x0 value:"
	Read(*,*) x0
	
	Write(*,*)"Enter y0 value:"
	Read(*,*) y0
	 
	
    Write(*,*)' Enter seed for random numbers '
    read(*,*) seed
	
	Write(*,*)'Enter sample(Lattice points for Monte Carlo) number '
	read(*,*) msamp
	
	Write(*,*)' Enter the radius value:'
	read(*,*) R
	
	
	
	pi = acos(-1.)
	

	

 	
	
	
	Call THREELOOPSADV(Np,Height,W,D,x0,y0,N,Nz,Ny,Nx,h,QuadBoole,QuadFlux,R)
	Call MonteCarloADV(W,Height,D,x0,y0,MCQuad,seed,msamp,stndev,R)
	
	Standard_deviation = stndev ! the standard deviation 
	sum = MCQuad
	QuadZYX = QuadFlux !this is on line 137  
	ApproxFlux = (W*Height*D)/(4.*pi*((Height/2)**2 + (x0 + D/2)**2 + (y0 - W/2)**2))

!These write statements will graph my approximated flux and my flux calculated via the 
!Boole method of integration. For my values of x0, I will manually input these values.
!the goal is to see what happens when x0 gets large as the Boole value should approach the 
!approximated value. I will start from x0 = 0 then go in small increments of .1 until I reach 1
!then I go to up 10 in steps of 1; I then continue going up until my Quadrature methods converge
!which you can see if you were to put in big values of x0
  

	


!these are our new values for the hollow sphere 

	Write(*,101) QuadFlux, ApproxFlux, MCQuad, stndev
	101 format("Boole's Quadrature = ", F20.5,/, "Approximated Flux =", F20.5,/, "MC =",F20.5&
            /,"standard deviation =",F20.5)  
	
	
	End Program Quadrature
	



Subroutine THREELOOPSADV(Np,Height,W,D,x0,y0,N,Nz,Ny,Nx,h,QuadBoole,QuadZYX,R)

!	Input:
!		Np: declared dimension of array
!		N: used step for Boole (multiple of 4)
!		array(0:Np): array of points of function to be integrated 
!		h: will be our dx, dy, dz
!		Nx: steps in the x-direction
!		Ny: steps in the y-direction
!		Nz: steps in the z-direction
!	Intermediate:
!		
!
!
!	Output:
!		QuadBoole : integral
	
	
	
	
	
	Implicit None
	Integer:: Nz, Ny, Nx, i, j, k, Np, N 
	Real:: array_x(0:Np)
	Real:: array_y(0:Np)
	Real:: array_z(0:Np)
	Real:: array(0:Np)
	Real:: Func_xyz, f_xyz
	Real:: x, y, z, dx, dy, dz, h, QuadBoole, QuadFlux, pi, QuadBoolex, QuadBooley, QuadBoolez, SumBoole
	Real:: x0, y0, f, D, Height, W, QuadZYX, R 

pi = acos(-1.)
	


! - - - - - - - - - - - - - - -  Boole's Quadrature - - - - - - - - - - - - - - !

Do i = 1, N ! loop over x
	
	dx = (D)/float(N)
	x = float(i)*dx
	
	Do j = 1, N ! loop over y
		
		dy = (W)/float(N)
		y =  float(j)*dy
			
			Do k = 1, N ! loop over z
			

			
			dz = (Height)/float(N)
			
			z = float(k)*dz
		
		   		If((x**2 + y**2 + z**2) .GT. R) then
			
			
					f_xyz = 1./((4.*pi) * (z**2 + (y-y0)**2 +(x+x0)**2))
	
			 
			
					Array_z(k) = f_xyz !filling our array for z
			
				End If		
			End do !end of loop over z
		
	
	Call Booles(SumBoole,Np,N,array_z,dz,QuadBoolez)
		
	 Array_y(j) = QuadBoolez  !filling our array for y integral 
		
	End Do ! end of loop over y

Call Booles(SumBoole,Np,N,array_y,dy,QuadBooley)

	

  Array_x(i) = QuadBooley !filling our array for x integral 

End Do !end of loop over x


Call Booles(SumBoole,Np,N,array_x,dx,QuadBoolex) ! We now solve our flux integral 

QuadBoole = QuadBoolex


QuadZYX = QuadBoolex		

	

Return 

End Subroutine THREELOOPSADV 


!!!!!!!!!!!!!! MONTE CARLO !!!!!!!!!!!!!!!

Subroutine MonteCarloADV(W,Height,D,x0,y0,sum,seed,msamp,standard_deviation,R)
	
	Implicit None
	
! - - - - - - -- -- - - - - - - - -Variables for random - - - - - - - - - - -- - - !
	Real:: ran3, R
	Integer:: seed 
	
	!Variables for integration
	
	Real:: sum,sum2,W,Height,D,x0,y0,f_xyz,x,y,z,standard_deviation
	Integer::i,msamp,pi
! - - - - - - - - - - - - - - - -- - - - - - - - -- - - - - - -- - - - - - - - !
	
	pi=acos(-1.)

	

	sum = 0.0
	sum2= 0.0 
     
      Do i = 1,msamp
        
         x = D*ran3(i*seed)
         y = W*ran3(2*i*seed)
         z = Height*ran3(3*i*seed)
    
    If((x**2 + y**2 + z**2) .GT. R) Then
    
        f_xyz = 1./((4.*pi) * (z**2 + (y-y0)**2 +(x+x0)**2))
        
        sum = sum + f_xyz
		sum2 = sum2+ (f_xyz)**2
		
		End If 
		End Do 
      
      sum = (sum/msamp)*(W*Height*D)  !Monte Carlo Quad, which will be returned to the main program 
      
      sum2 = (sum2/msamp) - (sum/(W*Height*D))**2
    
    Standard_deviation = W*Height*D*sqrt(sum2/msamp) ! the standard deviation or the error, it gets smaller quickly
    												 ! as I increase my x0 values, which is what I wanted because there
    												 ! is less spread in the values, meaning that my methods are acceptable. 
	
	Return 
	
	End Subroutine MonteCarloADV 





Subroutine Booles(SumBoole, Np, N, array_z, h, QuadBoole) 

!	Input:
!		Np: declared dimension of array
!		N: used dimension of array (multiple of 4)
!		array(0:Np): array of points of function to be integrated 
!		h: will be our dx, dy, dz
!
!	Output:
!		QuadBoole : integral

	Implicit None
	Integer:: Np, N !Dimensions
	Real array_z(0:Np)
	Integer:: i
	Real:: SumBoole, h, QuadBoole
	

!............. ERROR TRAP...........! 
!..... check that N is multiple of 4.....!



      if( (N/4)*4 .ne. N)then
        print*,' N is not a multiple of 4 ',N
        stop
      endif
!.............End ERROR TRAP..........!
	
	SumBoole = 0. !initializing our summation
	
	Do i = 4, N, 4
	
	
	SumBoole = SumBoole+((7.*array_z(i))+(32.*array_z(i-1))+(12.*array_z(i-2))+(32.*array_z(i-3))+(7.*array_z(i-4))) 
	
	End Do
	
	QuadBoole = ((2.*h)/(45.))*(SumBoole)
	
	Return
	
	End Subroutine Booles
	


	