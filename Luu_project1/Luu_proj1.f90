!THIS IS COMPUTATIONAL PROJECT 1
!Written by: Ken Luu
! 	 
! Purpose: The purpose of this program is to derive the given function using the 3-point
!			formula, the 5-point formula, and to get a feel of what is happening to this 
!			derivative with different methods 
!
!	Programmer: 		Date:			Description:
!	- - - - - -			- - - - -		- - - - - - -
!	Ken Luu				09/07/17		Created from the depths of Mordor
!
!	Quantity:			Units:			Description:
!	Input:
!	- - - - - 			- - - - -		- - - - - - - 
!	dx					unit-less		this is the limit for the derivative or our "h"
!	x					unit-less		this is the assigned x-value according to our last name
!
!	Output:
!	- - - - -
!	threept_dffdxx (our 3-point derivative of f(x)) 
!	fivept_dffdxx (our 5-point derivative of f(x))
!	dffdxx_analytic (this is the analytical derivative of the assigned f(x))
!
!	EXPERT PROJECT ANSWER: For the power law dependence I got 1.8905 for my 3 point formula
!							and 3.707 for my 5pt formula. This is what I expected because it
!							matches with the formula given where for the three point we have O(h^2)
!							and for the five point formula we have O(h^4) and my numbers were pretty close
!							to the expected values. 
!
! - - - - - - - - - - Coding Begins	- - - - - - - - - - - - -
!
	Program Derivatives
	Implicit None
	Real:: x, dx, threept_Doubleprime, fivept_Doubleprime, threept_dffdxx, fivept_dffdxx 
	Real:: dffdxx_analytic, difference_3pt, difference_5pt, fx, f_x
		
	
	Open(unit=10, file='3point.dat', status='unknown', position='append')
	Open(unit=11, file='5point.dat', status='unknown', position='append') 
	Open(unit=12, file='error3pt.dat', status='unknown', position='append')
	Open(unit=13, file='error5pt.dat', status='unknown', position='append')  
	
	
	Write(*,*) "Enter the value x" !since my last name starts with the letter L, my x value = 12
	Read(*,*) x
	
	
	!for the range of dx I started from the value 0.0001 up to .1 in small increments of 1 in
	!accordance with the respective decimal number, for example from 0.0001 I went to 0.0002 and then
	!0.0003. From 0.001 I went to 0.002 and then 0.003, etc. until I got to my final dx: 10
	!I had to do this many times because when I went to lower dx's my plot looked really really strange
	!so I had to choose something that does not blow up the axis too much 
	 
	Write(*,*) "Enter the value dx"
	Read(*,*) dx 
	
	fx = f_x(x)
	
	!here I call upon some subroutines so that I can solve my derivatives 
	Call threepoint(x, dx, threept_Doubleprime, f_x) !3doubleprime is the second deriv. of the 3pt form
	Call fivepoint(x, dx, fivept_Doubleprime, f_x)	!5doubleprime is the second deriv. of the 5pt form
	

	
	!here are the return values from my subroutine
	threept_dffdxx = threept_Doubleprime
	fivept_dffdxx = fivept_Doubleprime
	dffdxx_analytic = 2.0*cos(x) - x*sin(x)
	 
	 !these are are error values, we are checking if it is proportional to dx^2
	 
	 difference_3pt = abs(dffdxx_analytic - threept_Doubleprime) 
	 difference_5pt = abs(dffdxx_analytic - fivept_Doubleprime)
	 
  !here we write the results to our open statements so we can plot
	
	Write(10,*) dx, threept_Doubleprime
	Write(11,*) dx, fivept_doubleprime 
	Write(12,*) dx, difference_3pt
	Write(13,*) dx, difference_5pt
	
	
	Close(10)
	Close(11)
	Close(12)
	Close(13)
	
	!this write statement and formatted statement below is a way for me to check if my derivative formula matches 
	!my analytical answer and how much the answer varies for different quantities of dx's
		
	Write(*,101) dx, dffdxx_analytic, threept_Doubleprime, fivept_dffdxx
	101 format("dx = ", F8.5, /, "Analytic Derivative = ",F20.8, /,  "3pt derivative = ",F20.8, &
	          /, "5pt derivative = ",F20.8)
	 
	
	
	End Program Derivatives 			

! these subroutines are called upon to solve the 3point formula and the 5point (below)
	
	Subroutine threepoint(x, dx, threept_dffdxx, f_x)
	Implicit None
	Real:: x, dx, threept_dffdxx, f_x
	
	threept_dffdxx = (f_x(x+dx) - 2.*f_x(x) + f_x(x-dx))/(dx**2)
	
	
	Return
	End Subroutine threepoint
	
	Subroutine fivepoint(x, dx, fivept_dffdxx, f_x)
	Implicit None
	Real:: x, dx, fivept_dffdxx, f_x
	
	fivept_dffdxx=(1./(12.*dx**2))*(((-f_x(x-2.*dx)))+(16.*f_x(x-dx))-(30.*f_x(x))&
	               + (16.*f_x(x+dx))-(f_x(x+2.*dx)))
			
	
	Return
	End Subroutine fivepoint 
	
! I am doing a subprogram for the given function to plug it into my subroutine.
	
 Real Function f_x(x)
	Implicit None
	Real:: x
	
	f_x = x*sin(x)
	
	Return
	
	End Function f_x

	
	 