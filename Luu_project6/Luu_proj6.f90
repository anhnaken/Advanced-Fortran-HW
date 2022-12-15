
Program Gimli 


!
!	Purpose: The purpose of this project is to solve the time dependent Schrodinger equation
! 			 via the Crank Nicholson method. 
!
!	Programmer: 			Date:			Version:			Description:
! - - - - - - - - 		- - - - - - -      - - - - - - - 	   - - - - - - - - - - - - - - - - - - 
!	Ken Luu				12/20/17			1.0					Created from the depths of Mordor	
!
!	CODING BEGINS!
!

 
  
 Implicit None 															! dt: time step, k: spring constant, time: time iterations
 Real :: L,dx,x,k,dt,t													! L:size of box, dx:discretized steps,x:discretized x-direction
 Integer :: n,i,j,tsteps 												! n:lattice points,i:iteration point = j, tsteps:steps for dt 	
 real :: sigma 															! width of the Gaussian wave function
 real :: x0																! center of the Gaussian wave function 
 complex :: gaussian													! the given wave function
 real :: probdensity 													! the value of the probability density 
 real :: norm,expect,expect2,width 										! normalization, expectation value, spread 
 real :: awidth															! analytic result to the width of the Gaussian 
 real :: E																! energy for the infinite square well
 real :: SHOEN															! energy for the SHO 
 Real :: hbar 															! Planck's constant 
 Real :: m 																! mass
 Real :: pi 															! 3.14159...



 Real, dimension(:,:), allocatable :: V 								! the potential function
 Real, dimension(:,:), allocatable :: KE								! kinetic energy
 Real, dimension(:), allocatable :: psiRe								! real part of the wave function, psi. 
 Real, dimension(:,:), allocatable :: H 								! Hamiltonian matrix 
 Real, dimension(:,:), allocatable :: SuperMat1							! The "Super Matrix" 
 Real, dimension(:,:), allocatable :: supermat2 						! second super matrix 
 Real, dimension(:,:), allocatable :: supermatINV						! the inverted matrix 
 Real, dimension(:), allocatable :: index								! index for the separate matinv.f

 
 
 Print*,"~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  "
 print*," This is Project #6. "
 Print*," Written by Ken Luu. "
 Print*," "
 Print*,"~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  "
 
 print*,"Please enter the box size, L:"
 read*, L
 
 print*,"Please enter the # of lattice points, N:"
 read*, n
 
 print*,"Please enter the # of time steps:"
 read*, tsteps 
 
 print*,"Please enter the time step, dt:"
 read*, dt
 
 print*,"Please enter the width of the Gaussian:"
 read*, sigma
 
 print*,"Please enter the center of the Gaussian:"
 read*, x0 

 print*,"Please enter the spring constant value:"
 read*, k

 
! Some memory allocations: !
! V: Potential energy
! T: Kinetic energy 
! psi: wave function 
! e: eigenvalues 
! work: a one dimensional vector 
! H: the hamiltonian  
! SuperMat: the "Super matrix" 
 
 allocate( V(1:n,1:n) ) 
 allocate( KE(n+1,n+1) )
 allocate( psiRe(1:2*n) )
 allocate( H(1:n,1:n) )
 allocate( SuperMat1(1:2*n,1:2*n) )
 allocate( SuperMat2(1:2*n,1:2*n) )
 allocate( index(1:2*n) ) 
 allocate( supermatINV(1:2*n,1:2*n) )


 Open(unit=10, file='normalization.dat', status='unknown')
 Open(unit=11, file='expectation.dat', status='unknown') 
 Open(unit=12, file='width.dat', status='unknown')
 Open(unit=13, file='xvsProbx.dat', status='unknown')
 Open(unit=14, file='analytic.dat', status='unknown') 
 



dx = (2.0*L)/float(n)

pi = acos(-1.0)
hbar = 1.0
m = 1.0

probdensity = 0.
psiRe = 0.
supermatINV = 0.
supermat1 = 0.
supermat2 = 0. 
KE = 0.
V = 0.
H = 0. 


! filling out the matrix for the Woods Saxon potential 
do i = 1, n 

x = -L + dx*(i-1)

V(i,i) = 0.5*k*(x**2) 											! potential of a SHO 


KE(i,i) = (hbar**2)/(m*dx**2) 
KE(i,i+1) = -(0.5*hbar**2)/(m*dx**2)
KE(i+1,i) = -(0.5*hbar**2)/(m*dx**2) 

	do j = 1, n 

	H(i,j) = V(i,j) + KE(i,j)
	
	enddo

enddo

 
 call firstmatrix(dt,hbar,n,SuperMat1,H)									! calls the subroutine to fill in the supermatrix
 
 call secondmatrix(dt,hbar,n,SuperMat2,H)									! calls the subroutine to fill in the supermatrix
 
 call wavefunctions(pi,x0,sigma,dx,L,n,psiRe)								! this subroutine the psi vector and the gaussian wave function
 
 call matinv(SuperMat1,SuperMatINV,2*n,2*n,index) 							! matrix inversion
 
t = 0. !initialization of the time steps 

 do j = 1,tsteps ! this do loop is looping over the number of time steps 
 
 norm = 0. ! initializing the summation of the normalization
 expect = 0. ! initializing the summation of the expectation value: <x>
 expect2 = 0. ! initializing the summation of the expectation value: <x^2>
 probdensity = 0. ! initializing the summation of the probability density: rho 
 
 
 do i = 1, n ! the start of the x-loop
 
 call amplitudesqrd(n,i,psiRe,probdensity)
 
 x = -L + dx*(i-1)
 
 call normalization(n,probdensity,dx,norm)
 
 expect = expect + (x*probdensity*dx) 								    ! this is the expectation value: <x>
 
 expect2 = expect2 + ((x**2)*probdensity*dx)								! this is the expectation value: <x^2>
 
 enddo ! the end of the x-loop

 expect = expect/norm 													! the expectation value,<x>, divided by the normalization
 expect2 = expect2/norm													! <x^2> being divided by the normalization 
 
 width = sqrt(expect2 - (expect**2))									! this is the width of the gaussian 


 write(10,*) t, norm 									! this write-statement writes the time and the normalization to a file
 write(11,*) t, expect 									! this write-statement writes the time and expectation value to a file
 write(12,*) t, width									! this write-statement writes the time and spread to a file 

 call squarewell(m,L,hbar,pi,i,E)						! calling subroutine to check energy conservation for the infinite square well
 call SHO(i,K,m,hbar,SHOEN) 							! calling subroutine to check for energy conservation for SHO 

! this format statement is to put on the user's screen to see if the energy does indeed stay constant with time 
 
write(*,*) t
writ(*,*) E
write(*,*) SHOEN 



 psiRe = matmul(superMat2,psiRe) 						! this is the matrix multiplication of the "supermat2" and the psi vector
 psiRe = matmul(SuperMatINV,psiRe)						! this is the matrix multiplication of the inverse matrix and the psi vector
 
 
 t = t + dt 										! time iterations
 
 enddo ! the end of the time loop
 
! this do loop is for the probability at the last time step  

do i = 1,n

x = -L +dx*(i-1)

call amplitudesqrd(n,i,psiRe,probdensity)
 
 write(13,*) x, probdensity 
 
enddo 


t = 0. ! initialize t 
 
 do i = 1, tsteps 
 x = -L+dx*(i-1)
 
 call analyticspread(pi,sigma,hbar,t,m,awidth)
 
 write(14,*) t, awidth  
 
 t = t + dt
 
enddo
 
close(10)
close(11)
close(12)
close(13)
close(14)
 
 
 end program Gimli 


Subroutine firstmatrix(dt,hbar,n,SuperMat,H)
! this subroutine is part of the whole super matrix, i simply broke it into two parts. 
implicit none
integer :: n,i,j 
real :: H(1:n,1:n)
real :: SuperMat(1:2*n,1:2*n)
real :: dt,hbar 


do i = 1, n

 do j = 1, n 
 SuperMat(i+n,j) = -dt*H(i,j)/(2.0*hbar)
 SuperMat(i,j+n) = dt*H(i,j)/(2.0*hbar)
 SuperMat(i,i) = 1.0
 SuperMat(j+n,j+n) = 1.0

 
! this subroutine is for the first "super matrix"  

  enddo
enddo

return 
end subroutine firstmatrix 


Subroutine secondmatrix(dt,hbar,n,SuperMat,H)
 ! this subroutine fills out the super matrix notice that it has a different value for the 
 ! "SuperMat(i,j+n)" different from the previous filled out matrix 
implicit none
integer :: n,j,i
real :: H(1:n,1:n)
real :: SuperMat(1:2*n,1:2*n)
real :: dt,hbar 


do i = 1,n

 do j = 1,n


 SuperMat(i,j+n) = -dt*H(i,j)/(2.0*hbar)
 SuperMat(i+n,j) = dt*H(i,j)/(2.0*hbar)
 SuperMat(i,i) = 1.0
 SuperMat(j+n,j+n) = 1.0

 enddo
enddo

return 
end subroutine secondmatrix 

 Subroutine wavefunctions(pi,x0,sigma,dx,L,n,psiRe)
 ! this subroutine fills the one dimensional psi vectors 
 ! Also, the given wave function "Gaussian" is given here. 
 implicit none
 real :: x,dx,L,x0,sigma,pi 
 integer :: i,n 
 real :: psiRe(1:2*n)
 complex :: gaussiann 


 do i = 1,n

 x = -L + dx*float(i-1)
 
 gaussiann = ((2.0*pi*(sigma**2.))**(-1.0/4.0))*exp(-((x-x0)**2.)/(4.0*(sigma**2.)))
 
 psiRe(i) = real(real(gaussiann)) 
 psiRe(i+n) = real(aimag(gaussiann)) 					     ! aimag is an intrinsic command 
 
 enddo 

 return
 end Subroutine wavefunctions 


Subroutine amplitudesqrd(n,i,psireal,probdensity)
! This subroutine calculates the probability density  
 implicit none
 integer :: n,i
 real :: probdensity 
 real :: psiReal(1:n)
 
 
 probdensity = (psiReal(i))**2 + (psiReal(i+n))**2

 return 
 end subroutine amplitudesqrd
 
 
 Subroutine normalization(n,probdensity,dx,norm)
 ! This subroutine calculates the normalization 
 ! it will be returned to the main program
 implicit none 
 real :: norm, probdensity,dx 
 integer :: i, n 

 norm = norm + (probdensity*dx)

return
end subroutine normalization 

Subroutine SquareWell(m,L,hbar,pi,i,E)
implicit none 
real :: E														! energy of a square well 
real :: hbar,m,L 												! Planck's constant, mass, size of the box
integer :: i 
real :: pi

! this is the energy of the infinite square well, but since our box size is 2L and squared
! we have an 8 on the denominator 

E = ((i**2)*(hbar**2)*(pi**2))/(8.0*m*L**2) 
 
return 

end subroutine SquareWell 

Subroutine SHO(i,K,m,hbar,SHOEN)
! this subroutine is to check for energy conservation of the SHO
implicit none 
real :: SHOEN 										! energy of an harmonic oscillator  
real :: hbar,K,m									! Planck's constant, spring constant, mass
integer :: i 

SHOEN = (0.5+i)*hbar*sqrt(k/m)

return 

end subroutine SHO 

Subroutine analyticspread(pi,sigma,hbar,t,m,awidth)
implicit none
real :: pi,sigma,hbar
real :: t,m,awidth

awidth = sqrt(sigma**2+((hbar**4)*(t**2))/(4.**m**2)*(sigma**2))

return

end subroutine analyticspread






	