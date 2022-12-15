!	This is Project # 4 
!
!
!	Purpose: The purpose of this program is to calculate numerical results for the Schrodinger 
!			 equation to the Woods-Saxon potential, our ultimate goal is to find eigenenergies 
!			 and eigenfunctions for this specific potential for values that approximate the nuclear
!			 potential alongside with excited states. 
!
!	Programmer: 		Date: 			Version:		Description:		
! - - - - - - - -	   - - - - - -	   - - - - - -		- - - - - - - 
!	Ken Luu			   11/07/2017		1.0				Created from the depths of Mordor
!
! Variables used:		Units:				Description:
! - - - - - - - -      - - - - - -			- - - - - - - - 
! 	h 					MeV*fm		    	reduced Planck's constant 
!	fm					m 					1 fermi is 1.0x10^-15 meters 
!	m					MeV					the mass, technically it is MeV/c^2 but c = 1
!	V0					MeV					the "depth"
!	Psi					unit-less			the wave function
!	dx					unit-less			discretized step
!	R					m					radius of the well 
!	L					m					the size of the box (at least twice R)
!	N					unit-less			the number of integration points 
!	V					MeV					Our potential function 
!	a					m					surface thickness 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
! Note:
! - - - - - - - - 
! Other variables will be used, probably, but will label what they are where I used them.	
!
!
! ~ ~ ~ ~ ~ ~ ~ CODING BEGINS ~ ~ ~ ~ ~ ~ ~ !

 Module SomeValues 
 
 Implicit None 
 Real, parameter :: a = .2
 Real, parameter :: V0 = 50.0 !MeV
 Real, parameter :: hbar = 197.3 ! MeV*fm
 Real, parameter :: m = 939.0 ! MeV
 Real, parameter :: pi = acos(-1.) 
 
 

 
 End Module SomeValues
 

 Program Melkor 
 
 Use SomeValues 
 
 Implicit None 
 Real :: R,L,dx,x,k
 Integer :: n,np,i,j 
 
 
! Here I dynamically allocating my kinetic and potential energy to make them matrices !

 Real, dimension(:,:), allocatable :: V
 Real, dimension(:,:), allocatable :: T
 Real, dimension(:,:), allocatable :: psi
 Real, dimension(:), allocatable :: e
 Real, dimension(:), allocatable :: work
 Real, dimension(:,:), allocatable :: H 
 

 Open(unit=10, file='gs.dat', status='unknown', position='append')
 Open(unit=11, file='firstexcited.dat', status='unknown', position='append') 
 Open(unit=12, file='secondexcited.dat', status='unknown', position='append')
 Open(unit=13, file='thirdexcited.dat', status='unknown', position='append')

 Print*,"~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  "
 print*," This is Project #4. "
 Print*," Written by Ken Luu. "
 Print*," "
 Print*," About this project: "
 Print*," "
 Print*," This program deals with the Woods Saxon potential, we are solving for "
 Print*,"  eigenenergies and plotting them. "
 Print*,"~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  "
 
 
 
 
 !What I inputed: 
 ! I started from R: 2-10 in steps of 1
 ! for my L values, I varied it each time to see the different variations
 ! i.e L= 2, 4 ,8, 20, 25, 31, 40, 41, 43, 50
 ! for my N values I did: 100- 200 in steps of 10, however from 180 I went to 200
 
 
 write(*,*)"Enter the number of points on lattice:"
 read(*,*)N

Print*,"-------------------------------------------------------------------------------"
 
 write(*,*)"Enter the radius:"
 read(*,*)R
 
Print*,"-------------------------------------------------------------------------------"
 
 
 Write(*,*)"Enter size of the box:"
 read(*,*)L 

Print*,"-------------------------------------------------------------------------------"
 
 
 
! Some memory allocations: !
! V: Potential energy
! T: Kinetic energy 
! psi: wave function 
! e: eigenvalues 
! work: a one dimensional vector 
! H: the hamiltonian 
 
 
 allocate( V(n,n) )
 allocate( T(n,n) )
 allocate( psi(n,n) )
 allocate( e(n) )
 allocate( work(n) )
 allocate( H(n,n) )




dx = (2.0*L)/float(N-1)


! filling out the matrix for the Woods Saxon potential 
do i = 1, n 

x = -L + dx*float(i-1)

V(i,i) = -V0/(1.0+exp( (abs(x) - R)/a ) )


enddo

! filling the kinetic energy matrix 
do i = 1, n 

x = -L + dx*float(i-1)
 
T(i,i) = (hbar**2/(m*dx**2)) 
T(i,i+1) = -0.5*hbar**2/(m*dx**2)
T(i+1,i) = -0.5*hbar**2/(m*dx**2) 

	do j = 1, n 

	H(i,j) = V(i,j) + T(i,j)
	
	enddo

enddo

 
! our Hamiltonian Matrix ! 

!calling the subroutines to solve for the eigenenergies  
 call eig(H,n,n,e,psi,work)

 call eigsrt(e,psi,n,n)
 


! this will write to the output file for the radius and the excited states to graph
 Write(10,*) R, e(1)
 Write(11,*) R, e(2)
 Write(12,*) R, e(3)
 Write(13,*) R, e(4) 


 close(10)
 close(11)
 close(12)
 close(13)
 
 End program Melkor