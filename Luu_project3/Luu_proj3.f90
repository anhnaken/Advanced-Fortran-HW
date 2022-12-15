!Written by Ken Luu 
!date : 11/04/ 2017
! description: created from the depths of Mordor
!
!
! Purpose: The purpose of this program is to get use to dealing with matrices
! in this program we deal with the semi empirical mass formula, then we deal with the drip lines
! we also deal with the coefficients of the formula and we calculate the uncertainties in
! the coefficients and the binding energy 


	PROGRAM Nymeria
	
	Implicit None
	
	Character(len=70) file_name 
	
	
!Input values: NN stands for the input neutron number
!              ZZ stands for the input proton number
!			   K stands for the input matrix dimension
!			   i,j are used as running index for do loop 	
!			   M is the number of data points from the file, in this case it is 1842
!				delta is for the kronecker delta 	
	Integer :: NN, ZZ
	Integer :: i, j, h, K, M
	real :: Binden, BindenUncertainty
	Real :: BindEnMin
	Real :: delta 
	Integer :: NNN, ZZZ
	

! - - - - - - - - - - - - - - - - - - !


	

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !





! Now Allocating my arrays , these are my arrays for the Proton number and Neutron number 
! that is being read in from the file 

	Integer, dimension(:),  allocatable :: Z
	Integer, dimension(:),  allocatable :: N
! - - - - - - - - - - - - - - - - - - - - - - - !


! This section of memory allocations is for the matrices that will be used to solve the 
! problem. 

! A: this is the A(alpha,beta) matrix that we have to acquire
! AINV: this is the inverse matrix of our A matrix
! G: this is functions of Z and N that relates to solving our coefficients, which is a matrix
! GT: the transpose of G
! C: this will be used to solve for the coefficients 
! BE: filling in the array for the binding energy from the given data points
! BEerror: filling in the array for the error bars in our binding energy
! tempmatrix: this is the temporary matrix that will be used for the matrix inversion 
! INDX: something from the matinv that keeps track of the matrix inversion algorithm 

	Real, dimension(:,:), allocatable :: A
	Real, dimension(:,:), allocatable :: AINV
	Real, dimension(:,:), allocatable :: G
	Real, dimension(:,:), allocatable :: Gerror 
	Real, dimension(:), allocatable :: C
	Real, dimension(:), allocatable :: BE
	Real, dimension(:), allocatable :: BEerror
	Real, dimension(:,:), allocatable :: ATEMP
	Integer, dimension(:), allocatable :: INDX
	Real, dimension(:), allocatable :: tempproduct
	Real, dimension(:,:), allocatable :: Coef
	Real, dimension(:,:), allocatable :: GT
	Integer, dimension(:), allocatable :: index
	Real, dimension(:), allocatable :: UncertainC

 Open(Unit=1000, file="Zdrip.dat", status="unknown")
 Open(Unit=2000, file="Ndrip.dat", status="unknown")


! - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 	
	
	Write(*, "(A)", ADVANCE = "NO") "Enter the name of the file:"
	Read(*, "(A)") file_name
 
	
	Write(*,*)"Enter the dimension of the Matrix:"
	Read(*,*) K
	
	
	Write(*,*)"Enter the Number of Protons:"
	Read(*,*) ZZ
	
	Write(*,*)"Enter the Number of Neutrons:"
	Read(*,*) NN
		
	
!opening the data file to read in the data 

Open (UNIT=100, FILE=file_name, STATUS="old")

Read (100, *) M

	allocate( Z(1:M) )
	allocate( N(1:M) )
	allocate( BE(1:M) )
	allocate( BEerror(1:M) )


! here I am reading in the Data file BE.dat 
Do i = 1, M
	
	Read (100,*) Z(i), N(i), BE(i), BEerror(i)

End Do

Close(100)
	
! here I am allocating my arrays 	
	allocate( A(K,K) )
	allocate( AINV(K,K) )
	allocate( G(K,M) )
	allocate( Gerror(K,M) )
	allocate( C(1:K) )
	allocate( ATEMP(K,K) )
	allocate( INDEX(K) )
	allocate( tempproduct(1:K) )
	allocate( coef(K,1) )
	allocate( GT(M,K) )
	allocate( uncertainC(1:K) )

! Calling my subroutines ! 
Call Getmatrix(Tempproduct,BE,BEerror,G,Z,N,M,K,A)

Call matinv(A,AINV,K,K,INDEX)

Atemp = matmul(A,Ainv)
print*, Atemp

! This is our coefficients, I know these are incorrect but I've been trying different methods 
! many times, I even compared this with the "Matmul" intrinsic operation and surprisingly got 
! the same answer, I compared it with peers and we checked our codes, every allocation seems correct
! and it smells correct but of course as we know just because it SEEMS correct does not mean 
! that it is correct. I'd like to see where I went wrong, I will discuss it with you. Thank you

print*,"Coefficients"
Print*," "
Print*," "
C = matmul(tempproduct,Ainv)

Print*,"Cvol=", C(1), "KeV"
Print*,"Csurf=", C(2), "KeV"
Print*,"Csym=", C(3), "KeV"
Print*,"Ccoul=", C(4), "KeV"
Print*,"Cpair=", C(5), "KeV"

Print*," "
Print*," "
Print*,"Uncertainty in coefficients"

UncertainC = 0.0

Do i = 1, K
UncertainC(i) =  UncertainC(i) + sqrt(Ainv(i,i))
End Do

Print*," "
Print*," "

Print*,"Cvol Uncertainty:", UncertainC(1)
Print*,"CsurfF Uncertainty:", UncertainC(2)
Print*,"Csyml Uncertainty:", UncertainC(3)
Print*,"Ccoul Uncertainty:", UncertainC(4)
Print*,"Cpair Uncertainty:", UncertainC(5)


if( ( mod(nn,2).eq.0 ) .and. ( mod(zz,2) .eq. 0) ) then
	delta = 1.
	else if( (mod(nn,2) .ne. 0) .and. (mod(zz,2) .ne. 0) ) then
	delta = -1.
	else 
	delta = 0.
	end if 

BindEn = C(1)*(ZZ + NN) + (C(2)*(ZZ + NN)**(2./3.)) + (C(3)*((NN - ZZ)**2)/(FLOAT((ZZ + NN)))&
+C(4)*ZZ*(ZZ - 1.)/((ZZ + NN)**(1./3.))) + delta*(C(5)*(ZZ + NN)**(-3./4.))

Print*," "
Print*," "
Print*,"Binding Energy: ", BindEn, "KeV"


Print*,"Uncertainty in Binding Energy "
Print*," "
Print*," "

BindenUncertainty = 0.0

Do i = 1, K
	BindenUncertainty = BindenUncertainty +((BindEn/c(i)))**2 *(UncertainC(i))**2
enddo
BindenUncertainty = sqrt(BindenUncertainty) ! we do the square root to get our uncertainty of our BE

Print*,"Uncertainty in BE:", BindenUncertainty, "KeV" 


! Advanced Project: Part I !
! Drip Lines !


! - - - Neutron Drip Line - - - !

! NOTE!!!!!!!!!! MY PROTON DRIP LINE SHOWS A LINEAR FIT WITHOUT LOG SCALE
! HOWEVER MY NEUTRON DRIP LINE DOES SOMETHING WEIRD, I could not see the error, it drops off suddenly



Do ZZ = 2, 160


!The reason why I put NN = ZZ is because I want to start in the positive zone of the
! binding energy 
	
	NN = ZZ 

619 continue 

if( ( mod(nn,2).eq.0 ) .and. ( mod(zz,2) .eq. 0) ) then
	delta = 1.
	else if( (mod(nn,2) .ne. 0) .and. (mod(zz,2) .ne. 0) ) then
	delta = -1.
	else 
	delta = 0.
	end if 



! the binding energy we want via inputs 	
BindEn = C(1)*(ZZ + NN) + (C(2)*(ZZ + NN)**(2./3.)) + (C(3)*((NN - ZZ)**2)/(FLOAT((ZZ + NN)))&
+C(4)*ZZ*(ZZ - 1.)/((ZZ + NN)**(1./3.))) + delta*(C(5)*(ZZ + NN)**(-3./4.))

! the minimum energy for this binding energy !

BindEnMin = C(1)*(ZZ + NN-1)+(C(2)*(ZZ + NN-1)**(2./3.))+(C(3)*(((NN-1) - ZZ)**2)/(FLOAT((ZZ + NN-1)))&
+C(4)*ZZ*(ZZ - 1.)/((ZZ + NN-1)**(1./3.))) + delta*(C(5)*(ZZ + NN-1)**(-3./4.))


 If( Binden - BindEnMin >= 0 ) then 
 	
 	NN = NN +1 
 	
 	goto 619
 endif 
 
 Write(1000,*) NN-1, ZZ
 
Enddo

Close(1000)

! - - - Proton Drip Line - - - !

Do NN = 2, 160

ZZ = NN

310 Continue 


if( ( mod(nn,2).eq.0 ) .and. ( mod(zz,2) .eq. 0) ) then
	delta = 1.
	else if( (mod(nn,2) .ne. 0) .and. (mod(zz,2) .ne. 0) ) then
	delta = -1.
	else 
	delta = 0.
	end if 



! the original binding energy !
BindEn = C(1)*(ZZ + NN) + (C(2)*(ZZ + NN)**(2./3.)) + (C(3)*((NN - ZZ)**2)/(FLOAT((ZZ + NN)))&
+C(4)*ZZ*(ZZ - 1.)/((ZZ + NN)**(1./3.))) + delta*(C(5)*(ZZ + NN)**(-3./4.))


! the minimum energy for the proton drip line !

BindEnMin = C(1)*(ZZ-1.+NN)+(C(2)*(ZZ-1.+NN)**(2./3.))+(C(3)*((NN-ZZ-1.)**2)/(FLOAT((ZZ-1 + NN)))&
+C(4)*(ZZ-1)*(ZZ-1 - 1)/((ZZ-1. + NN)**(1./3.))) + delta*(C(5)*(ZZ-1 + NN)**(-3./4.))

 If( BindEn - BindEnMin >= 0 ) then
 
 	ZZ = ZZ + 1 
 	 
Goto 310
 
 endif
 
 
 
 Write(2000,*) NN, ZZ-1
 

enddo	

Close(2000)
	




End program Nymeria


 
 !this subroutine is used to solve for our A(ALPHA,BETA) matrix and our coefficients
 ! in this subroutine i calculate my a matrix and other things that will help with the code 
 Subroutine Getmatrix(Tempproduct,BE,BEerror,G,Z,N,M,K,Amat)
 Implicit none
 Integer :: i,j,K,M,h
 Integer :: Z(1:M)
 Integer :: N(1:M)
 Real :: G(K,M)
 Real :: BEerror(1:M)
 Real :: delta
 Real :: Amat(1:K,1:K)
 Real :: ATEMP(K,K)
 Real :: Tempproduct(1:K)
 Real :: BE(1:M)
 
!initializing my summations 
 Amat = 0.0 
 Tempproduct = 0.0

 Do i = 1, M	

 
 ! this if statement is for the kronecker delta 
 	if( ( mod(n(i),2).eq.0 ) .and. ( mod(z(i),2) .eq. 0) ) then
	delta = 1.
	else if( (mod(n(i),2) .ne. 0) .and. (mod(z(i),2) .ne. 0) ) then
	delta = -1.
	else 
	delta = 0.
	end if 
 
!this is the g functions that are analogous to the A number 
	G(1,i) = float((Z(i) + N(i)))

	G(2,i) = float((Z(i) + N(i)))**(2./3.)

	G(3,i) = float(((N(i) - Z(i)))**2)/float((Z(i) + N(i)))

	G(4,i) = (float(Z(i))*(float((Z(i) - 1))))/(float(Z(i) + N(i)))**(1./3.)

	G(5,i) = (float((Z(i))) + float(N(i)))**(-3./4.)*delta
 
	
	
!this is for our A matrix	
	
	Do j = 1, K
	 Do h = 1, K
		Amat(j,h) = Amat(j,h) + G(h,i)*(1./(BEerror(i)**2))*G(j,i)
		 		
	 		
	 enddo
	
	
	! tempproduct is the product between my binding energy and my G matrix
	
		tempproduct(j) = tempproduct(j) + G(j,i)*(1./(BEerror(i)**2))*BE(i)
	
	enddo
 
 
 Enddo

 
 Return
 
 End Subroutine Getmatrix
 



	
	
	
	
	 
	 
	 

	