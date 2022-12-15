PROGRAM Binding

IMPLICIT NONE
REAL :: BE_exp, BE_ex, BE_ex_1, pro, Var, BE, x
INTEGER :: A, Z, N, J, max, I, K, M, step, N_e, Z_e, N_a, Z_a

REAL, ALLOCATABLE :: g(:,:)
REAL, ALLOCATABLE :: INDX(:)
INTEGER, ALLOCATABLE :: Z_array(:)
INTEGER, ALLOCATABLE :: N_array(:)
REAL, ALLOCATABLE :: BE_array(:)
REAL, ALLOCATABLE :: BE_error(:)
REAL, ALLOCATABLE :: Matrix(:,:)
REAL, ALLOCATABLE :: Matrix_inv(:,:)
REAL, ALLOCATABLE :: C_v(:)
REAL, ALLOCATABLE :: C(:)
REAL, ALLOCATABLE :: delta_C(:)
REAL, ALLOCATABLE :: delta_f(:)
CHARACTER*15 file_name

WRITE(*,*) "Enter the dimensions of the matrix."
READ(*,*) K
WRITE(*,*) K

WRITE(*,*) "Enter your number of data points."
READ(*,*) M
WRITE(*,*) M

WRITE(*,*) "Enter your number of protons."
READ(*,*) Z
WRITE(*,*) Z

WRITE(*,*) "Enter your number of neutrons."
READ(*,*) N
WRITE(*,*) N

WRITE(*, '(A)', ADVANCE = 'NO') 'Enter name of the data file: '
READ(*, '(A)') file_name

OPEN(unit = 100, file=file_name, status='old')

ALLOCATE (g(1:K, 1:M))
ALLOCATE (INDX(1:K))
ALLOCATE (Matrix(1:K, 1:K))
ALLOCATE (Z_array(1:M))
ALLOCATE (N_array(1:M))
ALLOCATE (BE_array(1:M))
ALLOCATE (BE_error(1:M))
ALLOCATE (Matrix_inv(1:K, 1:K))
ALLOCATE (C_v(1:K))
ALLOCATE (C(1:K))
ALLOCATE (delta_C(1:K))
ALLOCATE (delta_f(1:K))

DO I = 1, M

READ(100,*) Z_array(I), N_array(I), BE_array(I), BE_error(I)

END DO

CLOSE (unit = 100)

CALL CO(Z_array, N_array, BE_error, M, K, g, Matrix, step)

CALL MatInv(Matrix, Matrix_inv, K, K, INDX)

CALL f(Z_array, N_array, BE_array, BE_error, K, C_v, M)

C = 0.

DO I = 1, K

DO J = 1, K

C(J) = C(J) + Matrix_inv(J,I)*C_v(I)

END DO

END DO

PRINT*, "The C_vol coefficient =", C(1)
PRINT*, "The C_surf coefficient =", C(2)
PRINT*, "The C_sym coefficient =", C(3)
PRINT*, "The C_coul coefficient =", C(4)
PRINT*, "The C_pair coefficient =", C(5)
PRINT*, "The 6th coefficient = ", C(6)

DO I = 1, K

delta_C(I) = SQRT(ABS(Matrix_inv(I,I)))

END DO

PRINT*, "The uncertainty in C_vol = ", delta_C(1)
PRINT*, "The uncertainty in C_surf = ", delta_C(2)
PRINT*, "The uncertainty in C_sym =", delta_C(3)
PRINT*, "The uncertainty in C_coul =", delta_C(4)
PRINT*, "The uncertainty in C_par =", delta_C(5)

DO I = 1, K

!delta_f(I) = g(I)*(delta_C(I))**2

END DO

Z_e = 1
N_e = 1

CALL B(C, N, Z, BE_exp, step, K, M)

PRINT*, "The Binding energy given Z and N =", BE_exp


OPEN(unit = 200, file = "Neutron_drip.dat", status = "unknown")
OPEN(unit = 300, file = "Proton_drip.dat", status = "unknown")


DO WHILE (N_e <= 157 .and. Z_e <= 108)

CALL S_n(C, BE_ex, BE_ex_1, step, K, M, N_e, Z_e)

IF (BE_ex_1 < BE_ex) THEN

WRITE(200,*) N_e, Z_e
Z_e = Z_e + 1
N_e = 1

ELSE
Z_e = Z_e
N_e = N_e + 1

END IF

END DO

CLOSE(unit = 200)

Z_e = 1
N_e = 1

DO WHILE (N_e <= 157 .and. Z_e <= 108)

CALL S_p(C, BE_ex, BE_ex_1, step, K, M, N_e, Z_e)

IF (BE_ex_1 < BE_ex) THEN

WRITE(300,*) N_e, Z_e
N_e = N_e + 1
Z_e = 1

ELSE
N_e = N_e
Z_e = Z_e + 1

END IF

END DO

CLOSE(unit = 300)


OPEN(unit = 400, file='Var_N.dat', status = 'unknown')
OPEN(unit = 500, file='Var_P.dat', status = 'unknown')



DO J = 1, M

N_a = N_array(J)
Z_a = Z_array(J)

CALL B(C, N_a, Z_a, BE, step, K, M)

Var = BE - BE_array(J)

WRITE(400,*) N_a, Var

WRITE(500,*) Z_a, Var

END DO

CLOSE(unit = 400)
CLOSE(unit = 500)

x = 0.

DO I = 1, M

N_a = N_array(I)
Z_a = Z_array(I)

CALL B(C, N_a, Z_a, BE, step, K, M)


x = x + (1./BE_error(I)**2)*(BE_array(I) - BE)**2

END DO

x = (1./(M - K))*(x)**2

WRITE(*,*) "Chi = ", x

END






SUBROUTINE CO(Z_array, N_array, BE_error, M, K, g, pro, step)

IMPLICIT NONE
INTEGER :: Z_array(1:M)
INTEGER :: N_array(1:M)
REAL :: BE_error(1:M), pro(1:K, 1:K)
REAL :: Matrix(1:K, 1:K)
INTEGER :: M, K, I, J, C, step

REAL :: g(1:K, 1:M)

g = 0.
pro = 0.

DO J = 1, M

g(1,J) = FLOAT((Z_array(J) + N_array(J)))/(BE_error(J))

g(2,J) = (FLOAT(Z_array(J) + N_array(J)))**(2./3.)/(BE_error(J))


g(3,J) = (FLOAT(N_array(J) - Z_array(J))**2)/(FLOAT(N_array(J) + Z_array(J)))/(BE_error(J))


g(4,J) = (FLOAT(Z_array(J))*(Z_array(J) - 1))/((N_array(J) + Z_array(J))**(1./3.))/(BE_error(J))

IF((Z_array(J)/2)*2 == Z_array(J) .and. (N_array(J)/2)*2 == N_array(J)) THEN
step = 1
ElSE IF ((Z_array(J)/2)*2 /= Z_array(J) .and. (N_array(J)/2)*2 /= N_array(J)) THEN
step = -1
ELSE
step = 0
END IF

g(5,J) = step*(Z_array(J) + N_array(J))**(-3./4.)/(BE_error(J))

g(6,J) = ((Z_array(J) + N_array(J))**(3))*(Z_array(J) - N_array(J))

DO I = 1, K

DO C = 1, K

pro(C,I) = pro(C,I) + g(C,J)*g(I,J)

END DO

END DO

END DO

RETURN

END SUBROUTINE






SUBROUTINE f(Z_array, N_array, BE_array, BE_error, K, pro, M)

IMPLICIT NONE
INTEGER :: Z_array(1:M)
INTEGER :: N_array(1:M)
REAL :: BE_array(1:M)
REAL :: BE_error(1:M)
REAL :: C_v(1:K, 1:M), pro(1:K)
INTEGER :: M, K, I, J, C, step

DO J = 1, M

C_v(1,J) = (Z_array(J) + N_array(J))*(BE_array(J))/((BE_error(J))**2)
C_v(2,J) = ((Z_array(J) + N_array(J))**(2./3.))*BE_array(J)/((BE_error(J))**2)
C_v(3,J) = ((N_array(J) - Z_array(J))**2)/&
(FLOAT(N_array(J) + Z_array(J)))*BE_array(J)/((BE_error(J))**2)

C_v(4,J) = ((Z_array(J))*(Z_array(J) - 1))&
/((N_array(J) + Z_array(J))**(1./3.))*BE_array(J)/((BE_error(J))**2)

IF((Z_array(J)/2)*2 == Z_array(J) .and. (N_array(J)/2)*2 == N_array(J)) THEN
step = 1
ElSE IF ((Z_array(J)/2)*2 /= Z_array(J) .and. (N_array(J)/2)*2 /= N_array(J)) THEN
step = -1
ELSE
step = 0
END IF

C_v(5,J) = (step*(Z_array(J) + N_array(J))**(-3./4.))*BE_array(J)/((BE_error(J))**2)

C_v(6,J) = ((Z_array(J) + N_array(J))**(3))*(Z_array(J) - N_array(J))*BE_array(J)/((BE_error(J))**2)

END DO

pro = 0.

DO J = 1, M

DO C = 1, K

pro(C) = pro(C) + C_v(C,J)

END DO

END DO

RETURN

END SUBROUTINE




SUBROUTINE B(C, N, Z, BE_exp, step, K, M)

IMPLICIT NONE
REAL :: BE_exp
REAL :: C(1:K)
INTEGER :: N, Z, step, K, M


IF((Z/2)*2 == Z .and. (N/2)*2 == N) THEN
step = 1
ElSE IF ((Z/2)*2 /= Z .and. (N/2)*2 /= N) THEN
step = -1
ELSE
step = 0
END IF

BE_exp = C(1)*(Z + N) + (C(2)*(Z + N)**(2./3.)) + (C(3)*((N - Z)**2)/(FLOAT((Z + N)))&
+ C(4)*Z*(Z - 1.)/((Z + N)**(1./3.))) + step*(C(5)*(Z + N)**(-3./4.))


RETURN

END SUBROUTINE






SUBROUTINE S_n(C, BE_exp, BE_exp_1, step, K, M, N, Z)

IMPLICIT NONE
REAL :: BE_exp, BE_exp_1
REAL :: C(1:K)
INTEGER :: step, K, M, I, J, N, Z

IF((Z/2)*2 == Z .and. ((N+1)/2)*2 == (N+1)) THEN
step = 1
ElSE IF ((Z/2)*2 /= Z .and. ((N+1)/2)*2 /= (N+1)) THEN
step = -1
ELSE
step = 0
END IF

BE_exp_1 = C(1)*(Z + N + 1) + C(2)*(Z + N + 1)&
**(2./3.) + C(3)*(((N + 1) - Z)**2)/(FLOAT((Z + N + 1)))+ &
(C(4)*Z*(Z - 1.))/((Z + N + 1)**(1./3.))&
+ step*C(5)*((Z + N + 1)**(-3./4.))

IF((Z/2)*2 == Z .and. (N/2)*2 == N) THEN
step = 1
ElSE IF ((Z/2)*2 /= Z .and. (N/2)*2 /= N) THEN
step = -1
ELSE
step = 0
END IF

BE_exp = C(1)*(Z + N) + C(2)*(Z + N)**&
(2./3.) + C(3)*((N - Z)**2)/(FLOAT((Z + N)))&
+ C(4)*Z*(Z - 1.)/((Z + N)**(1./3.))&
+ step*C(5)*(Z + N)**(-3./4.)

RETURN

END SUBROUTINE




SUBROUTINE S_p(C, BE_exp, BE_exp_1, step, K, M, N, Z)

IMPLICIT NONE
REAL :: BE_exp, BE_exp_1
REAL :: C(1:K)
INTEGER :: step, K, M, I, J, N, Z

IF(((Z+1)/2)*2 == Z+1 .and. (N/2)*2 == N) THEN
step = 1
ElSE IF (((Z+1)/2)*2 /= Z .and. (N/2)*2 /= N) THEN
step = -1
ELSE
step = 0
END IF

BE_exp_1 = C(1)*(Z + N + 1) + C(2)*(Z + N + 1)&
**(2./3.) + C(3)*((N  - (Z+1))**2)/(FLOAT((Z + N + 1)))+ &
(C(4)*Z*(Z+1))/((Z + N + 1)**(1./3.))&
+ step*C(5)*((Z + N + 1)**(-3./4.))

IF((Z/2)*2 == Z .and. (N/2)*2 == N) THEN
step = 1
ElSE IF ((Z/2)*2 /= Z .and. (N/2)*2 /= N) THEN
step = -1
ELSE
step = 0
END IF

BE_exp = C(1)*(Z + N) + C(2)*(Z + N)**&
(2./3.) + C(3)*((N - Z)**2)/(FLOAT((Z + N)))&
+ C(4)*Z*(Z - 1.)/((Z + N)**(1./3.))&
+ step*C(5)*(Z + N)**(-3./4.)

RETURN

END SUBROUTINE












