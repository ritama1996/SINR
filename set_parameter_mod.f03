MODULE set_parameter
    use read_file

!It gives the value of parameters used in SINR thermostat !!

IMPLICIT NONE

REAL*8 :: M, K_B_T, pi, omega, g_f
INTEGER :: Nrespa1
REAL*8, ALLOCATABLE, DIMENSION(:) :: w


!!..........................................................!!

CONTAINS

!!.........................................................!!
SUBROUTINE real_parameters

IMPLICIT NONE
REAL*8 :: w1,w2,w3,w4

pi=4.d0*ATAN(1.d0)
omega = 3.d0
g_f=0.1d0
M=1.d0
K_B_T=1.d0
del_t = del_t * pi/omega
Nrespa1 = INT(del_t/dt)  !! Keeping dt fixed

IF(N_sy==3)THEN
ALLOCATE(w(N_sy))

    w1=1.351207
    w2=-1.7024
    w3=1.351207

    w(1)=w1/DFLOAT(Nrespa2)
    w(2)=w2/DFLOAT(Nrespa2)
    w(3)=w3/DFLOAT(Nrespa2)

    ELSE IF(N_sy==7)THEN
    ALLOCATE(w(N_sy))

        w1=0.7845136d0
        w2=0.2355732d0
        w3=-1.1776799d0
        w4=1.d0-2.d0*(w1+w2+w3)

        w(1)=w1/DFLOAT(Nrespa2)
        w(2)=w2/DFLOAT(Nrespa2)
        w(3)=w3/DFLOAT(Nrespa2)
        w(4)=w4/DFLOAT(Nrespa2)
        w(5)=w3/DFLOAT(Nrespa2)
        w(6)=w2/DFLOAT(Nrespa2)
        w(7)=w1/DFLOAT(Nrespa2)
        
ENDIF

END SUBROUTINE real_parameters
!!.........................................................!!
END MODULE set_parameter
