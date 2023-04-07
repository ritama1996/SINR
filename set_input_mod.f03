MODULE set_input
    use read_file
    use set_parameter

!!This is to initialize the variables!!

IMPLICIT NONE

REAL(8), ALLOCATABLE :: q(:,:),v(:,:),lambda(:,:),v1(:,:,:),v2(:,:,:),cons(:,:)
REAL(8) :: lbylp1
!!......................................................................!!
CONTAINS
SUBROUTINE input_allocation

IMPLICIT NONE
INTEGER :: i,j

ALLOCATE(q(2,N_atom))
ALLOCATE(v(2,N_atom))
ALLOCATE(lambda(2,N_atom))
ALLOCATE(v1(2,N_atom,L))
ALLOCATE(v2(2,N_atom,L))
ALLOCATE(cons(2,N_atom))

lbylp1=DFLOAT(L)/DFLOAT(L+1)


DO i=1,2
    q(i,:)=0.25d0
    v(i,:)=DSQRT(lbylp1)/DSQRT(M)
    lambda(i,:)=M*v(i,:)**2.d0
    DO j=1,L
        v1(i,:,j)=1.d0/DSQRT(Q1)
        v2(i,:,j)=1.d0
        lambda(i,:)=lambda(i,:)+lbylp1*Q1*v1(i,:,j)**2.d0
    ENDDO
ENDDO

END SUBROUTINE input_allocation

END MODULE set_input
