MODULE get_force
    use set_parameter
    use set_input

!!......This is to calculate force......!!
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: force
LOGICAL :: fast_force

CONTAINS
SUBROUTINE force_calculation(q,force)

IMPLICIT NONE
REAL(8),DIMENSION(2,N_atom), INTENT(in) :: q
REAL(8),DIMENSION(2,N_atom), INTENT(out) :: force

IF(fast_force)THEN

    force=-omega**2.d0*q
    
    ELSE
    !force=-omega**2.d0*q-DFLOAT(Nrespa1)*g_f*q**3.d0
    force(1,1)=-16.d0*q(1,1)*(q(1,1)**2-1)
    force(2,1)=-8.d0*q(2,1)
    
ENDIF

END SUBROUTINE force_calculation
    
END MODULE get_force
