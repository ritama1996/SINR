PROGRAM sinr_algorithm
    use read_file
    use set_parameter
    use set_input
    use get_force
    use mts_method
    use molecular_dynamics

IMPLICIT NONE
REAL*8 :: t1,t2

CALL cpu_time(t1)

OPEN(10, FILE="sinr.in")

ALLOCATE(force(3,N_atom))

CALL read_input

RESPA_method : SELECT CASE(mts)

    CASE("XO")
        use_mts%XO_RESPA = .true.
    CASE("XI")
        use_mts%XI_RESPA = .true.
  
END SELECT RESPA_method


CALL real_parameters
CALL input_allocation

IF(use_mts%XO_RESPA)THEN
        
    OPEN(3,FILE="XYZ_XO.out")
    WRITE(3,"(6X,A,I1,/)")"No. of thermostat used=", L 
    WRITE(3,"(16X,A,15X,A,15X,A)")"X","Y","Z"
    
    OPEN(5,FILE="Conserved_XO.out")
    WRITE(5,"(6X,A,I1,/)")"No. of thermostat used=", L
    
    OPEN(20,FILE="Trajectory_XO_atom_1.out")

ELSEIF (use_mts%XI_RESPA)THEN
    
    OPEN(4,FILE="XYZ_XI.out")
    WRITE(4,"(6X,A,I1,/)")"No. of thermostat used=", L
    WRITE(4,"(16X,A,15X,A,15X,A)")"X","Y","Z"

    OPEN(7,FILE="Conserved_XI.out")
    WRITE(7,"(6X,A,I1,/)")"No. of thermostat used=", L

    OPEN(30,FILE="Trajectory_XI_atom_1.out")

ENDIF


CALL md_loop

CLOSE(3)
CLOSE(4)
CLOSE(5)
CLOSE(7)
CLOSE(10)
CLOSE(20)
CLOSE(30)

CALL cpu_time(t2)

WRITE(*,"(A,F16.6,A)")"Time required =",t2-t1,"s"

END PROGRAM sinr_algorithm
