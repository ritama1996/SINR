MODULE read_file

IMPLICIT NONE
INTEGER :: N_MD,L,N_sy,Nrespa2,N_atom
REAL(8) :: gamma,Q1,Q2,del_t,dt
CHARACTER (LEN=5) :: mts  

CONTAINS
SUBROUTINE read_input


!READ(1,"(A)")
READ(10,"(/,I10)")N_MD

READ(10,"(/,I10)")L 

READ(10,"(/,I10)")N_sy

READ(10,"(/,F16.6)")dt

READ(10,"(/,I10)")N_atom

READ(10,"(/,I10)")Nrespa2

READ(10,"(/,F16.6)")Q1
READ(10,"(/,F16.6)")Q2

READ(10,"(/,F16.6)")gamma

READ(10,"(/,A5)")mts 

READ(10,"(/,F16.6)")del_t

END SUBROUTINE read_input

END MODULE read_file
