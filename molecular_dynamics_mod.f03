MODULE molecular_dynamics
    use set_parameter
    use set_input
    use get_force
    use update_variables
    use mts_method

    IMPLICIT NONE
    INTEGER :: i, j, ia, io, ji 
    REAL(8) :: t,KE,vsq,q_2,q_4,PE

    CONTAINS

    SUBROUTINE md_loop
    IMPLICIT NONE

    CALL RANDOM_SEED()

    fast_force=.false.
    CALL force_calculation(q,force)
    
    io = 1
    t = 0.d0 !Time
    outer_loop : DO
         
        IF(use_mts%XO_RESPA)CALL nose_update(w,del_t,M,v1,v2,v)

        ji = 1

        inner_loop : DO

            IF(use_mts%XI_RESPA)CALL nose_update(w,dt,M,v1,v2,v)

            CALL velocity_update(dt,M,lambda,force,v,v1)

            CALL pos_ou_update(dt,M,v,q,v2)
        
            fast_force=.true.
            IF(ji==Nrespa1)fast_force=.false.

            CALL force_calculation(q,force)

            CALL velocity_update(dt,M,lambda,force,v,v1)

            IF(use_mts%XI_RESPA)CALL nose_update(w,dt,M,v1,v2,v)

            ji = ji + 1

            IF(ji>NINT(del_t/dt))EXIT inner_loop

        ENDDO inner_loop

        IF(use_mts%XO_RESPA)CALL nose_update(w,del_t,M,v1,v2,v)
        KE=0.d0
        PE=0.d0
        DO ia=1,N_atom
          vsq=DOT_PRODUCT(v(1:2,ia),v(1:2,ia))  !v=sqrt(v_x^2 + v_y^2 + v_z^2)
          KE=KE+0.5d0 * M * DSQRT(vsq)
         !q_2=DOT_PRODUCT(q(1:2,ia),q(1:2,ia))
         !q_4=SUM(q(ia,1:3)**4)
          PE=PE+4.d0*(q(1,1)**2-1.d0)**2+4.d0*(q(2,1)**2-1)
        ENDDO
        
        WRITE(11,*)t,q(1:2,1),v(1:2,1),KE,PE


   IF(use_mts%XO_RESPA)THEN
     WRITE(3,*)"No. of steps =",io
                
     DO ia=1,N_atom
       WRITE(3,"(I5,3F16.10)")ia,q(1,ia),q(2,ia)!(3,ia)
       IF(ia==1)WRITE(20,*)q(1,ia),q(2,ia)!q(3,ia)
     ENDDO

   ELSEIF(use_mts%XI_RESPA)THEN
     WRITE(4,*)"No. of steps =",io
     DO ia=1,N_atom
       WRITE(4,"(I5,3F16.10)")ia,q(1,ia),q(2,ia)!q(3,ia)
       IF(ia==1)WRITE(30,*)q(1,ia),q(2,ia)!q(3,ia)
     ENDDO

   ENDIF
        
        

        DO i=1,3
            DO ia=1,N_atom
                cons(i,ia) = M * v(i,ia) * v(i,ia)
                DO j=1,L 
                    cons(i,ia)=cons(i,ia)+lbylp1 * Q1 * v1(i,ia,j)**2.d0
                ENDDO
            ENDDO
        ENDDO

        IF(use_mts%XO_RESPA)THEN
            IF(mod(io,5000)==0)THEN
                !WRITE(5,"(6X,I10)")io 
                !DO ia=1,N_atom
                    WRITE(5,*)io,cons(1,1),cons(2,1),cons(3,1)
                !ENDDO
            ENDIF
        ENDIF 

        IF(use_mts%XI_RESPA)THEN
            IF(mod(io,5000)==0)THEN
                WRITE(7,"(6X,I10)")io 
                DO ia=1,N_atom
                    WRITE(7,*)ia,cons(1,ia),cons(2,ia),cons(3,ia) !To write conserved quantity for individual atoms in individual file
                ENDDO
            ENDIF
        ENDIF 

        io = io + 1
        IF(io>N_MD)EXIT outer_loop
     t = t+del_t  !Increment of time by del_t
    ENDDO outer_loop


    END SUBROUTINE md_loop

END MODULE molecular_dynamics
