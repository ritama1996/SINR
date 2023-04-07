MODULE update_variables
    use set_parameter
    use set_input
    use get_force
!!.............This is to to update variables................!!
    IMPLICIT NONE
    
    CONTAINS
!!*************Update nose variables and velocity****************!!
!!*************Time step for this step will be according to RESPA used**************!!
    SUBROUTINE nose_update(w,tx,M,v1,v2,v)

!!w = weight of Suzuki-Yoshida order of integration........!!
!!tx = Time step ; del_t for XO-RESPA & dt for XI-RESPA....!!

    IMPLICIT NONE
        REAL(8), INTENT(in) :: w(N_sy)
        REAL(8), INTENT(IN) :: tx, M
        REAL(8), INTENT(inout) :: v1(2,N_atom,L), v2(2,N_atom,L), v(2,N_atom)
        REAL(8), DIMENSION(2,N_atom) :: H, sum_H
        REAL(8) :: dti, G(2,N_atom,L)
        INTEGER :: irespa2, isy, i, j, ia
!!#################################################################
    DO irespa2=1,Nrespa2                                          !
                                                                  !
        DO isy=1,N_sy                                             !
                                                                  !
        dti=w(isy)*tx                                             !
        !! Upto this part in one subroutine#######################!
            DO i=1,2
                DO ia=1,N_atom

                    DO j=1,L
                        G(i,ia,j)=Q1*v1(i,ia,j)**2.d0 - K_B_T
                        v2(i,ia,j)=v2(i,ia,j)+(G(i,ia,j)/Q2) * (dti/4.d0)
                    ENDDO

                CALL summation(i,ia,dti,v1,v2,sum_H)
                H(i,ia)=DSQRT(lambda(i,ia)/(M*v(i,ia)*v(i,ia)+lbylp1*sum_H(i,ia)))
                v(i,ia)=v(i,ia)*H(i,ia)

                    DO j=1,L
                        v1(i,ia,j)=v1(i,ia,j)*H(i,ia)*DEXP(-(v2(i,ia,j)*dti/2.d0))
                        G(i,ia,j)=Q1*v1(i,ia,j)**2.d0 - K_B_T
                        v2(i,ia,j)=v2(i,ia,j)+(G(i,ia,j)/Q2) * (dti/4.d0)
                    ENDDO
                ENDDO
            ENDDO
        
        ENDDO
    ENDDO
!    WRITE(*,*)w(3)*tx
    END SUBROUTINE nose_update

    SUBROUTINE summation(i,ia,dti,v1,v2,sum_H)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: v1(2,N_atom,L), v2(2,N_atom,L),dti
    REAL(8), INTENT(out) :: sum_H(2,N_atom)
    INTEGER, INTENT(in) :: i,ia
    INTEGER :: j
    

    sum_H(i,ia) = 0.d0
    DO j=1,L
        sum_H(i,ia)=sum_H(i,ia)+Q1*v1(i,ia,j)**2.d0*DEXP(-(v2(i,ia,j)*dti))
    ENDDO

    END SUBROUTINE summation
!!....................................................................................!!
    SUBROUTINE velocity_update(ty,M,lambda,force,v,v1)

    IMPLICIT NONE
    
    REAL(8), INTENT(IN) :: ty, M, force(2,N_atom),lambda(2,N_atom)
    REAL(8), INTENT(INOUT) :: v(2,N_atom), v1(2,N_atom,L)
    REAL(8) :: a,b,root_b,s,sdot,arg
    INTEGER :: i,j,ia

    DO i=1,2
        DO ia = 1,N_atom

            a=(force(i,ia)*v(i,ia))/lambda(i,ia)
            b=force(i,ia)**2.d0/M

            root_b=DSQRT(b)
            arg=ty*root_b/2.d0

            IF (arg<(10.D0**(-5)))THEN
                s=(1.d0/root_b) * sinh_limit(arg)+(a/b) * (cosh_limit(arg)-1.d0)
                sdot=cosh_limit(arg)+(a/root_b) * sinh_limit(arg) 
            ELSE
                s=(1.d0/root_b) * sinh(arg) + (a/b) * (cosh(arg)-1.d0)
                sdot=cosh(arg) + (a/root_b) * sinh(arg)
            ENDIF

            v(i,ia) = (v(i,ia)+(force(i,ia)/M)*s)/sdot
    
            DO j=1,L
                v1(i,ia,j)=v1(i,ia,j)/sdot
            ENDDO

        ENDDO
    ENDDO

    END SUBROUTINE velocity_update 
!!...................................................................!!
    FUNCTION sinh_limit(x)result(sinh_l)

    REAL(8) :: x
    REAL(8) :: sinh_l

    sinh_l=x+(x**3.D0)/6.D0

    RETURN
    END FUNCTION sinh_limit
!!...................................................................!!
    FUNCTION cosh_limit(x)result(cosh_l)

    REAL(8) :: x
    REAL(8) :: cosh_l

    cosh_l=1.d0+(x**2.D0)/2.D0+(x**4.d0)/24.d0
    RETURN
    END FUNCTION cosh_limit

!!*******************************************************************************************************!!
    SUBROUTINE pos_ou_update(ty,M,v,q,v2)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: ty, v(2,N_atom), M
    REAL(8), INTENT(INOUT) :: q(2,N_atom), v2(2,N_atom,L)
    REAL(8) :: RND, sigma, GAUSS_DIST1(1)
    INTEGER :: i, ia, j

    sigma = DSQRT((2.d0 * gamma * K_B_T)/Q2)
    CALL box_mullar(M,RND)

    DO i=1,2
        DO ia=1,N_atom
            q(i,ia) = q(i,ia) + v(i,ia) * (ty/2.d0)
            DO j=1,L 
                v2(i,ia,j)=v2(i,ia,j)*DEXP(-(gamma*ty)) + sigma*RND*DSQRT((1.d0-DEXP(-2.d0*gamma*ty))/2.d0*gamma)
            ENDDO
            q(i,ia) = q(i,ia) + v(i,ia) * (ty/2.d0)
        ENDDO
    ENDDO
   

    END SUBROUTINE pos_ou_update
!.........................................................................................................................!
    SUBROUTINE box_mullar(M,R)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: M 
    REAL(8), INTENT(OUT) :: R 
    REAL(8) :: ran1, ran2

    CALL RANDOM_NUMBER(ran1)
    CALL RANDOM_NUMBER(ran2)

    ran1=1.d0-ran1 

    R =  DSQRT(-(2.D0*DLOG(ran1)))*DCOS(2.D0*M*pi*ran2)

    END SUBROUTINE box_mullar

SUBROUTINE GAUSS_DIST(MU,SIGMA_R,DIM,GAUSS_DIST1)
!
IMPLICIT NONE
!
REAL*8, INTENT(IN)    :: mu
REAL*8, INTENT(IN)    :: sigma_r
INTEGER,  INTENT(IN)        :: dim
REAL*8                :: gauss_dist1( dim )
!local variables
REAL*8                :: x1, x2, w, y1 , y2
INTEGER  :: i
!
!
DO i = 1, dim,2
  !
  gaussian_loop: DO
     !
     !$omp critical
     call random_number(y1)
     call random_number(y2)
     !$omp end critical
     !
     x1 = 2.0D0 * y1 - 1.0D0
     x2 = 2.0D0 * y2 - 1.0D0
     !
     w = x1 * x1 + x2 * x2
     !
     IF ( w < 1.0D0 ) EXIT gaussian_loop
    !
  END DO gaussian_loop
w = dSQRT( ( - 2.0D0 * dLOG( w ) ) / w )
 gauss_dist1(i) = x1 * w * sigma_r
  !
  IF ( i >= dim ) EXIT
  !
  gauss_dist1(i+1) = x2 * w * sigma_r
  !
END DO
!
gauss_dist1(:) = gauss_dist1(:) + mu
!
RETURN
!
END subroutine gauss_dist



END MODULE update_variables
