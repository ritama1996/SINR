PROGRAM histogram
IMPLICIT NONE
REAL(8) :: qmax, qmin,x,y,z,pint,xd,dq, pan, sum  
REAL(8),ALLOCATABLE, DIMENSION(:) :: prob_x,prob_y,prob_z
INTEGER :: nbins,i,ia,ios,ndata ,idata,ibin 
REAL(8), PARAMETER :: pi=3.14d0, g=0.1d0, omega=3.d0 
CHARACTER (LEN=50) :: trajectory_file

nbins=80
qmax=2.d0
qmin=-2.d0
dq=(qmax-qmin)/DFLOAT(nbins)

ALLOCATE(prob_x(nbins))
ALLOCATE(prob_y(nbins))
ALLOCATE(prob_z(nbins))

DO ibin=1,nbins
    prob_x(ibin)=0.d0
    prob_y(ibin)=0.d0
    prob_z(ibin)=0.d0
ENDDO

OPEN(UNIT=1,FILE="Trajectory.in",STATUS="old")
OPEN(2,FILE="Probability_x.out")

ndata=0
read_loop : DO 
    READ(1,*,IOSTAT=ios)
    IF(ios/=0)EXIT read_loop
    ndata=ndata+1
ENDDO read_loop

REWIND(1)

DO idata=1,ndata
    READ(1,*)x,y,z 
    ibin=INT((x-qmin)/dq)+1
    IF(ibin<=1)ibin=1
    IF(ibin>nbins)ibin=nbins 
    prob_x(ibin)=prob_x(ibin)+1.d0
    
    IF(mod(idata,ndata)==0)THEN
        
        pint=0.d0
    
        DO ibin=1,nbins
            pint=pint+prob_x(ibin)*dq
        ENDDO
        DO ibin=1,nbins
            xd=qmin+dfloat(ibin-1)*dq+0.5d0*dq
            pan=(1.d0/0.83477d0)*DEXP(-(0.5d0*(omega*xd)**2+0.25d0*g*xd**4))
            WRITE(2,*)xd,prob_x(ibin)/pint,pan  
        ENDDO
    ENDIF

ENDDO
CLOSE(1)
CLOSE(2)


END PROGRAM histogram
