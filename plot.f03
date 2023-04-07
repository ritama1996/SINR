PROGRAM file_plot

  IMPLICIT NONE
  REAL(8) :: del_t, lg
  INTEGER :: i,n,ios

  OPEN(1,FILE="plot.dat")
  OPEN(2,FILE="plot_lg_lt.dat")

  read_file : DO
        READ(1,*,IOSTAT=ios)
        IF(ios/=0)EXIT read_file
        n=n+1
  ENDDO read_file

  REWIND(1)

  DO i=1,n
    READ(1,*)del_t,lg
    WRITE(2,*)dlog10(del_t),lg
  ENDDO

  CLOSE(1)
  CLOSE(2)

END PROGRAM file_plot
