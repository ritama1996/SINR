MODULE mts_method

IMPLICIT NONE

TYPE method 

LOGICAL :: XI_RESPA, XO_RESPA 

END TYPE method 
 TYPE(method) :: use_mts

END MODULE mts_method