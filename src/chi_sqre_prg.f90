PROGRAM CHI2_U3_BENDING
  !
  !
  ! by Currix TM
  !
  USE nrtype
  USE defparam_2DVM
  !
  USE u3_2dvm_mod
  !
  USE FIT_2DVM
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP) :: CHI2
  REAL(KIND = DP) :: P11
  REAL(KIND = DP) :: P21, P22, P23
  REAL(KIND = DP) :: P31, P32, P33
  REAL(KIND = DP) :: P41, P42, P43, P44, P45, P46, P47
  CHARACTER(LEN = 64) :: input_file_name
  !
  !     NAMELIST DEFINITIONS (read in main program)
  !
  NAMELIST/INP0/ BENT, exp_data_file
  NAMELIST/INP1/ N_VAL, LMAX, VMAX, EMINL
  NAMELIST/INP2/ IPRINT, DIS_RES   
  NAMELIST/INP1b/ P11
  NAMELIST/INP2b/ P21, P22, P23
  NAMELIST/INP3b/ P31, P32, P33
  NAMELIST/INP4b/ P41, P42, P43, P44, P45, P46, P47
  !
  ! NAMELIST FILE NAME FROM STDIN
  READ(5,*) input_file_name
  !
  !!
  
  !     OPEN INPUT FILE
  OPEN(UNIT=10,FILE=TRIM(input_file_name),STATUS='OLD')
  !
  !     READING INPUT
  !
  READ(10,INP0)
  READ(10,INP1)
  READ(10,INP2)
  READ(10,INP1b)
  READ(10,INP2b)
  READ(10,INP3b)
  READ(10,INP4b)
  !
  CLOSE(10)
  !
  ! HAMILTONIAN PARAMETER VECTOR ALLOCATION AND INITIALIZATION
  ALLOCATE(H_pars(1:NPMAX), STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR ALLOCATING H_pars MATRIX'
  H_pars = 0.0_DP
  !
  H_pars(1) = P11
  H_pars(2) = P21
  H_pars(3) = P22
  H_pars(4) = P23
  H_pars(5) = P31
  H_pars(6) = P32
  H_pars(7) = P33
  H_pars(8) = P41
  H_pars(9) = P42
  H_pars(10) = P43
  H_pars(11) = P44
  H_pars(12) = P45
  H_pars(13) = P46
  H_pars(14) = P47
  !
  ! 
  CHI2 = CHI_SQRE(H_pars)
  !
  IF (.NOT. DIS_RES .OR. IPRINT > 0) WRITE(*,*) CHI2
  !
END PROGRAM CHI2_U3_BENDING
