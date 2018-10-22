program Minimization
  use nrtype
  use defparam_2DVM
  use u3_2dvm_mod
  use FIT_2DVM
  use needed_functions
  implicit none
  integer, parameter:: ird=5, iwe=6, isav=7 !read, write and save units
  character (len=65):: input_minuit, output_file, output_file_0
  character (len=65):: namelist_file, fixed_par

  REAL(KIND = DP) :: P11
  REAL(KIND = DP) :: P21, P22, P23
  REAL(KIND = DP) :: P31, P32, P33
  REAL(KIND = DP) :: P41, P42, P43, P44, P45, P46, P47
  !
  
  NAMELIST/INP0/ BENT, exp_data_file, output_file_0
  NAMELIST/INP1/ N_VAL, LMAX, VMAX, EMINL
  NAMELIST/INP2/ IPRINT, DIS_RES   
  NAMELIST/INP1b/ P11
  NAMELIST/INP2b/ P21, P22, P23
  NAMELIST/INP3b/ P31, P32, P33
  NAMELIST/INP4b/ P41, P42, P43, P44, P45, P46, P47
  NAMELIST/fix_par/ fixed_par
  read (*,*) namelist_file

  open(unit=11,file=trim(namelist_file),status="old",action="read")

  read(11,INP0)
  read(11,INP1)
  read(11,INP2)
  read(11,INP1b)
  read(11,INP2b)
  read(11,INP3b)
  read(11,INP4b)
  read(11,fix_par)
  close(11)

  !Now we need to build the Minuit-input file
  open(unit=111,file="minuit_input",status="replace")
 inquire(111)


 write(111,'(A9)') "SET TITLE"
 write(111,'(A21)') "'Minuit minimization'"
 write(111,'(A10)') "PARAMETERS"
 write(111,20) P11
 write(111,21) P21
 write(111,22) P22
 write(111,23) P23
 write(111,24) P31
 write(111,25) P32
 write(111,26) P33
 write(111,27) P41
 write(111,28) P42
 write(111,29) P43
 write(111,30) P44
 write(111,31) P45
 write(111,32) P46
 write(111,33) P47
 write(111,*)
 write(111,'(a)') adjustl(fixed_par)
 write(111,'(A10)') 'set stra 2'
 write(111,'(A14)') 'minimize 10000'
 write(111,'(A6)') 'call 3'
 write(111,'(A4)') 'exit'
 write(111,*)
 write(111,*)

20 format("1     'P11 ' ",D20.10,"        0.1D-02")
21 format("2     'P21 ' ",D20.10,"        0.1D-02")
22 format("3     'P22 ' ",D20.10,"        0.1D-02")
23 format("4     'P23 ' ",D20.10,"        0.1D-02")
24 format("5     'P31 ' ",D20.10,"        0.1D-02")
25 format("6     'P32 ' ",D20.10,"        0.1D-02")
26 format("7     'P33 ' ",D20.10,"        0.1D-02")
27 format("8     'P41 ' ",D20.10,"        0.1D-02")
28 format("9     'P42 ' ",D20.10,"        0.1D-02")
29 format("10    'P43 ' ",D20.10,"        0.1D-02")
30 format("11    'P44 ' ",D20.10,"        0.1D-02")
31 format("12    'P45 ' ",D20.10,"        0.1D-02")
32 format("13    'P46 ' ",D20.10,"        0.1D-02")
33 format("14    'P47 ' ",D20.10,"        0.1D-02")

 close(111)
 open(unit=ird,file="minuit_input",status='old')
 open(unit=iwe,file=trim(output_file_0),status='replace')
 !
 !Allocate H_4b_pars:
 allocate (H_4b_pars(1:npmax),stat=IERR)
 IF (IERR /= 0) STOP 'ERROR ALLOCATING H_4b_pars MATRIX'
 ! 
 !Call Minuit:
 call mintio(ird,iwe,isav)
 call minuit(FCN,PRE_CHI_SQRE)
 !
 close(ird,status='delete')
 close(iwe)

end program Minimization

