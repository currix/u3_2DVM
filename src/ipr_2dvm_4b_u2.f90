PROGRAM ipr_4b_2DVM_u2
  !
  ! Program to compute Energy and Participation Ratio (PR) of the
  ! U(3) 2DVM 4-Body Hamiltonian in the U(2) Dynamical Symmetry Basis
  !
  ! by Currix TM.
  !
  !
  USE nrtype
  !
  USE defparam_2DVM
  !
  USE u3_2dvm_mod
  !
  ! Lapack 95
#ifdef  __GFORTRAN__
  ! gfortran
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR
#ifdef _OPENMP
  USE OMP_LIB
#endif
#else
  !ifort
  USE F95_PRECISION, ONLY: WP => DP
  USE LAPACK95, ONLY: SYEVR
#endif 
  !
  !
  IMPLICIT NONE
  !
  !
  INTEGER(KIND = I4B) :: state_index, state_index_2
  !
  ! Hamiltonian Parameters
  REAL(KIND = DP) :: P11
  REAL(KIND = DP) :: P21, P22, P23
  REAL(KIND = DP) :: P31, P32, P33
  REAL(KIND = DP) :: P41, P42, P43, P44, P45, P46, P47
  !
  INTEGER(KIND = I4B) :: L, L_min, L_max
  !
  INTEGER(KIND = I4B) :: out_unit
  !
  CHARACTER(LEN=65) :: output_file, input_file_name
  !
#ifdef _OPENMP
  INTEGER(KIND = I4B) :: threads
#endif
  !
  !     NAMELIST DEFINITIONS
  !
  NAMELIST/par_aux/  IPRINT, Eigenvec_Log, Excitation_Log, Save_avec_Log
  NAMELIST/par_0/  N_val, L_min, L_max
  NAMELIST/INP1b/ P11
  NAMELIST/INP2b/ P21, P22, P23
  NAMELIST/INP3b/ P31, P32, P33
  NAMELIST/INP4b/ P41, P42, P43, P44, P45, P46, P47
  !
  ! NAMELIST FILE NAME FROM STDIN
  READ(5,*) input_file_name
  !
  !
  !
  !     OPEN INPUT FILE
  OPEN(UNIT=10,FILE=TRIM(input_file_name),STATUS='OLD')
  !
  !     READING INPUT
  !
  READ(10,par_aux)
  READ(10,par_0)
  READ(10,INP1b)
  READ(10,INP2b)
  READ(10,INP3b)
  READ(10,INP4b)
  !
  CLOSE(10)
  !
  ! HAMILTONIAN PARAMETER VECTOR ALLOCATION AND INITIALIZATION
  ALLOCATE(H_4b_pars(1:NPMAX), STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR ALLOCATING H_4b_pars MATRIX'
  H_4b_pars = 0.0_DP
  !
  H_4b_pars(1) = P11
  H_4b_pars(2) = P21
  H_4b_pars(3) = P22
  H_4b_pars(4) = P23
  H_4b_pars(5) = P31
  H_4b_pars(6) = P32
  H_4b_pars(7) = P33
  H_4b_pars(8) = P41
  H_4b_pars(9) = P42
  H_4b_pars(10) = P43
  H_4b_pars(11) = P44
  H_4b_pars(12) = P45
  H_4b_pars(13) = P46
  H_4b_pars(14) = P47
  !
  !
  !
  ! TESTS
  IF ((N_val < L_min) .or. (N_val < L_max)) &
       STOP 'ERROR :: N_VAL < L_VAL. SAYONARA BABY'
  !
  !
  IF (Excitation_Log) THEN
     !
     ! Compute L = 0 Ground state
     !
     !
     ! Initialize time
     CALL CPU_TIME(time_check_ref)
     !
     ! L_VAL = 0 BLOCK DIMENSION
     dim_block = DIM_L_BLOCK(N_val, 0)
     !
     ! Build U(3) > U(2) > SO(2) Basis
     ! ALLOCATE BASIS
     ALLOCATE(U2_Basis(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U2_Basis allocation request denied."
        STOP
     ENDIF
     !
     CALL U2_BASIS_VIBRON(N_val, 0, U2_Basis) ! Build U2 basis
     !
     !
     ! Allocate and Build Hamiltonian and Hamiltonian Operator Matrices
     !
     CALL Build_U2DS_Operator_Matrices(N_val, 0, dim_block, U2_Basis)
     !
     !
     CALL Build_Ham_4Body_U2(N_val, 0, dim_block, U2_Basis)  
     !
     ! Hamiltonian Diagonalization
     !
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = 0 ", &
          " Time (1) :: L = 0 Model Hamiltonian Building ", time_check - time_check_ref
     time_check_ref = time_check
     !
     Eigenval_vector = 0.0_DP
     !      
     ! Diagonalize Hamiltonian matrix (LAPACK95)
#ifdef  __GFORTRAN__
     !gfortran
     CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='N', UPLO='U', IL=1, IU=1)
#else
     !ifort
     CALL SYEVR(A=Ham_matrix, W=Eigenval_vector, UPLO='U', IL=1, IU=1)
#endif 
     !
     GS_energy = Eigenval_vector(1)
     !
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
          " Time (4) :: U3 Model Hamiltonian diagonalization ", time_check - time_check_ref
     time_check_ref = time_check
     !    
  ENDIF
  !
#ifdef _OPENMP
  threads = omp_get_max_threads()
  !
  write(*,6) threads, abs(L_max-L_min)+1
  !
6 format("Available threads: ", I3, "; Loops iterations: ", I4)
  !
  IF (threads .gt. abs(L_max-L_min)+1) THEN
     !
     threads = abs(L_max-L_min)+1
     !
  ENDIF
  !
  ! Threads needed: auto-set
  !
  !$OMP PARALLEL DO DEFAULT(private), SHARED(N_val,GS_energy, P11, P21, P22, &
  !$OMP P23, P31, P32, P33, P41, P42, P43, P44, P45, P46, P47, &
  !$OMP Save_avec_Log, Excitation_Log, Eigenvec_Log, &
  !$OMP Iprint, L_max, L_min), num_threads(threads)
  !
  !
#endif
#ifndef __GFORTRAN__
  call mkl_set_dynamic(0)
#endif
  !
  DO L = L_min, L_max
     !
     out_unit = L + 100
     !
     L_val = L
     !
     ! Build filename
     IF (L_val < 10) THEN
        IF ( N_val < 10) THEN !to avoid spaces
           WRITE(output_file, '("autoval_4b_u2_N",I1,"_L",I1,".dat")')  N_val, L_val
        ELSE IF ( N_val < 100) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I2,"_L",I1,".dat")')  N_val, L_val
        ELSE IF ( N_val < 1000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I3,"_L",I1,".dat")')  N_val, L_val
        ELSE IF ( N_val < 10000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I4,"_L",I1,".dat")')  N_val, L_val
        ELSE
           WRITE(output_file, '("autoval_4b_u2_N",I6,"_L",I1,".dat")')  N_val, L_val
        ENDIF
     ELSE IF (L_val >= 10 .AND. L_val<100) THEN
        IF ( N_val < 10) THEN !to avoid spaces
           WRITE(output_file, '("autoval_4b_u2_N",I1,"_L",I2,".dat")')  N_val, L_val
        ELSE IF ( N_val < 100) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I2,"_L",I2,".dat")')  N_val, L_val
        ELSE IF ( N_val < 1000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I3,"_L",I2,".dat")')  N_val, L_val
        ELSE IF ( N_val < 10000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I4,"_L",I2,".dat")')  N_val, L_val
        ELSE
           WRITE(output_file, '("autoval_4b_u2_N",I6,"_L",I2,".dat")')  N_val, L_val
        ENDIF
     ELSE
        IF ( N_val < 10) THEN !to avoid spaces
           WRITE(output_file, '("autoval_4b_u2_N",I1,"_L",I3,".dat")')  N_val, L_val
        ELSE IF ( N_val < 100) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I2,"_L",I3,".dat")')  N_val, L_val
        ELSE IF ( N_val < 1000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I3,"_L",I3,".dat")')  N_val, L_val
        ELSE IF ( N_val < 10000) THEN 
           WRITE(output_file, '("autoval_4b_u2_N",I4,"_L",I3,".dat")')  N_val, L_val
        ELSE
           WRITE(output_file, '("autoval_4b_u2_N",I6,"_L",I3,".dat")')  N_val, L_val
        ENDIF
     ENDIF
     !
     OPEN(UNIT=out_unit, FILE=output_file, STATUS='unknown', ACTION='write') 
     !
     !
     ! L_VAL BLOCK DIMENSION
     dim_block = DIM_L_BLOCK(N_val, l_val)
     !
     ! Build U(3) > U(2) > SO(2) Basis
     ! ALLOCATE BASIS
     IF (ALLOCATED(U2_Basis)) THEN    
        DEALLOCATE(U2_Basis, STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "U2_Basis deallocation request denied."
           STOP
        ENDIF
     ENDIF
     !
     ALLOCATE(U2_Basis(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U2_Basis allocation request denied."
        STOP
     ENDIF
     !
     CALL U2_BASIS_VIBRON(N_val, L_val, U2_Basis) ! Build U2 basis
     !
     ! Allocate and Build Hamiltonian and Hamiltonian Operator Matrices
     !
     CALL Build_U2DS_Operator_Matrices(N_val, L_val, dim_block, U2_Basis)
     !
     !
     CALL Build_Ham_4Body_U2(N_val, L_val, dim_block, U2_Basis)  
     !
     IF (Save_avec_Log) THEN
        !
        ALLOCATE(Diagonal_vector(1:dim_block), STAT = IERR)    
        IF (IERR /= 0) STOP 'Diagonal_vector allocation request denied.'
        !
        forall (state_index=1:dim_block) Diagonal_vector(state_index) = Ham_matrix(state_index, state_index)
        !
     ENDIF
     !
     IF (Iprint > 3) THEN
        WRITE(UNIT = out_unit, FMT = *) ' '
        WRITE(UNIT = out_unit, FMT = *) ' Hamiltonian Matrix'
        DO state_index = 1, dim_block                         
           WRITE(UNIT = out_unit, FMT = *) (Ham_matrix(state_index, state_index_2), state_index_2 = 1, dim_block)
        ENDDO
     ENDIF
     !
     ! Hamiltonian Diagonalization
     !
     !
#ifndef _OPENMP      
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
          " Time (3) :: U3 Hamiltonian Building ", time_check - time_check_ref
     time_check_ref = time_check
#endif
     !
     Eigenval_vector = 0.0_DP
     !      
     ! Diagonalize Hamiltonian matrix (LAPACK95)
#ifdef  __GFORTRAN__
     !gfortran
     IF (Eigenvec_Log .OR. Save_avec_Log) THEN
        CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='V', UPLO='U')
     ELSE
        CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='N', UPLO='U')
     ENDIF
#else
     !ifort
     IF (Eigenvec_Log .OR. Save_avec_Log) THEN
        !
        ALLOCATE(Eigenvec_array(1:dim_block, 1:dim_block), STAT = IERR)    
        IF (IERR /= 0) THEN
           WRITE(UNIT = *, FMT = *) "Eigenvec_array allocation request denied."
           STOP
        ENDIF
        !
        CALL SYEVR(A=Ham_matrix, W=Eigenval_vector, UPLO='U', Z = Eigenvec_array)
     ELSE
        CALL SYEVR(A=Ham_matrix, W=Eigenval_vector, UPLO='U')
     ENDIF
     Ham_matrix = Eigenvec_array
     !
     DEALLOCATE(Eigenvec_array, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Eigenvec_array deallocation request denied."
        STOP
     ENDIF
#endif  
     !
     IF (Excitation_Log) THEN
        !
        eigenval_vector = eigenval_vector - GS_energy
        !
     ENDIF
     !
#ifndef _OPENMP
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
          " Time (4) :: U3 Model Hamiltonian diagonalization ", time_check - time_check_ref
     time_check_ref = time_check
#endif
     !
     !
     ! CALCULATE IPR
     !
     IF (Iprint > 0) WRITE(UNIT = out_unit, FMT = *) "L_val = ", L_val
     !
     IF (Eigenvec_Log) THEN 
        !
        DO state_index = 1, dim_block
           !
           WRITE(UNIT = out_unit, FMT = *) state_index, &
                Eigenval_vector(state_index), &
                Inv_Part_Ratio(Ham_matrix(1:dim_block, state_index))
           !
        ENDDO
        ! 
     ELSE
        !
        DO state_index = 1, dim_block
           !
           WRITE(UNIT = out_unit, FMT = *) state_index, &
                Eigenval_vector(state_index)
           !
        ENDDO
        !
     ENDIF
     !
     ! Save eigenvector components
     IF (Save_avec_Log) THEN
        !
        IF (Excitation_Log) THEN
           !
           CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
                Eigenval_vector, Diagonal_vector - GS_energy, "u2", Ham_matrix, &
                unit=L+20)
        ELSE
           CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
                Eigenval_vector, Diagonal_vector, "u2", Ham_matrix, &
                unit=L+20)
        ENDIF
        !
        DEALLOCATE(Diagonal_vector, STAT = IERR)    
        IF (IERR /= 0) STOP 'Diagonal_vector deallocation request denied.'
        !
     ENDIF
     !
     CLOSE(out_unit,IOSTAT=IERR)
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Error closing the unit of eigenvalues."
        STOP
     ENDIF
     !
  ENDDO
  !
  !$OMP END PARALLEL DO
  !
END PROGRAM ipr_4b_2DVM_u2
