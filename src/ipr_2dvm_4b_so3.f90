PROGRAM ipr_4b_2DVM_so3
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
  CHARACTER(LEN = 64) :: input_file_name
  !
  !
  !     NAMELIST DEFINITIONS
  !
  NAMELIST/INP1/ N_VAL, L_VAL, IPRINT, Eigenvec_Log, Excitation_Log, Save_avec_Log
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
  READ(10,INP1)
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
  IF (N_val < L_val) STOP 'ERROR :: N_VAL < L_VAL. SAYONARA BABY'
  !
  !
  IF (Excitation_Log .AND. l_val /= 0) THEN
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
     ! Build U(3) > SO(3) > SO(2) Basis
     ! ALLOCATE BASIS
     ALLOCATE(SO3_Basis(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "SO3_Basis allocation request denied."
        STOP
     ENDIF
     !
     CALL SO3_BASIS_VIBRON(N_val, 0, SO3_Basis) ! Build SO3 basis
     !
     !
     ! Allocate and Build Hamiltonian and Hamiltonian Operator Matrices
     !
     CALL Build_SO3DS_Operator_Matrices(N_val, 0, dim_block, SO3_Basis)
     !
     !
     CALL Build_Ham_4Body_SO3(N_val, 0, dim_block, SO3_Basis)  
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
  ! L_VAL BLOCK DIMENSION
  dim_block = DIM_L_BLOCK(N_val, l_val)
  !
  ! Build U(3) > SO(3) > SO(2) Basis
  ! ALLOCATE BASIS
  IF (ALLOCATED(SO3_Basis)) THEN    
     DEALLOCATE(SO3_Basis, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "SO3_Basis deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  ALLOCATE(SO3_Basis(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "SO3_Basis allocation request denied."
     STOP
  ENDIF
  !
  CALL SO3_BASIS_VIBRON(N_val, L_val, SO3_Basis) ! Build SO3 basis
  !
  ! Allocate and Build Hamiltonian and Hamiltonian Operator Matrices
  !
  CALL Build_SO3DS_Operator_Matrices(N_val, L_val, dim_block, SO3_Basis)
  !
  !
  CALL Build_Ham_4Body_SO3(N_val, L_val, dim_block, SO3_Basis)  
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
     WRITE(*,*) ' '
     WRITE(*,*) ' Hamiltonian Matrix'
     DO state_index = 1, dim_block                         
        WRITE(*,*) (Ham_matrix(state_index, state_index_2), state_index_2 = 1, dim_block)
     ENDDO
  ENDIF
  !
  ! Hamiltonian Diagonalization
  !
  !      
  ! Check time
  CALL CPU_TIME(time_check)
  IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
       " Time (3) :: U3 Hamiltonian Building ", time_check - time_check_ref
  time_check_ref = time_check
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
     IF (L_val == 0_I4B) THEN
        GS_energy = eigenval_vector(1)
        IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "GS_energy = ", GS_energy
     ENDIF
     !
     eigenval_vector = eigenval_vector - GS_energy
     !
  ENDIF
  !
  ! Check time
  CALL CPU_TIME(time_check)
  IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
       " Time (4) :: U3 Model Hamiltonian diagonalization ", time_check - time_check_ref
  time_check_ref = time_check
  !
  !
  ! CALCULATE IPR
  !
  IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "L_val = ", L_val
  !
  IF (Eigenvec_Log) THEN 
     !
     DO state_index = 1, dim_block
        !
        WRITE(UNIT = *, FMT = *) state_index, Eigenval_vector(state_index), &
             Inv_Part_Ratio(Ham_matrix(1:dim_block, state_index))
        !
     ENDDO
     ! 
  ELSE
     !
     DO state_index = 1, dim_block
        !
        WRITE(UNIT = *, FMT = *) state_index, Eigenval_vector(state_index)
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
             Eigenval_vector, Diagonal_vector - GS_energy, "so3", Ham_matrix)
     ELSE
        CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
             Eigenval_vector, Diagonal_vector, "so3", Ham_matrix)
     ENDIF
     !
     DEALLOCATE(Diagonal_vector, STAT = IERR)    
     IF (IERR /= 0) STOP 'Diagonal_vector deallocation request denied.'
     !
  ENDIF
  !
  !    
  !
END PROGRAM ipr_4b_2DVM_so3
