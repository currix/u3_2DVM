PROGRAM ipr_2DVM_so3
  !
  ! Program to compute Energy or Inverse Participation Ratio 
  ! of the U(3) 2DVM 2 Body Hamiltonian using displaced oscillator basis
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
  REAL(KIND = DP) :: epsilon_p, alpha_p, beta_p, A_p ! Model Hamiltonian Parameters
  !
  INTEGER(KIND = I4B) :: state_index, state_index_2
  !
  !
  ! NAMELISTS
  NAMELIST/par_aux/ Iprint, Eigenvec_Log, Excitation_Log, Save_avec_Log
  NAMELIST/par_0/ N_val, l_val
  NAMELIST/par_1/ epsilon_p, alpha_p, beta_p, A_p
  !
  !
  ! Read parameters
  !
  READ(UNIT = *, NML = par_aux)
  !
  IF (Iprint > 1) PRINT 10
  READ(UNIT = *, NML = par_0)
  !
  IF (Iprint > 1) PRINT 20
  READ(UNIT = *, NML = par_1)
  !
  !
  IF (Iprint > 1) THEN
     WRITE(UNIT = *, FMT = 5) Iprint, Eigenvec_Log, Excitation_Log, Save_avec_Log
     WRITE(UNIT = *, FMT = 15) N_val, L_val
     WRITE(UNIT = *, FMT = 25) epsilon_p, alpha_p, beta_p, A_p
  ENDIF
  !
  ! TESTS
  IF (N_val < L_val) STOP 'ERROR :: N_VAL < L_VAL. SAYONARA BABY'
  !
  !
  ! HAMILTONIAN PARAMETERS
  !ham_param(1) => n
  Ham_parameter(1) = epsilon_p
  !ham_param(2) => n(n+1)
  Ham_parameter(2) = alpha_p
  !ham_param(3) => l^2
  Ham_parameter(3) = beta_p
  !ham_param(4) =>       => P = N(N+1) - W^2
  Ham_parameter(4) = A_p
  !
  IF (Excitation_Log .AND. l_val /= 0) THEN
     ! Compute L = 0 Ground state
     !
     ! Initialize time
     CALL CPU_TIME(time_check_ref)
     !
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
     ! Build Hamiltonian Matrix
     ! ALLOCATE Hamiltonian Matrix
     ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
        STOP
     ENDIF
     !
     !
     Ham_matrix = 0.0_DP
     CALL Build_Ham_SO3(N_val, 0, dim_block, SO3_Basis, Ham_matrix) 
     !
     ! Hamiltonian Diagonalization
     !
     !
     ! ALLOCATE EIGENVALUES VECTOR
     ALLOCATE(Eigenval_vector(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Eigenval_vector allocation request denied."
        STOP
     ENDIF
     !
     !      
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = 0 ", &
          " Time (1) :: L = 0 Model Hamiltonian Building ", time_check - time_check_ref
     time_check_ref = time_check
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
     ! DEALLOCATE EIGENVALUES VECTOR
     DEALLOCATE(Eigenval_vector, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Eigenval_vector deallocation request denied."
        STOP
     ENDIF
     ! DEALLOCATE Hamiltonian Matrix
     DEALLOCATE(Ham_matrix, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
        STOP
     ENDIF
     ! DEALLOCATE BASIS
     DEALLOCATE(SO3_Basis, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "SO3_Basis deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  ! L_VAL BLOCK DIMENSION
  dim_block = DIM_L_BLOCK(N_val, l_val)
  !
  ! Build U(3) > SO(3) > SO(2) Basis
  ! ALLOCATE BASIS
  ALLOCATE(SO3_Basis(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "SO3_Basis allocation request denied."
     STOP
  ENDIF
  !
  CALL SO3_BASIS_VIBRON(N_val, L_val, SO3_Basis) ! Build SO3 basis
  !
  ! Build Hamiltonian Matrix
  ! ALLOCATE Hamiltonian Matrix 
  ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
     STOP
  ENDIF
  !
  !
  Ham_matrix = 0.0_DP
  CALL Build_Ham_SO3(N_val, L_val, dim_block, SO3_Basis, Ham_matrix) 
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
  IF (Iprint > 2) THEN
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
  ! ALLOCATE EIGENVALUES VECTOR
  ALLOCATE(Eigenval_vector(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Eigenval_vector allocation request denied."
     STOP
  ENDIF
  !
  !      
  ! Check time
  CALL CPU_TIME(time_check)
  IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
       " Time (3) :: U3 Hamiltonian Building ", time_check - time_check_ref
  time_check_ref = time_check
  !
  ! Diagonalize Hamiltonian matrix (LAPACK95)
#ifdef  __GFORTRAN__
  !gfortran
  IF (Eigenvec_Log) THEN
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
  ! DEALLOCATE EIGENVALUES VECTOR
  DEALLOCATE(Eigenval_vector, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Eigenval_vector deallocation request denied."
     STOP
  ENDIF
  ! DEALLOCATE Hamiltonian Matrix
  DEALLOCATE(Ham_matrix, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
     STOP
  ENDIF
  !
  ! DEALLOCATE BASIS
  DEALLOCATE(SO3_Basis, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "SO3_Basis deallocation request denied."
     STOP
  ENDIF
  !
  !    
5 FORMAT(1X, " Iprint = ", I2, "; Eigenvec_LOG = ", L2, "; Excitation_Log = ", L2, "; Save_Avec_Log = ", L2)
10 FORMAT(1X, "Reading  N_val, L_val")
15 FORMAT(1X, "N_val = ", I6, "; L_val = ", I6)
20 FORMAT(1X, "Reading  Hamiltonian Paramenters")
25 FORMAT(1X, "epsilon = ", ES14.7, "; alpha = ", ES14.7, " ; beta = ", ES14.7, "; A = ", ES14.7)
  !
  !
  !
END PROGRAM ipr_2DVM_so3
