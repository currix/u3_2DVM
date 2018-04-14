PROGRAM ipr_4b_2DVM_u2
  !
  ! Program to compute Energy and Participation Ratio (PR)
  ! of the U(3) 2DVM 2 Body Hamiltonian
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
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR
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
  !     COMPUTED W2 AND SQUARED W2 (SO(3) CASIMIR) MATRICES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2MAT, W4MAT
  !
  !     COMPUTED W2·WBAR2+WBAR2·W2 (SO(3) CASIMIR x barSO(3) CASIMIR) MATRIX
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2W2BARMAT
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
     ! Build Hamiltonian Matrix
     ! ALLOCATE Hamiltonian Matrix and W^2 Casimir matrix
     ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
        STOP
     ENDIF
     Ham_matrix = 0.0_DP
     !
     ALLOCATE(W2MAT(1:dim_block,1:dim_block), STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR ALLOCATING W2MAT MATRIX'
     W2MAT = 0.0D0
     !     
     ALLOCATE(W4MAT(1:dim_block,1:dim_block), STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR ALLOCATING W4MAT MATRIX'
     W4MAT = 0.0D0
     !     
     ALLOCATE(W2W2BARMAT(1:dim_block,1:dim_block), STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR ALLOCATING W2W2BARMAT MATRIX'
     W2W2BARMAT = 0.0D0
     !
     CALL HBLDU3GEN(N_val, 0, HAM_matrix, W2MAT, W4MAT, W2W2BARMAT, IPRINT)
     !Build_Ham(N_val, 0, dim_block, U2_Basis, W2_matrix, Ham_matrix) 
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
     CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='N', UPLO='U', IL=1, IU=1)
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
     ! DEALLOCATE Hamiltonian Matrix and W^2 Casimir matrix
     DEALLOCATE(Ham_matrix, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
        STOP
     ENDIF
     DEALLOCATE(W2MAT, STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2MAT MATRIX'
     !     
     DEALLOCATE(W4MAT, STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR DEALLOCATING W4MAT MATRIX'
     !     
     DEALLOCATE(W2W2BARMAT, STAT = IERR)
     IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2W2BARMAT MATRIX'
     ! DEALLOCATE BASIS
     DEALLOCATE(U2_Basis, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U2_Basis deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  ! L_VAL BLOCK DIMENSION
  dim_block = DIM_L_BLOCK(N_val, l_val)
  !
  ! Build U(3) > U(2) > SO(2) Basis
  ! ALLOCATE BASIS
  ALLOCATE(U2_Basis(1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "U2_Basis allocation request denied."
     STOP
  ENDIF
  !
  CALL U2_BASIS_VIBRON(N_val, L_val, U2_Basis) ! Build U2 basis
  !
  ! Build Hamiltonian Matrix
  ! ALLOCATE Hamiltonian Matrix and W^2 Casimir matrix
  ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
     STOP
  ENDIF
  !
  !
  Ham_matrix = 0.0_DP
  !
  ALLOCATE(W2MAT(1:dim_block,1:dim_block), STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR ALLOCATING W2MAT MATRIX'
  W2MAT = 0.0D0
  !     
  ALLOCATE(W4MAT(1:dim_block,1:dim_block), STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR ALLOCATING W4MAT MATRIX'
  W4MAT = 0.0D0
  !     
  ALLOCATE(W2W2BARMAT(1:dim_block,1:dim_block), STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR ALLOCATING W2W2BARMAT MATRIX'
  W2W2BARMAT = 0.0D0
  !
  !
  CALL HBLDU3GEN(N_val, L_val, HAM_matrix, W2MAT, W4MAT, W2W2BARMAT, IPRINT)
  !  CALL Build_Ham(N_val, L_val, dim_block, U2_Basis, W2_matrix, Ham_matrix) 
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
  ! Allocate Eigenvalues vector
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
  IF (Eigenvec_Log .OR. Save_avec_Log) THEN
     CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='V', UPLO='U')
  ELSE
     CALL LA_SYEVR(A=Ham_matrix, W=Eigenval_vector, JOBZ='N', UPLO='U')
  ENDIF
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
             Eigenval_vector, Diagonal_vector - GS_energy, "u2", Ham_matrix)
     ELSE
        CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
             Eigenval_vector, Diagonal_vector, "u2", Ham_matrix)
     ENDIF
     !
     DEALLOCATE(Diagonal_vector, STAT = IERR)    
     IF (IERR /= 0) STOP 'Diagonal_vector deallocation request denied.'
     !
  ENDIF
  !
  !    
  ! DEALLOCATE EIGENVALUES VECTOR
  DEALLOCATE(Eigenval_vector, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Eigenval_vector deallocation request denied."
     STOP
  ENDIF
  ! DEALLOCATE Hamiltonian Matrix and W^2 Casimir matrix
  DEALLOCATE(Ham_matrix, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
     STOP
  ENDIF
  DEALLOCATE(W2MAT, STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2MAT MATRIX'
  !     
  DEALLOCATE(W4MAT, STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR DEALLOCATING W4MAT MATRIX'
  !     
  DEALLOCATE(W2W2BARMAT, STAT = IERR)
  IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2W2BARMAT MATRIX'
  !
  ! DEALLOCATE BASIS
  DEALLOCATE(U2_Basis, STAT = IERR)    
  IF (IERR /= 0) THEN
     WRITE(UNIT = *, FMT = *) "U2_Basis deallocation request denied."
     STOP
  ENDIF
  !
  !
! 5 FORMAT(1X, " Iprint = ", I2, "; Eigenvec_LOG = ", L2, "; Excitation_Log = ", L2, "; Save_Avec_Log = ", L2)
! 10 FORMAT(1X, "Reading  N_val, L_val")
! 15 FORMAT(1X, "N_val = ", I6, "; L_val = ", I6)
! 20 FORMAT(1X, "Reading  Hamiltonian Paramenters")
! 25 FORMAT(1X, "epsilon = ", ES14.7, "; alpha = ", ES14.7, " ; beta = ", ES14.7, "; A = ", ES14.7)
  !
  !
END PROGRAM ipr_4b_2DVM_u2
