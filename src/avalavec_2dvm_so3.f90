PROGRAM avalavec_2DVM_so3
  !
  ! Program to compute eigenvalues and eigenvectors of the U(3) 2DVM 
  ! Two-Body Hamiltonian with the Chain II basis (displaced oscillator)
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
  REAL(KIND = DP) :: epsilon_p, alpha_p, beta_p, A_p ! Model Hamiltonian Parameters
  !
  INTEGER(KIND = I4B) :: state_index, state_index_2, omega
  !
  INTEGER(KIND = I4B) :: L, L_min, L_max
  !
  INTEGER(KIND = I4B) :: out_unit
  !
  CHARACTER(LEN=65) :: output_file
  !
#ifdef _OPENMP
  INTEGER(KIND = I4B) :: threads
#endif
  !
  !
  ! NAMELISTS
  NAMELIST/par_aux/ Iprint, Eigenvec_Log, Excitation_Log, Save_avec_Log
  NAMELIST/par_0/ N_val, L_val, L_min, L_max
  NAMELIST/par_1/ epsilon_p, alpha_p, beta_p, A_p
  !
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
     WRITE(UNIT = *, FMT = 5) Iprint, Eigenvec_Log, Excitation_Log
     WRITE(UNIT = *, FMT = 15) N_val, L_min, L_max
     WRITE(UNIT = *, FMT = 25) epsilon_p, alpha_p, beta_p, A_p
  ENDIF
  !
  ! TESTS
  IF ((N_val < L_max) .OR. (N_val < L_min)) &
       STOP 'ERROR :: N_VAL < L_VAL. SAYONARA BABY'
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
  IF (Excitation_Log) THEN
     ! Compute L = 0 Ground state Energy
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
  !$OMP PARALLEL DO DEFAULT(private), SHARED(N_val,GS_energy, epsilon_p, &
  !$OMP alpha_p, beta_p, A_p, Save_avec_Log, Excitation_Log, Eigenvec_Log, &
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
          WRITE(output_file, '("eigval_so3_N",I1,"_L",I1,".dat")')  N_val, L_val
       ELSE IF ( N_val < 100) THEN 
          WRITE(output_file, '("eigval_so3_N",I2,"_L",I1,".dat")')  N_val, L_val
       ELSE IF ( N_val < 1000) THEN 
          WRITE(output_file, '("eigval_so3_N",I3,"_L",I1,".dat")')  N_val, L_val
       ELSE IF ( N_val < 10000) THEN 
          WRITE(output_file, '("eigval_so3_N",I4,"_L",I1,".dat")')  N_val, L_val
       ELSE
          WRITE(output_file, '("eigval_so3_N",I6,"_L",I1,".dat")')  N_val, L_val
       ENDIF
    ELSE IF (L_val >= 10 .AND. L_val<100) THEN
       IF ( N_val < 10) THEN !to avoid spaces
          WRITE(output_file, '("eigval_so3_N",I1,"_L",I2,".dat")')  N_val, L_val
       ELSE IF ( N_val < 100) THEN 
          WRITE(output_file, '("eigval_so3_N",I2,"_L",I2,".dat")')  N_val, L_val
       ELSE IF ( N_val < 1000) THEN 
          WRITE(output_file, '("eigval_so3_N",I3,"_L",I2,".dat")')  N_val, L_val
       ELSE IF ( N_val < 10000) THEN 
          WRITE(output_file, '("eigval_so3_N",I4,"_L",I2,".dat")')  N_val, L_val
       ELSE
          WRITE(output_file, '("eigval_so3_N",I6,"_L",I2,".dat")')  N_val, L_val
       ENDIF
    ELSE
       IF ( N_val < 10) THEN !to avoid spaces
          WRITE(output_file, '("eigval_so3_N",I1,"_L",I3,".dat")')  N_val, L_val
       ELSE IF ( N_val < 100) THEN 
          WRITE(output_file, '("eigval_so3_N",I2,"_L",I3,".dat")')  N_val, L_val
       ELSE IF ( N_val < 1000) THEN 
          WRITE(output_file, '("eigval_so3_N",I3,"_L",I3,".dat")')  N_val, L_val
       ELSE IF ( N_val < 10000) THEN 
          WRITE(output_file, '("eigval_so3_N",I4,"_L",I3,".dat")')  N_val, L_val
       ELSE
          WRITE(output_file, '("eigval_so3_N",I6,"_L",I3,".dat")')  N_val, L_val
       ENDIF
    ENDIF
    !
  OPEN(UNIT=out_unit, FILE=output_file, STATUS='unknown', ACTION='write') 
 
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
     WRITE(UNIT=out_unit, FMT = *) ' '
     WRITE(UNIT=out_unit,FMT = *) ' Hamiltonian Matrix'
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
#ifndef _OPENMP 
  ! Check time
  CALL CPU_TIME(time_check)
  IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
       " Time (3) :: U3 Hamiltonian Building ", time_check - time_check_ref
  time_check_ref = time_check
#endif
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
        IF (Iprint > 0) WRITE(UNIT = out_unit, FMT = *) "GS_energy = ", GS_energy
     ENDIF
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
  IF (Iprint > 0) WRITE(UNIT = out_unit, FMT = *) "L_val = ", L_val
  !
  DO state_index = 1, dim_block
     !
     omega = SO3_Basis(state_index)%omega_SO3_val
     !
     WRITE(UNIT = out_unit, FMT = *) omega, (N_val-omega)/2, Eigenval_vector(state_index)
     !
     ! Display eigenvectors
     IF (Eigenvec_Log .AND. Iprint > 0) THEN
        !
        DO state_index_2 = 1, dim_block
           !
           omega = SO3_Basis(state_index_2)%omega_SO3_val
           !
           WRITE(UNIT = out_unit, FMT = *) Ham_matrix(state_index_2, state_index), "|", omega, ">"
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  ! Save eigenvector components
  IF (Save_avec_Log) THEN
     !
     IF (Excitation_Log) THEN
        !
        CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
             Eigenval_vector, Diagonal_vector - GS_energy, "so3", Ham_matrix, &
             unit=L_val+20)
     ELSE
        CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, &
             Eigenval_vector, Diagonal_vector, "so3", Ham_matrix, &
             unit=L_val+20)
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
  ENDDO
  !
  !$OMP END PARALLEL DO
  !
5 FORMAT(1X, " Iprint = ", I2, "; Eigenvec_LOG = ", L2, "; Excitation_Log = ", L2, "; Save_Avec_Log = ", L2)
10 FORMAT(1X, "Reading  N_val, L_val")
15 FORMAT(1X, "N_val = ", I6, "; L_val = ", I6)
20 FORMAT(1X, "Reading  Hamiltonian Paramenters")
25 FORMAT(1X, "epsilon = ", ES14.7, "; alpha = ", ES14.7, " ; beta = ", ES14.7, "; A = ", ES14.7)
  !
  !
END PROGRAM avalavec_2DVM_so3
