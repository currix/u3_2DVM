MODULE u3_2dvm_mod
  !
  USE nrtype
  !
  USE defparam_2DVM
  !
  ! by Currix TM
  !
  IMPLICIT NONE
  !
  ! Hamiltonian Matrix
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_matrix 
  ! W^2 Casimir Matrix in U(2) DS
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2_matrix ! Atomic M23 in 4 body Hamiltonian with U(2) basis
  ! n Casimir Matrix in SO(3) DS
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: n_matrix ! Atomic M11 in 4 body Hamiltonian with SO(3) basis
  !
  ! Eigenvalue Vector
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Diagonal_vector ! Hamiltonian Diagonal Elements
  ! Hamiltonian Diagonal Vector
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Eigenval_vector ! Hamiltonian Eigenvalues
  ! Cylindrical Oscillator Basis
  TYPE(u2_bas), DIMENSION(:), ALLOCATABLE :: U2_Basis !  U(3) > U(2) > SO(2) basis
  ! Displaced Oscillator Basis
  TYPE(so3_bas), DIMENSION(:), ALLOCATABLE :: SO3_Basis !  U(3) > SO(3) > SO(2) basis
  !
  ! Model Hamiltonian Parameters
  INTEGER(KIND = I4B), PARAMETER :: n_modham_param = 2
  REAL(KIND = DP), DIMENSION(1:n_modham_param) :: ModHam_parameter 
  !
  ! General 2_body Hamiltonian Parameters
  INTEGER(KIND = I4B), PARAMETER :: n_ham_param = 4
  REAL(KIND = DP), DIMENSION(1:n_ham_param) :: Ham_parameter 
  !
  ! Number of parameters in the general 4-body hamiltonian
  INTEGER(KIND = I4B) :: NPMAX = 14
  ! Hamiltonian 4-Body Hamiltonian Parameters
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: H_4b_pars
  !
  ! 4 Body Hamiltonian case U(2) Basis
  ! W^4  P46 in 4 body Hamiltonian 
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M46_P13_matrix
  !
  ! 4 Body Hamiltonian case SO(3) Basis
  ! n^2 P21 in 4 body Hamiltonian
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M21_P2_matrix 
  ! n^3 P31 in 4 body Hamiltonian 
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M31_P5_matrix
  ! n l^2 P32 in 4 body Hamiltonian
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M32_P6_matrix
  ! n^4 P41 in 4 body Hamiltonian
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M41_P8_matrix
  ! n^2 l^2 P42 in 4 body Hamiltonian
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M42_P9_matrix
  !
  ! Operators common to both Dynamical Symmetries
  ! Wbar^2 Casimir Matrix
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2bar_matrix ! Atomic M23bar in 4 body Hamiltonian
  ! n x W^2 + W^2 x n  P33 in 4 body Hamiltonian 
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M33_P7_matrix
  ! n^2 x W^2 + W^2 x n^2  P45 in 4 body Hamiltonian 
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M45_P12_matrix
  ! W^2 x Wbar^2 + Wbar^2 x W^2  P47 in 4 body Hamiltonian 
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: M47_P14_matrix
  INTEGER(KIND = I4B) :: Ierr, Iprint
  !
  LOGICAL :: Eigenvec_Log     ! If .T. compute eigenvalues and eigenvectors
  LOGICAL :: Excitation_Log   ! If .T. compute excitation energy with respect to L = 0 g.s.
  LOGICAL :: Save_avec_Log   ! If .T. save eigenvector components.
  !
  !
  REAL(KIND = DP) :: GS_energy = 0.0_DP ! L = 0 Ground state energy 
  !
  !
  REAL(KIND = DP), PARAMETER :: Zero_Parameter = 1.0E-20_DP ! Zero limit to discard a parameter evaluation
  !
  REAL(KIND = SP) :: time_check, time_check_ref ! Auxiliary parameters for benchmarking
  !
CONTAINS
  !
  FUNCTION DIM_L_BLOCK(N_val, L_val)
    !
    IMPLICIT NONE
    !
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    !
    INTEGER(KIND = I4B) :: DIM_L_BLOCK 
    !
    !
    DIM_L_BLOCK = (N_val - MOD(N_val - l_val ,2_I4B) - l_val)/2_I4B + 1_I4B
    !
    !
  END FUNCTION DIM_L_BLOCK
  !
  SUBROUTINE U2_BASIS_VIBRON(N_val, L_val, U2_Basis) 
    !
    ! Subroutine to build the U(3) > U(2) > SO(2) basis in the 2DVM
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Vibrational Angular momentum
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(OUT) :: U2_Basis
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, np
    !
    index = 1_I4B
    !
    np = L_val
    !
    DO WHILE (np <= N_val) 
       !
       U2_Basis(index)%np_U2_val = np
       !
       index = index + 1_I4B
       np = np + 2_I4B
       !
    ENDDO
    !
    U2_Basis(1:Index-1)%N_U3_val = N_val
    U2_Basis(1:Index-1)%L_val = L_val
    !
    !
  END SUBROUTINE U2_BASIS_VIBRON
  !
  !
  SUBROUTINE SO3_BASIS_VIBRON(N_val, L_val, SO3_Basis) 
    !
    ! Subroutine to build the U(3) > SO(3) > SO(2) basis in the 2DVM
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Vibrational Angular momentum
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(OUT) :: SO3_Basis
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, omega
    !
    index = 1_I4B
    !
    omega = L_val + MOD(N_val-L_val,2)
    !
    DO WHILE (omega <= N_val) 
       !
       SO3_Basis(index)%omega_SO3_val = omega
       !
       index = index + 1_I4B
       omega = omega + 2_I4B
       !
    ENDDO
    !
    SO3_Basis(1:Index-1)%N_U3_val = N_val
    SO3_Basis(1:Index-1)%L_val = L_val
    !
    !
  END SUBROUTINE SO3_BASIS_VIBRON
  !
  SUBROUTINE SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2_matrix, W2bar_flag)
    !
    ! Subroutine to build the SO(3) W^2 Casimir operator in the 2DVM Chain I
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Vibrational Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Block dimension for N_val and L_val
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis ! Cylindrical oscillator basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2_matrix
    !
    ! Optional argument :: when equal to 1 compute Wbar2 Casimir
    INTEGER(KIND = I4B), OPTIONAL, INTENT(IN) :: W2bar_flag
    !
    ! Local variables
    INTEGER(KIND = I4B) :: index
    REAL(KIND = DP) :: Nvalue, Lvalue
    REAL(KIND = DP) :: npvalue, tempval, non_diag_phase
    !
    ! Initialize
    W2_matrix = 0.0_DP
    Nvalue = REAL(N_val,DP)
    Lvalue =  REAL(L_val,DP)
    !
    non_diag_phase = -1.0_DP
    IF (PRESENT(W2bar_flag) .AND. W2bar_flag == 1_I4B) non_diag_phase = 1.0_DP
    !
    ! NON-DIAGONAL PART
    DO index = 2, dim_block   
       npvalue = REAL(U2_Basis(index)%np_U2_val, DP)
       tempval = SQRT((Nvalue - npvalue + 2.0_DP)*(Nvalue - npvalue + 1.0_DP)*(npvalue + Lvalue)*(npvalue - Lvalue))
       W2_matrix(index-1,index) = non_diag_phase*tempval
       W2_matrix(index,index-1) =  W2_matrix(index-1,index)
    ENDDO
    !
    ! DIAGONAL PART
    DO index = 1, dim_block                     
       npvalue = REAL(U2_Basis(index)%np_U2_val, DP)
       W2_matrix(index,index) = (Nvalue - npvalue)*(npvalue + 2.0_DP) + (Nvalue - npvalue + 1.0_DP)*npvalue + Lvalue*Lvalue
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE SO3_Cas_Build
  !
  SUBROUTINE Build_Mod_Ham(N_val, L_val, dim_block, U2_Basis, W2_casimir, Ham_U3_mat) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Model Hamiltonian
    ! Cylindrical Oscillator Basis U(2)
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis   ! U2 basis   { | [N] np L > }
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2_casimir ! SO(3) W^2 casimir matrix 
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_U3_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: oper, U2_state
    REAL(KIND = DP) :: npvalue, Nvalue
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    ! Build Hamiltonian
    operator : DO oper = 1, n_modham_param
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(ModHam_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
       CASE (1) ! n ---> DIAGONAL
          !
          DO U2_state = 1, dim_block
             !
             npvalue = REAL(U2_Basis(U2_state)%np_U2_val, DP)
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  ModHam_parameter(1)*npvalue
             !
          ENDDO
          !
          !
       CASE (2) ! Pairing Operator
          !
          CALL SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2_casimir)
          !
          Ham_U3_mat = Ham_U3_mat -  ModHam_parameter(2)*W2_casimir
          !
          ! Diagonal Pairing contribution
          DO U2_state = 1, dim_block
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  ModHam_parameter(2)*Nvalue*(Nvalue + 1_I4B)
             !
          ENDDO
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
  END SUBROUTINE BUILD_MOD_HAM
  !
  SUBROUTINE Build_Ham(N_val, L_val, dim_block, U2_Basis, W2_casimir, Ham_U3_mat) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Two-Body Hamiltonian
    ! Cylindrical Oscillator Basis U(2)
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis   ! U2 basis   { | [N] np L > }
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2_casimir ! SO(3) W^2 casimir matrix 
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_U3_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: oper, U2_state
    REAL(KIND = DP) :: npvalue, Nvalue, lvalue2
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    ! Build Hamiltonian
    operator : DO oper = 1, n_ham_param
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(Ham_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
       CASE (1) ! n ---> DIAGONAL
          !
          DO U2_state = 1, dim_block
             !
             npvalue = REAL(U2_Basis(U2_state)%np_U2_val, DP)
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  Ham_parameter(1)*npvalue
             !
          ENDDO
          !
          !
       CASE (2) ! n(n+1) ---> DIAGONAL
          !
          DO U2_state = 1, dim_block
             !
             npvalue = REAL(U2_Basis(U2_state)%np_U2_val, DP)
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  Ham_parameter(2)*npvalue*(npvalue + 1.0_DP)
             !
          ENDDO
          !
          !
       CASE (3) ! l^2 ---> DIAGONAL
          !
          lvalue2 = REAL(U2_Basis(1)%L_val, DP)**2
          ! 
          DO U2_state = 1, dim_block
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  Ham_parameter(3)*lvalue2
             !
          ENDDO
          !
          !
       CASE (4) ! Pairing Operator
          !
          CALL SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2_casimir)
          !
          Ham_U3_mat = Ham_U3_mat - Ham_parameter(4)*W2_casimir ! notice the minus sign for P = N(N+1)-W^2
          !
          ! Diagonal Pairing contribution
          DO U2_state = 1, dim_block
             !
             Ham_U3_mat(U2_state, U2_state) = Ham_U3_mat(U2_state, U2_state) + &
                  Ham_parameter(4)*Nvalue*(Nvalue + 1_I4B)
             !
          ENDDO
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
  END SUBROUTINE BUILD_HAM
  !
  FUNCTION POCCHAMMER_S(a, b)
    !
    ! Pocchammer symbol as defined in Frank and Isacker 5.176 p.162
    ! (rising factorial)
    !
    ! by Currix TM
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=I4B), INTENT(IN) :: a, b
    INTEGER(KIND=I4B) :: POCCHAMMER_S
    !
    INTEGER(KIND=I4B) :: index  
    !
    POCCHAMMER_S = a
    !
    DO index = 1, b - 1
       POCCHAMMER_S = POCCHAMMER_S * (a + index)
    ENDDO
    !
    !
    !
  END FUNCTION POCCHAMMER_S
  !
  !
  SUBROUTINE U2_Cas_Build(N_val, L_val, dim_block, SO3_basis, n_matrix)
    !
    ! Subroutine to build the U(2) n Casimir operator in the 2DVM Chain II
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Vibrational Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Block dimension for N_val and L_val
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(IN) :: SO3_Basis ! Displaced oscillator basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: n_matrix
    !
    ! Local variables
    INTEGER(KIND = I4B) :: SO3_state
    REAL(KIND = DP) :: Nvalue, Lvalue
    REAL(KIND = DP) :: omega_val
    !
    ! Initialize
    n_matrix = 0.0_DP
    Nvalue = REAL(N_val,DP)
    Lvalue =  REAL(L_val,DP)
    !
    ! Number Operator n ---> NON DIAGONAL Eq. (35) PRA (Check Typos)
    !
    ! Diagonal contribution 
    DO SO3_state = 1, dim_block
       !
       omega_val = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
       !
       n_matrix(SO3_state, SO3_state) = n_matrix(SO3_state, SO3_state) + &
            n_elem_diag(Nvalue, omega_val, L_val)
       !
    ENDDO
    !
    !
    ! Non-Diagonal contribution w2 -> w1+2 (B_0)
    DO SO3_state = 1, dim_block - 1
       !
       omega_val = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
       !
       n_matrix(SO3_state, SO3_state + 1) = n_matrix(SO3_state, SO3_state + 1) + &
            B0(Nvalue,omega_val,L_val)
       n_matrix(SO3_state + 1, SO3_state) = n_matrix(SO3_state, SO3_state + 1) 
       !
    ENDDO
    !
    !
    RETURN
    !
  END SUBROUTINE U2_Cas_Build
  !
  SUBROUTINE Build_Mod_Ham_SO3(N_val, L_val, dim_block, SO3_Basis, Ham_U3_mat) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Model Hamiltonian
    ! Displaced Oscillator Basis SO(3)
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(IN) :: SO3_Basis   ! SO3 basis   { | [N] w L > }
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_U3_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: oper, SO3_state
    REAL(KIND = DP) :: omegaval, Nvalue
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    ! Build Hamiltonian
    operator : DO oper = 1, n_modham_param
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(ModHam_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
       CASE (1) ! Number Operator n ---> NON DIAGONAL Eq. (35) PRA (Check Typos)
          !
          ! Diagonal contribution 
          DO SO3_state = 1, dim_block
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  ModHam_parameter(1)*n_elem_diag(Nvalue,omegaval,L_val)
             !
          ENDDO
          !
          !
          ! Non-Diagonal contribution  w2 -> w1+2 (B_0)
          DO SO3_state = 1, dim_block - 1
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state + 1) = Ham_U3_mat(SO3_state, SO3_state + 1) + &
                  ModHam_parameter(1)*B0(Nvalue,omegaval,L_val)
             !
          ENDDO
          !
          !
       CASE (2) ! Pairing Operator P ---> DIAGONAL
          !
          ! Diagonal Pairing contribution
          DO SO3_state = 1, dim_block
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  ModHam_parameter(2)*( &
                  Nvalue*(Nvalue + 1.0_DP) - omegaval*(omegaval + 1.0_DP) )
             !
          ENDDO
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
  END SUBROUTINE BUILD_MOD_HAM_SO3
  !
  !
  SUBROUTINE Build_Ham_SO3(N_val, L_val, dim_block, SO3_Basis, Ham_U3_mat) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Two Body Hamiltonian
    ! Displaced Oscillator Basis SO(3)
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(IN) :: SO3_Basis   ! SO3 basis   { | [N] w L > }
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: Ham_U3_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: oper, SO3_state
    REAL(KIND = DP) :: omegaval, Nvalue, lval
    !
    !
    Nvalue = REAL(N_val, DP)
    !
    ! Build Hamiltonian
    operator : DO oper = 1, n_ham_param
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(Ham_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
       CASE (1) ! Number Operator n ---> NON DIAGONAL Eq. (35) PRA (Check Typos)
          !
          ! Diagonal contribution 
          DO SO3_state = 1, dim_block
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  Ham_parameter(1)*n_elem_diag(Nvalue,omegaval,L_val)
             !
          ENDDO
          !
          !
          ! Non-Diagonal contribution w2 -> w1+2 (B_0)
          DO SO3_state = 1, dim_block - 1
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state + 1) = Ham_U3_mat(SO3_state, SO3_state + 1) + &
                  Ham_parameter(1)*B0(Nvalue,omegaval,L_val)
             !
          ENDDO
          !
          !
       CASE (2) ! Anharmonicity Operator n(n+1) ---> NON DIAGONAL See notes
          !
          ! Diagonal contribution 
          DO SO3_state = 1, dim_block
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  Ham_parameter(2)*( &
                  (n_elem_diag(Nvalue,omegaval,L_val)**2 + &
                  B0(Nvalue,omegaval,L_val)**2 + &
                  B0(Nvalue,omegaval-2.0_DP,L_val)**2) +  &
                  n_elem_diag(Nvalue,omegaval,L_val) & ! + n
                  )
             !
          ENDDO
          !
          !
          ! Non-Diagonal contribution w2 -> w1+2
          DO SO3_state = 1, dim_block - 1
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state + 1) = Ham_U3_mat(SO3_state, SO3_state + 1) + &
                  Ham_parameter(2)*( &
                  (n_elem_diag(Nvalue,omegaval,L_val)*B0(Nvalue,omegaval,L_val) + &
                  n_elem_diag(Nvalue,omegaval+2,L_val)*B0(Nvalue,omegaval,L_val)) + & 
                  B0(Nvalue,omegaval,L_val) & ! + n
                  )
             !
          ENDDO
          !
          !
          ! Non-Diagonal contribution w2 -> w1+4
          DO SO3_state = 1, dim_block - 2
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state + 2) = Ham_U3_mat(SO3_state, SO3_state + 2) + &
                  Ham_parameter(2)*( &
                  (B0(Nvalue,omegaval,L_val)*B0(Nvalue,omegaval+2.0_DP,L_val)) &
                  )
             !
          ENDDO
          !
          !
       CASE (3) ! Angular Momentum Operator l^2 ---> DIAGONAL
          !
          ! Diagonal Pairing contribution
          DO SO3_state = 1, dim_block
             !
             lval = REAL(SO3_Basis(SO3_state)%L_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  Ham_parameter(3)*(lval**2)
             !
          ENDDO
          !
       CASE (4) ! Pairing Operator P ---> DIAGONAL
          !
          ! Diagonal Pairing contribution
          DO SO3_state = 1, dim_block
             !
             omegaval = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
             !
             Ham_U3_mat(SO3_state, SO3_state) = Ham_U3_mat(SO3_state, SO3_state) + &
                  Ham_parameter(4)*( &
                  Nvalue*(Nvalue + 1.0_DP) - omegaval*(omegaval + 1.0_DP) )
             !
          ENDDO
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
  END SUBROUTINE BUILD_HAM_SO3
  !
  FUNCTION n_elem_diag(Nv, wv, L_val)
    !
    REAL(KIND = DP), INTENT(IN) :: Nv, wv ! REAL N and omega
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    REAL(KIND = DP) :: n_elem_diag
    !
    !  Local variables
    REAL(KIND = DP) :: lv 
    !
    lv = REAL(L_val,DP)
    !
    IF (wv < 0.0_DP .OR. wv > Nv) THEN ! Test if wv is a valid value
       !
       n_elem_diag = 0.0_DP
       !
    ELSE
       !
       n_elem_diag = (Nv - wv) * &
            ((wv-lv+2.0_DP)*(wv-lv+1.0_DP) + (wv+lv+2.0_DP)*(wv+lv+1.0_DP)) / &
            (2.0_DP*(2.0_DP*wv+1.0_DP)*(2.0_DP*wv+3.0_DP)) + &
            (Nv + wv + 1.0_DP) * &
            ((wv+lv)*(wv+lv-1.0_DP) + (wv-lv)*(wv-lv-1.0_DP)) / &
            (2.0_DP*(2.0_DP*wv-1.0_DP)*(2.0_DP*wv+1.0_DP))
       !
    ENDIF
  END FUNCTION n_elem_diag
  !
  FUNCTION B0(Nv, wv, L_val)
    !
    ! <N w2 l | n | N w1 l > = A(N,w1,l) delta_w2,w1 +
    !                          B0(N,w1,l) delta_w2,w1+2 +
    !                          B1(N,w1,l) delta_w2,w1-2
    !
    ! And B1(N,w1,l) = B0(N,w1-2,l)
    !
    REAL(KIND = DP), INTENT(IN) :: Nv, wv ! REAL N and omega
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    REAL(KIND = DP) :: B0
    !
    !  Local variables
    REAL(KIND = DP) :: lv 
    !
    lv = REAL(L_val,DP)
    !
    ! 
    !
    IF (wv < 0.0_DP .OR. wv > Nv) THEN ! Test if wv is a valid value
       B0 = 0.0_DP
    ELSE
       B0 = SQRT( &
            (Nv - wv)*(Nv + wv + 3.0_DP) * &
            (wv-lv+2.0_DP)*(wv+lv+2.0_DP)*(wv+lv+1.0_DP)*(wv-lv+1.0_DP) / &
            ((2.0_DP*wv+1.0_DP)*(2.0_DP*wv+3.0_DP)**2*(2.0_DP*wv+5.0_DP)) &
            )
    ENDIF
    !
  END FUNCTION B0
  !
  !
  FUNCTION H_Inv_Part_Ratio(N_val, L_val, dim_block, U2_Basis, eigenvector)
    !
    ! Subroutine to compute the Husimi IPR for a 2DVM eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis   ! U2 basis   { | [N] np L > }
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (U(2) basis)
    !
    REAL(KIND = DP) :: H_Inv_Part_Ratio
    !
    ! External function
    REAL(KIND = DP), EXTERNAL :: loggamma
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: n1, n2, n3, n4, i1, i2, i3, i4
    REAL(KIND = DP) :: prodval
    !
    H_Inv_Part_Ratio = 0.0_DP
    !
    DO i1 = 1, dim_block
       n1 = U2_Basis(i1)%np_U2_val
       !
       DO i2 = 1, dim_block
          n2 = U2_Basis(i2)%np_U2_val
          !
          DO i3 = 1, dim_block
             n3 = U2_Basis(i3)%np_U2_val
             !
             DO i4 = 1, dim_block
                n4 = U2_Basis(i4)%np_U2_val
                !
                IF (n1 + n2 /= n3 + n4) CYCLE
                !
                prodval = EXP( &
                     loggamma(REAL(N_val + 1,DP)) + loggamma(REAL(N_val + 3,DP)) - loggamma(REAL(2_I4B*N_val + 3_I4B,DP)) &
                     + loggamma(REAL( (n1+n2)/2 - L_val + 1, DP)) + loggamma(REAL( (n1+n2)/2 + L_val + 1, DP)) &
                     + loggamma(REAL(2_I4B*N_val - n1 - n2 + 1_I4B,DP)) &
                     - 0.5_DP * &
                     (loggamma(REAL(N_val-n1+1,DP)) + loggamma(REAL(N_val-n2+1,DP)) &
                     + loggamma(REAL(N_val-n3+1,DP)) + loggamma(REAL(N_val-n4+1,DP)) &
                     + loggamma(REAL((n1 + l_val)/2 + 1,DP)) +  loggamma(REAL((n2 + l_val)/2 + 1,DP)) & 
                     + loggamma(REAL((n3 + l_val)/2 + 1,DP)) +  loggamma(REAL((n4 + l_val)/2 + 1,DP)) & 
                     + loggamma(REAL((n1 - l_val)/2 + 1,DP)) +  loggamma(REAL((n2 - l_val)/2 + 1,DP)) & 
                     + loggamma(REAL((n3 - l_val)/2 + 1,DP)) +  loggamma(REAL((n4 - l_val)/2 + 1,DP)) &
                     ) &
                     )
                !
                prodval = eigenvector(i1)* eigenvector(i2)* eigenvector(i3)* eigenvector(i4)*prodval
                !
                H_Inv_Part_Ratio = H_Inv_Part_Ratio + prodval
                !
             ENDDO
             !
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    !
  END FUNCTION H_Inv_Part_Ratio
  !
  FUNCTION Inv_Part_Ratio(eigenvector)
    !
    ! Subroutine to compute the IPR for a 2DVM eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (U(2) basis)
    !
    REAL(KIND = DP) :: Inv_Part_Ratio
    !
    !
    Inv_Part_Ratio = 1.0_DP/SUM(eigenvector**4)
    !
  END FUNCTION Inv_Part_Ratio
  !
  !
  SUBROUTINE SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, Eigvals, Diagonal, Basis, Ham_U3_mat)
    !
    ! Subroutine to save the energies, the diagonal of the Hamiltonian and the
    ! components of the eigenvectors of the U(3) 2DVM Model Hamiltonian.
    !
    ! Eigenvectors are saved column-wise.
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum    
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Eigvals ! Eigenvalues
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Diagonal ! Hamiltonian Diagonal 
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    CHARACTER(LEN=*), INTENT(IN) :: Basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: Ham_U3_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: basis_index
    !
    CHARACTER(LEN=65) :: output_filename
    !
    !
    ! Build filename
    IF ( N_val < 10) THEN !to avoid spaces
       WRITE(output_filename, '("eigvec_",A,"_N",I1,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 100) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I2,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 1000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I3,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 10000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I4,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE
       WRITE(output_filename, '("eigvec_",A,"_N",I6,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ENDIF
    !
    OPEN(UNIT = 76, FILE = output_filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !
    WRITE(UNIT=76, FMT=*) "# N = ", N_val, "; L = ", L_val, " ; ", Basis,  " basis"
    !
    WRITE(UNIT=76, FMT=*) "# Eigenvalues "
    WRITE(UNIT=76, FMT=*) Eigvals(1:dim_block)
    !
    WRITE(UNIT=76, FMT=*) "# Hamiltonian Diagonal "
    WRITE(UNIT=76, FMT=*) Diagonal(1:dim_block)
    !
    DO basis_index = 1, dim_block
       WRITE(UNIT=76, FMT=*) Ham_U3_mat(basis_index, 1:dim_block)
    ENDDO
    !
  END SUBROUTINE SAVE_EIGENV_COMPONENTS
  !
  SUBROUTINE Build_U2DS_Operator_Matrices(N_val, L_val, dim_block, U2_Basis)
    !
    ! Subroutine to allocate and build the Operator Matrices for the U(3) 2DVM Four-Body Hamiltonian
    !
    ! Cylindrical Oscillator Basis U(2)
    !
    !     H = P11 n + 
    !         P21 n^2 + P22 l^2 + P23 W^2 +  
    !         P31 n^3 + P32 n·l^2 + P33 (n·W^2 + W^2·n) +
    !         P41 n^4 + P42 n^2·l^2 + P43 l^4 + P44 l^2·W^2 + 
    !         P45 (n^2·W^2 + W^2·n^2) + P46 W^4 + P47 (W^2·Wbar^2 + Wbar^2·W^2)/2
    !
    !     HPAR(1)   P11
    !     HPAR(2)   P21
    !     HPAR(3)   P22
    !     HPAR(4)   P23
    !     HPAR(5)   P31
    !     HPAR(6)   P32
    !     HPAR(7)   P33
    !     HPAR(8)   P41
    !     HPAR(9)   P42
    !     HPAR(10)  P43
    !     HPAR(11)  P44
    !     HPAR(12)  P45
    !     HPAR(13)  P46
    !     HPAR(14)  P47
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis   ! U2 basis   { | [N] np L > }
    !
    !
    INTEGER(KIND = I4B) :: U2_state
    REAL(KIND = DP), DIMENSION(1:dim_block) :: np_values ! n values
    REAL(KIND = DP) :: npvalue, npvaluep1, sqnpvalue, sqnpvaluep1, Nvalue, lvalue2
    !
    IF (IPRINT >= 2) WRITE(*,*) 'Subroutine Build_U2DS_Operator_Matrices'
    !
    Nvalue = REAL(N_val, DP)
    !
    lvalue2 = REAL(L_val**2, DP)
    !
    !
    DO U2_state = 1, dim_block
       !
       np_values(U2_state) = REAL(U2_Basis(U2_state)%np_U2_val, DP)
       !
    ENDDO
    ! Allocate Hamiltonian Matrix
    IF (ALLOCATED(Ham_matrix)) THEN
       DEALLOCATE(Ham_matrix, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Ham_matrix  DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
       STOP
    ENDIF
    !
    ! Allocate Eigenvalues Matrix
    IF (ALLOCATED(Eigenval_vector)) THEN
       DEALLOCATE(Eigenval_vector, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Eigenval_vector DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(Eigenval_vector(1:dim_block), STAT = Ierr)
    IF (Ierr /= 0) STOP 'Eigenval_vector ALLOCATION ERROR'
    !
    ! W^2 Casimir matrix
    IF (ALLOCATED(W2_matrix)) THEN
       DEALLOCATE(W2_matrix, STAT = Ierr)
       IF (Ierr /= 0) STOP 'W2_matrix  DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(W2_matrix(1:dim_block,1:dim_block), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "W2_casimir allocation request denied."
       STOP
    ENDIF
    !
    ! Build W^2 Casimir matrix Parameter 4 and atomic operator 
    CALL SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2_matrix)
    !
    ! n x W^2 + W^2 x n matrix Parameter 7
    IF (H_4b_pars(7) /= 0) THEN
       !
       IF (ALLOCATED(M33_P7_matrix)) THEN
          DEALLOCATE(M33_P7_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'Ham_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M33_P7_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M33_P7_matrix allocation request denied."
          STOP
       ENDIF
       M33_P7_matrix = 0.0_DP
       !
       DO U2_state = 1, dim_block-1
          !
          npvalue = np_values(U2_state)
          npvaluep1 = np_values(U2_state+1)
          !
          M33_P7_matrix(U2_state, U2_state) = 2.0_DP*npvalue*W2_matrix(U2_state, U2_state)
          M33_P7_matrix(U2_state, U2_state+1) = (npvaluep1+npvalue)*W2_matrix(U2_state, U2_state+1)
          M33_P7_matrix(U2_state+1, U2_state) = M33_P7_matrix(U2_state, U2_state+1)
          !
       ENDDO
       !
       ! U2_state = dim_block case
       npvalue = np_values(dim_block)
       M33_P7_matrix(dim_block, dim_block) = 2.0_DP*npvalue*W2_matrix(dim_block, dim_block)
       !
    ENDIF
    !
    !
    ! n^2 x W^2 + W^2 x n^2 matrix Parameter 12
    IF (H_4b_pars(12) /= 0) THEN
       IF (ALLOCATED(M45_P12_matrix)) THEN
          DEALLOCATE(M45_P12_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M45_P12_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M45_P12_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M45_P12_matrix allocation request denied."
          STOP
       ENDIF
       M45_P12_matrix = 0.0_DP
       !
       DO U2_state = 1, dim_block-1
          !
          sqnpvalue = np_values(U2_state)**2
          sqnpvaluep1 = np_values(U2_state+1)**2
          !
          M45_P12_matrix(U2_state, U2_state) = 2.0_DP*sqnpvalue*W2_matrix(U2_state, U2_state)
          M45_P12_matrix(U2_state, U2_state+1) = (sqnpvaluep1+sqnpvalue)*W2_matrix(U2_state, U2_state+1)
          M45_P12_matrix(U2_state+1, U2_state) = M45_P12_matrix(U2_state, U2_state+1)
          !
       ENDDO
       !
       ! U2_state = dim_block case
       sqnpvalue = np_values(dim_block)**2
       M45_P12_matrix(dim_block, dim_block) = 2.0_DP*sqnpvalue*W2_matrix(dim_block, dim_block)
       !
    ENDIF
    !
    ! W^4 matrix Parameter 13
    IF (H_4b_pars(13) /= 0) THEN
       !
       IF (ALLOCATED(M46_P13_matrix)) THEN
          DEALLOCATE(M46_P13_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M46_P13_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M46_P13_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M46_P13_matrix allocation request denied."
          STOP
       ENDIF
       M46_P13_matrix = 0.0_DP
       !
       M46_P13_matrix = MATMUL(W2_matrix,W2_matrix)
       !
    ENDIF
    !
    ! Wbar^2 and W^2 x Wbar^2 + Wbar^2 x W^2 matrix Parameter 14
    IF (H_4b_pars(14) /= 0) THEN
       IF (ALLOCATED(W2bar_matrix)) THEN
          DEALLOCATE(W2bar_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'W2bar_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(W2bar_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "W2bar_matrix allocation request denied."
          STOP
       ENDIF
       !
       ! Build Wbar^2 Casimir matrix
       CALL SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2bar_matrix, 1_I4b)
       !
       IF (ALLOCATED(M47_P14_matrix)) THEN
          DEALLOCATE(M47_P14_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M47_P14_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M47_P14_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M47_P14_matrix allocation request denied."
          STOP
       ENDIF
       M47_P14_matrix = 0.0_DP
       !
       M47_P14_matrix = 0.5_DP*&
            (MATMUL(W2bar_matrix, W2_matrix) + MATMUL(W2_matrix, W2bar_matrix) )
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE Build_U2DS_Operator_Matrices
  !
  SUBROUTINE Build_Ham_4Body_U2(N_val, L_val, dim_block, U2_Basis) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Four-Body Hamiltonian
    !
    ! Cylindrical Oscillator Basis U(2)
    !
    !     H = P11 n + 
    !         P21 n^2 + P22 l^2 + P23 W^2 +  
    !         P31 n^3 + P32 n·l^2 + P33 (n·W^2 + W^2·n) +
    !         P41 n^4 + P42 n^2·l^2 + P43 l^4 + P44 l^2·W^2 + 
    !         P45 (n^2·W^2 + W^2·n^2) + P46 W^4 + P47 (W^2·Wbar^2 + Wbar^2·W^2)/2
    !
    !     HPAR(1)   P11
    !     HPAR(2)   P21
    !     HPAR(3)   P22
    !     HPAR(4)   P23
    !     HPAR(5)   P31
    !     HPAR(6)   P32
    !     HPAR(7)   P33
    !     HPAR(8)   P41
    !     HPAR(9)   P42
    !     HPAR(10)  P43
    !     HPAR(11)  P44
    !     HPAR(12)  P45
    !     HPAR(13)  P46
    !     HPAR(14)  P47
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(u2_bas), DIMENSION(:), INTENT(IN) :: U2_Basis   ! U2 basis   { | [N] np L > }
    !
    !
    INTEGER(KIND = I4B) :: U2_state
    REAL(KIND = DP), DIMENSION(1:dim_block) :: np_values ! n values
    REAL(KIND = DP) :: npvalue, Nvalue, lvalue2
    !
    !
    IF (IPRINT >= 2) WRITE(*,*) 'Subroutine Build_Ham_4Body_U2'
    !
    Nvalue = REAL(N_val, DP)
    !
    lvalue2 = REAL(L_val**2, DP)
    !
    !
    DO U2_state = 1, dim_block
       !
       np_values(U2_state) = REAL(U2_Basis(U2_state)%np_U2_val, DP)
       !
    ENDDO
    !
    ! Build Hamiltonian
    !
    Ham_matrix = 0.0_DP
    !
    ! Diagonal contributions in the U(2) dynamical symmetry
    !
    DO U2_state = 1, dim_block
       !
       npvalue = np_values(U2_state)
       !
       Ham_matrix(U2_state, U2_state) = Ham_matrix(U2_state, U2_state) + &
            ! npvalue polynomial term
            npvalue * ( &
            H_4b_pars(1) + & ! n -- P11 -- 1
            npvalue * ( &
            H_4b_pars(2) + & ! n^2 -- P21 -- 2
            npvalue * ( &
            H_4b_pars(5) + & ! n^3 -- P31 -- 5            !
            npvalue * &
            H_4b_pars(8) & ! n^4 -- P41 -- 8
            ) &
            ) &
            ) 
       !
       ! lvalue2 polynomial term
       Ham_matrix(U2_state, U2_state) = Ham_matrix(U2_state, U2_state) + &
            lvalue2 * ( &
            H_4b_pars(3) + &  ! l^2 -- P22 -- 3
            lvalue2 *  &
            H_4b_pars(10) &   ! l^4 -- P43 -- 10
            ) 
       !
       Ham_matrix(U2_state, U2_state) = Ham_matrix(U2_state, U2_state) + &
            H_4b_pars(6)*npvalue*lvalue2  ! n l^2 -- P32 -- 6
       !                    
       Ham_matrix(U2_state, U2_state) = Ham_matrix(U2_state, U2_state) + &
            H_4b_pars(9)*npvalue*npvalue*lvalue2  ! n^2 l^2 -- P42 -- 9
       !
    ENDDO
    !
    ! Non-diagonal U(2) Dynamical Symmetry contributions
    !
    ! W^2 matrix Parameter 4
    IF (H_4b_pars(4) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(4)*W2_matrix
       !
    ENDIF
    !
    ! n x W^2 + W^2 x n matrix Parameter 7
    IF (H_4b_pars(7) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(7)*M33_P7_matrix
       !
    ENDIF
    !
    ! l^2 x W^2 matrix Parameter 11
    IF (H_4b_pars(11) /= 0 .AND. L_val > 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(11)*lvalue2*W2_matrix
       !       
    ENDIF
    !
    ! n^2 x W^2 + W^2 x n^2 matrix Parameter 12
    IF (H_4b_pars(12) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(12)*M45_P12_matrix
       !
    ENDIF
    !
    ! W^4 matrix Parameter 13
    IF (H_4b_pars(13) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(13)*M46_P13_matrix
       !
    ENDIF
    !
    ! Wbar^2 and W^2 x Wbar^2 + Wbar^2 x W^2 matrix Parameter 14
    IF (H_4b_pars(14) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(14)*M47_P14_matrix
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE BUILD_HAM_4BODY_U2
  !
  !
  SUBROUTINE Build_SO3DS_Operator_Matrices(N_val, L_val, dim_block, SO3_Basis)
    !
    ! Subroutine to allocate and build the Operator Matrices for the U(3) 2DVM Four-Body Hamiltonian
    !
    ! Displaced Oscillator Basis SO(3)
    !
    !     H = P11 n + 
    !         P21 n^2 + P22 l^2 + P23 W^2 +  
    !         P31 n^3 + P32 n·l^2 + P33 (n·W^2 + W^2·n) +
    !         P41 n^4 + P42 n^2·l^2 + P43 l^4 + P44 l^2·W^2 + 
    !         P45 (n^2·W^2 + W^2·n^2) + P46 W^4 + P47 (W^2·Wbar^2 + Wbar^2·W^2)/2
    !
    !     HPAR(1)   P11
    !     HPAR(2)   P21
    !     HPAR(3)   P22
    !     HPAR(4)   P23
    !     HPAR(5)   P31
    !     HPAR(6)   P32
    !     HPAR(7)   P33
    !     HPAR(8)   P41
    !     HPAR(9)   P42
    !     HPAR(10)  P43
    !     HPAR(11)  P44
    !     HPAR(12)  P45
    !     HPAR(13)  P46
    !     HPAR(14)  P47
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(IN) :: SO3_Basis   ! SO(3) DS basis   { | [N] w L > }
    !
    !
    INTEGER(KIND = I4B) :: SO3_state, omega_v
    REAL(KIND = DP), DIMENSION(1:dim_block) :: omega_values ! omega values
    REAL(KIND = DP) :: omega, omegap, Nvalue, lvalue2
    !
    IF (IPRINT >= 2) WRITE(*,*) 'Subroutine Build_SO3DS_Operator_Matrices'
    !
    Nvalue = REAL(N_val, DP)
    !
    lvalue2 = REAL(L_val**2, DP)
    !
    !
    DO SO3_state = 1, dim_block
       !
       omega_values(SO3_state) = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
       !
    ENDDO
    ! Allocate Hamiltonian Matrix
    IF (ALLOCATED(Ham_matrix)) THEN
       DEALLOCATE(Ham_matrix, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Ham_matrix  DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "Ham_matrix allocation request denied."
       STOP
    ENDIF
    !
    ! Allocate Eigenvalues Matrix
    IF (ALLOCATED(Eigenval_vector)) THEN
       DEALLOCATE(Eigenval_vector, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Eigenval_vector DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(Eigenval_vector(1:dim_block), STAT = Ierr)
    IF (Ierr /= 0) STOP 'Eigenval_vector ALLOCATION ERROR'
    !
    ! n Casimir matrix
    IF (ALLOCATED(n_matrix)) THEN
       DEALLOCATE(n_matrix, STAT = Ierr)
       IF (Ierr /= 0) STOP 'n_matrix  DEALLOCATION ERROR'
    ENDIF
    ALLOCATE(n_matrix(1:dim_block,1:dim_block), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "n_casimir allocation request denied."
       STOP
    ENDIF
    !
    ! Build n Casimir matrix Parameter 1 and atomic operator 
    CALL U2_Cas_Build(N_val, L_val, dim_block, SO3_basis, n_matrix)
    !
    ! n^2 matrix Parameter 2
    IF (H_4b_pars(2) /= 0 .OR. H_4b_pars(9) /= 0 .OR. H_4b_pars(12) /= 0) THEN
       !
       IF (ALLOCATED(M21_P2_matrix)) THEN
          DEALLOCATE(M21_P2_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M21_P2_matrix DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M21_P2_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M21_P2_matrix allocation request denied."
          STOP
       ENDIF
       M21_P2_matrix = 0.0_DP
       !
       M21_P2_matrix = MATMUL(n_matrix, n_matrix)
       !
       !
    ENDIF
    !
    ! n^3 matrix Parameter 5
    IF (H_4b_pars(5) /= 0) THEN
       !
       IF (ALLOCATED(M31_P5_matrix)) THEN
          DEALLOCATE(M31_P5_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M31_P5_matrix DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M31_P5_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M31_P5_matrix allocation request denied."
          STOP
       ENDIF
       M31_P5_matrix = 0.0_DP
       !
       M31_P5_matrix = MATMUL(MATMUL(n_matrix, n_matrix), n_matrix)
       !
       !
    ENDIF
    !
    ! n^4 matrix Parameter 8
    IF (H_4b_pars(8) /= 0) THEN
       !
       IF (ALLOCATED(M41_P8_matrix)) THEN
          DEALLOCATE(M41_P8_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M41_P8_matrix DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M41_P8_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M41_P8_matrix allocation request denied."
          STOP
       ENDIF
       M41_P8_matrix = 0.0_DP
       !
       M41_P8_matrix = MATMUL(MATMUL(n_matrix, n_matrix), MATMUL(n_matrix, n_matrix))
       !
       !
    ENDIF
    !
    ! n x W^2 + W^2 x n matrix Parameter 7
    IF (H_4b_pars(7) /= 0) THEN
       IF (ALLOCATED(M33_P7_matrix)) THEN
          DEALLOCATE(M33_P7_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M33_P7_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M33_P7_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M33_P7_matrix allocation request denied."
          STOP
       ENDIF
       M33_P7_matrix = 0.0_DP
       !
       DO SO3_state = 1, dim_block-1
          ! 
          omega = omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP)
          omegap = omega_values(SO3_state+1)*(omega_values(SO3_state+1) + 1.0_DP)
          !
          M33_P7_matrix(SO3_state, SO3_state) = 2.0_DP*omega*n_matrix(SO3_state, SO3_state)
          M33_P7_matrix(SO3_state, SO3_state+1) = (omegap+omega)*n_matrix(SO3_state, SO3_state+1)
          M33_P7_matrix(SO3_state+1, SO3_state) = M33_P7_matrix(SO3_state, SO3_state+1)
          !
       ENDDO
       !
       ! SO3_state = dim_block case
       omega = omega_values(dim_block)*(omega_values(dim_block) + 1.0_DP)
       M33_P7_matrix(dim_block, dim_block) = 2.0_DP*omega*n_matrix(dim_block, dim_block)
       !
    ENDIF
    !
    !
    ! n^2 x W^2 + W^2 x n^2 matrix Parameter 12
    ! Double banded operator n^2 --> M21_P2_matrix
    IF (H_4b_pars(12) /= 0) THEN
       IF (ALLOCATED(M45_P12_matrix)) THEN
          DEALLOCATE(M45_P12_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M45_P12_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M45_P12_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M45_P12_matrix allocation request denied."
          STOP
       ENDIF
       M45_P12_matrix = 0.0_DP
       !
       !
       DO SO3_state = 1, dim_block
          ! 
          omega = omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP)
          !
          M45_P12_matrix(SO3_state, SO3_state) = 2.0_DP*omega*M21_P2_matrix(SO3_state, SO3_state)
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block-1
          ! 
          omega = omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP)
          omegap = omega_values(SO3_state+1)*(omega_values(SO3_state+1) + 1.0_DP)
          !
          M45_P12_matrix(SO3_state, SO3_state+1) = (omegap+omega)*M21_P2_matrix(SO3_state, SO3_state+1)
          M45_P12_matrix(SO3_state+1, SO3_state) = M45_P12_matrix(SO3_state, SO3_state+1)
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block-2
          ! 
          omega = omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP)
          omegap = omega_values(SO3_state+2)*(omega_values(SO3_state+2) + 1.0_DP)
          !
          M45_P12_matrix(SO3_state, SO3_state+2) = (omegap+omega)*M21_P2_matrix(SO3_state, SO3_state+2)
          M45_P12_matrix(SO3_state+2, SO3_state) = M45_P12_matrix(SO3_state, SO3_state+2)
          !
       ENDDO
       !
    ENDIF
    !
    ! Wbar^2 and W^2 x Wbar^2 + Wbar^2 x W^2 matrix Parameter 14
    IF (H_4b_pars(14) /= 0) THEN
       IF (ALLOCATED(W2bar_matrix)) THEN
          DEALLOCATE(W2bar_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'W2bar_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(W2bar_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "W2bar_matrix allocation request denied."
          STOP
       ENDIF
       !
       DO SO3_state = 1, dim_block ! Diagonal contribution
          ! 
          omega_v = SO3_Basis(SO3_state)%omega_SO3_val
          !
          W2bar_matrix(SO3_state, SO3_state) = A_SO4BAR(N_val, omega_v, L_val)**2 + &
               B_SO4BAR(N_val, omega_v, L_val)**2 + &
               C_SO4BAR(N_val, omega_v, L_val)**2 + lvalue2
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block - 1 ! First Diagonal contribution
          ! 
          omega_v = SO3_Basis(SO3_state)%omega_SO3_val
          !
          W2bar_matrix(SO3_state, SO3_state + 1) = A_SO4BAR(N_val, omega_v, L_val)*B_SO4BAR(N_val, omega_v + 2_I4B, L_val) + &
               C_SO4BAR(N_val, omega_v, L_val)*A_SO4BAR(N_val, omega_v + 2_I4B, L_val)
          !
          W2bar_matrix(SO3_state + 1, SO3_state) = W2bar_matrix(SO3_state, SO3_state + 1) 
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block - 2 ! Second Diagonal contribution
          ! 
          omega_v = SO3_Basis(SO3_state)%omega_SO3_val
          !
          W2bar_matrix(SO3_state, SO3_state + 2) = B_SO4BAR(N_val, omega_v + 4_I4B, L_val) * C_SO4BAR(N_val, omega_v, L_val)
          !
          W2bar_matrix(SO3_state + 2, SO3_state) = W2bar_matrix(SO3_state, SO3_state + 2) 
          !
       ENDDO
       !
       IF (ALLOCATED(M47_P14_matrix)) THEN
          DEALLOCATE(M47_P14_matrix, STAT = Ierr)
          IF (Ierr /= 0) STOP 'M47_P14_matrix  DEALLOCATION ERROR'
       ENDIF
       ALLOCATE(M47_P14_matrix(1:dim_block,1:dim_block), STAT = IERR)    
       IF (IERR /= 0) THEN
          WRITE(UNIT = *, FMT = *) "M47_P14_matrix allocation request denied."
          STOP
       ENDIF
       M47_P14_matrix = 0.0_DP
       !
       DO SO3_state = 1, dim_block ! Diagonal contribution
          !
          M47_P14_matrix(SO3_state, SO3_state) = W2bar_matrix(SO3_state, SO3_state)* &
               2.0_DP*omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP)
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block - 1 ! First Diagonal contribution
          !
          M47_P14_matrix(SO3_state, SO3_state+1) = W2bar_matrix(SO3_state, SO3_state+1)* &
               (omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP) + &
               omega_values(SO3_state+1)*(omega_values(SO3_state+1) + 1.0_DP))
          !
          M47_P14_matrix(SO3_state+1, SO3_state) =M47_P14_matrix(SO3_state, SO3_state+1)
          !
       ENDDO
       !
       DO SO3_state = 1, dim_block - 2 ! Second Diagonal contribution
          ! 
          M47_P14_matrix(SO3_state, SO3_state+2) = W2bar_matrix(SO3_state, SO3_state+2)* &
               (omega_values(SO3_state)*(omega_values(SO3_state) + 1.0_DP) + &
               omega_values(SO3_state+2)*(omega_values(SO3_state+2) + 1.0_DP))
          !
          M47_P14_matrix(SO3_state+2, SO3_state) =M47_P14_matrix(SO3_state, SO3_state+2)
          !
       ENDDO
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE Build_SO3DS_Operator_Matrices
  !
  !
  SUBROUTINE Build_Ham_4Body_SO3(N_val, L_val, dim_block, SO3_Basis) 
    !
    !
    ! Subroutine to build the U(3) 2DVM Four-Body Hamiltonian
    !
    ! Displaced Oscillator Basis SO(3)
    !
    !     H = P11 n + 
    !         P21 n^2 + P22 l^2 + P23 W^2 +  
    !         P31 n^3 + P32 n·l^2 + P33 (n·W^2 + W^2·n) +
    !         P41 n^4 + P42 n^2·l^2 + P43 l^4 + P44 l^2·W^2 + 
    !         P45 (n^2·W^2 + W^2·n^2) + P46 W^4 + P47 (W^2·Wbar^2 + Wbar^2·W^2)/2
    !
    !     HPAR(1)   P11
    !     HPAR(2)   P21
    !     HPAR(3)   P22
    !     HPAR(4)   P23
    !     HPAR(5)   P31
    !     HPAR(6)   P32
    !     HPAR(7)   P33
    !     HPAR(8)   P41
    !     HPAR(9)   P42
    !     HPAR(10)  P43
    !     HPAR(11)  P44
    !     HPAR(12)  P45
    !     HPAR(13)  P46
    !     HPAR(14)  P47
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    TYPE(so3_bas), DIMENSION(:), INTENT(IN) :: SO3_Basis   ! SO3 basis   { | [N] w L > }
    !
    ! Local variables
    INTEGER(KIND = I4B) :: SO3_state
    REAL(KIND = DP), DIMENSION(1:dim_block) :: omega_values ! w values
    REAL(KIND = DP) :: omega_val, Nvalue, lvalue2
    !
    !
    Nvalue = REAL(N_val, DP)
    lvalue2 = REAL(L_val**2, DP)
    !
    !
    DO SO3_state = 1, dim_block
       !
       omega_values(SO3_state) = REAL(SO3_Basis(SO3_state)%omega_SO3_val, DP)
       !
    ENDDO
    !
    ! Build Hamiltonian
    !
    Ham_matrix = 0.0_DP
    !
    ! Diagonal contributions in the SO(3) dynamical symmetry
    !
    DO SO3_state = 1, dim_block
       !
       omega_val = omega_values(SO3_state)
       !
       Ham_matrix(SO3_state, SO3_state) = Ham_matrix(SO3_state, SO3_state) + &
            !
            ! lvalue2 polynomial term
            lvalue2 * ( &
            H_4b_pars(3) + &  ! l^2 -- P22 -- 3
            lvalue2 *  &
            H_4b_pars(10) &   ! l^4 -- P43 -- 10
            ) 
       !
       Ham_matrix(SO3_state, SO3_state) = Ham_matrix(SO3_state, SO3_state) + &
            !
            ! W^2 polynomial term
            omega_val*(omega_val + 1.0_DP) * ( &
            H_4b_pars(4) + &  ! W^2 -- P23 -- 4
            omega_val*(omega_val + 1.0_DP) *  &
            H_4b_pars(13) &   ! W^4 -- P46 -- 13
            )
            !
       Ham_matrix(SO3_state, SO3_state) = Ham_matrix(SO3_state, SO3_state) + &
            !
            H_4b_pars(11)*omega_val*(omega_val + 1.0_DP) * lvalue2  ! W^2 l^2 -- P44 -- 11
       !            
    ENDDO
    !
    ! Non-diagonal SO(3) Dynamical Symmetry contributions
    !
    ! n matrix Parameter 1
    IF (H_4b_pars(1) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(1)*n_matrix
       !
    ENDIF
    !
    ! n^2 matrix Parameter 2
    IF (H_4b_pars(2) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(2)*M21_P2_matrix
       !
    ENDIF
    !
    ! n^3 matrix Parameter 5
    IF (H_4b_pars(5) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(5)*M31_P5_matrix
       !
    ENDIF
    !
    ! n l^2 matrix Parameter 6
    IF (H_4b_pars(6) /= 0 .AND. L_val /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(6)*lvalue2*n_matrix
       !
    ENDIF
    !
    ! n x W^2 + W^2 x n matrix Parameter 7
    IF (H_4b_pars(7) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(7)*M33_P7_matrix
       !
    ENDIF
    !
    ! n^4 matrix Parameter 8
    IF (H_4b_pars(8) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(8)*M41_P8_matrix
       !
    ENDIF
    !
    ! n^2 l^2 matrix Parameter 9
    IF (H_4b_pars(9) /= 0 .AND. L_val /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(9)*lvalue2*M21_P2_matrix
       !
    ENDIF
    !
    ! n^2 x W^2 + W^2 x n^2 matrix Parameter 12
    IF (H_4b_pars(12) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(12)*M45_P12_matrix
       !
    ENDIF
    !
    ! Wbar^2 and W^2 x Wbar^2 + Wbar^2 x W^2 matrix Parameter 14
    IF (H_4b_pars(14) /= 0) THEN
       !
       Ham_matrix = Ham_matrix + H_4b_pars(14)*M47_P14_matrix
       !
    ENDIF
    !
    RETURN
   !
    !
  END SUBROUTINE BUILD_HAM_4BODY_SO3
  !
  !
  FUNCTION A_SO4BAR(N_val, omega_val, l_val)
    !
    ! <N w l-1 | R- | N w l > 
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, omega_val, l_val ! N, omega and angular momentum
    !
    REAL(KIND = DP) :: A_SO4BAR
    !
    !  Local variables
    REAL(KIND = DP) :: Nv, wv, lv 
    !
    Nv = REAL(N_val,DP)
    wv = REAL(omega_val,DP)
    lv = REAL(L_val,DP)
    !
    A_SO4BAR = (2.0_DP*Nv+3.0_DP)*(2.0_DP*lv + 1.0_DP)* &
         SQRT(0.5_DP*(wv-lv+1.0_DP)*(wv+lv))/ &
         ((2.0_DP*wv-1.0_DP)*(2.0_DP*wv+3.0_DP))
    !
    RETURN
    !
  END FUNCTION A_SO4BAR
  !
  FUNCTION B_SO4BAR(N_val, omega_val, l_val)
    !
    ! <N w-2 l-1 | R- | N w l > 
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, omega_val, l_val ! N, omega and angular momentum
    !
    REAL(KIND = DP) :: B_SO4BAR
    !
    !  Local variables
    REAL(KIND = DP) :: Nv, wv, lv 
    !
    Nv = REAL(N_val,DP)
    wv = REAL(omega_val,DP)
    lv = REAL(L_val,DP)
    !
    B_SO4BAR = -SQRT(2.0_DP* &
         (Nv+wv+1.0_DP)*(Nv-wv+2.0_DP)* &
         (wv+lv)*(wv-lv)*(wv+lv-1.0_DP)*(wv+lv-2.0_DP)/ &
         ((2.0_DP*wv+1.0_DP)*(2.0_DP*wv-1.0_DP)**2*(2.0_DP*wv-3.0_DP)) &
         )
    !
    RETURN
    !
  END FUNCTION B_SO4BAR
  !
  FUNCTION C_SO4BAR(N_val, omega_val, l_val)
    !
    ! <N w+2 l-1 | R- | N w l > 
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, omega_val, l_val ! N, omega and angular momentum
    !
    REAL(KIND = DP) :: C_SO4BAR
    !
    !  Local variables
    REAL(KIND = DP) :: Nv, wv, lv 
    !
    Nv = REAL(N_val,DP)
    wv = REAL(omega_val,DP)
    lv = REAL(L_val,DP)
    !
    C_SO4BAR = SQRT(2.0_DP* &
         (Nv+wv+3.0_DP)*(Nv-wv)* &
         (wv-lv+1.0_DP)*(wv+lv+1.0_DP)*(wv-lv+2.0_DP)*(wv-lv+3.0_DP)/ &
         ((2.0_DP*wv+1.0_DP)*(2.0_DP*wv+3.0_DP)**2*(2.0_DP*wv+5.0_DP)) &
         )
    !
    RETURN
    !
  END FUNCTION C_SO4BAR
  !
END MODULE u3_2dvm_mod
