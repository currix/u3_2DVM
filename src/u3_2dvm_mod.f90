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
  ! W^2 Casimir Matrix
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2_matrix
  ! Eigenvalue Vector
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Diagonal_vector ! Hamiltonian Eigenvalues
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
  SUBROUTINE SO3_Cas_Build(N_val, L_val, dim_block, U2_basis, W2_matrix)
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
    ! Local variables
    INTEGER(KIND = I4B) :: index
    REAL(KIND = DP) :: Nvalue, Lvalue
    REAL(KIND = DP) :: npvalue, tempval
    !
    ! Initialize
    W2_matrix = 0.0_DP
    Nvalue = REAL(N_val,DP)
    Lvalue =  REAL(L_val,DP)
    !
    ! NON-DIAGONAL PART
    DO index = 2, dim_block  
       npvalue = REAL(U2_Basis(index)%np_U2_val, DP)
       tempval = SQRT((Nvalue - npvalue + 2.0_DP)*(Nvalue - npvalue + 1.0_DP)*(npvalue + Lvalue)*(npvalue - Lvalue))
       W2_matrix(index-1,index) = -tempval
       W2_matrix(index,index-1) = -tempval
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
          lvalue2 = REAL(U2_Basis(U2_state)%L_val, DP)**2
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HBLDU3GEN(N_val, L, HAM, W2MAT, W4MAT, W2W2BARMAT, IPRINT)
    !     
    !     SUBROUTINE THAT BUILDS THE GENERAL 1-,2-,3-, AND 4-BODY 
    !     HAMILTONIAN MATRIX IN THE U(3) MODEL FOR BENDING VIBRATIONS
    !
    !     INPUT
    !     N_val : U(3) IRREP LABEL (BENDING)
    !     L     : VIBRATIONAL ANGULAR MOMENTUM LABEL
    !
    !     OUTPUT
    !     HAM     : HAMILTONIAN MATRIX
    !     W2MAT   : ARRAY WITH THE SO(3) CASIMIR (W^2) BLOCK)
    !     W4MAT   : ARRAY WITH THE SQUARED SO(3) CASIMIR (W^4 BLOCK) 
    !     W2W2BARMAT : ARRAY WITH THE OPERATOR W^2·WBAR^2+WBAR^2·W^2 
    !     
    !     BENDING HAMILTONIAN
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
    !     
    !     by Currix TM
    !     
    IMPLICIT NONE
    !
    !     
    !     DEFINITION OF VARIABLES
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, L
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: HAM
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2MAT, W4MAT, W2W2BARMAT
    !
    INTEGER(KIND = I4B), INTENT(IN) :: IPRINT
    !
    !     LOCAL VARIABLES
    INTEGER(KIND = I4B) :: I, I1
    !     TEMPORAL VALUES TO STORE INTEGERS AS DOUBLE PREC. REALS
    REAL(KIND = DP) :: V2, VL
    !     
    !     
    IF (IPRINT.GT.2) WRITE(*,*) 'HAMILTONIAN BUILDING SUBROUTINE STARTS HERE'
    !
    !     MATRIX DIMENSIONS
    dim_block = (N_val-MOD(N_val-L,2)-L)/2+1
    !
    !     
    VL = REAL(L,DP)
    !
    !     NON-DIAGONAL INTERACTIONS
    !
    !     BUILDING W2 BLOCK
    !
    CALL SO3CASBUILD(W2MAT, N_val, L, Iprint)
    !
    !     HAM MULTIPLYING W2 TIMES P23 PARAMETER
    HAM = H_4b_pars(4)*W2MAT
    !
    !     BUILDING W^2·n + n·W^2 BLOCK
    !
    IF (H_4b_pars(7) /= 0) THEN
       !
       !     HAM MULTIPLYING (W^2·n + n·W^2) TIMES P33 PARAMETER
       DO I1 = 1, dim_block
          V2 = REAL(N_val - (2*I1 - 2 + MOD(N_val - L, 2)), DP)
          HAM(I1,I1) = HAM(I1,I1) + H_4b_pars(7)*2.0D0*V2*W2MAT(I1,I1)
       ENDDO
       DO I1 = 1, dim_block - 1
          V2 = REAL(N_val - (2*I1 - 2 + MOD(N_val - L, 2)), DP)
          HAM(I1+1,I1) = HAM(I1+1,I1) + &
               H_4b_pars(7)*2.0_DP*(V2-1.0_DP)*W2MAT(I1+1,I1)
          HAM(I1,I1+1) = HAM(I1,I1+1) + &
               H_4b_pars(7)*2.0_DP*(V2-1.0_DP)*W2MAT(I1,I1+1)
       ENDDO
    ENDIF
    !
    !   HAM MULTIPLYING W^2·l^2 TIMES P44 PARAMETER
    IF (H_4b_pars(11) /= 0 .AND. VL > 0) THEN
       DO I1 = 1, dim_block
          HAM(I1,I1) = HAM(I1,I1) + H_4b_pars(11)*VL*VL*W2MAT(I1,I1)
       ENDDO
       DO I1 = 1, dim_block-1
          HAM(I1+1,I1)=HAM(I1+1,I1) + H_4b_pars(11)*VL*VL*W2MAT(I1+1,I1)
          HAM(I1,I1+1)=HAM(I1,I1+1) + H_4b_pars(11)*VL*VL*W2MAT(I1,I1+1)
       ENDDO
    ENDIF
    !
    !     BUILDING W^2·n^2 + n^2·W^2 BLOCK
    !
    IF (H_4b_pars(12) /= 0) THEN
       !
       !     HAM MULTIPLYING (W^2·n^2 + n^2·W^2) TIMES P45 PARAMETER
       DO I1 = 1, dim_block
          V2 = REAL(N_val - (2*I1-2+MOD(N_val-L,2)), DP)
          HAM(I1,I1) = HAM(I1,I1) + H_4b_pars(12)*2.0_DP*V2*V2*W2MAT(I1,I1)
       ENDDO
       DO I1 = 1, dim_block - 1
          V2 = REAL(N_val - (2*I1-2+MOD(N_val-L,2)), DP)
          HAM(I1+1,I1) = HAM(I1+1,I1) + &
               H_4b_pars(12)*2.0_DP*(V2*V2-2.0_DP*V2+2.0_DP)*W2MAT(I1+1,I1)
          HAM(I1,I1+1) = HAM(I1,I1+1) + &
               H_4b_pars(12)*2.0_DP*(V2*V2-2.0_DP*V2+2.0_DP)*W2MAT(I1,I1+1)
       ENDDO
    ENDIF
    !
    !     BUILDING W2·WBAR2 + WBAR2·W2 BLOCK (HAS TO BE EVALUATED BEFORE W4 BLOCK)
    !
    IF (H_4b_pars(14) /= 0) THEN
       !
       CALL SO3SO3BARBUILD(W2MAT, W4MAT, W2W2BARMAT, N_val, L, Iprint)
       !
       !     HAM MULTIPLYING W2·WBAR2 + WBAR2·W2 TIMES P47 PARAMETER
       DO I1 = 1, dim_block
          HAM(I1,I1) = HAM(I1,I1) + H_4b_pars(14)*W2W2BARMAT(I1,I1)
       ENDDO
       DO I1 = 1, dim_block - 1
          HAM(I1+1,I1) = HAM(I1+1,I1) + H_4b_pars(14)*W2W2BARMAT(I1+1,I1)
          HAM(I1,I1+1) = HAM(I1,I1+1) + H_4b_pars(14)*W2W2BARMAT(I1,I1+1)
       ENDDO
       DO I1 = 1, dim_block - 2
          HAM(I1+2,I1) = HAM(I1+2,I1) + H_4b_pars(14)*W2W2BARMAT(I1+2,I1)
          HAM(I1,I1+2) = HAM(I1,I1+2) + H_4b_pars(14)*W2W2BARMAT(I1,I1+2)
       ENDDO
    ENDIF
    !
    !     BUILDING W4 BLOCK
    !
    IF (H_4b_pars(13) /= 0) THEN
       CALL SO32CASBUILD(W2MAT,W4MAT,N_val,L, Iprint)
       !
       !     HAM MULTIPLYING P^2 TIMES P46 PARAMETER
       DO I1 = 1, dim_block
          HAM(I1,I1) = HAM(I1,I1) + H_4b_pars(13)*W4MAT(I1,I1)
       ENDDO
       DO I1 = 1, dim_block-1
          HAM(I1+1,I1) = HAM(I1+1,I1) + H_4b_pars(13)*W4MAT(I1+1,I1)
          HAM(I1,I1+1) = HAM(I1,I1+1) + H_4b_pars(13)*W4MAT(I1,I1+1)
       ENDDO
       DO I1 = 1, dim_block-2
          HAM(I1+2,I1) = HAM(I1+2,I1) + H_4b_pars(13)*W4MAT(I1+2,I1)
          HAM(I1,I1+2) = HAM(I1,I1+2) + H_4b_pars(13)*W4MAT(I1,I1+2)
       ENDDO
    ENDIF
    !
    !     DIAGONAL INTERACTIONS
    DO I = 1, dim_block
       !
       V2 = REAL(N_val - (2*I-2+MOD(N_val-L,2)), DP)
       !     
       HAM(I,I) = HAM(I,I) +  H_4b_pars(1)*V2 + &
                                !     
            H_4b_pars(2)*V2*V2 + &
                                !     
            H_4b_pars(3)*VL*VL + &
                                !     
            H_4b_pars(5)*V2*V2*V2 + &
                                !     
            H_4b_pars(6)*V2*VL*VL + &
                                !     
            H_4b_pars(8)*V2*V2*V2*V2 + &
                                !     
            H_4b_pars(9)*V2*V2*VL*VL + &
                                !     
            H_4b_pars(10)*VL*VL*VL*VL 
       !     
    ENDDO
    !
    !
    IF (IPRINT > 2) WRITE(*,*)'HAMILTONIAN BUILDING SUBROUTINE ENDS HERE'
    !
    RETURN
  END SUBROUTINE HBLDU3GEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  SUBROUTINE SO3CASBUILD(W2BLOCK, N_val, L, Iprint)
    !
    !     BUILDING THE SO(3) CASIMIR W^2 MATRIX IN
    !     THE CYLINDRICAL OSCILLATOR BASIS. 
    !     SO(3) = {L,D+,D-}
    !
    !     PROGRAM OF THE U(3) ALGEBRAIC MODEL FOR BENDING DYNAMICS
    !     
    !     INPUT 
    !     N_val        : U(3) REPRESENTATION (TOTAL NUMBER OF U(3) BOSONS)
    !     L            : VIBRATIONAL ANGULAR MOMENTUM
    !
    !     OUTPUT
    !     W2BLOCK : LOCAL BASIS MATRIX
    !     
    !
    !     
    !     by Currix TM
    !
    IMPLICIT NONE
    !     
    !    DEFINITION OF VARIABLES
    INTEGER(KIND = I4B), INTENT(IN)  ::  N_val, L
    INTEGER(KIND = I4B), INTENT(IN)  ::  Iprint
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2BLOCK
    !
    !
    !     LOCAL VARIABLES                                     
    INTEGER(KIND = I4B) :: I1, DIMB
    REAL(KIND = DP) :: VN2, VNN, VL, VTMP
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SO3CASBUILD STARTS HERE'
    !
    DIMB = (N_val - MOD(N_val - L, 2) - L)/2 + 1
    !
    VNN = REAL(N_val, DP)
    VL = REAL(L, DP)
    !
    !     NON-DIAGONAL PART
    DO I1 = 1, DIMB-1
       VN2 = REAL(N_val - MOD(N_val - L,2) - 2*(I1-1), DP)
       VTMP = DSQRT((VNN-VN2+2.0D0)*(VNN-VN2+1.0D0)*(VN2+VL)*(VN2-VL))
       W2BLOCK(I1+1,I1) = -VTMP
       W2BLOCK(I1,I1+1) = -VTMP
    ENDDO
    !
    !     DIAGONAL PART
    DO I1 = 1, DIMB                     
       VN2 = REAL(N_val - MOD(N_val - L,2) - 2*(I1-1), DP)
       W2BLOCK(I1,I1)=(VNN-VN2)*(VN2+2.0D0)+(VNN-VN2+1.0D0)*VN2+VL*VL
    ENDDO
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SO3CASBUILD ENDS HERE'
    !
    RETURN
    !
  END SUBROUTINE SO3CASBUILD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SO3SO3BARBUILD(W2BLOCK,WB2BLOCK,W2WB2BLOCK,NN,L, Iprint)
    !
    !     SUBROUTINE THAT BUILD THE SO(3) CASIMIR (W^2·Wbar^2 + Wbar^2·W2)/2 MATRIX IN
    !     THE CYLINDRICAL OSCILLATOR BASIS. 
    !     SO(3) = {L,D+,D-}
    !     SObar(3) = {L,R+,R-}
    !
    !     PROGRAM OF THE U(3) ALGEBRAIC MODEL FOR BENDING DYNAMICS
    !     
    !     INPUT 
    !     W2BLOCK   ; LOCAL BASIS MATRIX OF W2
    !     WB2BLOCK  : ARRAY WHERE TO BUILD THE SObar(3) CASIMIR
    !     NN        : U(3) REPRESENTATION (TOTAL NUMBER OF U(3) BOSONS)
    !     L         : VIBRATIONAL ANGULAR MOMENTUM
    !
    !     OUTPUT
    !     W2WB2BLOCK : LOCAL BASIS MATRIX OF (W^2·Wbar^2 + Wbar^2·W2)/2
    !     
    !
    !     
    !     by Currix TM
    !
    !     
    IMPLICIT NONE
    !     
    !  DEFINITION OF VARIABLES     
    INTEGER(KIND = I4B), INTENT(IN)  :: NN, L
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2BLOCK, WB2BLOCK
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2WB2BLOCK
    INTEGER(KIND = I4B), INTENT(IN)  :: Iprint
    !
    !     TEMPORAL VARIABLES                                     
    INTEGER(KIND = I4B) :: I1, DIMB 
    REAL(KIND = DP) :: VN2, VNN, VL, VTMP
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SO3SO3BARBUILD STARTS HERE'
    !
    DIMB = (NN - MOD(NN-L,2) - L)/2 + 1
    !
    !     INITIALIZING
    W2BLOCK = 0.0_DP
    WB2BLOCK = 0.0_DP
    W2WB2BLOCK = 0.0_DP
    !
    VNN = REAL(NN, DP)
    VL = REAL(L, DP)
    !
    !    NON-DIAGONAL PART
    DO I1 = 1, DIMB-1
       VN2 = DFLOAT(NN - MOD(NN-L,2) - 2*(I1-1))
       VTMP = DSQRT((VNN-VN2+2.0_DP)*(VNN-VN2+1.0_DP)*(VN2+VL)*(VN2-VL))
       W2BLOCK(I1+1,I1) = -VTMP
       W2BLOCK(I1,I1+1) = -VTMP
       WB2BLOCK(I1+1,I1) = VTMP
       WB2BLOCK(I1,I1+1) = VTMP
    ENDDO
    !
    !     DIAGONAL PART
    DO I1 = 1, DIMB                     
       VN2 = DFLOAT(NN - MOD(NN-L,2) - 2*(I1-1))
       W2BLOCK(I1,I1)=(VNN-VN2)*(VN2+2.0_DP)+(VNN-VN2+1.0_DP)*VN2+VL*VL
       WB2BLOCK(I1,I1) = W2BLOCK(I1,I1)
    ENDDO
    !
    !
    W2WB2BLOCK = 0.5_DP*(MATMUL(W2BLOCK, WB2BLOCK) + MATMUL(WB2BLOCK, W2BLOCK))
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SO3SO3BARBUILD ENDS HERE'
    !
    RETURN
    !
  END SUBROUTINE SO3SO3BARBUILD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SO32CASBUILD(W2BLOCK, W4BLOCK, NN, L, Iprint)
    !
    !     SUBROUTINE THAT BUILD THE SO(3) CASIMIR (W^2)^2 MATRIX IN
    !     THE CYLINDRICAL OSCILLATOR BASIS. 
    !     SO(3) = {L,D+,D-}
    !
    !
    !     PROGRAM OF THE U(3) ALGEBRAIC MODEL FOR BENDING DYNAMICS
    !     
    !     INPUT 
    !     W2BLOCK   : LOCAL BASIS MATRIX OF THE PAIRING OPERATOR
    !     NN        : U(3) REPRESENTATION (TOTAL NUMBER OF U(3) BOSONS)
    !     L         : VIBRATIONAL ANGULAR MOMENTUM
    !     Iprint    : Control level of output
    !
    !     OUTPUT
    !     W4BLOCK   : LOCAL BASIS MATRIX OF W^4
    !     
    !
    !     
    !     by Currix TM
    !
    !     
    IMPLICIT NONE
    !     
    !     DEFINITION OF VARIABLES     
    INTEGER(KIND = I4B), INTENT(IN)  :: NN, L, Iprint
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: W2BLOCK
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W4BLOCK
    !
    !     TEMPORAL VARIABLES                                     
    INTEGER(KIND = I4B) ::  DIMB
    !
    IF (Iprint > 2) WRITE(*,*) 'SUBROUTINE SO32CASBUILD STARTS HERE'
    !
    DIMB = (NN - MOD(NN-L,2) - L)/2 + 1
    !
    !     INITIALIZING
    W4BLOCK = 0.0_DP
    !
    !     MATRIX PRODUCT
    W4BLOCK = MATMUL(W2BLOCK, W2BLOCK)
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SO32CASBUILD ENDS HERE'
    !
    RETURN
  END SUBROUTINE SO32CASBUILD

END MODULE u3_2dvm_mod
