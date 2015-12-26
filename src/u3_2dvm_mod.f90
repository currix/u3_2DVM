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
  ! Eigenvalue Matrix
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
  ! General Hamiltonian Parameters
  INTEGER(KIND = I4B), PARAMETER :: n_ham_param = 4
  REAL(KIND = DP), DIMENSION(1:n_ham_param) :: Ham_parameter 
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
          ! Non-Diagonal contribution (B_0)
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
    !
    n_elem_diag = (Nv - wv) * &
         ((wv-lv+2.0_DP)*(wv-lv+1.0_DP) + (wv+lv+2.0_DP)*(wv+lv+1.0_DP)) / &
         (2.0_DP*(2.0_DP*wv+1.0_DP)*(2.0_DP*wv+3.0_DP)) + &
         (Nv + wv + 1.0_DP) * &
         ((wv+lv)*(wv+lv-1.0_DP) + (wv-lv)*(wv-lv-1.0_DP)) / &
         (2.0_DP*(2.0_DP*wv-1.0_DP)*(2.0_DP*wv+1.0_DP))
    !
  END FUNCTION n_elem_diag
  !
  FUNCTION B0(Nv, wv, L_val)
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
    B0 = SQRT( &
         (Nv - wv)*(Nv + wv + 3.0_DP) * &
         (wv-lv+2.0_DP)*(wv+lv+2.0_DP)*(wv+lv+1.0_DP)*(wv-lv+1.0_DP) / &
         ((2.0_DP*wv+1.0_DP)*(2.0_DP*wv+3.0_DP)**2*(2.0_DP*wv+5.0_DP)) &
         )
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
  SUBROUTINE SAVE_EIGENV_COMPONENTS(N_val, L_val, xi, dim_block, Basis, Ham_U3_mat)
    !
    ! Subroutine to save the components of the eigenvectors 
    ! of the U(3) 2DVM Model Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(3) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum    
    !
    REAL(KIND = DP),  INTENT(IN) :: xi ! Control Parameter
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
    ELSE IF ( dim_block < 100) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I2,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( dim_block < 1000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I3,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( dim_block < 10000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I4,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE
       WRITE(output_filename, '("eigvec_",A,"_N",I6,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ENDIF
    !
    OPEN(UNIT = 76, FILE = output_filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !
    WRITE(UNIT=76, FMT=*) "# N = ", N_val, "; L = ", L_val, ";  xi = ", xi, " ; ", Basis,  " basis"
    !
    DO basis_index = 1, dim_block
       WRITE(UNIT=76, FMT=*) Ham_U3_mat(basis_index, 1:dim_block)
    ENDDO
    !
  END SUBROUTINE SAVE_EIGENV_COMPONENTS
  !
END MODULE u3_2dvm_mod
