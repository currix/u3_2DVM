MODULE FIT_2DVM
  !
  USE nrtype
  USE defparam_2DVM
  !
  USE u3_2dvm_mod
  !
  ! Lapack 95
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR
  !
  IMPLICIT NONE
  !     
  !     IF .T. BENT ELSE LINEAR CASE
  LOGICAL :: BENT
  !
  !
  !     IF .T. THEN ENERGIES REFERRED THE MINIMUM OF EACH L BLOCK, 
  !     ELSE REFERRED TO THE GROUND STATE (MINIMUM L = 0)
  LOGICAL :: EMINL
  !
  !     IF .T. DISPLAY ONLY RESIDUALS FOR PYTHON INTEGRATION
  LOGICAL :: DIS_RES
  !
  !     NAME OF EXPERIMENTAL DATA FILE
  CHARACTER(LEN = 64) ::  exp_data_file
  !
  !     MAXIMUM VIBRATIONAL ANGULAR MOMENTUM AND EXPERIMENTAL NUMBER OF QUANTA 
  !     TO BE CONSIDERED
  INTEGER(KIND = I4B) :: LMAX, VMAX
  !
  !     EXPERIMENTAL ENERGIES AND ERRORS MATRICES
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: EXPEN, EXPERR
  !     
  !     EXPERIMENTAL ASSIGNMENTS AND POINTER:  EXPAS(n,1) v_n EXPAS(n,2) l_n 
  INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE :: EXPAS
  INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE ::  LXP
  !
  !     POINTER TO COMPARE WITH EXPERIMENTAL RESULTS
  INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE :: SPNTEX
  !
  !     STATISTICAL VARIABLES
  REAL(KIND = DP) ::  CHSQP, CHSQT
  !
  !
  !     TOTAL NUMBER OF EXPERIMENTAL DATA
  INTEGER(KIND = I4B) :: TOTDAT
  !
  !     VECTOR WITH n's(v's) ASSIGNED TO COMPUTED EIGENVECTORS 
  INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE :: BLAS
  !
  !
CONTAINS
  !
  !
  FUNCTION CHI_SQRE(H_pars)
    !
    ! H_pars VECTOR:
    ! P11,P21,P22,P23,P31,P32,P33,P41,P42,P43,P44,P45,P46,P47
    !
    IMPLICIT NONE
    !     DEFINITION OF VARIABLES
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: H_pars
    REAL(KIND = DP) :: CHI_SQRE
    !
    !     LOCAL VARIABLES
    !
    INTEGER(KIND = I4B) :: I, JMIN, L, NDAT
    !
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Eigenval_vector ! Hamiltonian Eigenvalues
    !
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_matrix 
    !
    !     COMPUTED W2 AND SQUARED W2 (SO(3) CASIMIR) MATRICES
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2MAT, W4MAT
    !
    !     COMPUTED W2·WBAR2+WBAR2·W2 (SO(3) CASIMIR x barSO(3) CASIMIR) MATRIX
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: W2W2BARMAT
    !
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: AVEZERO
    !
    REAL(KIND = DP) :: EMIN
    !
    !
    INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE :: VEXPAS
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: VEXPEN, VEXPERR
    !     LOGICAL VARIABLE THAT CONTROL IF ALL EIGENVECTORS ARE ASSIGNED
    !     HAS TO BE MIFLAG = 3 AND IPRINT > 0
    LOGICAL :: ASSIGNALL 
    !     LEVELS IN THE EXPERIMENTAL ENERGY LIST BUT NOT INCLUDED IN CHI^2
    INTEGER(KIND = I4B) :: NOT_FIT_LEVELS
    !
    !
    IF (IPRINT > 2) THEN
       WRITE(*,*) 'STARTING CHISQRE CALCULATION'
    ENDIF
    !
    !
    !     TESTS     
    IF (LMAX > N_val) STOP 'L_max > N, SAYONARA BABY'
    !
    ! ONE BODY
    !     n
    !   HPAR(1) = P11
    !
    !  TWO BODY
    !      n^2
    !       HPAR(2) = P21
    !      L^2
    !       HPAR(3) = P22
    !      W^2
    !       HPAR(4) = P23
    ! 
    !  THREE BODY
    !      n^3
    !       HPAR(5) = P31
    !      n·l^2
    !       HPAR(6) = P32
    !      n·W^2 + W^2·n
    !       HPAR(7) = P33
    ! 
    !  FOUR BODY
    !      n^4
    !       HPAR(8) = P41
    !      n^2·l^2
    !       HPAR(9) = P42
    !      l^4
    !       HPAR(10) = P43
    !      W^2·l^2
    !       HPAR(11) = P44
    !      n^2·W^2 + W^2·n^2
    !       HPAR(12) = P45
    !      W^4 
    !       HPAR(13) = P46
    !      (W^2·Wb^2 + Wb^2·W^2)/2
    !       HPAR(14) = P47
    !     
    !     
    !     READING EXPERIMENTAL DATA FILE      
    CALL RDXENRU3(BENT, exp_data_file, LMAX, VMAX, EXPEN, EXPERR, EXPAS, LXP, Iprint)
    !
    !     
    !     INITIALIZE CHI SQUARE
    CHSQT = 0.0D0
    !     
    !     INITIALIZE TOTAL NUMBER OF DATA
    TOTDAT = 0
    !     
    !     INITIALIZE EMIN
    EMIN = 0.0D0
    !     
    !     MAIN LOOP FOR CHISQRE CALCULATION 
    DO L = 0, LMAX
       !     
       IF (IPRINT > 3) THEN
          WRITE(*,*)
          WRITE(*,*) 'N_val = ',N_val,' L = ',L
          WRITE(*,*)
       ENDIF
       !     SELECT EXPERIMENTAL DATA FOR VIBRATIONAL ANGULAR MOMENTUM L
       CALL SLCTLU3(L,EXPEN,EXPERR,EXPAS,LXP,VEXPAS,VEXPEN,VEXPERR,NDAT, iprint)
       !
       !"!"!"!"!"!"!"!"!"!"PRINT*, "NDAT", NDAT
       !"!"!"!"!"!"!"!"!"!"PRINT*, "EXP_EN", VEXPEN
       !"!"!"!"!"!"!"!"!"!"PRINT*, VEXPAS
       !
       ! Block dimension
       dim_block = (N_val - MOD(N_val - L, 2) - L)/2 + 1
       !
       ! MATRIX ALLOCATION AND INITIALIZATION
       ALLOCATE(Ham_matrix(1:dim_block,1:dim_block), STAT = IERR)
       IF (IERR /= 0) STOP 'ERROR ALLOCATING HAM MATRIX'
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
       ALLOCATE(Eigenval_vector(1:dim_block), STAT = Ierr)
       IF (Ierr /= 0) STOP 'Eigenval_vector ALLOCATION ERROR'
       !
       ALLOCATE(AVEZERO(1:dim_block), STAT = Ierr)
       IF (Ierr /= 0) STOP 'AVEZERO ALLOCATION ERROR'
       !
       !     BUILD HAMILTONIAN MATRIX OF ALGEBRAIC HAMILTONIAN
       CALL HBLDU3GEN(N_val, L, HAM_matrix, W2MAT, W4MAT, W2W2BARMAT, IPRINT)  
       !     
       !     DIAGONALIZE HAMILTONIAN MATRIX
       CALL LA_SYEVR(A=HAM_matrix, W=Eigenval_vector, JOBZ='V', UPLO='U')
       !
       !    ASSIGN  COMPUTED DATA TO LOCAL BASIS STATES
       IF (IPRINT > 0) THEN 
          ASSIGNALL = .TRUE.
       ELSE
          ASSIGNALL = .FALSE.
       ENDIF
       !
       CALL ASSGNU3(BENT, ASSIGNALL, N_val, L, dim_block, W2MAT, W4MAT, W2W2BARMAT, &
            HAM_matrix, Eigenval_vector, VEXPAS, NDAT, BLAS, Iprint)
       !
       !     REFER ENERGIES TO GROUND L=0 STATE
       !
       IF (EMINL) THEN
          JMIN = 0
          DO I = 1, dim_block
             IF (BENT) THEN 
                IF (BLAS(I).EQ.0) JMIN = I
             ELSE
                IF (BLAS(I).EQ.L) JMIN = I
             ENDIF
          ENDDO
          EMIN = Eigenval_vector(JMIN)
          DO I = 1, dim_block
             AVEZERO(I) = HAM_matrix(I,JMIN)
          ENDDO
          !
          IF (IPRINT > 2) WRITE(*,*) 'J(E_MIN) = ', JMIN
          !
       ELSE IF (L == 0) THEN
          JMIN = 0
          DO I = 1, dim_block
             IF (BLAS(I) == 0) JMIN = I
          ENDDO
          EMIN = Eigenval_vector(JMIN)
          AVEZERO = Ham_matrix(:,JMIN)
          !    
          IF (IPRINT >= 2) WRITE(*,*) 'J(E_MIN) = ', JMIN
          !     
       ENDIF
       !
       Eigenval_vector = Eigenval_vector - EMIN
       !
       !     
       !     COMPARE WITH EXPERIMENT AND COMPUTE STATISTICAL PARAMETERS
       CALL CPXPDTU3(N_val,L,VEXPAS,Eigenval_vector,BLAS,VEXPEN,VEXPERR,NDAT,SPNTEX,CHSQP,IPRINT)
       !     
       !     DISPLAY RESULTS
       IF (IPRINT > 0) CALL DSDTU3(BENT, N_val, L, Eigenval_vector, &
            VEXPEN, VEXPERR, NDAT, VEXPAS, HAM_matrix, BLAS, SPNTEX, CHSQP, IPRINT)
       !    
       !     EXCLUDE NOT FIT LEVELS
       NOT_FIT_LEVELS = 0
       DO I = 1, NDAT
          IF (VEXPERR(I) == 0) NOT_FIT_LEVELS = NOT_FIT_LEVELS + 1
       ENDDO
       !
       !     
       TOTDAT = TOTDAT + NDAT - NOT_FIT_LEVELS
       !     
       CHSQT = CHSQT + CHSQP
       !
       ! DEALLOCATE MATRICES
       DEALLOCATE(Ham_matrix, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Ham_matrix DEALLOCATION ERROR'
       !  
       ! DEALLOCATE(Ham2, STAT = Ierr)
       ! IF (Ierr /= 0) STOP 'Ham2 DEALLOCATION ERROR'
       !  
       DEALLOCATE(W2MAT, STAT = IERR)
       IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2MAT MATRIX'
       !     
       DEALLOCATE(W4MAT, STAT = IERR)
       IF (IERR /= 0) STOP 'ERROR DEALLOCATING W4MAT MATRIX'
       !     
       DEALLOCATE(W2W2BARMAT, STAT = IERR)
       IF (IERR /= 0) STOP 'ERROR DEALLOCATING W2W2BARMAT MATRIX'
       !
       DEALLOCATE(Eigenval_vector, STAT = Ierr)
       IF (Ierr /= 0) STOP 'Eigenval_vector DEALLOCATION ERROR'
       !
       DEALLOCATE(AVEZERO, STAT = Ierr)
       IF (Ierr /= 0) STOP 'AVEZERO DEALLOCATION ERROR'
       !
    ENDDO
    !
    CHI_SQRE = CHSQT
    !
    IF (IPRINT > 0) THEN
       WRITE(*,*) "NUMBER OF DATA = ", TOTDAT, "CHI2 = ", CHI_SQRE, &
            " SDEV = ", SQRT(CHI_SQRE/TOTDAT)
    ENDIF
    !
    IF (IPRINT > 0)  THEN
       WRITE(*,10) H_PARS(1)
       WRITE(*,20) H_PARS(2),H_PARS(3), H_PARS(4)
       WRITE(*,30) H_PARS(5),H_PARS(6),H_PARS(7)
       WRITE(*,40) H_PARS(8),H_PARS(9),H_PARS(10)
       WRITE(*,41) H_PARS(11),H_PARS(12),H_PARS(13)
       WRITE(*,42) H_PARS(14)
    ENDIF
    !
10  FORMAT (1X,"ONE BODY :: ", " P11 = ", G15.5)
20  FORMAT (1X,"TWO BODY :: ", " P21 = ", G15.6, &
         " P22 = ", G15.6, " P23 = ", G15.6)
30  FORMAT (1X,"THREE BODY :: ", " P31 = ", G15.6, &
         " P32 = ", G15.6, " P33 = ", G15.6)
40  FORMAT (1X,"FOUR BODY (i):: ", " P41 = ", G15.6,  &
         " P42 = ", G15.6, " P43 = ", G15.6)
41  FORMAT (1X,"FOUR BODY (ii):: ", " P44 = ", G15.6, &
         " P45 = ", G15.6, " P46 = ", G15.6)
42  FORMAT (1X,"FOUR BODY (iii):: ", " P47 = ", G15.6)
    !
    RETURN 
    !
  END FUNCTION CHI_SQRE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RDXENRU3(BENT,DTFL,LMAX,VMAX,EXPEN,EXPERR,EXPAS,LXPOIN, IPRINT)
    !
    !     SUBROUTINE THAT READS EXPERIMENTAL DATA (ENERGIES, ERRORS 
    !     ASSIGNMENTS IN THE U(3) MODEL FOR BENDING VIBRATIONS.
    !
    !     FORMAT OF THE INPUT FILE:
    !     -------------------------
    !     N   # NUMBER OF DATA
    !     E1   ERR1   V1  L1
    !     E2   ERR2   V2  L2
    !             .................
    !
    !     EN   ERRN   VN2 LN
    !
    !     THE VALUES HAVE NOT TO BE IN ANY DETERMINED ORDER 
    !     AND IN CASE ERRORS ARE UNKNOWN ALL MUST BE GIVEN 
    !     A 1.0 VALUE.
    !
    !     INPUT
    !     BENT : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !     DTFL : FILE WITH EXPERIMENTAL DATA SET
    !     LMAX : MAXIMUM VIBRATIONAL ANGULAR MOMENTUM LABEL
    !     VMAX : MAXIMUM NUMBER OF QUANTA OF VIBRATION
    !
    !     OUTPUT
    !     EXPEN  : VECTOR WITH EXPERIMENTAL ENERGIES
    !     EXPERR : VECTOR WITH EXPERIMENTAL ENERGY ERRORS
    !     EXPAS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS
    !     LXPOIN : POINTER LOCATING L VALUES
    !     
    !
    !     by Currix TM
    !     
    IMPLICIT NONE
    !     
    !     DEFINITION OF VARIABLES
    !     
    LOGICAL, INTENT(IN) :: BENT
    CHARACTER(LEN = 64), INTENT(IN) :: DTFL
    INTEGER(KIND = I4B), INTENT(IN) ::  LMAX, VMAX, IPRINT
    !
    INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: EXPAS
    INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: LXPOIN
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: EXPEN, EXPERR
    !
    !     TEMPORAL VARIABLES 
    INTEGER(KIND = I4B) ::  I, J, M, INDX, JT, N_total, COUNT, LTEMP, Ierr
    REAL(KIND = DP) :: TMP1, TMP2
    INTEGER(KIND = I4B), DIMENSION(1:2) :: TMP
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE READENERG'
    !
    !
    !    READ ALL EXPERIMENTAL DATA SELECTING RIGHT V'S AND L
    OPEN(UNIT=30,FILE=DTFL,STATUS='OLD')
    !
    ! Total number of experimental data
    READ(30,*) N_total
    !
    ! ALLOCATE AND INITIALIZE MATRICES
    !
    ALLOCATE(EXPAS(1:N_total, 1:2), STAT = Ierr)
    IF (Ierr /= 0) STOP 'EXPAS ALLOCATION ERROR'
    EXPAS = 0
    ALLOCATE(LXPOIN(0:LMAX+1), STAT = Ierr)
    IF (Ierr /= 0) STOP 'LXPOIN ALLOCATION ERROR'
    LXPOIN = 0
    ALLOCATE(EXPEN(1:N_total), STAT = Ierr)
    IF (Ierr /= 0) STOP 'EXPEN ALLOCATION ERROR'
    EXPEN = 0.0_DP
    ALLOCATE(EXPERR(1:N_total), STAT = Ierr)
    IF (Ierr /= 0) STOP 'EXPERR ALLOCATION ERROR'
    EXPERR = 0.0_DP
    !
    COUNT = 0
    DO I = 1, N_total
       READ(30,*) TMP1, TMP2, (TMP(J), J = 1, 2)
       !  
       !     TEST LINEAR LABELS
       IF (.NOT.BENT .AND. MOD(TMP(1),2) /= MOD(TMP(2),2)) THEN
          WRITE(*,*) "Wrong parity: LINEAR CASE WITH n = ",TMP(1), " l = ", TMP(2)
          STOP 'ERROR READING EXPERIMENTAL ENERGIES'
       ENDIF
       !
       IF (TMP(2) <= LMAX .AND. TMP(1) <= VMAX) THEN
          COUNT = COUNT+1
          EXPEN(COUNT) = TMP1
          EXPERR(COUNT) = TMP2
          EXPAS(COUNT,:) = TMP
          !
          LXPOIN(TMP(2)+1) = LXPOIN(TMP(2)+1) + 1
          !"!"!"!"!"!"!"!"!"!"print*, "tmp", tmp1, tmp(2), "lxpoin", lxpoin
       ENDIF
    ENDDO
    !
    !"!"!"!"!"!"!"!"!"!"print*, "lxpoin", lxpoin
    !
    ! WHERE DO L EXPERIMENTAL VALUES START IN THE POINTER
!!!!!!    LXPOIN(0) = 1
!!!!!!    LXPOIN(LMAX + 1) = 0
    DO I = 1, LMAX+1
       LXPOIN(I) = LXPOIN(I)+LXPOIN(I-1)
    ENDDO
    !
    !"!"!"!"!"!"!"!"!"!"print*, "lxpoin", lxpoin
    !
    !    ORDERING DATA ACCORDING TO VIBRATIONAL ANGULAR MOMENTUM L
    !     
    DO J = 2, COUNT
       LTEMP = EXPAS(J,2)
       !     
       TMP1 = EXPEN(J)
       TMP2 = EXPERR(J)
       TMP = EXPAS(J,:)
       !
       JT = J
       !
       ido: DO I = J-1, 1, -1
          IF (EXPAS(I,2) > LTEMP) THEN
             EXPEN(JT) = EXPEN(I)
             EXPERR(JT) = EXPERR(I)
             EXPAS(JT,:) = EXPAS(I,:)
             !
             JT = I
          ELSE
             EXIT ido
          ENDIF
       ENDDO ido
       EXPEN(JT) = TMP1
       EXPERR(JT) = TMP2
       EXPAS(JT,:) = TMP
    ENDDO
    !     
    !     ORDERING IN INCREASING ENERGY FOR EACH L 
    !     
    DO I = 0, LMAX
       !     INDEX-1 OF THE FIRST ELEMENT ELEMENT OF THE I-TH VIBRATIONAL ANGULAR MOMENTUM
       INDX = LXPOIN(I)
       !
       DO J = 2, LXPOIN(I+1)-LXPOIN(I)
          !     
          TMP1 = EXPEN(INDX+J)
          TMP2 = EXPERR(INDX+J)
          TMP = EXPAS(INDX+J,:)
          !     
          JT = J
          !
          !"!"!"!"!"!"!"!"!"!"print*, "j do", jt, tmp1, tmp2, tmp
          !     
          mdo: DO M = J-1, 1, -1
             !"!"!"!"!"!"!"!"!"!"print*, m
             IF (EXPEN(INDX+M) > TMP1) THEN
                EXPEN(INDX+JT) = EXPEN(INDX+M)
                EXPERR(INDX+JT) = EXPERR(INDX+M)
                EXPAS(INDX+JT,:) = EXPAS(INDX+M,:)
                !
                JT = M
             ELSE
                EXIT mdo
             ENDIF
          ENDDO mdo
          EXPEN(INDX+JT) = TMP1
          EXPERR(INDX+JT) = TMP2
          EXPAS(INDX+JT,:) = TMP
       ENDDO
    ENDDO
    !
    IF (IPRINT > 2) THEN
       WRITE(*,*)
       WRITE(*,*) 'EXPERIMENTAL ENERGIES'
       WRITE(*,*)
       DO I = 0, LMAX
          WRITE(*,*) 'VIBRATIONAL ANGULAR MOMENTUM', I
          WRITE(*,*) 'LXPOIN', LXPOIN(I+1), LXPOIN(I)
          DO J = 1, LXPOIN(I+1) - LXPOIN(I)
             WRITE(*,*) EXPAS(LXPOIN(I)+J,:), EXPEN(LXPOIN(I)+J), EXPERR(LXPOIN(I)+J)
          ENDDO
       ENDDO
       WRITE(*,*)
    ENDIF
    !
    CLOSE(30)
    !
    RETURN
    !
  END SUBROUTINE RDXENRU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SLCTLU3(L, EXPEN, EXPERR, EXPAS, LXPOIN, &
       VEXPAS, VEXPEN, VEXPERR, NDAT, Iprint)
    !      
    !      
    !      SUBROUTINE THAT SELECTS THE EXPERIMENTAL DATA FOR 
    !      VIBRATIONAL ANGULAR MOMENTUM L AMONG EXPERIMENTAL 
    !      ENERGIES, ERRORS AND ASSIGNMENTS IN THE U(3) MODEL.
    !      
    !      INPUT
    !      L      : VIBRATIONAL ANGULAR MOMENTUM
    !      EXPEN  : VECTOR WITH ALL EXPERIMENTAL ENERGIES
    !      EXPERR : VECTOR WITH ALL EXPERIMENTAL ENERGY ERRORS
    !      EXPAS  : MATRIX WITH ALL EXPERIMENTAL ASSIGNMENTS FORMAT:(V L)
    !      LXPOIN : POINTER TO LOCATE THE DIFFERENT VALUES OF L
    !      IPRINT : VERBOSITY CONTROL
    ! 
    !      OUTPUT
    !      VEXPAS : EXPERIMENTAL ASSIGNMENTS FOR POLYAD (V,L) FORMAT:(V L)
    !      VEXPEN : EXPERIMENTAL ENERGIES FOR POLYAD (V,L)
    !      VEXPERR: EXPERIMENTAL ERRORS FOR POLYAD (V,L)
    !      NDAT   : NUMBER OF EXPERIMENTAL DATA FOR VIBRATIONAL ANGULAR MOMENTUM L
    ! 
    !      by Currix TM
    !      
    IMPLICIT NONE
    !
    !     DEFINITION OF VARIABLES
    INTEGER(KIND = I4B), INTENT(IN) :: L
    REAL(KIND = DP), DIMENSION(:), INTENT(IN)  :: EXPEN, EXPERR
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: EXPAS
    INTEGER(KIND = I4B), DIMENSION(0:),  INTENT(IN) :: LXPOIN
    !
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: VEXPEN, VEXPERR
    INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: VEXPAS
    INTEGER(KIND = I4B), INTENT(OUT) :: NDAT
    !     
    INTEGER(KIND = I4B), INTENT(IN) :: IPRINT
    !
    !     TEMPORAL VARIABLES                                     
    INTEGER(KIND = I4B) :: I, J, INDX, NTOT, IERR
    !     
    !     
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SELECT VL'
    !     
    NDAT = 0
    !
    !"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"!"
    !!"!"!"!"!"!"!"print*, "l", l
    !"!"!"!"!"!"!"!"print*, "lxpoin", lxpoin
    !
    IF (ALLOCATED(VEXPEN)) THEN
       DEALLOCATE(VEXPEN, STAT = IERR)
       IF (IERR /= 0) STOP "ERROR DEALLOCATING VEXPEN"
       DEALLOCATE(VEXPERR, STAT = IERR)
       IF (IERR /= 0) STOP "ERROR DEALLOCATING VEXPERR"
       DEALLOCATE(VEXPAS, STAT = IERR)
       IF (IERR /= 0) STOP "ERROR DEALLOCATING VEXPAS"
    ENDIF
    !
    !  NUMBER OF EXPERIMENTAL DATA
    NTOT = LXPOIN(L+1) - LXPOIN(L)
    !
    ALLOCATE(VEXPEN(1:NTOT), STAT = IERR)
    IF (IERR /= 0) STOP "ERROR ALLOCATING VEXPEN"
    !
    ALLOCATE(VEXPERR(1:NTOT), STAT = IERR)
    IF (IERR /= 0) STOP "ERROR ALLOCATING VEXPERR"
    !
    ALLOCATE(VEXPAS(1:NTOT,1:2), STAT = IERR)
    IF (IERR /= 0) STOP "ERROR ALLOCATING VEXPAS"
    !
    INDX = LXPOIN(L)
    DO I = 1, NTOT
       !"!"!"!"!"!"!"!"!"!"print*, "EXPAS", I, EXPAS(INDX+I, :)
       IF (EXPAS(INDX+I,2) == L) THEN
          NDAT = NDAT+1
          VEXPEN(NDAT) = EXPEN(INDX+I)
          VEXPERR(NDAT) = EXPERR(INDX+I)
          VEXPAS(NDAT,:) = EXPAS(INDX+I,:)
       ENDIF
    ENDDO
    !     
    !
    IF (IPRINT > 2) THEN
       WRITE(*,*)
       WRITE(*,*) 'EXPERIMENTAL ENERGIES FOR L =',L
       WRITE(*,*)
       WRITE(*,*) 'NDAT', NDAT
       DO I = 1, NDAT
          WRITE(*,*)(VEXPAS(I,J),J=1,2), VEXPEN(I), VEXPERR(I)
       ENDDO
       WRITE(*,*)
    ENDIF
    !
    !     
    RETURN 
  END SUBROUTINE SLCTLU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ASSGNU3(BENT, ASSIGNALL, N2, L, DIM, W2, W4, W2W2B, HAM, EIGEN, &
       VEXPAS, NDAT, BLAS, IPRINT)
    !
    !     SUBROUTINE THAT LOOKS FOR THE BEST ASSIGNATIONS  
    !     FOR THE CALCULATED LEVELS TO
    !     LOCAL BASIS STATES IN THE U(3) MODEL.
    !
    !     IF LINEAR IT ASSIGNS TO A n,l BASIS IF BENT TO A v,K BASIS
    !
    !     IN CASE OF HIGH MIXING THE ASSIGNMENT COULD BE AMBIGUOUS AND 
    !     THEN IT IS FOLLOWED BY SCALING DOWN THE MIXING PARAMETERS UNTIL 
    !     AN UNAMBIGUOUS ASSIGNMENT CAN BE FOUND AND THEN SCALE THEM UP 
    !     AGAIN TO THE ORIGINAL VALUE, PROJECTING IN EACH STEP ONTO THE 
    !     PREVIOUS ONES (PROCEDURE USED TO RESOLVE  STARK MIXED ROTATIONAL
    !     LEVELS, PROGRAMMED BY THOMAS MULLER). 
    !     
    !     INPUT
    !     BENT     : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !     ASSIGNALL: IF .F. ASSIGN UNAMBIGUOUSLY ONLY DATA WITH EXPERIMENTAL INFO  
    !     N2       : U(3) IRREP LABEL (BENDING)
    !     L        : VIBRATIONAL ANGULAR MOMENTUM LABEL
    !     DIM      : BLOCK DIMENSION
    !     W2       : ARRAY WITH THE SO(3) CASIMIR BLOCK 
    !     W4       : ARRAY WITH THE SO(3) CASIMIR BLOCK 
    !     W2W2B    : ARRAY WITH THE SO(3) W2·W2B + W2B·W2 CASIMIR BLOCK 
    !     HAM      : HAMILTONIAN MATRIX (EIGENVECTOR MATRIX IN INPUT)
    !     EIGEN    : EIGENVALUES VECTOR
    !     VEXPAS   : EXPERIMENTAL ASSIGNMENTS FOR POLYAD (L) FORMAT:(V l) OR (n l)
    !     NDAT     : NUMBER OF EXPERIMENTAL DATA FOR POLYAD (L)
    !
    !     OUTPUT
    !     BLAS  :  LOCAL BASIS (BENT OR LINEAR) VECTOR WITH ASSIGNMENTS
    !
    !     by Currix TM
    !     
    IMPLICIT NONE
    !
    !     DEFINITION OF VARIABLES
    !
    LOGICAL, INTENT(IN) :: BENT, ASSIGNALL
    INTEGER(KIND = I4B), INTENT(IN) :: N2
    INTEGER(KIND = I4B), INTENT(IN) :: L, DIM, NDAT, IPRINT
    REAL(KIND = DP), DIMENSION(:,:), INTENT(INOUT)  :: W2, W4, W2W2B, HAM
    REAL(KIND = DP), DIMENSION(:), INTENT(INOUT)  :: EIGEN
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: VEXPAS
    !
    INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: BLAS
    !
    !
    ! LOCAL VARIABLES                                     
    LOGICAL :: FLAG
    INTEGER(KIND = I4B) :: IERR
    REAL(KIND = DP), DIMENSION(1:NPMAX)  :: HPART
    INTEGER(KIND = I4B) :: MIX
    INTEGER(KIND = I4B), DIMENSION(1:DIM) :: PTEMP
    REAL(KIND = DP) :: FAC, MIXSTEP
    REAL(KIND = DP) :: VMAX, VT
    REAL(KIND = DP), DIMENSION(1:DIM,1:DIM) :: HAM2
    REAL(KIND = DP), DIMENSION(1:DIM) :: EIGEN2
    REAL(KIND = DP), DIMENSION(1:DIM,1:DIM) :: TBRACK
    INTEGER(KIND = I4B) :: I, J, J2, JMAX, ICOUNT
    !
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE ASSIGN STARTS HERE'
    !
    ! !     POLYAD DIMENSION
    ! DIM = (N2 - MOD(N2 - L,2) - L)/2 + 1
    !
    !      INITIALIZE COUNTER
    !
    ICOUNT = 0
    !
    !
    ! ALLOCATE(HAM2(1:DIM, 1:DIM), STAT = IERR)
    ! IF (IERR /= 0) STOP 'ERROR ALLOCATING HAM2 - ASSGNU3 MATRIX'
    ! HAM2 = 0.0_DP
    ! !
    ! ALLOCATE(EIGEN2(1:DIM), STAT = IERR)
    ! IF (IERR /= 0) STOP 'ERROR ALLOCATING EIGEN2 - ASSGNU3 MATRIX'
    ! EIGEN2 = 0.0_DP
    ! !
    ! ALLOCATE(PTEMP(1:DIM), STAT = IERR)
    ! IF (IERR /= 0) STOP 'ERROR ALLOCATING PTEMP - ASSGNU3 MATRIX'
    ! PTEMP = 0      
    !
    ALLOCATE(BLAS(1:DIM), STAT = IERR)
    IF (IERR /= 0) STOP 'ERROR ALLOCATING BLAS - ASSGNU3 MATRIX'
    PTEMP = 0      
    !
    IF (BENT) THEN
       !
       !     BENT CASE
       !
       !     PROJECTION OF THE EIGENVECTORS TO THE BENT BASIS
       !     (I) DIAGONALIZATION OF THE PAIRING OPERATOR TO GET TRANSF. BRACK.
       ! ALLOCATE(TBRACK(1:DIM, 1:DIM), STAT = IERR)
       ! IF (IERR /= 0) STOP 'ERROR ALLOCATING TBRACK - ASSGNU3 MATRIX'
       ! TBRACK = 0.0_DP
       !
       CALL SO3CASBUILD(W2, N2, L, Iprint)
       !
       !     CHANGE SIGN TO W2 TO GET THE CORRECT GROUND STATE WF
       !       CALL LA_SYEVR(A=-W2, W=HAM2(:,1), JOBZ='V', UPLO='U')
       W2 = -W2
       CALL LA_SYEVR(A=W2, W=EIGEN2, JOBZ='V', UPLO='U')
       TBRACK = W2
       !
       W2 = MATMUL(TBRACK, HAM)
       !
       !      LOOK FOR MAXIMUM COMPONENT
       IF (ASSIGNALL) THEN
          CALL MAXCU3(BENT, N2, L, W2, DIM, BLAS, FLAG, IPRINT)
       ELSE
          CALL MAXCU3EXP(BENT, VEXPAS, NDAT, N2, L, W2, DIM, BLAS, FLAG, IPRINT)
       ENDIF
       !
       IF (IPRINT >= 3) THEN
          WRITE(*,*) 'MAXC', FLAG
          WRITE(*,*) (BLAS(I),I=1,DIM)
       ENDIF
       !
       !     SAVE COPY OF INITIAL PARAMETERS
       HPART = H_4b_pars
       !     
       !     MIX: NUMBER OF ITERATIONS UP AND DOWN
       !     
       MIX = 0
       !
       !     SCALE DOWN
       !     
       FAC = 1.0_DP
       !     
       !     IF FLAG = .T. PROBLEMS WITH ASSIGNMENT     
       DO WHILE (FLAG)
          !     
          FAC = FAC/2.0D0
          !     
          MIX = MIX + 1
          !     
          !     RECOVER COPY OF INITIAL PARAMETERS
          H_4b_pars = HPART
          !     
          IF (IPRINT >= 2) WRITE(*,*)'SCALING DOWN MIXING BY ', FAC
          !
          ICOUNT = ICOUNT + 1
          !     
          CALL SCLHAM(FAC, BENT, N2, L, HAM, W2, W4, W2W2B, IPRINT)
          !     
          !     DIAGONALIZE HAMILTONIAN
          !
          EIGEN = 0.0D0
          !     
          CALL LA_SYEVR(A=HAM, W=EIGEN, JOBZ='V', UPLO='U')
          !     
          !     LOOK FOR MAXIMUM COMPONENT
          !     
          !     (II) MATRIX PRODUCT
          W2 = MATMUL(TRANSPOSE(TBRACK), HAM)
          !
          IF (ASSIGNALL) THEN
             CALL MAXCU3(BENT, N2, L, W2, DIM, BLAS, FLAG, IPRINT)
          ELSE
             CALL MAXCU3EXP(BENT, VEXPAS, NDAT, N2, L, W2, DIM, BLAS, FLAG, IPRINT)
          ENDIF
          !     
          IF (IPRINT >= 4) THEN
             write(*,*) 'maxc', flag
             WRITE(*,*) (BLAS(I),I=1,DIM)
          ENDIF
          !
       END DO
       !
       IF (IPRINT >= 1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
       !            
       IF (ICOUNT == 0) RETURN ! WELL DEFINED EIGENVECTORS
       !
       !     
       !     SCALE UP
       MIXSTEP = 1.0_DP/(1.5_DP*2.0_DP**MIX)
       FAC = 1.0_DP/(2.0_DP**MIX)
       !
       HAM = W2 
       !
       scaleup_so3: DO
          !
          !     RECOVER COPY OF INITIAL PARAMETERS
          H_4b_pars = HPART
          !     
          MIXSTEP = 1.5_DP*MIXSTEP
          FAC = FAC + MIXSTEP
          !
          IF (FAC > 1.0_DP) FAC = 1.0_DP
          !     
          IF (IPRINT >= 2) WRITE(*,*) 'SCALING MIXING UP BY ', FAC 
          !
          ICOUNT = ICOUNT + 1
          !
          CALL SCLHAM(FAC,BENT,N2,L,HAM2,W2,W4,W2W2B, iprint)
          !     
          !     DIAGONALIZE HAMILTONIAN
          EIGEN = 0.0_DP
          !
          CALL LA_SYEVR(A=HAM2, W=EIGEN, JOBZ='V', UPLO='U')
          !
          !     (II) MATRIX PRODUCT
          W2 = MATMUL(TRANSPOSE(TBRACK), HAM2)
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          !   
          !   PROJECTING ONTO FORMER EIGENVECTORS
          !   THE EIGENVECTORS ARE PROJECTED ONTO THE PREVIOUS ONES AND 
          !   IF EITHER SOME OF THE PROJECTIONS SQUARED IS LESS THAN 0.5
          !   OR TWO HAVE THE SAME MAXIMAL PROJECTIONS THE PARAMETERS ARE 
          !   RESCALED DOWN AGAIN.
          !   
          IF (ASSIGNALL) THEN
             !     SCALAR PRODUCT
             DO I = 1, DIM 
                VMAX = 0.0_DP
                PTEMP(I) = 0
                DO J = 1, DIM
                   VT = DOT_PRODUCT(HAM(:,J), W2(:,I))
                   !     
                   !     MAXIMAL COMPONENT
                   VT = VT * VT
                   IF (VT > VMAX) THEN
                      VMAX = VT
                      JMAX = BLAS(J)
                   ENDIF
                ENDDO
                PTEMP(I) = JMAX 
                !     
                !     CHECK FOR AMBIGUITIES 
                IF (VMAX < 0.5_DP) THEN
                   FAC = FAC - MIXSTEP
                   MIXSTEP = MIXSTEP/3.0_DP
                   CYCLE scaleup_so3 
                ELSE
                   DO J2 = 1, I-1
                      IF (PTEMP(J2).EQ.JMAX) THEN
                         FAC = FAC - MIXSTEP
                         MIXSTEP = MIXSTEP/3.0_DP
                         CYCLE scaleup_so3
                      ENDIF
                   ENDDO
                ENDIF
                !
             ENDDO
             !
          ELSE
             !
             !     COMPARING ONLY TO STATES WITH BLAS(I)/= -1
             !     SCALAR PRODUCT
             PTEMP = -1
             !
             !
             DO I = 1, DIM
                IF (BLAS(I).NE.-1) THEN
                   VMAX = 0.0_DP
                   !
                   DO J = 1, DIM
                      !
                      VT = DOT_PRODUCT(HAM(:,I),W2(:,J))
                      !
                      !     MAXIMAL COMPONENT
                      VT = VT * VT
                      IF (VT > VMAX) THEN
                         VMAX = VT
                         JMAX = J
                      ENDIF
                   ENDDO
                   !     
                   PTEMP(JMAX) = BLAS(I)
                   !
                   IF (IPRINT > 2) WRITE(*,*) I, BLAS(I), JMAX, VMAX
                   !
                   !     CHECK FOR AMBIGUITIES 
                   IF (VMAX < 0.5_DP) THEN
                      FAC = FAC - MIXSTEP
                      MIXSTEP = MIXSTEP/3.0_DP
                      CYCLE scaleup_so3
                   ENDIF
                   DO J2 = 1, DIM
                      IF (PTEMP(J2) == BLAS(I).AND. J2 /= JMAX) THEN
                         FAC = FAC - MIXSTEP
                         MIXSTEP = MIXSTEP/3.0_DP
                         CYCLE scaleup_so3
                      ENDIF
                   ENDDO
                   !                    
                ENDIF
                !
             ENDDO
             !
          ENDIF
          !
          HAM = W2
          !
          BLAS = PTEMP
          !
          IF (FAC == 1.0_DP) THEN
             !
             CALL SCLHAM(FAC,BENT,N2,L,HAM,W2,W4,W2W2B,iprint)
             !
             EIGEN = 0.0_DP
             !
             CALL LA_SYEVR(A=HAM, W=EIGEN, JOBZ='V', UPLO='U')
             !     
             IF (IPRINT >= 1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
             !
             RETURN
             !    
          ENDIF
          !
       ENDDO scaleup_so3
       !
       !
    ELSE
       !     LINEAR CASE
       !     
       !     LOOK FOR MAXIMUM COMPONENT
       IF (ASSIGNALL) THEN
          CALL MAXCU3(BENT, N2, L, HAM, DIM, BLAS, FLAG, IPRINT)
       ELSE
          CALL MAXCU3EXP(BENT, VEXPAS, NDAT, N2, L, HAM, DIM, BLAS, FLAG, IPRINT)
       ENDIF
       !     
       IF (IPRINT >= 3) THEN
          WRITE(*,*) 'MAXC', FLAG
          WRITE(*,*) (BLAS(I),I=1,DIM)
       ENDIF
       !      
       !  SAVE COPY OF INITIAL PARAMETERS
       HPART = H_4b_pars
       !     
       !     MIX: NUMBER OF ITERATIONS UP AND DOWN
       MIX = 0
       !     
       !     SCALE DOWN
       !     
       FAC = 1.0D0
       !
       !      IF FLAG = .T. : PROBLEMS WITH ASSIGNMENT
       DO WHILE (FLAG)
          !     
          FAC = FAC/2.0_DP
          !
          MIX = MIX + 1
          !     
          !     RECOVER COPY OF INITIAL PARAMETERS
          H_4b_pars = HPART
          !
          IF (IPRINT >= 2) WRITE(*,*)'SCALING DOWN MIXING BY ', FAC
          !     
          ICOUNT = ICOUNT + 1
          !     
          CALL SCLHAM(FAC, BENT, N2, L, HAM, W2, W4, W2W2B, IPRINT)
          !
          !     DIAGONALIZE HAMILTONIAN
          !
          EIGEN = 0.0_DP
          !
          CALL LA_SYEVR(A=HAM, W=EIGEN, JOBZ='V', UPLO='U')
          !     
          !     LOOK FOR MAXIMUM COMPONENT
          !
          ! print*, "E1", HAM(:,1)**2
          ! print*, ""
          ! print*, "E2", HAM(:,2)**2
          ! print*, ""
          ! print*, "E3", HAM(:,3)**2
          !
          IF (ASSIGNALL) THEN
             CALL MAXCU3(BENT, N2, L, HAM, DIM, BLAS, FLAG, IPRINT)
          ELSE
             CALL MAXCU3EXP(BENT, VEXPAS, NDAT, N2, L, HAM, DIM, BLAS, FLAG, IPRINT)
          ENDIF
          !
          IF (IPRINT >= 3) THEN
             WRITE(*,*) 'MAXC DOWN', FLAG
             WRITE(*,*) (BLAS(I),I=1,DIM)
          ENDIF
          !
       ENDDO
       !
       IF (IPRINT >= 1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
       !            
       IF (ICOUNT == 0) RETURN ! WELL DEFINED EIGENVECTORS
       !
       !     SCALE UP
       !     
       MIXSTEP = 1.0_DP/(1.5_DP*2.0_DP**MIX)
       FAC = 1.0_DP/(2.0_DP**MIX)      
       !
       scaleup_u2: DO
          !     RECOVER COPY OF INITIAL PARAMETERS
          H_4b_pars = HPART
          !     
          MIXSTEP = 1.5_DP*MIXSTEP
          FAC = FAC + MIXSTEP
          !
          IF (FAC > 1.0_DP) FAC = 1.0_DP
          !
          IF (IPRINT >= 2) WRITE(*,*) 'SCALING MIXING UP BY ',FAC 
          !     
          ICOUNT = ICOUNT + 1
          !          
          CALL SCLHAM(FAC, BENT, N2, L, HAM2, W2, W4, W2W2B, iprint)
          !     
          !     DIAGONALIZE HAMILTONIAN
          !
          EIGEN = 0.0_DP
          !
          CALL LA_SYEVR(A=HAM2, W=EIGEN, JOBZ='V', UPLO='U')
          !
          !
          !     PROJECTING ONTO FORMER EIGENVECTORS
          !     THE EIGENVECTORS ARE PROJECTED ONTO THE PREVIOUS ONES AND 
          !     IF EITHER SOME OF THE PROJECTIONS SQUARED IS LESS THAN 0.5
          !     OR TWO HAVE THE SAME MAXIMAL PROJECTIONS THE PARAMETERS ARE 
          !     RESCALED DOWN AGAIN.
          !     
          IF (ASSIGNALL) THEN
             !     SCALAR PRODUCT
             DO I = 1, DIM 
                VMAX = 0.0D0
                PTEMP(I) = 0
                DO J = 1, DIM
                   !
                   VT = DOT_PRODUCT(HAM(:,J), HAM2(:,I))
                   !
                   !     MAXIMAL COMPONENT
                   VT = VT * VT
                   IF (VT > VMAX) THEN
                      VMAX = VT
                      JMAX = BLAS(J)
                   ENDIF
                ENDDO
                PTEMP(I) = JMAX
                print*, i, vmax, jmax, ptemp(i)
                !     
                !     CHECK FOR AMBIGUITIES
                IF (VMAX < 0.5_DP) THEN
                   FAC = FAC - MIXSTEP
                   MIXSTEP = MIXSTEP/3.0_DP
                   CYCLE scaleup_u2 
                ELSE
                   DO J2 = 1, I - 1
                      IF (PTEMP(J2).EQ.JMAX) THEN
                         FAC = FAC - MIXSTEP
                         MIXSTEP = MIXSTEP/3.0_DP
                         CYCLE scaleup_u2
                         !
                      ENDIF
                   ENDDO
                ENDIF
                !
             ENDDO
             !
          ELSE
             !
             !     COMPARING ONLY TO STATES WITH BLAS(I).NE.-1
             !     SCALAR PRODUCT
             PTEMP = -1
             !
             DO I = 1, DIM
                IF (BLAS(I) /= -1) THEN
                   VMAX = 0.0D0
                   !
                   DO J = 1, DIM
                      !
                      VT = DOT_PRODUCT(HAM(:,I),HAM2(:,J))
                      !
                      !     MAXIMAL COMPONENT
                      VT = VT * VT
                      IF (VT > VMAX) THEN
                         VMAX = VT
                         JMAX = J
                      ENDIF
                   ENDDO
                   !
                   PTEMP(JMAX) = BLAS(I)
                   !
                   IF (IPRINT > 2) WRITE(*,*) I, BLAS(I), JMAX, VMAX
                   !
                   !     CHECK FOR AMBIGUITIES 
                   IF (VMAX < 0.5_DP) THEN
                      FAC = FAC - MIXSTEP
                      MIXSTEP = MIXSTEP/3.0_DP
                      CYCLE scaleup_u2
                   ENDIF
                   DO J2 = 1, DIM
                      IF (PTEMP(J2) == BLAS(I) .AND. J2 /= JMAX) THEN
                         FAC = FAC - MIXSTEP
                         MIXSTEP = MIXSTEP/3.0D0
                         CYCLE scaleup_u2
                      ENDIF
                   ENDDO
                   !                    
                ENDIF
                !
             ENDDO
             !
          ENDIF
          !
          HAM = HAM2
          !
          BLAS = PTEMP
          !
          print*, "FAC = ", FAC
          print*, "BLAS = ", BLAS
          !
          IF (FAC == 1.0_DP) THEN
             !
             CALL SCLHAM(FAC,BENT,N2,L,HAM,W2,W4,W2W2B, iprint)
             !
             EIGEN = 0.0_DP
             !
             CALL LA_SYEVR(A=HAM, W=EIGEN, JOBZ='V', UPLO='U')
             !     
             IF (IPRINT >= 1) WRITE(*,*) ICOUNT, ' ITERATIONS TO PROJECT'
             !
             IF (IPRINT > 2) WRITE(*,*) 'FINAL BLAS ', BLAS
             !
             RETURN
             !    
          ENDIF
          !
       END DO scaleup_u2
       !
    ENDIF
    !
  END SUBROUTINE ASSGNU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CPXPDTU3(N2, L, EXPAS, EIGEN, BLAS, EXPEN, EXPERR, NDAT, SPNTEX, CHSQP, IPRINT)
    !
    !    SUBROUTINE THAT BASED IN EXPERIMENTAL ASSIGNMENTS
    !    COMPARE THE DIAGONALIZATION RESULTS WITH THE
    !    EXPERIMENTAL ENERGIES IN THE BENDING U(3) MODEL.
    !    IN ADDITION IT COMPUTES STATISTICAL PARAMETERS.
    !
    !    
    !    INPUT
    !    N2      : SU(3) IRREP LABEL (BENDING)
    !    L       : VIBRATIONAL ANGULAR MOMENTUM
    !    EIGEN   : EIGENVALUES VECTOR
    !    BLAS    : LOCAL BASIS ASSIGNMENTS
    !    EXPEN   : VECTOR WITH EXPERIMENTAL ENERGIES FOR POLYAD V
    !    EXPERR  : VECTOR WITH EXPERIMENTAL ENERGY ERRORS FOR POLYAD V
    !    EXPASS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS FOR POLYAD V
    !    NDAT    : NUMBER OF EXPERIMENTAL DATA IN POLYAD V
    !
    !    OUTPUT
    !    SPNTX   : POINTER WITH ASSIGNMENTS TO EXPERIMENTAL DATA
    !    CHSQP   : PARTIAL CHI-SQUARE
    !
    !    by Currix TM
    !    
    IMPLICIT NONE
    !
    !     DEFINITION OF VARIABLES
    !     
    INTEGER(KIND = I4B), INTENT(IN) :: N2, L,  NDAT, IPRINT
    INTEGER(KIND = I4B), DIMENSION(:), INTENT(IN) :: BLAS
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: EXPAS 
    REAL(KIND = DP), DIMENSION(:), INTENT(IN)  :: EIGEN, EXPEN, EXPERR
    !
    INTEGER(KIND = I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SPNTEX
    REAL(KIND = DP), INTENT(OUT)  :: CHSQP
    !
    !     TEMPORAL VARIABLES                                     
    INTEGER(KIND = I4B) :: DIM, I, J, NOT_FIT_LEV
    !
    !
    IF (IPRINT > 2) WRITE(*,*) 'STARTING SUBROUTINE CPXPDTU3'
    !"!"!"!"!"!"!"!"!"!"print*, ndat
    !"!"!"!"!"!"!"!"!"!"print*, blas
    !"!"!"!"!"!"!"!"!"!"print*, expen
    !"!"!"!"!"!"!"!"!"!"print*, expas
    !
    IF (ALLOCATED(SPNTEX)) THEN
       DEALLOCATE(SPNTEX, STAT = IERR)
       IF (IERR /= 0) STOP "ERROR DEALLOCATING SPNTEX"
    ENDIF
    ALLOCATE(SPNTEX(1:NDAT), STAT = IERR)
    IF (IERR /= 0) STOP "ERROR ALLOCATING SPNTEX"
    SPNTEX = 0
    !
    !"!"!"!"!"!"!"!"!"!"!"print*, spntex              !
    !
    CHSQP = 0.0_DP
    NOT_FIT_LEV = 0
    !
    !     POLYAD DIMENSION
    DIM = (N2-MOD(N2-L,2)-L)/2+1
    !
    DO I = 1, NDAT
       DO J = 1, DIM
          !    
          !     COMPARE EXPERIMENTAL ASSIGNMENTS AND LOCAL BASIS ASSIGN,
          IF (BLAS(J) == EXPAS(I,1)) THEN
             !     POINTER
             SPNTEX(I) = J
             !     CHI SQUARE               
             IF (EXPERR(I) /= 0) THEN
                CHSQP = CHSQP + &
                     (EIGEN(J)-EXPEN(I))*(EIGEN(J)-EXPEN(I))/(EXPERR(I)*EXPERR(I))
                IF (DIS_RES) WRITE(*,*) EIGEN(J)-EXPEN(I) 
             ELSE
                NOT_FIT_LEV = NOT_FIT_LEV + 1
             ENDIF
             !
          ENDIF
          !  
       ENDDO
    ENDDO
    !
    IF (IPRINT >= 2) WRITE(*,*) 'EXPER. LEVELS', NDAT, "NOT FITTED LEVELS", NOT_FIT_LEV
    IF (IPRINT > 2) WRITE(*,*) 'LEAVING SUBROUTINE CPXPDTU3'
    !
    RETURN 
    !     
  END SUBROUTINE CPXPDTU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DSDTU3(BENT, N2, L, EIGEN, EXPEN, EXPERR, NDAT, EXPAS, AVEC, BLAS, SPNTX, CHSQP, IPRINT)
    !
    !    DISPLAYS THE RESULTS OF THE HAMILTONIAN 
    !   DIAGONALIZATION IN THE U(3) ALGEBRAIC MODEL.
    !
    !   
    !   INPUT
    !   BENT   : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !   L      : VIBRATIONAL ANGULAR MOMENTUM
    !   EIGEN  : EIGENVALUES VECTOR
    !   EXPEN  : VECTOR WITH EXPERIMENTAL ENERGIES
    !   EXPERR : VECTOR WITH EXPERIMENTAL ERRORS
    !   NDAT   : NUMBER OF EXPERIMENTAL DATA TO CONSIDER                
    !   EXPAS  : MATRIX WITH EXPERIMENTAL ASSIGNMENTS
    !   AVEC   : EIGENVECTORS
    !   BLAS   : POINTER TO LOCAL BASIS
    !   SPNTX  : POINTER WITH ASSIGNMENTS TO EXPERIMENTAL DATA
    !   CHSQP  : PARTIAL CHI-SQUARE
    !
    !   by Currix TM
    !   
    IMPLICIT NONE
    !
    !  DEFINITION OF VARIABLES
    !     
    LOGICAL, INTENT(IN) :: BENT
    INTEGER(KIND = I4B), INTENT(IN) :: N2, L, NDAT, IPRINT
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: EXPAS
    INTEGER(KIND = I4B), DIMENSION(:), INTENT(IN) :: BLAS, SPNTX
    REAL(KIND = DP), INTENT(IN)  :: CHSQP
    REAL(KIND = DP), DIMENSION(:), INTENT(IN)  :: EXPEN, EXPERR, EIGEN
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN)  :: AVEC
    !
    !     TEMPORAL VARIABLES                                     
    INTEGER(KIND = I4B) :: I, J, DIM
    !
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE DSDT'
    !
    !     POLYAD DIMENSION
    DIM = (N2 - MOD(N2-L, 2) - L)/2 + 1
    !
    !     
    IF (BENT) THEN 
       WRITE(*,1061) N2, L
    ELSE
       WRITE(*,1060) N2, L 
    ENDIF
    !
    WRITE(*,*)
    WRITE(*,1070) NDAT, CHSQP
    WRITE(*,*)
    !    
    WRITE(*,1080)
    !     
    IF (BENT) THEN
       !     
       !     BENT CASE
       !     
       IF (IPRINT >= 0 .AND. NDAT /= 0) THEN
          WRITE(*,1086)
          DO I = 1, NDAT
             WRITE(*,1090) (EXPAS(I,J),J=1,2), EIGEN(SPNTX(I)), &
                  EXPEN(I),EXPERR(I), EIGEN(SPNTX(I))-EXPEN(I)
          ENDDO
       ENDIF
       !
       IF (IPRINT >= 1) THEN
          !
          WRITE(*,1080)
          DO I = 1, DIM
             WRITE(*,1110) EIGEN(I), BLAS(I), L
          ENDDO
       ENDIF
       !     
       IF (IPRINT >= 2) THEN
          !
          WRITE(*,1080)
          DO I = 1, DIM
             WRITE(*,1110) EIGEN(I), BLAS(I), L
             WRITE(*,*)
             DO J = 1, DIM
                WRITE(*,1120) AVEC(J,I), N2-(2*J-2+MOD(N2-L,2)), L
             ENDDO
             WRITE(*,*)
          ENDDO
       ENDIF
       ! 
    ELSE
       !     
       !     LINEAR CASE
       !     
       IF (IPRINT>=0 .AND. NDAT/=0) THEN
          WRITE(*,1085)
          DO I = 1, NDAT
             WRITE(*,1090) (EXPAS(I,J),J=1,2), EIGEN(SPNTX(I)), &
                  EXPEN(I),EXPERR(I),EIGEN(SPNTX(I))-EXPEN(I)
          ENDDO
       ENDIF
       !
       IF (IPRINT>=1) THEN
          !
          WRITE(*,1080)
          DO I = 1, DIM
             WRITE(*,1110) EIGEN(I), BLAS(I), L
          ENDDO
       ENDIF
       !     
       IF (IPRINT>=2) THEN
          !
          WRITE(*,1080)
          DO I = 1, DIM
             WRITE(*,1110) EIGEN(I), BLAS(I), L
             WRITE(*,*)
             DO J = 1, DIM
                WRITE(*,1120) AVEC(J,I), N2-(2*J-2+MOD(N2-L,2)), L
             ENDDO
             WRITE(*,*)
          ENDDO
       ENDIF
       ! 
    ENDIF
    !  
    WRITE(*,1080)
    WRITE(*,*)
    WRITE(*,*)
    !
    RETURN
    !      
1060 FORMAT(4X,'LINEAR CASE',6X,'N2 = ',I4,4X,'VIBR. ANG. MOM. = ',I4)
1061 FORMAT(4X,'BENT CASE',6X,'N2 = ',I4,4X,'VIBR. ANG. MOM. = ',I4)
1070 FORMAT(4X,'NUMBER OF EXP. DATA = ', I4, 6X, 'PARTIAL CHI SQ.',D15.5)
1080 FORMAT(4X,'***********************************************************',/, &
         4X,'***********************************************************')
1085 FORMAT(5X,'EXP.(n  L)',4X,'Calc.Energy',6X,'Exp.Energy',3X,'Calc.-Exp.')
1086 FORMAT(5X,'EXP.(v  K)',4X,'Calc.Energy',6X,'Exp.Energy',3X,'Calc.-Exp.')
1090 FORMAT(7X,'(',I3,I3,')',2X,F12.4,2X,F10.2,'(',F4.1,')',2X,F10.5)
1110 FORMAT(7X,D15.8,2X,'(',I3,I3,')')
1120 FORMAT(7X,D15.8,7X,'|',I3,I3,' >')
    !
  END SUBROUTINE DSDTU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MAXCU3(BENT, N2, L, EVEC, DIM, LMC, FLAG, IPRINT)
    !
    !     SUBROUTINE THAT LOOKS FOR MAXIMAL TYPE I BASIS (LINEAR)
    !     COMPONENTS IN THE EIGENVECTORS OF THE BENDING SU(3)
    !     MODEL HAMILTONIAN.
    !
    !     INPUT
    !     BENT  : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !     N2    : U(3) IRREP LABEL (BENDING)
    !     L     : VIBRATIONAL ANGULAR MOMENTUM
    !     EVEC  : EIGENVECTORS MATRIX                   
    !     DIM   : BLOCK DIMENSION               
    !
    !     OUTPUT
    !     FLAG  :  .TRUE. IF THERE IS AMBIGUITY IN THE ASSIGNMENT
    !     LMC   :  LOCAL BASIS ASSIGNMENTS 
    !
    !     by Currix TM
    !     
    IMPLICIT NONE
    !     
    !     DEFINITION OF VARIABLES
    LOGICAL, INTENT(IN) :: BENT
    INTEGER(KIND = I4B), INTENT(IN) :: N2, L, DIM, IPRINT
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: EVEC
    !     
    LOGICAL, INTENT(OUT) :: FLAG
    INTEGER(KIND = I4B), DIMENSION(:), INTENT(OUT) :: LMC
    !
    ! LOCAL VARIABLES
    INTEGER (KIND = I4B) :: I, I2, J, JMAX
    REAL(KIND = DP) :: VMAX,VT
    !     
    !     COMPUTATION OF MAXIMAL COMPONENTS
    !     
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE MAXC'
    !
    FLAG = .FALSE.
    !
    DO I=1, DIM
       VMAX = 0.0_DP
       LMC(I) = -1
       DO J=1, DIM
          VT = EVEC(J,I) * EVEC(J,I)
          IF (VT > VMAX) THEN
             JMAX = J
             VMAX = VT
          END IF
       END DO
       DO I2=1, I-1
          IF (BENT .AND. LMC(I2) == (JMAX - 1)) FLAG = .TRUE. ! NONLINEAR CASE
          IF (.NOT. BENT .AND. LMC(I2) == (N2 - 2*JMAX+2-MOD(N2-L,2))) FLAG = .TRUE. ! LINEAR CASE
       END DO
       IF (VMAX < 0.5D0) FLAG = .TRUE.
       IF (BENT) THEN
          LMC(I) = JMAX - 1
       ELSE
          LMC(I) =  N2 - 2*JMAX+2-MOD(N2-L,2)
       ENDIF
       !
    END DO
    !     
    RETURN
    !   
  END SUBROUTINE MAXCU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MAXCU3EXP(BENT, VEXPAS, NDAT, N2, L, EVEC, DIM, LMC, FLAG, IPRINT)
    !
    !
    !     INPUT
    !     BENT  : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !     LDEXP : LEADING DIMENSION OF MATRIX VEXPAS
    !     VEXPAS: EXPERIMENTAL ASSIGNMENTS FOR POLYAD (L) FORMAT:(V l) OR (n l)
    !     NDAT  : NUMBER OF EXPERIMENTAL DATA FOR POLYAD (L)
    !     N2    : SU(3) IRREP LABEL (BENDING)
    !     L     : VIBRATIONAL ANGULAR MOMENTUM
    !     EVEC  : EIGENVECTORS MATRIX                   
    !     DIM   : POLYAD DIMENSION               
    !
    !     OUTPUT
    !     FLAG  :  .TRUE. IF THERE IS AMBIGUITY IN THE ASSIGNMENT
    !     LMC   :  LOCAL BASIS ASSIGNMENTS 
    !
    !   by Currix TM
    !     
    IMPLICIT NONE
    !     
    !     DEFINITION OF VARIABLES
    !   
    LOGICAL, INTENT(IN) :: BENT
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: VEXPAS
    INTEGER(KIND = I4B), INTENT(IN) :: NDAT, N2, L, DIM, IPRINT
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: EVEC
    !
    LOGICAL, INTENT(OUT) :: FLAG
    INTEGER(KIND = I4B), DIMENSION(:), INTENT(OUT) :: LMC
    !     
    ! LOCAL VARIABLES
    INTEGER(KIND = I4B) :: I, I2, I3, J, JMAX, NASSIGNED, NUT
    REAL(KIND = DP) :: VMAX,VT
    !     
    !     COMPUTATION OF MAXIMAL COMPONENTS
    !    
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE MAXC EXP'
    !
    FLAG = .FALSE.
    !
    NASSIGNED = 0
    !
    LMC = -1
    !
    DO I=1, DIM
       VMAX = 0.0D0
       DO J=1, DIM
          VT = EVEC(J,I) * EVEC(J,I)
          IF (VT > VMAX) THEN
             JMAX = J
             VMAX = VT
          END IF
       END DO
       !     ASSIGN GROUND STATE OF EACH POLYAD
       IF (BENT .AND. (JMAX - 1) == 0) LMC(I) = 0
       !
       IF (.NOT. BENT .AND. (N2 - 2*JMAX+2-MOD(N2-L,2)) == L) LMC(I) = L
       !
       IF (VMAX < 0.5D0) FLAG = .TRUE.
       !     COMPARE WITH EXPERIMENTAL DATA
       IF (BENT) THEN
          NUT = (JMAX - 1) ! NONLINEAR CASE
       ELSE
          NUT = (N2 - 2*JMAX+2-MOD(N2-L,2)) ! LINEAR CASE
       ENDIF
       DO I2 = 1, NDAT
          IF (VEXPAS(I2,1) == NUT) THEN
             NASSIGNED = NASSIGNED + 1
             LMC(I) = VEXPAS(I2,1)
             !
             DO I3 = 1, I-1
                IF (LMC(I3) == NUT) FLAG = .TRUE.
             END DO
             !     
             IF (VMAX < 0.5D0) FLAG = .TRUE.
             !
             EXIT
             !
          ENDIF
       ENDDO
       !
       !
       IF (NASSIGNED == NDAT) EXIT
       !
    ENDDO
    ! FLAG = .TRUE. IF NOT ALL EXPERIMENTAL ENERGY LEVELS ARE ASSIGNED
    IF (NASSIGNED /= NDAT) FLAG = .TRUE.
    !
    RETURN
    !
  END SUBROUTINE MAXCU3EXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SCLHAM(FAC,BENT,N2,L,HAM,W2MAT,W4MAT,W2W2BARMAT,IPRINT)
    !
    !     SUBROUTINE THAT SCALES DOWN (x1/2) THE MIXING PARAMETERS AND
    !     RECOMPUTES HAMILTONIAN IN THE U(3) MODEL.
    !     
    !     INPUT
    !     FAC  : FACTOR THAT RESCALES THE MIXING PARAMETERS
    !     BENT : .T. BENT MOLECULE, .F. LINEAR MOLECULE
    !     N2   : U(3) IRREP LABEL (BENDING)
    !     L    : VIBRATIONAL ANGULAR MOMENTUM LABEL
    !     W2MAT : ARRAY WITH THE SO(3) CASIMIR P BLOCK 
    !     W4MAT : ARRAY WITH THE SO(3) CASIMIR P^2 BLOCK 
    !     W2W2BARMAT : ARRAY WITH THE OPERATOR W2·WBAR2 + WBAR2·W2 BLOCK 
    !     
    !
    !     OUTPUT
    !     HAM      : RECOMPUTED HAMILTONIAN MATRIX
    !
    !
    !     by Currix TM
    !     
    IMPLICIT NONE
    !
    !     
    !     DEFINITION OF VARIABLES     
    REAL(KIND = DP), INTENT(IN) :: FAC
    LOGICAL, INTENT(IN) :: BENT
    INTEGER(KIND = I4B), INTENT(IN) :: N2, L, IPRINT
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: W2MAT, W4MAT, W2W2BARMAT
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: HAM
    !
    !
    IF (IPRINT > 2) WRITE(*,*) 'SUBROUTINE SCALE MIXING STARTS HERE'
    !
    !     RESCALE PARAMETERS
    !
    IF (.NOT.BENT) THEN
       H_4b_pars(4) = H_4b_pars(4)*FAC         
       H_4b_pars(11) = H_4b_pars(11)*FAC
       H_4b_pars(13) = H_4b_pars(13)*FAC
    ELSE
       H_4b_pars(1) = H_4b_pars(1)*FAC
       H_4b_pars(2) = H_4b_pars(2)*FAC
       H_4b_pars(5) = H_4b_pars(5)*FAC
       H_4b_pars(6) = H_4b_pars(6)*FAC
       H_4b_pars(8) = H_4b_pars(8)*FAC
       H_4b_pars(9) = H_4b_pars(9)*FAC
    ENDIF
    H_4b_pars(7) = H_4b_pars(7)*FAC ! Both cases
    H_4b_pars(12) = H_4b_pars(12)*FAC ! Both cases
    H_4b_pars(14) = H_4b_pars(14)*FAC ! Both cases
    CALL HBLDU3GEN(N2, L, HAM, W2MAT, W4MAT, W2W2BARMAT, IPRINT)  
    !
    RETURN
  END SUBROUTINE SCLHAM
  !
END MODULE FIT_2DVM
      



