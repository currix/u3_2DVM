MODULE defparam_2DVM
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! Type definitions
  TYPE :: so3_bas
     INTEGER(KIND = I4B) :: N_U3_val
     INTEGER(KIND = I4B) :: omega_SO3_val
     INTEGER(KIND = I4B) :: L_val
  END TYPE so3_bas
  !
  TYPE :: u2_bas
     INTEGER(KIND = I4B) :: N_U3_val
     INTEGER(KIND = I4B) :: np_U2_val
     INTEGER(KIND = I4B) :: L_val
  END TYPE u2_bas
  !
  TYPE :: exp_level
     REAL(KIND = DP) :: exp_term_value
     REAL(KIND = DP) :: exp_err_value
     INTEGER(KIND = I4B) :: exp_v
     INTEGER(KIND = I4B) :: exp_L
     CHARACTER(LEN = 3) :: reference
  END TYPE exp_level
  !
  TYPE term_values
     TYPE(exp_level) :: exper_level
     TYPE(term_values), POINTER :: next
  END TYPE term_values
  !
  !
  INTEGER(KIND = I4B) :: N_val ! U(3) [N] value
  !
  INTEGER(KIND = I4B) :: l_val ! Vibrational Angular Momentum
  !
  INTEGER(KIND = I4B) :: dim_block ! L Block dimension
  !
END MODULE defparam_2DVM
