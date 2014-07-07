FUNCTION loggamma(z)
  !
  ! Value of ln[Gamma(z)]. For z > gamma_limit compute
  ! the asymptotic expansion of Ln(Gamma(z)) Abram. 6.1.41
  !
  USE nrtype
  !
  REAL(KIND = DP), intent(in) :: z
  !
  REAL(KIND = DP) :: loggamma
  !
  REAL(KIND = DP), PARAMETER :: gamma_limit = 60.0_DP
  !
  IF (z < gamma_limit) THEN
     !
     loggamma = LOG(GAMMA(z))
     !
  ELSE
     !
     loggamma = -z + LOG(z)*(z-0.5_DP) + 0.5_DP*LOG(2.0_DP*PI_D)  &
          + 1.0_DP/(12.0_DP*z) &
          - 1.0_DP/(360.0_DP*z**3) 
     !
  ENDIF
END FUNCTION loggamma
