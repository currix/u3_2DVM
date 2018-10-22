module needed_functions
  use nrtype
  use defparam_2DVM
  use FIT_2DVM
  use u3_2dvm_mod
  implicit none
contains
  subroutine FCN(npar,grad,fval,xval,iflag,PRE_CHI_SQRE)
    implicit none
    Real(kind=DP):: grad(*),xval(*),fval
    integer:: iflag,npar
    Real(kind=DP):: PRE_CHI_SQRE,rms
    
    if (iflag .eq. 1)  then
       if (iprint.gt.1) write(*,*) " FCN iflag = ", iflag
    endif
    fval=PRE_CHI_SQRE(XVAL(1),XVAL(2),XVAL(3),XVAL(4),XVAL(5),XVAL(6), &
         XVAL(7),XVAL(8),XVAL(9),XVAL(10),XVAL(11), &
         XVAL(12),XVAL(13),XVAL(14))
    rms=0.0_DP
    rms=sqrt(fval/(dble(totdat-npar)))
    write(*,*)" ___________________________________ "
    write(*,*)"| rms = ",rms," |"
    write(*,*)"|___________________________________|"
  end subroutine FCN
end module needed_functions
