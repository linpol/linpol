module ttm3f_mod
!----------------------------------------------------------------------------!
! Parameters of the TTM3F potential                                         !
!----------------------------------------------------------------------------!
! (12-10-6) inverse polynomial + exponential parameters for the vdw ineractions
double precision, parameter ::   &
                                 vdwC=-0.72298855D+03,vdwD=0.10211829D+06, vdwE=0.37170376D+01
! smearing factors for dipole-dipole(aDD), charge-charge/charge-dipole(aCCaCD)
double precision, parameter :: aDD=0.175d0, aCCaCD=0.175d0
!   ........... polarizabilities
double precision, parameter :: polarM=1.444d0
double precision, parameter :: polfacO=0.837d0, polfacH = 0.496d0, polfacM=0.837d0
double precision, parameter :: dms_param1 =0.5d0,dms_param2 = 0.9578d0, dms_param3=0.012d0
!   ...........   gammaM
double precision, parameter :: gammaM=0.46d0
!----------------------------------------------------------------------------!
! Parameters related to the accuracy of the energy/derivatives calculation   !
!----------------------------------------------------------------------------!
integer, parameter          :: MAXITER  = 400
double precision, parameter :: diptol   = 1.d-15 
double precision, parameter :: dmix     = 0.7d0
!----------------------------------------------------------------------------!
! CONSTANTS                                                                  !
!----------------------------------------------------------------------------!
double precision, parameter :: CHARGECON = 18.22234397655801030455d0
double precision, parameter :: DEBYE  = 4.8033324d0
!----------------------------------------------------------------------------!
! Variables and allocatable arrays needed by the "ttm3f" subroutine         !
!----------------------------------------------------------------------------!
integer :: fO, lO, fH, lH, fM, lM, fO3, lO3, fH3, lH3, fM3, lM3
integer :: Nw_old
integer :: Nats, Natsd, Natsq
double precision, dimension(:,:), allocatable :: RM
double precision, dimension(:,:), allocatable :: dRM
double precision, dimension(:,:), allocatable :: DDT
double precision, dimension( : ), allocatable :: dip
double precision, dimension( : ), allocatable :: pr_dip
double precision, dimension( : ), allocatable :: Efq
double precision, dimension( : ), allocatable :: Efd
double precision, dimension( : ), allocatable :: charge
double precision, dimension( : ), allocatable :: phi
double precision, dimension(:,:,:,:), allocatable :: grdq
end module ttm3f_mod
