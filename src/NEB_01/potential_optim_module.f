        module potential_optim_module
        use type_module

! Added by MFH:
        use opt_events, only: opt_event_next_pot_grad
! End of additions

*** choice 1
!      use potential_module, only: pot, grad,
!     &                      record_walker_position

*** choice 2
!      use ring_polymer_potential_module, only: ring_polymer_action,
!     &                      ring_polymer_grad
!      use ring_polymer_data_module, only: record_ring_position

*** choice 3
!      use neb_module, only: neb_gradient, record_neb_position

      implicit none

      private
      public :: pot_and_grad, set_min_pot, num_calls,
     &          record_position

      integer, save :: iopt, num_calls=0

      contains

*********************************************************************
      subroutine set_min_pot(ichoice)
      integer, intent(in) :: ichoice

      iopt=ichoice

      num_calls=0

      end subroutine set_min_pot

*********************************************************************
      subroutine record_position(x)
      real(kind=dp), dimension(:), intent(in) :: x

      select case(iopt)

      case(1)

!         call record_walker_position(x)

      case(2)

!         call record_ring_position(x)

      case(3)

!         call record_neb_position(x)

      case default
         print *,'naming: no case for iopt'
         stop

      end select

      end subroutine record_position

*********************************************************************
      subroutine pot_and_grad(x,f,g)
      real(kind=dp), dimension(:), intent(in) :: x
      real(kind=dp), intent(out) :: f
      real(kind=dp), dimension(:), intent(out) :: g

      select case(iopt)

      case(1)

!         f=pot(x)
!         call grad(x,g)

      case(2)

!         f=ring_polymer_action(x)
!         call ring_polymer_grad(x,g)

      case(3)

!         f=1.0_dp
!         call neb_gradient(x,g)

      case(4)
        ! Added by mfh:
        call opt_event_next_pot_grad(x,f,g)

      case default
         print *,'naming: no case for iopt'
         stop

      end select

      num_calls=num_calls+1

      end subroutine pot_and_grad

*********************************************************************
      end module potential_optim_module
