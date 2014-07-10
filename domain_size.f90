!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain_size

!BOP
! !MODULE: domain_size
!
! !DESCRIPTION:
!  This module contains parameters for the global model domain size
!  decomposition block size.  It is used by the domain and block
!  modules for decomposing the model domain across processors.
!
! !REVISION HISTORY:
!  CVS:$Id: domain.F90,v 1.16 2002/11/18 17:35:33 jfd Exp $
!  CVS:$Name:  $

! !USES:

   use kinds_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public ::  &  ! model size parameters
      nx_global =  320 ,&! extent of horizontal axis in i direction
      ny_global =  384 ,&! extent of horizontal axis in j direction


!EOP
!BOC
!EOC
!***********************************************************************

 end module domain_size

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
