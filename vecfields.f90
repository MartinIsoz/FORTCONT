!
!=====NOTES=============================================================
!
!     Module with sample vecfields for the set of subrooutines FORTCONT
!
!     Author:      Martin Isoz
!     Date  :      01/11/2013
!     System:      Linux 3.2.0-54-generic
!     Lic.  :      FreeBSD (viz. LICENSE)
!
!=====MODULE=WITH=SUBROUTINES===========================================
!
      module Vecfields
      implicit none
      private
      public      :: elastica1D1P
      contains
!
!=====SUBROUTINES=======================================================
!
            subroutine elastica1D1P( xVec,RSVec,pars)                     !vypocet pravych stran vektoroveho pole
!
!     vector field for one parameter continuation, copied from the
!     MATCONT propagation material.
!
!     it is a model of a column from imperfect material under a
!     concentric axial load exhibiting the characteristic deformation of
!     buckling. load is represented by parameter La, which is used as
!     the continuation parameter. the material imperfections and other
!     influences are lumped into parameter De
!
!     for more details on the problem, see for example wiki:
!     http://en.wikipedia.org/wiki/Buckling
!
!     xVec  ... input vector (2 original state variables + 1 parameter) 
!     RSVec ... right hand side vector - values returned by vecfield
!     pars  ... fixed parameter of the problem
            implicit none
!           
            double precision,intent(in),dimension(3) :: xVec            !input, state variable
            double precision,intent(in),dimension(3) :: pars            !input, fixed parameters
            double precision,intent(out),dimension(2) :: RSVec          !output
!     local variables (only for clarity of the code)
            double precision  :: MM,CC,De                               !parameters
            double precision  :: uu1,uu2,La                             !state variables
!     state variables extraction
            uu1         = xVec(1)
            uu2         = xVec(2)
            La          = xVec(3)
!     parameters extraction
            MM          = pars(1)
            CC          = pars(2)
            De          = pars(3)
!     output values calculation
            RSVec(1)    = uu2
            RSVec(2)    = 1.0d0/MM*(-(uu1-De) + &
     &2.0d0*La*sin(uu1)) - CC*uu2
            return
      end subroutine elastica1D1P
!
!=======================================================================
!
      end module Vecfields
