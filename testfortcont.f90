!
!=====NOTES=============================================================
!
!	program for testing FORTCONT subroutines on some vector fields
!	defined in the module VECFIELDS
!
!	At the time, FORTCONT can be used only for 1 parameter
!	continuation. More specificly, it is a set of predictor-corrector
!	(PC) methods.
!
!	The PC methods availible are corrector using the pseudo-arclength,
!	correctorAL and corrector based on finding the nearest point on
!	the solution curve, using the Moore-Penrose matrix pseudoinverse,
!	correctorPM
!
!	Author:	Martin Isoz
!	Date  :	01/11/2013
!	System:	Linux 3.2.0-54-generic
!	Lic.  :	FreeBSD (viz. LICENSE)
!
!-----------------------------------------------------------------------
!
!	in the following text:
!	nVar	... number of state variables
!	nParC	... number of continued parameters
!	nEqs	... number of vector field elements (number of equations)
!	nParF	... number of vector field fixed parameters
!
!-----------------------------------------------------------------------
!
!	vector field in VECFIELDS must be a subroutine callable with:
!	call VECFIELDS (xVec, RSVec, pars), where
!	xVec	... state variable, vector of points for which I am
!		    calculating the return values of vecfield
!		... double precision, dimension(nVar+nParC),intent(in)
!	RSVec	... vector of vector field right hand sides on output
!		... double precision, dimension(nEqs),intent(out)
!	pars	... fixed parameters of the subroutine
!		... double precision, dimension(nParF),intent(in)
!
!-----------------------------------------------------------------------
!
!
!     continuation itself is done by calling the subroutine SOLVER by:
!
!          call solver (VECFIELD,PREDICTOR,CORRECTOR,xVec,ip,&
!    & BMin,BMax,BTrg,pars,(/nDim,nVar/),nParF,g,h,rstep,nprstep,eps,&
!    & maxiter,wrtunt,retcode)
!
!     and input parameters have the following meaning
!
!     subroutines
!     VECFIELD    ... name of the subroutine from VECFIELDS module (or
!                     any other vector field with requierred structure)
!     PREDICTOR   ... the predictor used for solution continuation
!     CORRECTOR   ... corrector used for the predicted solution refining
!                 ... correctorAL or correctorPM
!
!     problem variables
!     xVec  ...   initial solution point, state variables + continued
!                 parameters (state variables in the next)
!           ...   double precision, dimension(nVar+nParC),intent(in)
!     ip    ...   index (position) of continued parameters in xVec
!           ...   integer(kind=4),dimension(nParC),intent(in)
!     BMin  ...   minimal admissible values of state variables
!           ...   double precision, dimension(nVar+nParC),intent(in)
!     BMax  ...   maximal admissible values of state variables
!           ...   double precision, dimension(nVar+nParC),intent(in)
!     BTrg  ...   target values of state variables on input, actually
!                 obtained values on the output
!           ...   double precision, dimension(nVar+nParC),intent(inou)
!     pars  ...   fixed parameters of the solved problem (parameters for
!                 the called VECFIELD)
!           ...   double precision,dimension(nParF),intent(in)
!     n     ...   dimension of the solved problem,
!           ...   n(1)=nEqs, n(2)=nVar+nParC
!           ...   integer(kind=4),dimension(2),intent(in)
!     m     ...   number of fixed parameters of the solved problem
!           ...   m = nParF
!           ...   integer(kind=4),dimension(1),intent(in)
!
!     auxiliary variables of the FORTCONT subroutines
!     g     ...   step for the numerical differentiation. it is not
!                 adaptable. At the time, the symmetric forward/bacwards
!                 formula of second order is used, but due to its cost 
!                 (number of operations needed for its evaluation),
!                 it may be changed to first order forwards formula in
!                 the following versions
!     h     ...   vector with steps for the predictor jump length
!           ...   (/h0,hMin,hMax/) - initial, minimal and maximal step
!                 the step length of predictor is adaptively changed in
!                 dependence of the speed of convergence of corrector
!                 methods
!     maxstep..   maximal number of steps of the main algorithm on input
!                 actual number of steps on output
!                 integer(kind=4),dimension(1),intent(inout)
!     nprstep..   for the predictor, there is a possibility to use
!                 higher order explicit numerical integration formula.
!                 if some multiple step method is used, nprstep is a
!                 number of steps needed to initiate the method. for
!                 example in case of Adams-Bashforts integration of 3rd
!                 order, nprstep = 3
!           ...   integer(kind=4),dimension(1),intent(in)
!     eps   ...   desired accuracy of corrector methods
!           ...   double precision,dimension(1),intent(in)
!     maxiter..   maximal number of iterations for the corrector method
!                 integer(kind=4),dimension(1),intent(in)
!
!     program runtime and IO variables
!     wrtunt...   unit to which the solution will be written.
!                 in accordance with fortran write statement, 6 is
!                 reserved for the terminal output.
!           ...   the unit must be openned and prepared for writing
!                 before the call to SOLVER subroutine
!           ...   integer(kind=4),dimension(1),intent(in)
!     retcode..   return code of the solver subroutine, basic
!                 diagnostics of the program runtime behavior
!                 -1- premature ending, probably due to a bug
!                 0 - some of the target values of state variables was
!                     reached, x(i) = BTrg(i)
!                 1 - algorithm stepped outside of region of interest,
!                     maximal(minimal) allowed value of state variable
!                     was exceeded
!                 2 - current step size is lower than minimal allowed
!                     step, hAct < hMin
!                 3 - corrector did not converge, the number of
!                     iterations in corrector method in current step
!                     exceeded maxiter parameter
!                 4 - maxstep was exceeded
!
!     the solution is written to wrtunt and is of the following form:
!
!-----------------------------------------------------------------------
!
!     before call to SOLVER subroutine, it is mandatory to prepare all
!     the needed variables as listed above as well as the unit for
!     writing the result in.
!
!     the results are plotted throught the GNUFOR2 interface between
!     FORTRAN and GNUPLOT written by Alexey Kuznetsov. For displaying
!     the graphics, one has to have GNUPlot installed and obey the rules
!     for GNUFOR2 use.
!
!=====TESTING=PROGRAM=ITSELF============================================
!
      program testfortcont
!
!=====HEAD==============================================================
!
!     Author:      Martin Isoz
!     Date  :      01/11/2013
!     System:      Linux 3.2.0-54-generic
!     Lic.  :      FreeBSD (viz. LICENSE)
!
!     Program for testing the FORTCONT subroutines
!     
!     compilation: make/make run/make clean
!     NOTE: the makefile is set up for my computer, so if you will wish
!           to compile the code for yourself, you will need to modify it
!     NOTE: as the LAPACK and BLAS subroutines are used in FORTCONT
!           module, one has to have them installed and properly linked
!           before the compilation
!
!=======================================================================
      use Vecfields                                                     !load the module vecfields
      use FORTCONT                                                      !load the FORTCONT module
      use gnufor2                                                       !load GNUFOR2 for graphics
!
      implicit none
!
!=====INTERFACE=DEFINITION=AND=MODULES=LOADING==========================
!
!
!=====VARIABLES=DECLARATION=============================================
!
!     state variables and parameters
      double precision,dimension(2):: uu1,uu2,La                        !state variables - initial points
      double precision            :: MM,CC,De                           !used solution parameters
      integer(kind=4),parameter   :: nDim=2,nParF=3,nVar=3              !number of equations, fixed pars and state vars (all)
      integer(kind=4),parameter   :: ip(nVar-nDim)=(/3/)                !position of the parameter in state variables vector
!     SOLVER subroutines variables
      double precision            :: xVec(nVar),pars(nVar)              !vectors of state variables and parameters
      double precision,parameter  :: g=1.0d-4,eps=1.0d-10               !numerical derivation step and desired sol. accuracy
      double precision,parameter  :: h(3)=(/1.0d-4,1.0d-9,3.0d-2/)      !steps for the predictor
      double precision            :: BMin(nVar),BMax(nVar),BTrg(nVar)   !boundaries of the region of interest
      integer(kind=4)             :: retcode,maxstep=1000               !return code, maximal/acual number of steps of the alg.
      integer(kind=4),parameter   :: maxiter=15                         !maximal number of iterations in corrector methods
      integer(kind=4),parameter   :: nprstep=5                          !number of kept tangent vectors
      integer(kind=4),parameter   :: wrtunt=6                           !write unit
!     MAIN program variables
      integer(kind=4)             :: cyccount,rstep                     !cycles counter, resulting number of steps
!     junk variables
      double precision            :: tVec(1,nVar)
!     strings
!
!
!=====EXECUTABLE=PART===================================================
!
!-----solution-initiation-----------------------------------------------
!     fixed parameters values
      MM          = 1.0d0                                               !1.0e0
      CC          = 1.0d0                                               !1.0e0
      De          = 5.0d-2                                              !0.05d0
!
!     starting points for different runs of the algorithm
      uu1(2)         = -1.883798570d0                                   !-1.883798570d0
      uu2(2)         = 0.0d0                                            !0.0d0
      La(2)          = 0.999d0                                          !0.990d0
      uu1(1)         = 0.0d0
      uu2(1)         = 0.0d0
      La(1)          = 0.0d0
!
!     preparation of the SOLVER variables
      pars = (/MM,CC,De/)                                               !fixed parameters of elastica1D1P
!      xVec = (/uu1,uu2,La/)                                             !state variables - call without cycle
      BMin  = (/-1.0d3, -1.0d3, -1.0d3/)
      BMax  = (/1.0d3, 1.0d3, 1.0d0/)
!     openning of the writing unit
!      open(1, file='dataES.dat', status='new')                          !open the file for writing
      
      
      tVec(1,:)  = (/1.0d0,1.0d0,1.0d0/)                                !initial tanget vector
!
!-----call-to-SOLVER-subroutine-----------------------------------------
!
      do cyccount=1,2,1
!     variables restart
      xVec = (/uu1(cyccount),uu2(cyccount),La(cyccount)/)               !state variables - call throught cycle
      BTrg  = BMax
      rstep = maxstep
!     call to SOLVER
      call solver (elastica1D1P,predictorPT,correctorPM,xVec,ip,&
     & BMin,BMax,BTrg,pars,(/nDim,nVar/),nParF,g,h,rstep,nprstep,eps,&
     & maxiter,wrtunt,retcode)
!
!     write retcode to terminal
      write(unit=*,fmt=*) 'Return code: ',retcode

!
      end do
!
!-----graphical-outputs-------------------------------------------------
!
!      call run_gnuplot('commands1.txt')                                 !Eq. points continuation
!
!     delete the file with results
!      close(1, status='delete')
!
!=======================================================================
!
      end program testfortcont
