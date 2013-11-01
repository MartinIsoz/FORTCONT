!
!=====POZNAMKY=K=PROGRAMU===============================================
!
!     Modul s ukazkovymi vektorovymi poli pro testovani a ukazku
!     souboru subroutin FORTCONT
!
!     Autor:      Martin Isoz
!     Datum:      31/10/2013
!     System:     Linux 3.2.0-54-generic
!     Licence:    FreeBSD (viz. README.txt)
!
!=====ZAVEDENI=MODULU=SE=SUBROUTINAMI===================================
!
      module Vecfields
      implicit none
      private
      public      :: vecField1P
      contains
!
!=====SUBROUTINES=======================================================
!
            subroutine vecfield1P( xVec,RSVec,pars)                     !vypocet pravych stran vektoroveho pole
!
!     vektorove pole pro kontinuaci podle jednoho parametru, prevzato z
!     textu "korupce v demokraticke spolecnosti" - propagacni material
!     programu MATCONT
!
!     jedna se o model vyduti nosniku pri zatizeni - hledame stacionarni
!     stavy tohoto systemu v zavislostni na zatizeni (lambda). dalsim
!     parametrem pole potom je prirozene nedokonalost materialu, delta.
!
!     xVec  ... vektor vstupu (stav. prom. 1 parametr) 
!     RSVec ... vektor pravych stran
!     pars  ... pevne volene hodnoty parametru        
            implicit none
!           
            double precision,intent(in),dimension(3) :: xVec            !vstupni, stavova promenna
            double precision,intent(in),dimension(3) :: pars            !vstupni promenna - parametry
            double precision,intent(out),dimension(2) :: RSVec          !vystupni promenna
!     lokalni promenne
            double precision  :: MM,CC,De                               !parametry
            double precision  :: uu1,uu2,La                             !stav. prom
!     rozbaleni stavovych promennych
            uu1         = xVec(1)
            uu2         = xVec(2)
            La          = xVec(3)
!     rozbaleni parametru
            MM          = pars(1)
            CC          = pars(2)
            De          = pars(3)
!     vypocetni vystupnich hodnot
            RSVec(1)    = uu2
            RSVec(2)    = 1.0d0/MM*(-(uu1-De) + &
     &2.0d0*La*sin(uu1)) - CC*uu2
            return
      end subroutine vecfield1P
!
!=======================================================================
!
      end module Vecfields
