!
!=====POZNAMKY=K=PROGRAMU===============================================
!
!     program pro testovani souboru subroutin FORTCONT na nekterych
!     vektorovych polich zavedenych v modulu VECFIELDS
!
!     FORTCONT zatim zvlada pouze kontinuaci podle jednoho parametru,
!     konkretne kontinuaci metodou predictor-corrector.
!
!     k dispozici je korektor podle pseudo-delky oblouku, korektorAL a
!     korektor pomoci Moore-Penrose pseudoinverzni matice (hledani
!     bodu na krivce reseni jehoz vzdalenost od predikovaneho bodu je
!     nejmensi mozna), korektorPM
!
!     Autor:      Martin Isoz
!     Datum:      31/10/2013
!     System:     Linux 3.2.0-54-generic
!     Licence:    FreeBSD (viz. README.txt)
!
!-----------------------------------------------------------------------
!
!     vektorove pole ve VECFIELDS musi byt subroutina volatelna:
!     call VECFIELDS ( xVec, RSVec, pars ), kde
!     xVec  ... vektor bodu ve kterych vycisluji vektorove pole
!           ... double precision, dimension(nVar+nParK),intent(in)
!     RSVec ... vektor pravych stran na vystupu
!           ... double precision, dimension(nEqs),intent(out)
!     pars  ... nekontinuovane parametry
!           ... double precision, dimension(nParF),intent(in)
!     a nVar je pocet stavovych promennych, nParK pocet parametru ktere
!     kontinuuji, nEqs je pocet prvku vektoroveho pole (pocet rovnic) a 
!     nParF je pocet fixnich parametru vektoroveho pole
!
!-----------------------------------------------------------------------
!
!     
!     reseni kontinuace je provadeno pomoci subroutiny SOLVER, ktera se
!     vola dne nasledujiciho klice:
!
!     call solver (VECFIELD,PREDIKTOR,KOREKTOR,xVec,ip,&
!    & BMin,BMax,BTrg,pars,(/nDim,nVar/),nParF,g,h,rstep,nprstep,eps,&
!    & maxiter,wrtunt,retcode) , kde
!
!     subroutiny
!     VECFIELD    ... jmeno subroutiny z modulu VECFIELDS
!     PREDIKTOR   ... vyuzivany prediktor pro kontinuaci
!     KOREKTOR    ... korektor vyuzivany pro zpresneni reseni
!                 ... korektorAL nebo korektorPM
!     stavove promenne
!     xVec  ...   pocatecni bod reseni, stavove promenne + kontinuovane
!                 parametry
!           ...   double precision, dimension(nVar+nParK),intent(in)
!     ip    ...   indexy (pozice) parametru v xVec
!           ...   integer(kind=4),dimension(nParK),intent(in)
!     BMin  ...   minimalni pripustne hodnoty stavovych promennych a
!                 kontinuovanych parametru
!           ...   double precision, dimension(nVar+nParK),intent(in)
!     BMax  ...   maximalni pripustne hodnoty stavovych promennych a
!                 kontinuovanych parametru
!           ...   double precision, dimension(nVar+nParK),intent(in)
!     BTrg  ...   cilove hodnoty stavovych promennych a kontinuovanych
!                 parametru na vstupu, skutecne dosazene hodnoty na
!                 vystupu
!           ...   double precision, dimension(nVar+nParK),intent(inout)
!     pars  ...   fixovane parametry reseneho problemu
!           ...   double precision,dimension(nParF),intent(in)
!     n     ...   rozmer reseneho problemu, n(1)=nEqs, n(2)=nVar+nParK
!           ...   integer(kind=4),dimension(2),intent(in)
!     m     ...   pocet fixovanych parametru reseneho problemu, m=nParF
!           ...   integer(kind=4),dimension(1),intent(in)
!
!     pomocne promenne subroutin FORTCONT
!     g     ...   krok pro numerickou derivaci, neupravovan, zatim se
!                 pro vycislovani derivaci pouziva symetricke formule 2.
!                 radu, coz mozna bude zmeneno, jelikoz je to drahe
!           ...   double precision,dimension(1),intent(in)
!     h     ...   vektor s kroky pro numerickou integraci (pohyb
!                 prediktoru), (/h0,hMin,hMax/) neboli pocatecni,
!                 minimalni a maximalni krok. velikost kroku je v
!                 algorimtu adaptivne upravovana
!           ...   double precision,dimension(3),intent(in)
!     maxstep..   maximalni pocet kroku algoritmu
!           ...   integer(kind=4),dimension(1),intent(inout)
!     nprstep..   pocet pripravnych kroku - kolikakrokovou integraci
!                 pouzivam v prediktoru - kolik musim ukladat tecnych
!                 vektoru a jak dlouho tedy trva, nez mohu zacit
!                 se samotnym behem algoritmu
!           ...   integer(kind=4),dimension(1),intent(in)
!     eps   ...   pozadovana presnost pro korektora (kdyz jiz upravuje
!                 reseni v norme o mene, nez tuto hodnotu, ukonci
!                 iteracni proces
!           ...   mozna bude zmeneno na vektor s pozadovanou presnosti
!                 pro upravu reseni i rezidual pravych stran
!           ...   double precision,dimension(1),intent(in)
!
!     informacni a IO promenne
!     wrtunt...   jednotka do ktere bude zapisovano reseni, 6 je pro
!                 terminal. tato jednotka musi byt pred volanim programu
!                 otevrena a pripravena k zapisovani
!           ...   integer(kind=4),dimension(1),intent(in)
!     retcode..   navratova hodnota funkce, zakladni diagnostika
!                 0 - dosazeni nektere cilove hodnoty stavove promenne,
!                     x(i) = BTrg(i)
!                 1 - prekroceni maximalni (minimalni) hodnoty stavove
!                     promenne, x(i) >(<) BMax(i)(BMin(i))
!                 2 - prilis maly krok, h < hMin
!                 3 - korektor nezkonvergoval, niter > maxiter
!                 4 - preroceni maximalniho povoleneho poctu kroku
!           ...   integer(kind=4),dimension(1),intent(out)
!
!-----------------------------------------------------------------------
!
!     pred volanim resitele je treba pripravit vsechny promenne a
!     otevrit jednotku pro zapisovani
!
!     pro graficke zobrazeni vysledku je vyuzivan modul pro spolupraci
!     FORTRANu s GNUPlotem - gnufor2, jehoz autorem je Alexey Kuznetsov
!     pro graficka zobrazeni je tedy treba mit GNUPlot nainstalovany v
!     systemu a dale se pak ridit pravidly pouzivani modulu gnufor2
!
!=====TESTOVACI=PROGRAM=================================================
!
      program testKorektor
!
!=====HLAVICKA==========================================================
!
!     Autor: Martin Isoz
!     Datum: 18.10.2013
!
!     Program pro testovani subroutiny korektor, zatim neni nijak reseno
!     ziskani tecneho vektoru ani spousta dalsich veci
!     
!     kompilace: make/make run/make clean
!
!=======================================================================
      use Vecfields                                                     !nactu vektorova pole
      use FORTCONT                                                      !nactu subroutiny pro numerickou kontinuaci
      use gnufor2
!
      implicit none
!
!=====DEFINICE=INTERFACE================================================
!
!     
!=====DEKLARACE=PROMENNYCH==============================================
!
!     stavove promenne a parametry
      double precision,dimension(2):: uu1,uu2,La                        !stav. promenne, 3 pocatecni body
      double precision            :: MM,CC,De                           !parametry
      integer(kind=4),parameter   :: nDim=2,nParF=3,nVar=3              !pocet eqs, pars, vars
      integer(kind=4),parameter   :: ip(nVar-nDim)=(/3/)                !umisteni par ve stav prom
!     promenne subroutiny SOLVER
      double precision            :: xVec(nVar),pars(nVar)              !stav. prom a vektor parametru
      double precision,parameter  :: g=1.0d-4,eps=1.0d-10               !krok pro num. derivaci a pozadovana presnost res.
      double precision,parameter  :: h(3)=(/1.0d-4,1.0d-9,3.0d-2/)      !kroky pro num. integraci, starovaci, min, max
      double precision            :: BMin(nVar),BMax(nVar),BTrg(nVar)   !hranice oblasti zajmu kontinuace
      integer(kind=4)             :: retcode,maxstep=1000               !navratovy kod solver,max/akt pocet kroku
      integer(kind=4),parameter   :: maxiter=15                         !max pocet kroku, max pocet iteraci korektora
      integer(kind=4),parameter   :: nprstep=5                          !pocet uchovavanych tecnych vektoru
      integer(kind=4),parameter   :: wrtunt=6                           !zapisovaci jednotka, jak vypsat na obrazovku?
!     promenne volatele
      integer(kind=4)             :: cyccount,rstep                     !pocitadlo cyklu, vysledne pocty kroku
!     junk
      double precision            :: tVec(1,nVar)
!     stringy
!
!
!=====EXECUTABLE=PART===================================================
!
!-----priprava-pro-reseni-ulohy-----------------------------------------
!     zadani hodnot pevnych parametru
      MM          = 1.0d0                                               !1.0e0
      CC          = 1.0d0                                               !1.0e0
      De          = 5.0d-2                                              !0.05d0
!
!     body ve kterych vyhodnocuji
      uu1(2)         = -1.883798570d0                                   !restart, dopocteni odvetveni, od spodu, -1.883798570d0
      uu2(2)         = 0.0d0                                            !nastrel byl dopocten predem num. res. v(x) = 0, 0.0d0
      La(2)          = 0.999d0                                          !nemuzu zacit primo z 1.0d0 kvuli podmince na pMax, 0.990d0
      uu1(1)         = 0.0d0
      uu2(1)         = 0.0d0
      La(1)          = 0.0d0
!
!     priprava promennych pro volani SOLVER
      pars = (/MM,CC,De/)                                               !fixovane parametry vecfield1P
!      xVec = (/uu1,uu2,La/)                                             !stavove promenne - jedno volani
      BMin  = (/-1.0d3, -1.0d3, -1.0d3/)
      BMax  = (/1.0d3, 1.0d3, 1.0d0/)
!     priprava vystupu
!      open(1, file='dataES.dat', status='new')                          !otevru si soubor
      
      
      tVec(1,:)  = (/1.0d0,1.0d0,1.0d0/)
!
!-----volani-subroutiny-SOLVER------------------------------------------
!
      do cyccount=1,2,1
!     priprava promennych
      xVec = (/uu1(cyccount),uu2(cyccount),La(cyccount)/)               !stavove promenne - vice volani
      BTrg  = BMax
      rstep = maxstep
!     volani
      call solver (vecfield1P,prediktorPT,korektorPM,xVec,ip,&
     & BMin,BMax,BTrg,pars,(/nDim,nVar/),nParF,g,h,rstep,nprstep,eps,&
     & maxiter,wrtunt,retcode)
!
!     vypis retcode
      write(unit=*,fmt=*) 'Return code: ',retcode

!
      end do
!
!-----graficka-zobrazeni------------------------------------------------
!
!      call run_gnuplot('commands1.txt')                                 !Rov. stavy, output do .ps
!
!     zahozeni souboru s vysledky (chci mit moznost jet rovnou dokola)
!      close(1, status='delete')
!
!=======================================================================
!
      end program testKorektor
