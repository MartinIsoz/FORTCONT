!
!=====POZNAMKY=K=PROGRAMU===============================================
!
!     Soubor subroutin pro numerickou kontinuaci s vyuzitim jazyka
!     FORTRAN, konkretne povetsinou standardu F95, i kdyz nektere
!     pouzivane funkce (MATMUL...) jsou z pozdejsich vydani.
!
!     Autor si je vedom, existence dalsich, mnohem pokrocilejsich baliku
!     pro numerickou kontinuaci, napriklad AUTO ci MATCONT.
!
!     Ucelem tohoto software neni poskytnout uceleny soubor nastroju pro
!     studium bifurkaci dynamickych systemu. Autor se snazi poskytnout
!     nazornou ukazku algoritmizace nekterych zakladnich numerickych
!     metod tohoto odvetvi vyuzitelnou pro demonstracni ucely, pochopeni
!     principu kontinuace reseni podle parametru, identifikace
!     bifurkacnich bodu a dalsich elementarnich partii tohoto rozsahleho
!     odvetvi matematiky.
!
!     Dale by mel byt software praktickou ukazkou vyuziti jazyka FORTRAN
!     pro reseni matematickych a inzenyrskych uloh.
!
!     Prezentovane algoritmy nejsou psany s durazem na rychly prubeh,
!     ale na nazornost a prehlednost.
!
!     V programech jsou vyuzivany algoritmy obsazene v knihovnach LAPACK
!     a BLAS, takze pri samostatne kompilaci software je nutne mit tyto
!     instalovane v systemu a linkovane v kompilatoru
!
!     Autor:      Martin Isoz
!     Datum:      31/10/2013
!     System:     Linux 3.2.0-54-generic
!     Licence:    FreeBSD (viz. README.txt)
!
!     Je nutno:
!     - vyresit ukonceni korektorPM i s ohledem na tecny vektor, je
!       totiz pouzivan i k nalezeni prvniho tecneho vektoru
!     - upravit + ozkouset korektorAL tak, aby vracel tecny vektor
!     - zkusit prediktorAB, je to zbytecne komplikovana subroutina, kdyz
!       budou korektory vracet odhady tecnych vektoru v dalsim bode
!       reseni a navic bude zrusena zavislost na GAUSE
!     - zkontrolovat argumenty poustene v programu do DSGESV
!
!=====ZAVEDENI=MODULU=SE=SUBROUTINAMI===================================
!
      module FORTCONT
      implicit none
      private
      public      :: jacU,detStab,initRES,&
     &               prediktorAB,prediktorPT,&
     &               korektorAL,korektorPM,&
     &               solver
      contains
!
!=====SUBROUTINES=======================================================
!
      subroutine jacU( subrF,xVec,pars,n,m,g,jacB,sk )                   !vraci diferencial zobrazeni (s i bez prom. parametru)
!     subrF ... subroutina vracejici hodnoty vektoru pravych stran
!     xVec  ... bod, ve kterem jacobiho mat. vycisluji (vektor)(n1+n2)
!     pars  ... vektor parametru pro subroutunu subrF
!     n     ... vektor s rozmery problemu, (/integer,integer/),
!               n1 - pocet rovnic
!               n2 - pocet stavovych promennych (vcetne parametru)
!     m     ... pocet (nekont.) parametru subroutiny subrF (integer)
!     g     ... krok pro numerickou derivaci (double precision)
!     jacB  ... vystupni matice, pro ulozeni hodnot (array)
!               rozmer: (n) x (n-size(sk))
!     sk    ... indexy vynechavanych sloupcu radku (/integer,sorted/),v
!               pripade, ze chci pracovat s celou matici, sk(0)
      implicit none
!
      external subrF
!     
!     vstupne-vystupni promenne
      integer(kind=4),intent(in)  :: n(2),m                             !vstupni promenne, rozmer problemu, pocet par. subroutiny
      integer(kind=4),intent(in)  :: sk(:)                              !vstupni promenna, indexy vynechavanych radku
      double precision,intent(in) :: xVec(n(2)),pars(m)                 !vstupni promenne, bod vycisleni, parametry subroutiny
      double precision,intent(out):: jacB(n(1),n(2)-size(sk))           !vystupni promenna, jacobiho matice pro fix. k prom.
      double precision,intent(in) :: g                                  !vstupni promenna, krok num. deriv.
!     lokalni promenne
      double precision  :: e_i(n(2)),f_r(n(1)),f_l(n(1))                !jednotkovy vektor a vektory vycisleni zprava a zleva
      integer(kind=4)   :: i,j,k,dk=0                                   !indexovaci promenne
      integer(kind=4)   :: skk(size(sk)+2)
!
!     vypocetni cyklus - mel by pokryt oba pripady
      skk = (/0,sk,n(2)+1/)                                             !rozsireny vektor pro vynechavani indexu
      dk  = 0                                                           !restartuji zpozdovaci promennou
      do k=1, size(skk)-1, 1
            do i=skk(k)+1, skk(k+1)-1, 1
                  do j=1, n(2), 1
                        e_i(j) = 0.0d0                                  !inicializace vektoru
                  end do
                  e_i(i) = 1.0d0                                        !e_i - jednotkovy vektor s 1 na i-tem miste
                  call subrF(xVec + g*e_i,f_r,pars,n(2),m)              !vyhodnoceni funkce zprava
                  call subrF(xVec - g*e_i,f_l,pars,n(2),m)              !vyhodnoceni funkce zleva
                  jacB(1:n(1),i+dk) = 1.0d0/(2.0d0*g)*(f_r-f_l)         !vypocet i-teho sloupce jacobianu, soumer. dvojbod. formule
            end do
            dk  = dk-1                                                  !musim si odecist umisteni, abych nepretekl z matice
      end do
!
      return
      end subroutine jacU
!
!=======================================================================
!
      logical function detStab (subrF,xVec,pars,n,m,g,fp)
!
!=====HLAVICKA==========================================================
!
!     Funkce na urceni stability rovnovazneho stavu. v xVec zavola
!     JACU s fixovanym parametrem (fp-ta stavova promenna), obecne muze
!     byt vektor, je-li F:R(n+m)->Rn, ale zatim vyuzivana pouze pro
!     F:R(n+1)->Rn. Takto je ziskana (n(1) x n(1)) matice, pro niz jsou
!     nasledne pomoci suboutiny DGEES vypoctena vlastni cisla. Nakonec
!     je z vlastnich cisel usouzeno na stabilitu RS (zatim pouze
!     stable/unstable, dalsi studium, napriklad je-li to bod HB ci jinak
!     bude mozna nasledovat.
!
!
!     Program vola subroutiny:
!     SUBRF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     JACU  ...   subroutina vracejici diferencial zobrazeni Rm->Rn,
!                 je mozne fixovat nektere promenne, coz vyusti v matici
!                 R(n x k), kde k je m-nFix. podrobneji v testjacU
!     DGEES ...   subroutina pro vypocet vlastnich cisel matice, mela by
!                 byt obsazena v knihovne lapack. slouzi k urceni
!                 stability nalezeneho rovnovazneho stavu. SNAD, podari
!                 li se ji spravne zavolat a vyuzit
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     xVec  ...   na vstupu vektor s predchozim resenim, xi
!           ...   na vystupu vektor s predikovanym resenim, X0
!     pars  ...   vektor parametru pro subroutinu subrF
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     g     ...   velikost kroku pro numerickou derivaci
!     fp    ...   index sledovaneho parametru ve stavove promenne -
!                 fixovana promenna pri volani JACU na vstupu
!     stab  ...   stabilita RS v xVec
!
!     Pozn: Bylo by idealni, kdyby se podarilo subroutinu DGEES upravit
!           tak, aby prestala pocitat hned ve chvili, kdy narazi na
!           prvni eigenvalue s RE > 0, ale to by zatim byla zbytecna
!           prace
!
!     Autor:      Martin Isoz
!     Datum:      21.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
      external subrF
!     vstupni promenne
      integer(kind=4),intent(in)    :: n(2),m,fp(:)                     !rozmer problemu, pocet pars, indexy pars ve stav prom
      double precision,intent(in)   :: g                                !velikost kroku pro derivaci
      double precision,intent(in)   :: pars(m)                          !pars pro subroutinu subrF
      double precision,intent(in)   :: xVec(n(2))                       !vektor s bodem ve kterem urcuji stabilitu
!     vystupni promenne
!     promenne pro JACU
      double precision              :: jac(n(1),n(1))                   !jacobiho matice systemu
!     promenne pro DGEES - zatim cert vi co, ale bude toho hodne
      integer(kind=4)               :: sdim                             !switch pro sortovani vl.c.
      double precision              :: wr(n(1)),wi(n(1))                !vektory pro ukladani RE a IM casti vl.c.
      double precision              :: vs(1,1)                          !ortogonalni matice (Z) shurrovych vekt. A=ZTZ*
      integer(kind=4)               :: ldvs=1                           !radkovy rozmer vs
      double precision              :: work(1,4*n(1))                   !pracovni prommena pro hledani vl. c.
      integer(kind=4)               :: lwork                            !pocet prvku work
      double precision              :: bwork(1,4*n(1))                  !pracovni promenna pro sortovani
      integer(kind=4)               :: info
!     jine interni promenne
      logical                       :: stab
!
!
!=====EXECUTABLE=PART===================================================
!
!-----konstrukce-matice-linearizovaneho-systemu-------------------------
      call jacU( subrF,xVec,pars,n,m,g,jac,fp )                         !vystup je matice linearizovaneho systemu
!
!-----nalezeni-vlastnich-cisel-matice-linearizovaneho-systemu-----------
!     priprava promennych
      lwork = 4*n(1)
      call dgees('N', 'S', SELECT, n(1), jac, n(1), sdim, wr, wi, vs,&
     & ldvs, work, lwork, bwork, info)
!
!-----vyhodnoceni-behu-a-vraceni-hodnoty--------------------------------
      if (sdim .gt. 0) then                                             !nalezeno alespon jedno vl.c. s RE > 0 -> UNSTAB
            stab = .false.
      else
            stab = .true.
      end if
      detStab = stab
      return
!-----modul-pro-dgees---------------------------------------------------
      contains
            logical function select(ar, ai)
!     .. Scalar Arguments ..
!
!     Logical function SELECT for use with DGEES
!
!     Returns the value .TRUE. if the imaginary part of the eigenvalue
!     (AR + AI*i) is zero, i.e. the eigenvalue is real
!
!     vstupni promenne
            double precision,intent(in) :: ai, ar
!     lokalni promenne
            logical d
            double precision junk
!     vypocetni cast
            junk = ai                                                   !ai me nezajima a potlacim tak varovani unused variable
            if (ar .ge. 0.0d0) then
                  d = .true.
            else
                  d = .false.
            end if
            
!
            select = d
!
            return
            end Function select
      end function detStab
!
!=======================================================================
!
      subroutine initRES (subrF,xVec,pars,tVec,n,m,g,eps,maxiter)
!
!=====HEAD==============================================================
!
!     subroutina pro vypocteni prvotniho tecneho vektoru. v pripade, ze
!     neni pouzivan KOREKTORPM, nebyl by tento vektor pri startu
!     algoritmu k dispozici -> je nutno jej dopocist
!
!     vypocet je zalozen na reseni soustavy rovnic
!     [Fx;tVec*]tVec = [0;1] iterativne pomoci neco-jako-newtona
!     projizdi souradny system podle bazi [1,0,...,0],[0,1,0,...],...
!     dopokavad nezkonverguje newton a nenabidne prvnotni reseni
!
!     Program vola subroutiny:
!     SUBRF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     JACU  ...   subroutina vracejici diferencial zobrazeni Rm->Rn,
!                 je mozne fixovat nektere promenne, coz vyusti v matici
!                 R(n x k), kde k je m-nFix. podrobneji v testjacU
!     KOREKTORPM  subroutina na zpresneni bodu v okoli krivky reseni,
!                 vraci i tecny vektor v danem bode
!
!     Vstupni argumenty programu:
!     SUBRF ...   subroutina vracejici vektor pravych stran/rezidui
!     xVec  ...   vektor s resenim v pocatecnim bode, vstup
!     pars  ...   vektor parametru pro subroutinu subrF
!     tVec  ...   tecny vektor ke krivce v x0, vystup
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     g     ...   velikost kroku pro numerickou derivaci
!     eps   ...   pozadovana presnost, ||xEnd - xEnd-1|| <= eps
!     maxiter..   maximalni pocet iteraci metody, navratova hodnota je
!                 skutecny pocet iteraci newtona   
!
!     Autor:      Martin Isoz
!     Datum:      10/29/2013
!     System:     Linux 3.2.0-54-generic
!     Licence:    FreeBSD (viz. README.txt) 
!     
      implicit none
!
!=====VARIABLES=DECLARATION=============================================
!
!     externi subroutiny
      external subrF                                                    !vektor. pole, vypocet diferencialu zobrazeni
!     vstupni parametry
      integer(kind=4),intent(in)    :: n(2),m                           !rozmer problemu, pocet pars, indexy pars ve stav prom
      double precision,intent(inout):: xVec(n(2))                       !bod ve kterem zjistuji vektor
      double precision,intent(in)   :: pars(m)                          !parametry pro subrF
      double precision,intent(in)   :: eps,g                            !pozadovana presnost pro newtona, krok pro numerickou derivaci
      double precision,intent(out)  :: tVec(n(2),1)                     !tecny vektor ke krivce
      integer(kind=4),intent(inout) :: maxiter                          !maximalni pocet iteraci pro newtona,skutecny pocet iteraci
!     lokalni promenne - jacU
      double precision  :: jac(n(1),n(2))                               !matice diferencialu zobrazeni
      integer(kind=4)   :: jSW(0)                                       !switch pro vynechavani v jacU
!     lokalni promenne - hledani tecneho vektoru
      double precision  :: sMat(n(2),n(2)),RHS(n(2),n(2)-n(1))          !matice resene soustavy, vektor pravych stran
      double precision  :: tVecO(n(2),n(2)-n(1))                        !aktualizovany tecny vektor
      integer(kind=4)   :: iterOut,iterIn                               !pocitadla iteraci
      
!
!=====EXECUTABLE=PART===================================================
!
!     vypocet diferencialu zobrazeni subrF v xVec
      call jacU( subrF,xVec,pars,n,m,g,jac,jSW )                        !vystup je matice linearizovaneho systemu
!
!     newton - iniciace
      do iterOut=1,n(2),1                                               !pro kazdy jednotkovy vektor - bazi
            tVecO                   = 0.0d0                             !iterOut-ty bazovy vektor
            tVecO(iterOut,n(2)-n(1))= 1.0d0                             !normovany vektor
            RHS                     = 0.0d0                             !iniciace vektoru pravych stran
            RHS(n(2),n(2)-n(1))     = 1.0d0
            sMat(1:n(1),:)          = jac                               !prvnich n(1) radku je matice jacobianu
            sMat(n(1)+1:n(2),:)     = transpose(tVecO)                  !zbytek je tecny vektor
            iterIn                  = maxiter                           !restartuji pocitadlo iteraci
            call korektorPM (subrF,xVec,tVecO,pars,n,m,g,eps,maxiter)
            if  (norm2(matmul(sMat,tVecO) - RHS) .lt. eps ) then
                  go to 666
            end if
      end do
666   tVec    = tVecO                                                   !ulozim vysledek
      maxiter = iterIn
      return
      end subroutine initRES
!
!=======================================================================
!
      subroutine prediktorAB (subrF,xVec,tMat,pars,n,m,nSol,g,h)
!
!=====HLAVICKA==========================================================
!
!     Prediktor pro metodu prediktor-korektor konstrukce evolucnich
!     diagramu. pro nalezeni tecneho vektoru a promenne vuci ktere je
!     reseni regularni je vyuzita kubickova subroutina gause
!
!     prediktor pro dany system vraci tecny vektor v bode x0
!
!     Program vola subroutiny:
!     gause ...   kubickova subroutina, vraci tecny vektor podle jedne
!                 ze souradnic, vuci ktere je reseni regularni a index
!                 teto souradnice.
!           ...   subroutina je psana v F77 a bylo by vhodne ji prepsat
!                 do F90 a zaradit do tohoto modulu
!     subrF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     jacU  ...   subroutina vracejici diferencial zobrazeni Rm->Rn,
!                 je mozne fixovat nektere promenne, coz vyusti v matici
!                 R(n x k), kde k je m-nFix. podrobneji v testjacU
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     xVec  ...   na vstupu vektor s predchozim resenim, xi
!           ...   na vystupu vektor s predikovanym resenim, X0
!     tMat  ...   na vstupu matice s predchozimi 5 tecnymi vektory
!                 [t(i-1);...;t(i-5)]
!           ...   na vystupu matice s poslednim dopoctenym tecnym,
!                 vektorem na prvnim radku,[t(i);...;t(i-4)]
!     pars  ...   vektor parametru pro subroutinu subrF
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     nSol  ...   pocet reseni pritomnych v xMat a tMat (pocet obecne
!                 nenulovych radku)
!     g     ...   velikost kroku pro numerickou derivaci
!     h     ...   velikost kroku pro numerickou integraci AB
!
!     Pozn: tVec0 je vracen kvuli penrose-mooreovu korektoru, korektor
!           zalozeni na pseudodelkce oblouku jej nevyuzije
!     Pozn: na to, ze od GAUSE potrebuju pouze tecny vektor (nadrovinu)
!           toho pocita zbytecne mnoho navic...
!     Pozn: Prace s promennymi je sice psana obecne pro R(n+m)->Rn, ale
!           kolkolem pridavam pouze jednu rovnici a tedy potrebuji
!           kontinuovat pouze podle jednoho parametru R(n+1)->Rn
!
!     Autor:      Martin Isoz
!     Datum:      19.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
      external subrF
!     vstupni promenne
      integer(kind=4),intent(in)    :: n(2),m,nSol                      !rozmer problemu, pocet pars, pocet res PO
      double precision,intent(in)   :: g,h                              !velikost kroku pro derivaci, integraci
      double precision,intent(in)   :: pars(m)                          !pars pro subroutinu subrF
!     vstupne-vystupni a vystupni promenne
      double precision,intent(inout):: xVec(n(2))                       !matice s predchozim/naslednym resenim
      double precision,intent(inout):: tMat((n(2)-n(1))*5,n(2))         !5 tecnych vektoru/nadrovin v xO
!     lokalni promenne
      double precision  :: jac(n(1),n(2))                               !matice pro vypocet diferencialu -> gause
      double precision  :: tVec(n(2)-n(1),n(2))                         !tecny vektor v x0 -> gause
!     promenne pro integraci AB
      double precision  :: AB(5,5)                                      !matice koeficientu
      double precision  :: PV(n(2))                                     !pomocna promenna
!     JACU
      integer(kind=4)   :: jSW(0)                                       !zero sized integer pro preskakovane prom. pri volani jacU
!     GAUSE
      integer(kind=4)   :: ierr=1,kfix=1                                !kontrola singularity, fix. promenna
      double precision  :: kor(n(2))                                    !nevyuzita promenna pro korekci vysledku
!
      integer(kind=4)   :: i                                            !pocitadla cyklu
!
!
!=====EXECUTABLE=PART===================================================
!
!-----koeficienty pouzitych metod---------------------------------------
!     koeficienty pro adams-bashforthovu integraci do 5. radu
      AB(1,:)=(/   1.0d0,    0.0d0,   0.0d0,    0.0d0,  0.0d0/)/  1.0d0 !koeficienty A-B, 1. rad
      AB(2,:)=(/   3.0d0,   -1.0d0,   0.0d0,    0.0d0,  0.0d0/)/  2.0d0 !koeficienty A-B, 2. rad
      AB(3,:)=(/  23.0d0,  -16.0d0,   5.0d0,    0.0d0,  0.0d0/)/ 12.0d0 !koeficienty A-B, 3. rad
      AB(4,:)=(/  55.0d0,  -59.0d0,  37.0d0,   -9.0d0,  0.0d0/)/ 24.0d0 !koeficienty A-B, 4. rad
      AB(5,:)=(/1901.0d0,-2774.0d0,2616.0d0,-1274.0d0,251.0d0/)/720.0d0 !koeficienty A-B, 5. rad
!-----samotna vypocetni cast--------------------------------------------
      kor = 0                                                           !nutno inicializovat, jinak bordel s detStab, ZBAVIT SE GAUSE
      call jacU( subrF,xVec,pars,n,m,g,jac,jSW )                        !ziskam soustavu n rovnic pro n+1 nezn.
!     volam kubickuv algoritmus gause - pocita toho zbytecne moc/malo
      call gause( n(1),n(2),n(1),jac,kor,tVec,kfix,ierr )               !zavolam kubickovu gause, korektor je zde nulovy vektor
!     integrace dalsiho kroku - pomoci AB najdu presnejsi odhad tVec
!      write (unit=*,fmt=*) 'tMat before'
!      do i=1,size(tMat,1)
!            write(*,'(20G12.4)') tMat(i,:)
!      end do
      do i=nSol-1,1, -1                                                 !dostupna res presunu o n(2)-n(1) souradnic niz
                  tMat((i)*(n(2)-n(1))+1:(i+1)*(n(2)-n(1)),:) = &
     &       tMat((i-1)*(n(2)-n(1))+1:(i)*(n(2)-n(1)),:)
      end do
      tMat(1:n(2)-n(1),:) = tVec                                        !nahoru ulozim posledni vysledek
!      write (unit=*,fmt=*) 'tMat after'
!      do i=1,size(tMat,1)
!            write(*,'(20G12.4)') tMat(i,:)
!      end do
      do i=1, n(2), 1
            PV(i) = sum(AB(nSol,:)*tMat(:,i))                           !naplnim pomocny vektor
      end do
      xVec = xVec + h*PV
      return
      end subroutine prediktorAB
!
!=======================================================================
!      
      subroutine prediktorPT (subrF,xVec,tMat,pars,n,m,nSol,g,h)
!
!=====HLAVICKA==========================================================
!
!     Prediktor pro metodu prediktor-korektor konstrukce evolucnich
!     diagramu. tento prediktor potrebuje jiz pripravenou tMat a posune
!     xVec na zaklade integrace AB formuli radu nSol
!
!     prediktor pro dany system vraci tecny vektor v bode x0
!
!     Program vola subroutiny:
!     ZADNE externi subroutiny nejsou volany
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     xVec  ...   na vstupu vektor s predchozim resenim, xi
!           ...   na vystupu vektor s predikovanym resenim, X0
!     tMat  ...   na vstupu matice s predchozimi 5 tecnymi vektory
!                 [t(i-1);...;t(i-5)]
!           ...   na vystupu matice s poslednim dopoctenym tecnym,
!                 vektorem na prvnim radku,[t(i);...;t(i-4)]
!           ...   oproti prediktorAB nemuze byt na vstupu nulovy vektor,
!                 je xVec(OUT) = xVec(IN)
!     pars  ...   vektor parametru pro subroutinu subrF, zde tento 
!                 argument neni vyuzit, ale je mezi volanymi, pro
!                 podrzeni stejnych volani pro vsechny korektory
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     nSol  ...   pocet reseni pritomnych v xMat a tMat (pocet obecne
!                 nenulovych radku)
!     g     ...   velikost kroku pro numerickou derivaci, zde tento 
!                 argument neni vyuzit, ale je mezi volanymi, pro
!                 podrzeni stejnych volani pro vsechny korektory
!     h     ...   velikost kroku pro numerickou integraci AB
!
!     Pozn: tVec0 je vracen kvuli penrose-mooreovu korektoru, korektor
!           zalozeni na pseudodelkce oblouku jej nevyuzije
!     Pozn: na to, ze od GAUSE potrebuju pouze tecny vektor (nadrovinu)
!           toho pocita zbytecne mnoho navic...
!     Pozn: Prace s promennymi je sice psana obecne pro R(n+m)->Rn, ale
!           kolkolem pridavam pouze jednu rovnici a tedy potrebuji
!           kontinuovat pouze podle jednoho parametru R(n+1)->Rn
!
!     Autor:      Martin Isoz
!     Datum:      22.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
      external subrF
!     vstupni promenne
      integer(kind=4),intent(in)    :: n(2),m,nSol                      !rozmer problemu, pocet pars, pocet res PO
      double precision,intent(in)   :: g,h                              !velikost kroku pro derivaci, integraci
      double precision,intent(in)   :: pars(m)                          !pars pro subroutinu subrF
!     vstupne-vystupni a vystupni promenne
      double precision,intent(inout):: xVec(n(2))                       !matice s predchozim/naslednym resenim
      double precision,intent(inout):: tMat((n(2)-n(1))*5,n(2))         !5 tecnych vektoru/nadrovin v xO
!     lokalni promenne
      double precision  :: gg,pp(m)                                     !pouze na zbaveni se warning unused dummy
!     promenne pro integraci AB
      double precision  :: AB(5,5)                                      !matice koeficientu
      double precision  :: PV(n(2))                                     !pomocna promenna
!
      integer(kind=4)   :: i                                            !pocitadla cyklu
!
!
!=====EXECUTABLE=PART===================================================
!
!-----koeficienty pouzitych metod---------------------------------------
!     koeficienty pro adams-bashforthovu integraci do 5. radu
      AB(1,:)=(/   1.0d0,    0.0d0,   0.0d0,    0.0d0,  0.0d0/)/  1.0d0 !koeficienty A-B, 1. rad
      AB(2,:)=(/   3.0d0,   -1.0d0,   0.0d0,    0.0d0,  0.0d0/)/  2.0d0 !koeficienty A-B, 2. rad
      AB(3,:)=(/  23.0d0,  -16.0d0,   5.0d0,    0.0d0,  0.0d0/)/ 12.0d0 !koeficienty A-B, 3. rad
      AB(4,:)=(/  55.0d0,  -59.0d0,  37.0d0,   -9.0d0,  0.0d0/)/ 24.0d0 !koeficienty A-B, 4. rad
      AB(5,:)=(/1901.0d0,-2774.0d0,2616.0d0,-1274.0d0,251.0d0/)/720.0d0 !koeficienty A-B, 5. rad
!     borden na zbaveni se warning unused dummy argument
      gg = g
      pp = pars
!-----samotna vypocetni cast--------------------------------------------
      do i=1, n(2), 1
            PV(i) = sum(AB(nSol,:)*tMat(:,i))                           !naplnim pomocny vektor
      end do
      xVec = xVec + h*PV
      return
      end subroutine prediktorPT
!
!=======================================================================
!      
      subroutine korektorAL (subrF,x0,tVec0,pars,n,m,h,eps,maxiter)
!
!=====HLAVICKA==========================================================
!
!     Korektor pro metodu prediktor-korektor konstrukce evolucnich
!     diagramu, zalozen na n-rozmerne newtonove metode, vyuzivajici
!     knihovny llapack k reseni soustav lin. rovnic
!     Jac*Deltax = -f
!
!     korektor je zalozen na kontinuaci podle pseudodelky oblouku,
!     pseudoarclength continuation, popsane v manualu programu MATCONT
!
!     Program vola subroutiny:
!     subrF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     DSGESV...   subroutina z knihovny BLAS, reseni SLAR, iterative
!                 solver (single prec + refinement .or. double prec
!                 factorization)
!     jacU  ...   subroutina vracejici diferencial zobrazeni Rm->Rn,
!                 je mozne fixovat nektere promenne, coz vyusti v matici
!                 R(n x k), kde k je m-nFix. podrobneji v testjacU
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     x0    ...   pocatecni nastrel
!                 v teto promenne je ulozen vysledny vektor reseni
!     tVec0 ...   tecny vektor ke krivce v bode x0
!     pars  ...   vektor parametru pro subroutinu subrF
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     h     ...   velikost kroku pro vycisleni num. derivaci
!     eps   ...   pozadovana presnost, ||xEnd - xEnd-1|| <= eps
!     maxiter..   maximalni pocet iteraci metody, navratova hodnota je
!                 skutecny pocet iteraci newtona
!
!     Pozn: Oproti predchozimu korektoru zalozenem na algoritmu DERPAR
!           (kubicek) jiz nevynechavam zadnou promennou, ale naopak,
!           pridavam rovnici (parametrizace podle delky oblouku)
!     Pozn: Prace s promennymi je sice psana obecne pro R(n+m)->Rn, ale
!           kolkolem pridavam pouze jednu rovnici a tedy potrebuji
!           kontinuovat pouze podle jednoho parametru R(n+1)->Rn
!
!     Autor:      Martin Isoz
!     Datum:      18.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
      external subrF
!     vstupni promenne
      integer(kind=4),intent(in)    :: n(2),m                           !rozmer problemu, ind. vyn. prom., pocet pars
      double precision,intent(inout):: x0(n(2))                         !vektor s poc. nastrelem a vyslednym resenim
      double precision,intent(inout):: tVec0(n(2)-n(1),n(2))            !tecny vektor/tecna nadrovina v xO
      double precision,intent(in)   :: h,eps,pars(m)                    !velikost kroku, presnost, vektor parametru do subrF
      integer(kind=4),intent(inout) :: maxiter                          !max. pocet iteraci
!     dalsi promenne pro NEWTONA
      double precision  :: xOld(n(2)),xNew(n(2))                        !pomocne pracovni promenne, DeltaX, vycisleni funkcnich hodnot
      double precision  :: Dx(n(2)),FVec(n(1))                          !pomocne pracovni promenne, DeltaX, vycisleni funkcnich hodnot
      double precision  :: delta                                        !sledovani prubezne presnosti
      integer(kind=4)   :: iter                                         !pocitadlo iteraci
!      integer(kind=4)   :: i                                            !indexovadlo pro vypisy
      double precision  :: jac(n(1),n(2))                               !matice pro vypocet jacobiho matice
      double precision  :: RHS(n(2)),dRHS(n(2),n(2))                    !prave strany metody a jejich diferencial
!     tlumici faktor Newtona
!      double precision  :: FN(n(1)),FO(n(1))                            !vycisleni fn hodnot pro xNew a xOld
      double precision  :: TLkoef                                       !tlumici faktor
!     DSGESV
      double precision  :: swork(n(2),1)
      double precision  :: work(n(2),n(2))
      integer(kind=4)   :: pivot(n(2),n(2)),iter2,info
!     jacU
      integer(kind=4)   :: jSW(0)                                       !zero sized integer pro preskakovane prom. pri volani jacU
!
!
!=====EXECUTABLE=PART===================================================
!
!     naplneni pomocnych promennych pro newtona
      iter        = 0
      delta       = eps + 1.0d0
      xOld        = x0
      TLkoef      = 1.0d0
!     
!     vypis nastrelu
!     write(unit=*,fmt=*) 'NASTREL====================================='
!     write(unit=*,fmt=*) 'k', iter, 'x(k)', xOld, 'delta', delta& 
!     &     , 'dumpF', TLkoef
!     write(unit=*,fmt=*) 'ITERACE====================================='
!     vypocetni cyklus newtonovy metody
      do while( eps .lt. delta .and. iter .le. maxiter )
!
            call jacU( subrF,xOld,pars,n,m,h,jac,jSW )                  !spocitam diferencial puvodniho zobrazeni
!
            call subrF(xOld,FVec,pars)                                  !vypocet pravych stran puvodniho zobrazeni v danem bode
!
            RHS         = 0.0d0                                         !restartovat right hand side vektor
            RHS(1:n(1)) = FVec                                          !naplnit prvnich n(1) pozic rezidui0
!
            dRHS                = 0.0d0                                 !restartovat diferencial RHS vektoru
            dRHS(1:n(1),:)      = jac                                   !prvnich n(1) radku je diferencial puvodniho zobrazeni
            dRHS(n(1)+1:n(2),:) = tVec0                                 !n(1)+1 - n(2) radek je tecna nadrovina
!
!     kontrolni vypis 1
!            write (unit=*,fmt=*) 'dRHS'
!            do i=1,size(dRHS,1)
!                  write(*,'(20G12.4)') dRHS(i,:)
!            end do
!            write (unit=*,fmt=*), '-RHS'
!            write(*,'(20G12.4)') -RHS
!
            call DSGESV(n,1,dRHS,n(2),pivot,-RHS,n,Dx,n,swork,work,&    !reseni soustavy lin. rovnic
     &                  iter2,info)
!     kontrolni vypis 2
!            write (unit=*,fmt=*), 'Dx'
!            write(*,'(20G12.4)') Dx
!
            xNew = xOld + TLkoef*Dx                                     !vypocet nove iterace
!
!!     uprava tlumiciho faktoru
!            call subrF(xNew,FN,pars)                                    !vycislim residualy funkci v xNew a xOld
!            call subrF(xOld,FO,pars)
!
!           if( sqrt(sum(FN**2.0d0)) .lt. sqrt(sum(FO**2.0d0)) ) then   !pri uspechu preulozim xNew, jinak zpulim tlumici faktor
!                 delta = dsqrt(sum(xNew**2.0d0 - xOld**2.0d0))         !vypocet rozdilu xNew a xOld (norma rozdilu)
!                 xOld = xNew
!           else
!                 TLkoef = TLkoef/2.0d0
!                 delta = delta
!                 xOld = xOld
!           end if
!
            delta = dsqrt(sum((xNew - xOld)**2.0d0))
            xOld = xNew
!
            iter = iter + 1                                             !zvyseni pocitadla iteraci
!
!     vypis iterace
!           write(unit=*,fmt='(20G12.4)') 'k', iter, 'x(k)', xNew , 'delta',delta&
!     &     , 'dumpF', TLkoef
!
      end do
      x0      = xOld                                                    !ulozim vysledek
      maxiter = iter
      RHS              = 0.0d0                                          !restart pravych stran
      RHS(n(1)+1:n(2)) = 1.0d0
      call DSGESV(n(2),n(2)-n(1),dRHS,n(2),pivot,RHS,&
     & n(2),tVec0,n(2),swork,work,iter2,info)                           !reseni soustavy lin. rovnic - dopocteni noveho tecneho vekt.
      return
      end subroutine korektorAL
!
!=======================================================================
!
      subroutine korektorPM (subrF,x0,tVec0,pars,n,m,h,eps,maxiter)
!
!=====HLAVICKA==========================================================
!
!     Korektor pro metodu prediktor-korektor konstrukce evolucnich
!     diagramu, zalozen na n-rozmerne newtonove metode, vyuzivajici
!     knihovny llapack k reseni soustav lin. rovnic
!     Jac*Deltax = -f
!
!     korektor je zalozen na penrose-moorove pseudoinverzni matici, je
!     upravovana nejenom hodnota nastrelu, ale i tecneho vektoru.
!     bylo by vhodne prozkoumat myslenku, nemohl-li by tento tencny
!     vektor slouzit primo jako tecny vektor v prediktoru
!
!     Program vola subroutiny:
!     subrF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     DSGESV...   subroutina z knihovny BLAS, reseni SLAR, iterative
!                 solver (single prec + refinement .or. double prec
!                 factorization)
!     jacU  ...   subroutina vracejici diferencial zobrazeni Rm->Rn,
!                 je mozne fixovat nektere promenne, coz vyusti v matici
!                 R(n x k), kde k je m-nFix. podrobneji v testjacU
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     x0    ...   pocatecni nastrel
!                 v teto promenne je ulozen vysledny vektor reseni
!     tVec0 ...   tecny vektor ke krivce v bode x0
!     pars  ...   vektor parametru pro subroutinu subrF
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     h     ...   velikost kroku pro vycisleni num. derivaci
!     eps   ...   pozadovana presnost, ||xEnd - xEnd-1|| <= eps
!     maxiter..   maximalni pocet iteraci metody, navratova hodnota je
!                 skutecny pocet iteraci newtona
!
!     Pozn: Oproti predchozimu korektoru zalozenem na algoritmu DERPAR
!           (kubicek) jiz nevynechavam zadnou promennou, ale naopak,
!           pridavam rovnici (parametrizace podle delky oblouku)
!     Pozn: Prace s promennymi je sice psana obecne pro R(n+m)->Rn, ale
!           kolkolem pridavam pouze jednu rovnici a tedy potrebuji
!           kontinuovat pouze podle jednoho parametru R(n+1)->Rn
!
!     Autor:      Martin Isoz
!     Datum:      21.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
      external subrF
!     vstupni promenne
      integer(kind=4),intent(in)    :: n(2),m                           !rozmer problemu, ind. vyn. prom., pocet pars
      double precision,intent(inout):: x0(n(2))                         !vektor s poc. nastrelem a vyslednym resenim
      double precision,intent(in)   :: tVec0(n(2)-n(1),n(2))            !tecny vektor/tecna nadrovina v xO
      double precision,intent(in)   :: h,eps,pars(m)                    !velikost kroku, presnost, vektor parametru do subrF
      integer(kind=4),intent(inout) :: maxiter                          !max. pocet iteraci
!     dalsi promenne pro NEWTONA
      double precision  :: xOld(n(2)),xNew(n(2))                        !pomocne pracovni promenne, DeltaX, vycisleni funkcnich hodnot
      double precision  :: tVec(n(2),n(2)-n(1))                         !aktualizovany tecny vektor
      double precision  :: Dx(n(2)),FVec(n(1))                          !pomocne pracovni promenne, DeltaX, vycisleni funkcnich hodnot
      double precision  :: Dv(n(2),1)                                   !pomocna pracovni promenna, DeltaV
      double precision  :: delta                                        !sledovani prubezne presnosti
      integer(kind=4)   :: iter                                         !pocitadlo iteraci
!      integer(kind=4)   :: i                                            !indexovadlo pro vypisy
      double precision  :: jac(n(1),n(2))                               !matice pro vypocet jacobiho matice
      double precision  :: RHSX(n(2)),dRHSX(n(2),n(2))                  !prave strany metody a jejich diferencial
      double precision  :: RHSV(n(2),1)                                 !prave strany metody pro tecny vektor
!     tlumici faktor Newtona
!      double precision  :: FN(n(1)),FO(n(1))                           !vycisleni fn hodnot pro xNew a xOld
      double precision  :: TLkoef                                       !tlumici faktor
!     DSGESV
      double precision  :: swork(n(2),1)
      double precision  :: work(n(2),n(2))
      integer(kind=4)   :: pivot(n(2),n(2)),iter2,info
!     jacU
      integer(kind=4)   :: jSW(0)                                       !zero sized integer pro preskakovane prom. pri volani jacU
!
!
!=====EXECUTABLE=PART===================================================
!
!     naplneni pomocnych promennych pro newtona
      iter        = 0
      delta       = eps + 1.0d0
      xOld        = x0
      TLkoef      = 1.0d0
      tVec        = transpose(tVec0)
!     
!     vypis nastrelu
!     write(unit=*,fmt=*) 'NASTREL====================================='
!     write(unit=*,fmt=*) 'k', iter, 'x(k)', xOld, 'delta', delta& 
!     &     , 'dumpF', TLkoef
!     write(unit=*,fmt=*) 'ITERACE====================================='
!     vypocetni cyklus newtonovy metody
      do while( eps .lt. delta .and. iter .le. maxiter )
!
            call jacU( subrF,xOld,pars,n,m,h,jac,jSW )                  !spocitam diferencial puvodniho zobrazeni
!
            call subrF(xOld,FVec,pars)                                  !vypocet pravych stran puvodniho zobrazeni v danem bode
!
            RHSX         = 0.0d0                                        !restartovat right hand side vektor
            RHSX(1:n(1)) = FVec                                         !naplnit prvnich n(1) pozic rezidui0
!
            RHSV          = 0.0d0                                       !pripravim si vektor pravych stran pro novy tecny vektor
!            write (unit=*,fmt=*) 'jac'
!            do i=1,size(jac,1)
!                  write(*,'(20G12.4)') jac(i,:)
!            end do
!            write (unit=*,fmt=*) 'tVec'
!            do i=1,size(tVec,1)
!                  write(*,'(20G12.4)') tVec(i,:)
!            end do
!            do i=1,n(1),1                                               !n(2)-n(1) == 1
!                  RHSV(i,n(2)-n(1)) = sum(jac(i,:)*tVec(:,n(2)-n(1)))
!            end do
!            write (unit=*,fmt=*) 'RHSV'
!            do i=1,size(RHSV,1)
!                  write(*,'(20G12.4)') RHSV(i,:)
!            end do
!
            dRHSX                = 0.0d0                                !restartovat diferencial RHS vektoru
            dRHSX(1:n(1),:)      = jac                                  !prvnich n(1) radku je diferencial puvodniho zobrazeni
            dRHSX(n(1)+1:n(2),:) = transpose(tVec)                      !n(1)+1 - n(2) radek je tecna nadrovina
!
!     kontrolni vypis 1
!            write (unit=*,fmt=*) 'dRHS'
!            do i=1,size(dRHS,1)
!                  write(*,'(20G12.4)') dRHS(i,:)
!            end do
!            write (unit=*,fmt=*), '-RHS'
!            write(*,'(20G12.4)') -RHS
!
            call DSGESV(n(2),n(2)-n(1),dRHSX,n(2),pivot,-RHSX,n(2),Dx,&
     & n(2),swork,work,iter2,info)                                      !reseni soustavy lin. rovnic - nalezeni Dx
            call DSGESV(n(2),n(2)-n(1),dRHSX,n(2),pivot,-RHSV,n(2),Dv,&
     & n(2),swork,work,iter2,info)                                      !reseni soustavy lin. rovnic - nalezeni Dv
!     kontrolni vypis 2
!            write (unit=*,fmt=*), 'Dx'
!            write(*,'(20G12.4)') Dx
!
            xNew = xOld + TLkoef*Dx                                     !vypocet nove iterace, xVec
            tVec = (tVec + TLkoef*Dv)/norm2(tVec)                       !vypocet nove iterace, tVec + normovani
!
!!     uprava tlumiciho faktoru
!            call subrF(xNew,FN,pars)                                    !vycislim residualy funkci v xNew a xOld
!            call subrF(xOld,FO,pars)
!
!           if( sqrt(sum(FN**2.0d0)) .lt. sqrt(sum(FO**2.0d0)) ) then   !pri uspechu preulozim xNew, jinak zpulim tlumici faktor
!                 delta = dsqrt(sum(xNew**2.0d0 - xOld**2.0d0))         !vypocet rozdilu xNew a xOld (norma rozdilu)
!                 xOld = xNew
!           else
!                 TLkoef = TLkoef/2.0d0
!                 delta = delta
!                 xOld = xOld
!           end if
!
            delta = dsqrt(sum((xNew - xOld)**2.0d0))
            xOld = xNew
!
            iter = iter + 1                                             !zvyseni pocitadla iteraci
!
!     vypis iterace
!           write(unit=*,fmt='(20G12.4)') 'k', iter, 'x(k)', xNew , 'delta',delta&
!     &     , 'dumpF', TLkoef
!
      end do
      x0 = xOld                                                         !ulozim vysledek
      maxiter = iter
      return
      end subroutine korektorPM
!
!=======================================================================
!
      subroutine solver (subrF,prediktor,korektor,x0,ip,BMin,BMax,BTrg,&
     & pars,n,m,g,h,maxstep,nprstep,eps,maxiter,wrtunt,retcode)
!
!=====HLAVICKA==========================================================
!
!     Solver pro kontinuacni ulohy, zatim pouze pro ulohy s kontinuaci
!     jednoho parametru. Solver vola stridave prediktor a korektor, az
!     do chvile, nez dosahne zastavovaciho kriteria nebo selze metoda
!
!     Vysledky jsou zapisovany po krocich do wrtunt, je-li toto soubor,
!     musi byt pri volani resitele otevreny a pripraveny k zapisovani
!
!     
!
!     Program vola subroutiny:
!     SUBRF ...   subroutina vracejici vektor residui /neboli vektorove
!                 pole z Rn, vracejici vektor pravych stran/
!                 subrF musi mit specialni strukturu argumentu,
!                 call subrF(xVec,RSVec,pars), kde xVec jsou vsechny
!                 stavove promenne (vcetne kontinuovanych parametru),
!                 RSVec je vektor pravych stran a pars jsou vsechny
!                 fixovane parametry
!     PREDIKTOR.  subroutina poskytujici nastrel reseni na zaklade
!                 predchoziho bodu. podrobneji v hlavicce prediktoru
!     KOREKTOR..  subroutina upresnujici nastrel vytvoreny prediktorem,
!                 presneji v hlavicce korektoru
!     DETSTAB...  subroutina na urceni stability RS pro vynechany pocet
!                 parametru
!
!     Pozn: Idea je takova, ze v dalsim bude mozne solver volat pro
!           ruzne korektory (a mozna i prediktory), zatim je napsan
!           korektor podle pseudodelky oblouku (AL), ale v dohledne dobe
!           by mel byt dopsan i penrose-moore (PM) korektor
!
!     Vstupni argumenty programu:
!     subrF ...   subroutina vracejici vektor residui (rovnice resene
!                 soustavy)
!                 vyzaduje subroutinu volatelnou func(xVec,FVec)
!     prediktor   subroutina vracejici nastrel dalsiho bodu na krivce
!                 reseni
!     korektor.   subroutina zpresnujici nastrel prediktora
!     x0    ...   pocatecni bod reseni
!     ip    ...   indexy (pozice) parametru v x0
!     BMin  ...   minimalni hodnoty stavovych promennych
!     BMax  ...   maximalni hodnoty stavovych promennych
!     BTrg  ...   cilove hodnoty stavovych promennych
!     pars  ...   vektor parametru pro subroutinu subrF
!     n     ...   vektor s rozmery problemu, (/integer,integer/),
!                 n1 - pocet rovnic
!                 n2 - pocet stavovych promennych (vcetne parametru)
!     m     ...   pocet parametru subroutiny subrF
!     g     ...   delka kroku pro numerickou derivaci
!     h     ...   velikost kroku pro kontinuaci, vektor (/h0,hMin,hMax/)
!                 h0 - pocatecni delka kroku
!                 hMin-minimalni delka kroku, h < hMin -> alg. failed
!                 hMax-maximalni delka kroku, adaptaci jiz nenatahovan
!     maxstep..   maximalni pocet kroku algoritmu
!     nprstep..   kolika-krokovy je prediktor - kolik radku bude mit
!                 matice tMat (pocet ukladanych tecnych vektoru)
!     eps   ...   pozadovana presnost, ||xEnd - xEnd-1|| <= eps
!     maxiter..   maximalni pocet iteraci korektora
!     wrtunt...   jednotka pro zapisovani vysledku, je-li to soubor,
!                 musi byt pripraven pro zapisovani
!     retcode..   navratova hodnota algoritmu
!                 0 - dosazeni nektere cilove hodnoty stavove promenne,
!                     x(i) = BTrg(i)
!                 1 - prekroceni maximalni (minimalni) hodnoty stavove
!                     promenne, x(i) >(<) BMax(i)(BMin(i))
!                 2 - prilis maly krok, h < hMin
!                 3 - korektor nezkonvergoval, niter > maxiter
!                 4 - preroceni maximalniho povoleneho poctu kroku
!
!     Pozn: Prace s promennymi je sice psana obecne pro R(n+m)->Rn, ale
!           kolkolem pridavam pouze jednu rovnici a tedy potrebuji
!           kontinuovat pouze podle jednoho parametru R(n+1)->Rn
!
!     Autor:      Martin Isoz
!     Datum:      19.10.2013
!     Licence:    FreeBSD (viz. README.txt)
!
!
!=======================================================================
!
      implicit none
!
!=====DEKLARACE=PROMENNYCH==============================================
!
!-----argumenty---------------------------------------------------------
!     externi subroutiny
      external subrF,prediktor,korektor
!
!     vstupni promenne
      integer(kind=4),intent(in)  :: n(2),m                             !n(1) rovnic, n(2) stav. prom, m pars
      integer(kind=4),intent(in)  :: ip(n(2)-n(1))                      !indexy parametru, delka je rozdil pocti rcic a prom
      integer(kind=4),intent(in)  :: nprstep                            !pocet ukladanych tecnych vektoru
      integer(kind=4),intent(in)  :: maxiter                            !maximalni pocet iteraci korektora
      integer(kind=4),intent(in)  :: wrtunt                             !jednotky pro zapisovani vysledku - stable a unstable EP
      double precision,intent(in) :: x0(n(2))                           !pocatecni nastrel reseni
      double precision,intent(in) :: BMin(n(2)),BMax(n(2))              !hranice reseni
      double precision,intent(in) :: pars(m)                            !parametry pro SUBRF
      double precision,intent(in) :: g,h(3)                             !kroky metod
      double precision,intent(in) :: eps                                !pozadovana presnost reseni
!     vstupne-vystupni promenne
      integer(kind=4),intent(inout) :: maxstep                          !maximalni pocet kroku/vysledny pocet kroku
      double precision,intent(inout):: BTrg(n(2))                       !cil pro stavove promenne/dosazene s.prom.

!     vystupni promenne
      integer(kind=4),intent(out) :: retcode                            !navratova hodnota
!
!-----lokalni promenne--------------------------------------------------
!     PREDIKTOR
      integer(kind=4)         :: nSol                                   !pocet predresenych tVec je 0
      double precision        :: tMat(nprstep,n(2))                     !matice tecnych vektoru
      double precision        :: hAct                                   !aktualni krok
!     KOREKTOR
      double precision        :: tVec(n(2))                             !tecny vektor (muze byt upravovan)
      integer(kind=4)         :: niter                                  !pocitadlo iteraci korektora
!     DETSTAB
      logical                 :: stab=.true.                                   !promenna pro rozhodnuti o stabilite RS
!     vseobecne
      integer(kind=4)         :: miter,i                                !pocitadla cyklu
      double precision        :: xVec(n(2))                             !aktualni reseni
!
!
!=====EXECUTABLE=PART===================================================
!     
!-----inicializace-promennych-------------------------------------------
      tMat  = 0
      tVec  = dsqrt(dble(n(2)))/dble(n(2))                              !normovany jednotkovy vektor?
      xVec  = x0
      hAct  = max(h(1)/1.0d2,1.0d-12)                                   !v pripravny cyklus->maly krok/double presnost
      miter = 0
      retcode=-1                                                        !zastavovaci hodnota
      nSol  = 1                                                         !pocet pristupnuch reseni v tMat
      niter = maxiter                                                   !restartuji pocitadlo iteraci newtona
!
!-----inicializace-tecneho-vektor-tVec----------------------------------
!
      call initRES (subrF,xVec,pars,tVec,n,m,g,eps,niter)               !inicializuji reseni - zpresneni x0 a prvni tecny vektor
      if (niter .eq. maxiter) then                                      !flag check
            retcode = 3                                                 !napln return code
            BTrg    = xVec                                              !posbirej vystupy
            maxstep = miter
            go to 666                                                   !vrat reseni
      end if
      niter = maxiter                                                   !restartuji pocitadlo iteraci newtona
      go to 666
!-----kontrola-a-zpresneni-nastrelu-------------------------------------
      call korektor(subrF,xVec,tVec,pars,n,m,g,eps,niter)
      if (niter .eq. maxiter) then                                      !flag check
            retcode = 3                                                 !napln return code
            BTrg    = xVec                                              !posbirej vystupy
            maxstep = miter
            go to 666                                                   !vrat reseni
      end if
      niter = maxiter                                                   !restartuji pocitadlo iteraci newtona
      tMat(1,:) = tVec                                                  !ulozim si zpresneny tVec
!     zjisteni stability reseni
      stab = detStab(subrF,xVec,pars,n,m,g,ip)
!     vypis vysledku zpresneni
      call wrt2unit(wrtunt,miter,xVec,hAct,niter,stab)                  !vypis vysledku
!-----vypocetni-cyklus---pripravny-(naplnim tMat)-----------------------
      do miter=1,nprstep-1,1                                            !zastaveni je reseno pres if/go to
!     volani prediktora
            call prediktor(subrF,xVec,tMat,pars,n,m,nSol,g,hAct)
            nSol  = nSol + 1
!     preulozeni tecneho vektoru
            tVec  = tMat(1,:)                                           !ulozim vysledny tecny vektor
!     volani korektora
            call korektor(subrF,xVec,tVec,pars,n,m,g,eps,niter)
            if (niter .eq. maxiter) then                                !flag check
                  retcode = 3                                           !napln return code
                  BTrg    = xVec                                        !posbirej vystupy
                  maxstep = miter
                  go to 666                                             !vrat reseni
            end if
!     uprava delky kroku - v pripravnem cyklu krok neupravuji
!     zjisteni stability reseni
            stab = detStab(subrF,xVec,pars,n,m,g,ip)
!     vypis vysledku kroku do wrtunt
            call wrt2unit(wrtunt,miter,xVec,hAct,niter,stab)            !vypis vysledku
!     kontrola vysledku a dalsi pomocne operace
            do i = 1,n(2),1
                  if (xVec(i) .lt. BMin(i) .or. xVec(i) .gt. BMax(i)) then
                        retcode = 1
                        BTrg    = xVec
                        maxstep = miter
                        go to 666
                  end if
                  if (dabs(xVec(i) - BTrg(i)) .le. hAct/1.0d1 ) then    !tohle je krajne nepekne
                        retcode = 0
                        BTrg    = xVec
                        maxstep = miter
                        go to 666
                  end if
            end do
            niter = maxiter                                             !restartuji pocitadlo iteraci newtona
      end do
!-----vypocetni-cyklus---hlavni-s-plnou-tMat----------------------------
      hAct = h(1)                                                       !restartuji velikost kroku
      do miter=nprstep,maxstep-1,1                                      !zastaveni je reseno pres if/go to
!     volani prediktora
            call prediktor(subrF,xVec,tMat,pars,n,m,nSol,g,hAct)
!     preulozeni tecneho vektoru
            tVec  = tMat(1,:)                                           !ulozim vysledny tecny vektor
!     volani korektora
            call korektor(subrF,xVec,tVec,pars,n,m,g,eps,niter)
!     uprava delky kroku
            if (niter .eq. maxiter) then                                !flag check
                  retcode = 3                                           !nezkonvergoval newton
                  BTrg    = xVec
                  maxstep = miter
                  go to 666
            end if
            if(dabs(hAct) .le. dabs(h(3))) then                         !hlidam maximalni delku kroku
                  if    (niter .lt. 2) then
                        hAct = hAct*2.0d0                               !pro maly pocet iteraci natahnu krok
                  elseif(niter .eq. 2) then
                        hAct = 1.5d0*hAct                               !mensi natahnuti kroku
                  elseif(niter .gt. 3) then
                        hAct = hAct/2.0d0                               !pro velky pocet iteraci pulim krok
                  elseif(niter .gt. 4) then
                        hAct = hAct/4.0d0
                  end if
            end if
            if(dabs(hAct) .le. dabs(h(2))) then                         !zkontrolovat delku kroku
                  retcode = 2                                           !prilis maly krok
                  BTrg    = xVec
                  maxstep = miter
                  go to 666
            end if
!     zjisteni stability reseni
            stab = detStab(subrF,xVec,pars,n,m,g,ip)
!     vypis vysledku kroku do wrtunt
            call wrt2unit(wrtunt,miter,xVec,hAct,niter,stab)            !vypis vysledku
!     kontrola vysledku a dalsi pomocne operace
            do i = 1,n(2),1
                  if (xVec(i) .lt. BMin(i).or.xVec(i) .gt. BMax(i)) then
                        retcode = 1                                     !vyboceni z regionu zajmu
                        BTrg    = xVec
                        maxstep = miter
                        go to 666
                  end if
                  if (dabs(xVec(i) - BTrg(i)) .le. hAct ) then          !tohle je krajne nepekne
                        retcode = 0                                     !dosazeni cilove hodnoty
                        BTrg    = xVec
                        maxstep = miter
                        go to 666
                  end if
            end do
            niter = maxiter                                             !restartuji pocitadlo iteraci pro newtona
      end do
!-----konecne-kontroly-a-upravy-----------------------------------------
      if (miter .ge. maxstep ) then
            retcode = 4                                                 !doslo k prekroceni poctu kroku
            BTrg    = xVec
            maxstep = miter
            go to 666
      end if
!     nic mezi touto radkou a return by nemelo byt vykonano
      write (unit=*,fmt=*) 'Unexpected program ending'                  !vypis na obrazovku, neco je spatne
666   return
      contains 
            subroutine wrt2unit(wrtunt,miter,xVec,hAct,niter,stab)
!           function for writing results to some unit
!           data declaration
            integer(kind=4),intent(in) ::wrtunt,miter,niter
            double precision,intent(in)::xVec(:),hAct
            logical,intent(in)         :: stab
!           writing itself
            write(unit=wrtunt,fmt='(20G12.4)') &
     &      'k', miter, &                                               !aktualni poradi kroku
     &      'x(k)', xVec , &                                            !vektor reseni
     &      'hAct', hAct, &                                             !pristi pouzity krok
     &      'niter', niter, &                                            !pocet iteraci korektora
     &      'stab', stab                                                !stabilita RS
            end subroutine wrt2unit
      end subroutine solver
!
!=======================================================================
!
      end module FORTCONT
!
!=======================================================================
!
