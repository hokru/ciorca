!********************************************
! Command Line Input Generator for ORCA 
! author: Holger Kruse (holger.kruse@uni-muenster.de)
! idea and concept (and some code..) based on 'cefine' for Turbomole by Stefan Grimme 
! Last change: 
! To Do:
!********************************************
Program comand_line_define_ORCA
implicit none
      
integer maxarg,io,i,nn,j,ii
integer charge,scfcv
integer nat,gcart,mode
integer iat(10000),giter,siter,grid
! parameter (maxarg=40)
character*80, allocatable :: arg(:)    
character*80 outfile
character*80 func ,fcore
character*80 bas,home,optim,coordsys,optstep
character*80 pwd            
character*80 atmp       
character*80 infile
character*80 solvent,citype,KCOPT,xyzfile
character*10 version
character*250 add
logical KEEP,HYBRID,GGA,SETSYM,ZORA,CCT
logical da,TS,COSMO,ECHO,RANGST,JOB,FREQ,RESTART,TE
logical VDW,PR,OPT,FXYZ,FTMOL,SMOOTH
logical QUICK,SMEAR,HOPT,NOFC,jbas,rijk,cbas
logical ANC,CNEW,SCS,MP2,T,MDCI,lgrid,Lbas,Lfunc,lsolvent,lfcore,ladd


real*8  xx(5),trustrad
real*8 xyz(3,10000)
character*80, allocatable :: TOT(:)

      
version='0.2b'



io=1

GGA=.true.   ! assuming always GGA
HYBRID=.false.
ECHO =.false.

!cccccccccccccccccccccccccccccccccc
!ccccccccc   defaults ccccccccccccc
!cccccccccccccccccccccccccccccccccc

! cccccccccc   logical for the default:
lgrid=.false.
Lfunc =.false.
Lbas    =.false.
lsolvent=.false.
lfcore=.false.
lsolvent=.false.
lfcore=.false.
ladd=.false.

VDW  =.false.
FTMOL=.false.
JOB=.false.
FREQ=.false.
TS=.false.
ZORA=.false.
RESTART=.true.
OPT=.false.
QUICK=.false.
SMOOTH=.false.
TE=.false.
SMEAR=.true.
ANC=.false.
 CNEW=.true.
 HOPT=.false.
setsym=.false.
 
 
 cbas=.false.
 jbas=.true.
 rijk=.false.


func='revPBE'
bas='def2-TZVP'
xyzfile='start.xyz'
outfile='orca.in'
infile='XYZ.in'

solvent='water'
charge =0
grid = 4
gcart=4
scfcv=7
 COSMO=.false. 


trustrad=0.3
optim='Delocal'

giter=250
siter=125

 coordsys='redundant'
 OptStep='rfo'
 NOFC=.false.
!ccccccccccccccccccccccccc
       write(*,*)'  CIORCA  V ', trim(version)

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccc    read and evaluate /$HOME/.ciorcarc  file c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
      call get_environment_variable('PWD',pwd)
      call get_environment_variable('HOME',atmp)
      home=trim(atmp)//'/.ciorcarc'
!      write(*,*) home
      inquire(file=home,exist=da)
!      inquire(file='~/.ciadfrc',exist=da) ! doesnt work with gfortran
      if(da)then
         open(unit=20,file=home)
!         open(unit=20,file='~/.ciadfrc')   ! doesnt work with gfortran
 842     read(20,'(a)',end=942)atmp
         if(index(atmp,'func').ne.0)then         
            call backstring(atmp,func,4)
         endif
         if(index(atmp,'basis').ne.0)then         
            call backstring(atmp,bas,5)
         endif
!          if(index(atmp,'fc ').ne.0)then         
!             call backstring(atmp,fcore,2)
!          endif
!          if(index(atmp,'smear').ne.0)then         
!             smear=.true.
!          endif
         if(index(atmp,'vdw').ne.0) then
            if(index(atmp,'on').ne.0)VDW=.true.   
         endif
         if(index(atmp,'trustrad').ne.0)then         
            call readl(atmp,xx,nn)
            trustrad=xx(nn)
         endif
!          if(index(atmp,'grid').ne.0)then         
!             call readl(atmp,xx,nn)
!             grid=xx(nn)
!          endif
         if(index(atmp,'maxcycle').ne.0)then         
            call readl(atmp,xx,nn)
            giter=xx(nn)
         endif
!          if(index(atmp,'scfiter').ne.0)then         
!             call readl(atmp,xx,nn)
!             siter=xx(nn)
!          endif
         if(index(atmp,'format').ne.0)then         
           if(index(atmp,'xyz').ne.0)then         
            FXYZ=.true.
           elseif(index(atmp,'tmol').ne.0)then         
            FTMOL=.true.
           endif
         endif
         if(index(atmp,'cosmo').ne.0)then         
            call backstring(atmp,solvent,5)
         endif
!          if(index(atmp,'restart').ne.0)then         
!             RESTART=.true.
!            if(index(atmp,'off').ne.0)then
!            RESTART=.false.
!            endif         
!          endif
         if(index(atmp,'ifile').ne.0)then         
            call backstring(atmp,infile,5)
         endif
         if(index(atmp,'ofile').ne.0)then         
            call backstring(atmp,outfile,5)
         endif
         goto 842
 942     close(20)
      endif

!cccccccccccccccccccccccccccccccccccccccccccc
!cccccccccc   get command line arguments   cc
!cccccccccccccccccccccccccccccccccccccccccccc
!c Test if an explicit option was set. Otherwise use defaults AFTER the argument evaluation
      arg=''
maxarg=iargc()
if(maxarg.gt.0) then
!maxarg=iargc()
      write(*,*) 'arguments(debug) :',maxarg
allocate(arg(maxarg),TOT(maxarg))
      do i=1,maxarg
         call getarg(i,arg(i))
      enddo
!ccc first print help
      if(arg(1).eq.'-h' .or. arg(1).eq.'?' .or. arg(1).eq.'-help')then
      pr=.true. 
      if(pr) then
        write(*,*) '                                       '
        write(*,'(10x,''*** Command line Input generator for ORCA ***'')') 
        write(*,'(10x,''          H. Kruse Dez. 2010      '')') 
        write(*,'(10x,''            | V '',a,'' |       '')')trim(version)
        write(*,'(10x,''          developer version     '')') 
      endif
          write(*,*)'OPTIONS:'
!          write(*,*)'   -v (print defaults, after reading .ciadfrc)'
          write(*,*)'   -h / -help / ? (this help output)'
          write(*,*)'   -H <int> <string1> ... <stringN>  (adds lines for simple input)'
!          write(*,*)' more information in the manual'
!          write(*,*)'   '
!          write(*,*)'* input/ouput:'
         write(*,*)'   -n (re-name output file )'
         write(*,*)'   -f   <string>  input file (def: XYZ.in)'
         write(*,*)'   -tm (TMOL input format)'
         write(*,*)'   -xyz (XYZ (xmol) input format (def))'
!          write(*,*)'  [ -job (sets JOB option=true)                ]'
!          write(*,*)'  [ -o    <string> actual JOB file)            ]' 
!          write(*,*)'  [ -wall <int>  (walltime in h for "jobfile") ]'
         write(*,*)' needs: <input> file in XYZ or TMOL format'
         write(*,*)' '
         write(*,*)'* calculation setup:'
         write(*,*)'   -func   <string>'
         write(*,*)'   -bas    <string>'
!          write(*,*)'   -fc     <string> (frozen core, def: none)'
!          write(*,*)'   -grid   <integer>'
         write(*,*)'   -vdw06 // -vdw10 (DFT-D) // -novdw'
!          write(*,*)'   -noopt / -nopt (do single point)'
         write(*,*)'   -opt  (do optimization)'
         write(*,*)'   -chrg   <integer>'
!          write(*,*)'   -scfcv  <integer> (scf convergence, def: 7)'
         write(*,*)'   -siter  <integer> (max scf iterations)'
!          write(*,*)'   -grad   <integer> (gradient convergence,def: 4)'
!          write(*,*)'   -c      <integer> (max cycle in GeoOpt)'
         write(*,*)'   -cart (cartesian coordinates, def: delocal)'
!          write(*,*)'   -intal (internal coordinates)'
!          write(*,*)'   -sym    <string> (symmetry, eg: "C\(4h\)" )'
         write(*,*)'   -cosmo  (COSMO with solvent=default=water)' 
         write(*,*)'   -cosmoS  <string> (COSMO with solvent=<string>)' 
         write(*,*)'   -freq (numerical freq)'
         write(*,*)'   -zora (turns scalar zora on)'
!          write(*,*)'   -ts <integer> (TS search, <int>=modes)'
!          write(*,*)'   '
!          write(*,*)'* special calls:'
!          write(*,*)'   -frag (handle self defined fragment files)'
!          write(*,*)'   -quick (pbe-d/dz grid=3,cart,Geocycle=30)'
!          write(*,*)'   -dftb (dftb input, quick pre-opt,-noopt for SP)'
!          write(*,*)'   -quild (quild input with deloc coordinates)'
!          write(*,*)'   -ancopt (prepares ancopt input)'
!          write(*,*)'   '
!          write(*,*)'   '
         stop
      endif
! now process the arguments. (A 'space' after the option is a good idea)
add=''
TOT=''
ii=0

      do i=1,maxarg
!       write(*,*) arg
!        TOT(i)=trim(add)
         if(arg(i).ne.'')then
            if(index(arg(i),'-H').ne.0)then !   '-add <int> <arg_1 ...arg_ii>
             ladd=.true.
!               call readl(arg(i+1),xx,nn)
             do j=i+1,maxarg
              if(index(arg(j),' -').ne.0)then  ! next arguments, goto end of loop
               goto 123
              endif
             TOT(j)=trim(arg(j))
!              write(*,*) TOT(j),i,j
             ii=j
             enddo
             exit
            endif
            if(index(arg(i),'-vdw').ne.0)then
            VDW=.false. ! DFT-D on VDW10 default
            TOT(i)='VDW10'
            endif
            if(index(arg(i),'-vdw10').ne.0)then
            VDW=.false. ! DFT-D on VDW10 default
            TOT(i)='VDW10'
            endif
            if(index(arg(i),'-vdw06').ne.0)then
            VDW=.false. ! DFT-D on VDW10 default
            TOT(i)='VDW06'
            endif
!             if(index(arg(i),'-ancopt ').ne.0) then  ! ancopt external optimizer
!             ANC=.true.  
!             endif
!             if(index(arg(i),'-v ').ne.0)   ECHO=.true. ! print some defaults
            if(index(arg(i),'-quick ').ne.0) quick=.true.
             if(index(arg(i),'-n ').ne.0) then       ! name of the output file
             outfile=arg(i+1) 
             endif
            if(index(arg(i),'-freq').ne.0) TOT(i)='NUMFREQ' ! numerical freq
!             if(index(arg(i),'-noopt').ne.0)   OPT=.false. ! SP
!             if(index(arg(i),'-nopt').ne.0)   OPT=.false. ! SP
            if(index(arg(i),'-opt').ne.0) then
            OPT=.true. ! GEO OPT
            TOT(i)='TightOpt'
            endif
!             if(index(arg(i),'-novdw').ne.0) VDW=.false.  ! DFT-D off
            if(index(arg(i),'-zora').ne.0) then
            ZORA=.true.  ! scalar zora
            TOT(i)='ZORA'
            endif
            if(index(arg(i),'-xyz ').ne.0) then ! xzy files as input
             FXYZ=.true.
             FTMOL=.false.
            endif
!             if(index(arg(i),'-norestart ').ne.0)then  !infile
!                restart=.false.           
!             add='NOAUTOSTART'
!             endif
             if(index(arg(i),'-tm ').ne.0) then ! turbomole format as input
             FTMOL=.true.
             FXYZ=.false.
             endif
             if(index(arg(i),'-f ').ne.0)then  !infile
                infile=arg(i+1)           
             endif
             if(index(arg(i),'-o ').ne.0)then ! OUTPUT FILE
                outfile=arg(i+1)           
             endif
             if(index(arg(i),'-chrg').ne.0)then !Charge
                call readl(arg(i+1),xx,nn)
                charge=idint(xx(1))
            endif
            if(index(arg(i),'-cosmoS').ne.0)then ! COSMO
            lsolvent=.true.
               COSMO=.true.
               solvent=arg(i+1)
            write(TOT(i),'(''COSMO('',a,'')'')') trim(solvent)
            endif
            if(index(arg(i),'-cosmo').ne.0)then ! COSMO
!             lsolvent=.true.
               COSMO=.true.
!                solvent=arg(i+1)
!             write(TOT(i),'(''COSMO('',a,'')'')') trim(solvent)
            endif
            if(index(arg(i),'-ts ').ne.0)then ! TS search
               call readl(arg(i+1),xx,nn)
               TS=.true.
               mode=idint(xx(1))
               TOT(i)='OptTS'
            endif
!             if(index(arg(i),'-scfcv').ne.0)then ! SCFCONV
!                call readl(arg(i+1),xx,nn)
!                scfcv=idint(xx(1))
!             endif
!             if(index(arg(i),'-grad').ne.0)then ! GRADCONV
!                call readl(arg(i+1),xx,nn)
!                gcart=idint(xx(1))
!             endif
            if(index(arg(i),'-siter ').ne.0)then ! SCFCONV
               call readl(arg(i+1),xx,nn)
               siter=idint(xx(1))
            endif
            if(index(arg(i),'-c ').ne.0)then ! geometry cycles
               call readl(arg(i+1),xx,nn)
               giter=idint(xx(1))
            endif
            if(index(arg(i),'-sym').ne.0)then  ! SET SYMMETRY
            TOT(i)='UseSym'
            endif
            if(index(arg(i),'-grid').ne.0) then ! NUM INT Grid
            lgrid=.true.
             call readl(arg(i+1),xx,nn)
             grid =xx(1)
!              TOT(i)='NOFINALGRID'
             write(TOT(i),'(''GRID'',I1,'' NOFINALGRID'')') grid
!             write(*,*) grid, xx(1)
            endif
            if(index(arg(i),'-bas').ne.0) then
!             Lbas=.true.
             bas  =arg(i+1)  ! basis set
!             write(TOT(i),'(a)') bas
            endif
            if(index(arg(i),'-cbas').ne.0) cbas=.true.
            if(index(arg(i),'-jbas').ne.0) jbas=.true.
            if(index(arg(i),'-rijk').ne.0) rijk=.true.
!             if(index(arg(i),'-cart').ne.0) then ! cartesian optimizer
!              optstep  ='Cartesian'
!             endif
!             if(index(arg(i),'-deloc').ne.0) then ! internal/ Zmatrix optimizer
!              optstep  ='deloc'
!             endif
!             if(index(arg(i),'-nofc').ne.0) then ! frozen core
!              NOFC=.true.
!              add='NOFROZENCORE'
!             endif
            if(index(arg(i),'-func').ne.0) then    ! functional
            Lfunc=.true.
               func=arg(i+1)           
              write(TOT(i),'(a)') func
            endif
            
!!!!!          MDCI         !!!!!!

            if(index(arg(i),'-ccsdt').ne.0) CCT=.true.
!             if(index(arg(i),'-angst').ne.0)RANGST=.true. ! change cartesian coords unit to
! angstrom
! !c keep outputs for debuging purposes  (no purpose yet)
!             if(index(arg(i),'-keep').ne.0)KEEP=.true. 
         endif
123 continue
      enddo
endif
!ccccccccccccccccccccccccc 


!cc hidden file
      da=.false.
!      inquire(file='.UHF',exist=da)
!      if(da)then
!        open(unit=21,file='.UHF')
!        read(21,*) nopen
!        write(*,'(''!!! UHF enforced by user in <.UHF> . Na-Nb:'',i4,'' !!!'')')nopen
!      endif
      inquire(file='.CHRG',exist=da)
      if(da)then
        open(unit=21,file='.CHRG')
        read(21,*) charge
        write(*,'(''!!! charge in <.CHRG> :'',i4,'' !!!'')')charge
      endif




!ccccccccccccccccccccccccccc
!cccc read coordinates     c
!ccccccccccccccccccccccccccc
!  more formats?
!       if(.not.Frag)
       if (FXYZ)call xyzrd(xyz,iat,nat,infile)
       if (FTMOL)call tmolrd(xyz,iat,nat,infile)
       open(unit=98,file=xyzfile)
       call wxyz(xyz,iat,nat,98)
       close(98)
!ccccccccccccccccccccccccccccccccccc 
!ccc prepare coord file for ANCOPT c
!ccccccccccccccccccccccccccccccccccc
        if(ANC) then
        open(unit=99,file='coord')
        call wtm(xyz,iat,nat,99)
        close(99)
        write(*,*)' ancopt needs: orca.in '
!c        call system('touch control')
!c        inquire(file='control',exist=da)
!c        if(.not.da) call system('echo ''$symmetry xx'' > control ')
        endif

!ccccccccccccccccccccccccccccccc
!ccccc   ECHO DEFAULTS        cc
!ccccccccccccccccccccccccccccccc
      if(ECHO) then
      write(*,*) ' '
      write(*,*) 'some default settings (after reading .ciorcarc) '
      write(*,*)'functional  ',trim(func), ', basis set ', trim(bas) 
      write(*,'('' num grid '',f4.1,'' grad conv'',i2)')grid,gcart
      write(*,*)'frozen core ',trim(fcore), ', vdw ', VDW
      write(*,*)'COSMO ',COSMO,', solvent ', trim(solvent)
      write(*,*)'inputfile ', trim(infile),', outputfile ',trim(outfile)
      write(*,'('' scfconv '',i2,'', optim '',a)')scfcv,trim(optim)
      stop 
      endif

!ccccccccccccccccccccccccccccccccccccc     
!cccccccccc   QUICK Option           c
!ccccccccccccccccccccccccccccccccccccc
       if(QUICK) then
      OPT=.true.
       grid=4
       func='revPBE'
       bas='SV(P)'
       scfcv=6
       gcart=3
       vdw=.true.
       optim='COPT'  ! cartesian for safety
       TOT(1)='revPBE SVP SVP/C LooseOpt COPT GRID3 NOFINALGRID'
       TOT(2:maxarg)=''
       endif
     


!ccccccccccccccccccccccccccccccccccccc
!ccccccc Print what you are doing   cc
!ccccccccccccccccccccccccccccccccccccc
          if(.not.ladd) then
          if(vdw) then
           write(*,'(5x,a,''/'',a,"-D")') trim(func),trim(bas)
          else
           write(*,'(5x,a,''/'',a)') trim(func),trim(bas)
          endif
          endif


!cccccccccccccccccccccccccccccccccccc
!cccc check the functional type  cccc 
!cccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccc
!cccc adjust aux basis  cccc 
!cccccccccccccccccccccccccccccccccccc
atmp=bas
if(rijk)  bas=trim(bas)//' '//trim(atmp)//'/JK '
if(jbas.and..not.rijk)  bas=trim(atmp)//'/J '
if(cbas) bas=trim(bas)//' '//trim(atmp)//'/C'


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccc   print output file                 cccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
if(ANC) outfile='orca.in'
   inquire(file=outfile,exist=da)
    if(da) then
      call system('cp '//outfile//' input.bak')
!       write(*,*) '  !!!  Warning input file already there  !!!'
!       write(*,*) '  old input now: input.bak'
     endif
if(io.ne.6)open(unit=io,file=outfile)
! START PARSING 
! write(io,'(''# One line for the "simple input" to ensure keyword depence'')')
! NOW WRITE THE SIMPLE INPUT LINE FOR ORCA: ONE LINE PLEASE
write(io,'($,''! '')')
do i=1,maxarg
write(io,'($,'' '',a)') trim(TOT(i))
write(6,'($,'' '',a)') trim(TOT(i))
enddo
write(6,'('' '')')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccc now process defaults of not specified in the arguments  cccc 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if(.not.ladd) then
      if(vdw) write(io,'($,'' VDW10'')')
      if(.not.lfunc) write(io,'($,'' '',a)') trim(func)
      if(.not.lbas) write(io,'($,'' '',a)') trim(bas)
      if(.not.lgrid) write(io,'($,'' GRID'',I1)') grid
      if(COSMO.and..not.lsolvent) write(io,'($,'' COSMO('',a,'')'')') trim(solvent)
endif
write(io,'($,'' PrintGap NoPop '')')
write(io,'(''#NOSOSCF AHSCF '')')
write(io,'(''#CHEAPINTS '')')
write(io,'(''#NOITER '')')
!********************************************
!**************************************** SCF
!********************************************
      write(io,'(''%scf Guess PModel # PAtom Hueckel'')')
      write(io,'(''#%scf Guess MORead'')')
      write(io,'(2x,''#MOInp "orca.gbw" # "orca.ges"'')')
      write(io,'(2x,''#Convergence Tight # loose  normal VeryTight Extreme '')')
      write(io,'(2x,''#LShift 0.25 #(default)'')')
      write(io,'(2x,''MaxIntMem 500 '')')
      write(io,'(2x,''MaxIter '',I5)') siter
      write(io,'(2x,''#Thresh 1e-8 #depends on "Convergence"'')')
      write(io,'(2x,''#TolE 1e-'',I1)')scfcv
write(io,'(''end '')')


!********************************************
!**************************************** MP2
!********************************************
if(MP2) then
      write(io,'(''$mp2 '')')
      write(io,'(''Maxcore 256 # '')')
      write(io,'(''RI true '')')
      write(io,'(''#CalcS2 true '')')
      if(SCS) write(io,'(''DoSCS true '')')
      write(io,'(''end '')')
endif
!********************************************
!*************************************** MDCI
!********************************************
if(MDCI) then
      write(io,'(''#%mdci'')')
      write(io,'(''#citype'',a)') trim(CITYPE)
! if(CEPA) write(io,'(''CEPA_1 '')')
! if(MP3)write(io,'(''MP3'')')
      IF(T)write(io,'(''# Triples 1'')')
      write(io,'(''#KCOpt  '',a)') trim(KCOPT)
      write(io,'(''#MaxIter 35 '')')
      write(io,'(''#TrafoType trafo_jk # trafo_ri '')')

      write(io,'(''#maxCoreWork 500 '')')
      write(io,'(''#MaxCoreIntAmp 500 '')')
endif
!********************************************
!********************************* GEOM BLOCK
!********************************************
if(OPT) then
write(io,'(''%geom '')')
write(io,'(2x,''MaxIter '',I4)') giter
write(io,'(2x,'' Trust '',F4.2,'' # Trust -0.4 # no radius update '')') trustrad
write(io,'(2x,''MaxStep 0.3 '')')
write(io,'(2x,''Inhess Lindh # Read unit Schlegel'')')
write(io,'(2x,'' #InHessName "orca.hess"'')')
if(TS) then
write(io,'(2x,''TS_search '')')
write(io,'(2x,'' TS_Mode {M '',I2,''} '')') mode
write(io,'(2x,'' end'')')
write(io,'(2x,'' #Recalc_Hess 5'')')
endif
if(HOPT) write(io,'(2x,''# Optimizehydrogens true # just relaxing the hydrogends'')')
write(io,'(2x,''Coordsys '',a,'' # deloc, cartesian'')') trim(coordsys)
write(io,'(2x,''Step '',a)') trim(optstep)
write(io,'(2x,''#TolE 5e-6 '')')
write(io,'(2x,''#TolRMSG 1e-'',I1)') gcart
write(io,'(2x,''#TolRMSD 2e-3 '')')
write(io,'(2x,''#TolMaxG 3e-3 '')')
write(io,'(2x,''#TolMaxD 4e-3 '')')
write(io,'(2x,'' '')')
write(io,'(''end '')')
endif



write(io,'(''xyzfile'',I2,'' 1 '',a)') charge,trim(xyzfile)
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')
write(io,'('' '')')


!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccc parsing done! Now additional or alternative files  cc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 If (.not.JOB) write(*,'(5x,''output : '',a)') trim(outfile)
if(io.ne.6)close(io)











end program



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCC           S. Grimme's helper subroutines                   CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!c returns the number of heavy (n>10) atoms in a coord file, and the number of diff atom types
!c sets cu_pd true for some special elements that need more input
!       subroutine atoms(n,nat,nt,cu_pd)
!       implicit none
!       integer n,i,j,nn,nat
!       logical cu_pd
!       character*80 a80
!       real*8 xx(10)
!       integer na(110),nt,na2(110)
! 
!       cu_pd=.false.
!       na = 0
!       na2=0
!       n=0
!       nat=0
!       j=0
!       open(unit=1,file='coord')
!       read(1,'(a)',end=100) a80
!  10   read(1,'(a)',end=100) a80
!       call readl(a80,xx,nn)
!       if(index(a80,'$').ne.0)goto 100
!       if(nn.eq.3)then
!          nat=nat+1
!          j=j+1
!          call elem(a80,i)
!          na2(i)=na2(i)+1
!          if(i.gt.10) na(i)=na(i)+1
! !cts check for Cu/Pd problem
!          if(i.eq.29.or.i.eq.46) cu_pd=.true.
!       endif
!       goto 10
! 100   close(1)
! 
!       n=0
!       do i=1,110
!          if(na(i).gt.0)n=n+1
!          if(na2(i).gt.0) nt=nt+1 ! hok   nt= number of different atom types
!       enddo
! 
!       end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ELEM(KEY1, NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) KEY1
      CHARACTER*2 ELEMNT(107),E

      DATA ELEMNT/'h ','he',&
      'li','be','b ','c ','n ','o ','f ','ne',&
      'na','mg','al','si','p ','s ','cl','ar',&
      'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
      'zn','ga','ge','as','se','br','kr',&
      'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
      'cd','in','sn','sb','te','i ','xe',&
      'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
      'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
      'au','hg','tl','pb','bi','po','at','rn',&
      'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',&
      'fm','md','cb','xx','xx','xx','xx','xx'/
     
      nat=0
      e='  '
      k=1
      DO J=1,len(key1)
         if (k.gt.2)exit
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            call lower(key1(j:j))
            N=ICHAR(key1(J:J))
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
      enddo

      DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

!C     *****************************************************************         

      SUBROUTINE lower(AS)
      CHARACTER*1 AS
      AS=CHAR(ICHAR(AS)-ICHAR('A')+ICHAR('a'))
      END

!C     *****************************************************************         

      SUBROUTINE backstring(A1,A2,lena2)
      CHARACTER*(*) A1
      CHARACTER*(*) A2
      integer n,lena2
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo
      DO J=1,lena2   
         a2(j:j)=' '
      enddo
      a1=a2
      a2='                                                            '
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo

      END


!C     *****************************************************************         
                                                                                
      SUBROUTINE READL(A1,X,N)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*(*) A1                                                      
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READAA(A1,IS,IB,IE)                                               
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READAA                                                             
      CHARACTER*(*) A                                                      
      NINE=ICHAR('9')                                                           
      IZERO=ICHAR('0')                                                          
      MINUS=ICHAR('-')                                                          
      IDOT=ICHAR('.')                                                           
      ND=ICHAR('D')                                                             
      NE=ICHAR('E')                                                             
      IBL=ICHAR(' ')                                                            
      IEND=0                                                                    
      IEND2=0                                                                   
      IDIG=0                                                                    
      C1=0                                                                      
      C2=0                                                                      
      ONE=1.D0                                                                  
      X = 1.D0                                                                  
      NL=LEN(A) 
      DO 10 J=ISTART,NL-1                                                       
         N=ICHAR(A(J:J))                                                          
         M=ICHAR(A(J+1:J+1)) 
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20                 
                               
   10 CONTINUE                                                                  
      READAA=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IDIG=0                                                                    
      DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   50 CONTINUE                                                                  
!C                                                                               
!C PUT THE PIECES TOGETHER                                                       
!C                                                                               
   60 CONTINUE                                                                  
      READAA= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READAA=READAA*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       

!C     *****************************************************************         

