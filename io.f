!  routines for in- and output
!  xzyrd: reading xyz files
!  wxyz: write xyz files
!  tmolrd: read turbomole (tmol) coord files
!  
!
      subroutine xyzrd(xyz,iat,nat,infile)
      implicit none
      character*2 elemnt(107)
      character*80 infile, outfile,atmp
      real*8 xyz(3,10000),xx(5)
      integer iat(10000),nat,nel,i,nn
      real*8 bohr
      logical da
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

      bohr=0.52917726
      nat=0
      inquire(file=infile,exist=da)

      if(da)then

      write(*,'(5x,''reading...'',$)')
! read XYZ file
      open(unit=3,file=infile)
       read(3,'(a)',end=100) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do i=1,10000   ! while loop would be better
            nat=nat+1
            read(3,'(a)',end=123) atmp  
           enddo 
          else
            nat=idint(xx(1))
           read(3,'(a)',end=100) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(3,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,iat(i))
            xyz(1:3,i)=xx(1:3)
!       write(*,'(a2,5x,3F18.12)') elemnt(iat(i)),xyz(1:3,i)
       enddo
 100  close(3)
      write(*,'(5x,''XYZ file : '',a)')  trim(infile)
      else
      write(*,*) ' no input file found !! '
      endif

      write(*,*) '    number of atoms:  ',nat
      end

      subroutine wxyz(xyz,iat,nat,io)
      implicit none
      character*2 elemnt(107)
      character*80 outfile,atmp
      real*8 xyz(3,10000),xx(5)
      integer iat(10000),nat,nel,i,nn,io
      real*8 bohr
      logical da
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

!       open(unit=4,file=outfile,access='append')
!       write(*,*) 'writing coords'
       write(io,*) nat
       write(io,*) '   '
       do i=1,nat
       write(io,'(a2,5x,3F18.12)') elemnt(iat(i)),xyz(1:3,i)
       enddo
!       close(4)
       end



      subroutine tmolrd(xyz,iat,nat,infile)
      implicit none
      character*2 elemnt(107)
      character*80 infile, outfile,atmp
      real*8 xyz(3,10000),xx(5)
      integer iat(10000),nat,nel,i,nn
      real*8 bohr
      logical da
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

      bohr=0.52917726
      i=0
      inquire(file=infile,exist=da)
      if(da)then
      write(*,'(''reading...'',$)')
! read TMOL file
      open(unit=3,file=infile)
      do while (da)
       read(3,'(a)',end=100) atmp ! $coord
        if(index(atmp,'$coord').ne.0) cycle     
        if(index(atmp,'$').ne.0) exit
        i=i+1     
        call readl(atmp,xx,nn)
        call elem(atmp,iat(i))
        xyz(1:3,i)=xx(1:3)*bohr
!        write(*,'(a2,5x,3F18.12)') elemnt(iat(i)),xyz(1:3,i)
      enddo
      nat=i
      write(*,*) 'turbomole file :  ', trim(infile)
 100  close(3)
      else
      write(*,*) ' no input file found !! '
      endif
      write(*,*) 'number of atoms:  ',nat
      end


! write XYZ from Fragment file, only slightly modified wxyz
      subroutine wfrag(xyz,iat,nat,io,frag)
      implicit none
      character*2 elemnt(107)
      character*80 outfile,atmp,frag
      real*8 xyz(3,10000),xx(5)
      integer iat(10000),nat,nel,i,nn,io
      real*8 bohr
      logical da
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

!       open(unit=4,file=outfile,access='append')
!       write(*,*) 'writing coords'
       do i=1,nat
       write(io,'(a2,3x,3F18.12,2x,2a)')
     .        elemnt(iat(i)),xyz(1:3,i),'f=',trim(frag)
       enddo
!       close(4)
       end

      subroutine wtm(xyz,iat,nat,io)
      implicit none
      character*2 elemnt(107)
      character*80 outfile,atmp
      real*8 xyz(3,10000),xx(5)
      integer iat(10000),nat,nel,i,nn,io
      real*8 bohr
      logical da
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

        bohr=0.52917726  
!       open(unit=4,file=outfile,access='append')
!       write(*,*) 'writing coords'
       write(io,'(a)')'$coord'
       do i=1,nat
       write(io,'(3F18.12,2x,a2)') xyz(1:3,i)/bohr , elemnt(iat(i))
       enddo
       write(io,'(a)')'$end'
!       close(4)
       end
