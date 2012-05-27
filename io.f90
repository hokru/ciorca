subroutine wxyz(xyz,iat,nat,io)
implicit none
character*2 elemnt(107)
character*80 outfile,atmp
real*8 xyz(3,10000),xx(5)
integer iat(10000),nat,nel,i,nn,io
real*8 bohr
logical da
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

      bohr=0.52917726
       write(io,*) nat
       write(io,*) '   '
       do i=1,nat
       write(io,'(a2,5x,3F18.12)') elemnt(iat(i)),xyz(1:3,i)
       enddo
       end



subroutine tmolrd(xyz,iat,nat,infile)
implicit none
character*2 elemnt(107)
character*80 infile, outfile,atmp
real*8 xyz(3,10000),xx(5)
integer iat(10000),nat,nel,i,nn
real*8 bohr
logical da
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
      enddo
      nat=i
      write(*,*) 'turbomole file :  ', trim(infile)
 100  close(3)
      else
      write(*,*) ' no input file found !! '
      endif
      write(*,*) 'number of atoms:  ',nat
      end


subroutine wtm(xyz,iat,nat,io)
implicit none
character*2 elemnt(107)
character*80 outfile,atmp
real*8 xyz(3,10000),xx(5)
integer iat(10000),nat,nel,i,nn,io
real*8 bohr
logical da
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



!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd2(maxat,xyz,iat,ifrez,nat,infile,echo)
implicit none
integer maxat
character*2 cc,ff
character*80  atmp
character*(*) infile
real(kind=8) xyz(3,maxat),xx(5)
real(kind=8) bohr
integer iat(maxat),nat,i,nn,j,ifrez(maxat),istat
logical da,echo

bohr=0.52917726d0
i=0
ifrez=0

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(*,'('' reading...'',$)')

open(unit=3,file=infile)
! test for tmol or xyz file

 read(3,'(a)') atmp ! $coord
rewind(3)
if(index(atmp,'$coord').ne.0) then
 ! count number of atoms
 do while (da)
  read(3,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue

 rewind(unit=3)

 ! read TMOL file
 read(3,*) atmp ! $coord
 do j=1,nat
    read(3,'(a)') atmp ! $coord
    backspace(3)
    if(index(atmp,' f ').ne.0) then
    read(3,*) xyz(1,j),xyz(2,j),xyz(3,j),cc,ff
     ifrez(j)=1
    else ! default
     read(3,*) xyz(1,j),xyz(2,j),xyz(3,j),cc
   endif
   call elem(cc,iat(j))
   xyz(1:3,j)=xyz(1:3,j)*bohr
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(3)
 else ! xyz file
       read(3,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming theyare coords
           do
            nat=nat+1
            read(3,'(a)',end=123) atmp
           enddo
          else
            nat=idint(xx(1))
           read(3,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(3,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,iat(i))
            xyz(1:3,i)=xx(1:3)
       enddo
 101  close(3)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif

if(maxval(ifrez,nat).eq.1) then
 if(echo) then
  print*,'  found frozen cart. coordinates'
  if(nat.lt.50) then ! dont spam the output to much ...
   write(*,'($,a)') '  '
   do i=1,nat
    write(*,'($,x,I1)') ifrez(i)
   enddo
   print*,''
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select
end subroutine


!********************************
!* convert a word to lower case *
!********************************
subroutine lower_case(word)
character (len=*) , intent(in out) :: word
integer :: i,ic,nlen
nlen = len(word)
do i=1,nlen
ic = ichar(word(i:i))
if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
end do
end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine

!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
return
end subroutine


