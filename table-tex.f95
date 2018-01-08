  program table
  implicit none
  integer::i,j,nlines,ncolumns,Z,io,wVxc
  double precision::a,b,c
  double precision, allocatable, dimension(:,:)::eigenvalues
  character*4,allocatable, dimension(:)::insert_cha
  character*15::filename,word(3)
  CHARACTER*200:: fileplace,my_fmt
  CHARACTER(len=255) :: cwd 

 do 11 wVxc=0,1
  nlines=0
  if(wVxc==0)then
    open(unit=11,file='table-tex-xc=0.txt')
  elseif(wVxc==1)then
    open(unit=11,file='table-tex-xc=1.txt')
  endif

  CALL getcwd(cwd)!! this intrinsic function Get current working directory
  WRITE(*,*)'current working directory is', TRIM(cwd)

  do 3 Z=1,54
  write(fileplace,'(a,a,i1,a,i2.2,a)')TRIM(cwd),"/Atom/xc=",wVxc,"/atom",Z,'/'
  write(filename,'(a,i2.2,a)')"Z=",Z,"-R-E-EV.txt"

  !print*,TRIM(ADJUSTL(fileplace)),filename


!!! determine number of lines in the file
 nlines=0
  open(unit=10,file=TRIM(ADJUSTL(fileplace))//filename, status='old',action='read')
10   continue   
     read(10,*,end=12) a
        nlines=nlines+1
goto 10
12 continue
     !print*,'number of lines',nlines
  close(10)

  allocate(eigenvalues(1,50))

!!! determine number of columns containing numbers in the file
 ncolumns=0
   do ncolumns=1,50
  open(unit=10,file=TRIM(ADJUSTL(fileplace))//filename, status='old',action='read')
     read(10,*,iostat=io) (eigenvalues(1,i),i=1,ncolumns),word(1)
     if(io.ne.0)goto 13
  close(10) 
   enddo 
13 continue
  close(10)
   deallocate(eigenvalues) 
   !print*,'number of columns', ncolumns


!!! save the eigenvalues in the array
  allocate(eigenvalues(nlines,ncolumns),insert_cha(ncolumns))
  open(unit=10,file=TRIM(ADJUSTL(fileplace))//filename, status='old',action='read')
   do j=1,nlines
     if(word(1)=='elec.')then
       read(10,*,iostat=io) (eigenvalues(j,i),i=1,ncolumns-1),(word(i),i=1,3),eigenvalues(j,ncolumns)
     else
       read(10,*,iostat=io) (eigenvalues(j,i),i=1,ncolumns-1)     
     endif    
   enddo
  close(10)
  !print*,(eigenvalues(nlines-1,i),i=1,ncolumns) 
  !print*,word,'io=',io

!!! number of  negative eigenvalues for xc=0
 if(wVxc==0)then
  j=0
  do i=3,ncolumns
     if(eigenvalues(nlines-1,i)<0)then
       j=j+1
     else
       goto 14
     endif
  enddo
14 continue
  !print*,'number of negative eigenvalues',j
  elseif(wVxc==1)then
    if(Z.le.2)j=1
    if(Z.ge.3.and.Z.le.10)j=3
    if(Z.ge.11.and.Z.le.12)j=4
    if(Z.ge.13.and.Z.le.18)j=5
    if(Z.ge.19.and.Z.le.20)j=6
    if(Z.ge.21.and.Z.le.30)j=7
    if(Z.ge.31.and.Z.le.36)j=8
    if(Z.ge.37.and.Z.le.38)j=9
    if(Z.ge.39.and.Z.le.44)j=10
    if(Z.ge.45.and.Z.le.46)j=9!!
    if(Z.ge.47.and.Z.le.48)j=10
    if(Z.ge.49.and.Z.le.54)j=11         
  endif

  if(j>1)then
    do i=1,j-1
     insert_cha(i)=' & '
    enddo
  endif  

 if(word(1)=='elec.')then
  insert_cha(j)=' & '
  write(my_fmt, '(a, i0, a)') '(i2,a4,',j, '(f15.6,a4),f15.5,a4)'
  write(11,my_fmt)Z,'&',(eigenvalues(nlines-1,i+2),insert_cha(i),i=1,j),eigenvalues(nlines-1,ncolumns),'\\'
 else
  insert_cha(j)=' \\ '
  write(my_fmt, '(a, i0, a)') '(i2,a4,',j, '(f15.7,a4))'
  write(11,my_fmt)Z,'&',(eigenvalues(nlines-1,i+2),insert_cha(i),i=1,j)
 endif
  write(11,*)'\hline'
  deallocate(eigenvalues,insert_cha)
3  enddo
   close(11)
11 enddo
  end program 

