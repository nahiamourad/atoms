Subroutine eigenfunctions(evec,mmax,Ls,N,Z,Ng,xg,hk,rk,ninter,FormMat,Robin)

 Implicit none
 integer::mmax,Ls,N,Z,Ng,ninter,Lambda,Robin
 double precision::evec(mmax+1,(Ls+1)*N,Z),xg(Ng),hk(ninter),rk(ninter),FormMat(5,Ng)

 integer::i,j,k,p,q,mn,nn,finish(mmax+1),start(mmax+1),l(mmax+1,Z/2+1),pl
 double precision ::temp(mmax+1,Z/2+1,Ng,ninter)
 character(len=30):: filename

 l=0!!! this l is designed to pick the part of the eigenfunction which is not zero 
 finish=0
 start=0
 if(Z.ge. 5 .and. Z.le. 10)then
  start(1)=1
  finish(1)=4
  start(2)=1
  finish(2)=1
  l(0+1,3)=1!the p part
  pl=1
 elseif(Z==11 .or. Z==12)then
  start(1)=1
  finish(1)=4
  start(2)=1
  finish(2)=1
  l(0+1,3)=1
  pl=1
 elseif(Z.ge. 23 .and. Z.le. 26)then
  start(1)=6
  finish(1)=8
  start(2)=3
  finish(2)=3
  start(3)=1
  finish(3)=1
  l(0+1,8)=2!the d part
  l(1+1,3)=1!the d part
  pl=2
 elseif(Z.ge. 41 .and. Z.le. 42)then
  start(1)=9
  finish(1)=11
  start(2)=5
  finish(2)=5
  start(3)=2
  finish(3)=2
  l(0+1,11)=2!the d part
  l(1+1,5)=1!the d part
  pl=2
 else
  print*, 'not defined see eigenfunctions.f95'
  stop
 endif

 write(filename,'(i2.2,a)')Z,"_eigenfunction*r.txt"
 open(unit=8,file=filename)

 do k=1,ninter
  do j=1,Ng 
   do mn=0,mmax
      if(finish(mn+1)==0)goto 1
     do nn=start(mn+1),finish(mn+1)
      p=4
      q=0
      if(k==1)p=3
      if(k==ninter .and. Robin==0)q=1
     temp(mn+1,nn,j,k)=0.d0
      do i=q,p
         temp(mn+1,nn,j,k)=temp(mn+1,nn,j,k)+evec(mn+1,(mn+l(mn+1,nn))*N+4*k-i,nn)*FormMat(Lambda(i),j)
      enddo 
    enddo
1 continue
   enddo
  write(8,*)xg(j)*hk(k)+rk(k),((temp(mn+1,p,j,k),p=start(mn+1),finish(mn+1)),mn=0,pl)!,&
                                             !exp(-0.5d0*(xg(j)*hk(k)+rk(k)))
  enddo
 enddo

 close(8)
End Subroutine
