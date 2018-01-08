!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Function CG(l1,l2,l,m,CGcoef,CGcoefint,nwigner)
  implicit none
  integer:: l1,l2,l,m,nwigner,i
  double precision:: CG,CGcoef(nwigner)
 integer, dimension(:,:)::CGcoefint(nwigner,4)
     
     CG=0.d0
      do i=1,nwigner
       if(l1==CGcoefint(i,1) .and. l2==CGcoefint(i,2) .and.&
          l==CGcoefint(i,3) .and. m==CGcoefint(i,4)) CG=CGcoef(i)
      enddo

 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Function constant(l1,l2,l,m,pi,CGcoef,CGcoefint,nwigner)
 implicit none 
 integer::l1,l2,l,m,nwigner
 double precision:: pi,constant,CG,CGcoef(nwigner)
 integer, dimension(:,:)::CGcoefint(nwigner,4)
 constant=0.d0
 if(abs(l1-l2).le.l.and.l.le.l1+l2.and.m.le.l2)then 
   constant=CG(l1,l2,l,0,CGcoef,CGcoefint,nwigner)
   if(constant.ne.0.) constant=sqrt(dble((2*l+1)*(2*l1+1)*(2*l2+1))/(4.d0*pi))*&
                                 CG(l1,l2,l,m,CGcoef,CGcoefint,nwigner)*constant

  if(l1<l2) print*,'You mistake here'
 endif
 end function  
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function F(i,j,n,h,r,wg,xg,Ng,FormMat,tFormMat)
 Implicit none
 integer::i,j,n,Ng,p
 double precision:: F,h,r
 double precision, dimension(:)::wg(Ng),xg(Ng)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 F=0.d0
 if(r .ne. 0) then
  do p=1,Ng
    F=F+wg(p)*FormMat(i,p)*FormMat(j,p)*FormMat(n,p)/(xg(p)*h+r)
  enddo
  F=h*F
 elseif(r==0) then
   if(n.ne.1) then 
   goto 1
   elseif(n==1.and.i.ne.1) then 
      n=i
      i=1
      goto 1
    elseif(n==1.and.j.ne.1) then
      n=j
      j=1
     else
      print*, 'There is something wrong with the F function'
   endif
1 continue
  do p=1,Ng
    F=F+wg(p)*FormMat(i,p)*formMat(j,p)*tFormMat(n-1,p)
  enddo
 endif
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Subroutine construct(nband,N,Banded,Matrix)
  Implicit none
  integer:: i,j,nband,N
  double precision:: Matrix(N,N),Banded(nband+1,N)
  Matrix=0.d0
   do j=1,N
    do i=max(1,j-nband),j
        Matrix(i,j)=Banded(nband+1+i-j,j)
    enddo
   enddo 
 
    do i=2,N
    do j=1,i-1
       Matrix(i,j)=Matrix(j,i)
    enddo
    enddo  
  end Subroutine

