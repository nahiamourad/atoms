 Subroutine  test(N,nband,A,Ng,xg,wg,ninter,hk,rk,Robin,Length,FormMat)

 Implicit none
 integer:: N,Ng,Robin,ninter,nband,Lambda
 double precision::A(nband+1,N),wg(Ng),xg(Ng),hk(ninter),rk(ninter),FormMat(5,Ng),Length
 integer::p,q,k,i,j
 double precision::rHS(N),right_integral,AB(nband+1,N),temp(Ng,ninter)


!!!!!The right hand side vector
 rHS=0.d0
   do k=1,ninter
     if(k>=2) then
      rHS(4*(k-1))=right_integral(Ng,xg,wg,hk(k-1),rk(k-1),2,FormMat)+&
                 right_integral(Ng,xg,wg,hk(k),rk(k),1,FormMat)
     endif
      do i=1,3
      rHS(4*(k-1)+i)=right_integral(Ng,xg,wg,hk(k),rk(k),i+2,FormMat)
      enddo
   enddo    
      
 if(Robin==1 .or. Robin==2)rHS(4*ninter)=right_integral(Ng,xg,wg,hk(ninter),rk(ninter),2,FormMat)
 if(Robin==2)then
      temp(1,1)=0.d0
    do i=1,Ng
       temp(1,1)=temp(1,1)+wg(i)*(3.d0*xg(i)**2-Length**2)/(xg(i)**2+Length**2)**3*(Length/xg(i)**2)
    enddo
    rHS(4*ninter)=rHS(4*ninter)+2.d0*temp(1,1)
 endif

  AB=A
  if(Robin==1)AB(nband+1,N)=AB(nband+1,N)+1.d0/Length
 call DPBSV('U', N,nband, 1, AB, nband+1,rHS, N, i)
 if(i.ne.0)print*,'in solving the equation...info=',i

 open(unit=2,file='test.txt')
 do k=1,ninter
  do j=1,Ng 
      p=4
      q=0
      if(k==1)p=3
      if(k==ninter .and. Robin==0)q=1
     temp(j,k)=0.d0
      do i=q,p
         temp(j,k)=temp(j,k)+rHS(4*k-i)*FormMat(Lambda(i),j)
      enddo 
    write(2,*)xg(j)*hk(k)+rk(k),temp(j,k), (xg(j)*hk(k)+rk(k))/(1.d0+(xg(j)*hk(k)+rk(k))**2)
  enddo
 enddo
 close(2)
 stop
 End Subroutine
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function  right_integral(Ng,xg,wg,h,r,i,FormMat)

  implicit none
    integer::Ng,nn,i
    double precision:: right_integral,wg(Ng),xg(Ng),FormMat(5,Ng),h,r 
  right_integral=0.d0
    do nn=1,Ng
       right_integral=right_integral+wg(nn)*FormMat(i,nn)*&
                      (xg(nn)*h+r)*(3.d0-(xg(nn)*h+r)**2)/(1.d0+(xg(nn)*h+r)**2)**3
    enddo
       right_integral=2.d0*right_integral*h 
 end function
