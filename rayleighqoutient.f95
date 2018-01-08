 Subroutine rayleighqoutient(lrho,mmax,Ls,ninter,nband,N,Z,pi,cef,evec,occn,&
                FF,FF1,A,B,G,KKC,EE,kinM,constantM,ev,Robin,Length,Res,cxc,&
                Ng,xg,wg,Na,wa,LegMat,FormMat,tFormMat,hr1,thr13,wVxc)
  Implicit none
  integer::ninter,nband,N,Z,mn,k,l,i,j,lrho,mmax,Ls,l1,l2,Robin
  double precision, dimension(lrho+1,N,N)::rho1,rho
  double precision, dimension(lrho+1,N)::Ql0,Ql1,rhoF0,rhoF1
  double precision, dimension(N,N)::KKC,EE
  double precision::Ekin,pi,cef,temp,temp2,Length,Res
  double precision::FF(5,5,5,ninter-1),FF1(4,4,4),A(nband+1,N),B(nband+1,N),G(N),GG(N)
  double precision::evec(mmax+1,(Ls+1)*N,Z),kinM(Ls+1,N,N),occn(mmax+1,Z),ev(Z,mmax+1)
  double precision::constantM(Ls+1,Ls+1,2*Ls+1,mmax+1),ev1(Z,mmax+1)
  !%%%special for exchange correlation
 integer::Ng,Na,wVxc
 double precision::cxc,AB(nband+1,N),AA(N,N)
 double precision::wg(Ng),xg(Ng),wa(Na),hr1(ninter),thr13(Ng,ninter)
 double precision::FormMat(5,Ng),tFormMat(4,Ng),LegMat(2*Ls+1,Na)
 double precision::Vxc(2*Ls+1,nband+1,N)
 
 rho=0.d0
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,mmax    !!will be changed later
       do k=1,Z
           if(occn(mn+1,k).ne.0) then
          do l1=mn,Ls
            do l2=max(mn,abs(l-l1)),min(Ls,l+l1)!mn,Ls!
               rho(l+1,i,j)=rho(l+1,i,j)+occn(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                                         evec(mn+1,l1*N+i,k)*evec(mn+1,l2*N+j,k)  
            enddo
          enddo 
           endif
      enddo
    enddo
   enddo
  enddo
 enddo  



 if(wVxc.eq.1)then
 call exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho,&
                           lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,wVxc)
  Vxc=-cxc*Vxc
 endif
 call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho,Ql0,rhoF0,Robin,Length)

 Res=0.d0
 rho1=0.d0
  do mn=0,mmax    !!will be changed later
   do k=1,Z
      if(occn(mn+1,k).ne.0) then
     rho1=0.d0
       do l=0,lrho
         do i=1,N
           do j=1,N
             do l1=mn,Ls
              do l2=max(mn,abs(l-l1)),min(Ls,l+l1)!mn,Ls!
              rho1(l+1,i,j)=rho1(l+1,i,j)+constantM(l1+1,l2+1,l+1,mn+1)*& 
                                evec(mn+1,l1*N+i,k)*evec(mn+1,l2*N+j,k)  
              enddo
            enddo 
          enddo
        enddo
      enddo

  GG=G/(dble(Z))   !!! just one eigenvector
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,GG,A,B,rho1,Ql1,rhoF1,Robin,Length) !!since\int\rho_1=1

 temp2=0.d0
   temp=0.d0
    do l=0,lrho
       do i=1,N
          temp=temp+rhoF1(l+1,i)*Ql0(l+1,i)
       enddo
    enddo
  temp2=temp
   temp=0.d0
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         temp=temp+rho1(1,i,j)*KKC(i,j) !!!no more C
      enddo
    enddo
   temp=-2.d0*sqrt(pi)*dble(Z)*temp
  temp2=temp2+temp

  if(cef.ne.0)then 
  temp=0.d0
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         temp=temp+rho1(2,i,j)*EE(i,j)
      enddo
    enddo
 temp=-cef*temp
 temp2=temp2+temp
 endif

 if(wVxc.eq.1)then
 temp=0.d0
 do l=0,lrho
 AB=Vxc(l+1,:,:)
 call construct(nband,N,AB,AA)                                                    
 do i=1,N
    do j=max(1,i-nband),min(i+nband,N)
        temp=temp+rho1(l+1,i,j)*AA(i,j)
    enddo
 enddo
 enddo
 temp2=temp2+temp
 endif

 !kinetic energy
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+evec(mn+1,l*N+i,k)*kinM(l+1,i,j)*evec(mn+1,l*N+j,k)
               enddo
             enddo
           enddo
 Ekin=0.5d0*temp

        ev1(k,mn+1)=Ekin+temp2
      Res=max(Res,abs(ev(k,mn+1)-ev1(k,mn+1)))
      print*,Res,ev(k,mn+1),ev1(k,mn+1)
       endif
    enddo
   enddo
 end subroutine
