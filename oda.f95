 Subroutine oda(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,rho,rho0,evec,occn,Ekin0,&
                FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,indicator,Robin,Length,hk,M,rk)
 
  Implicit none
  integer::ninter,nband,N,Z,mn,k,l,i,j,lrho,ms,Ls,wVxc,Robin
  double precision, dimension(lrho+1,N,N)::rho1,rho,rho0
  double precision, dimension(lrho+1,N)::Ql0,Ql1,rhoF0
  double precision, dimension(N,N)::KKC,EE,CC
  double precision::Ekin,pi,t0,cef,temp,Ekin0,Zp,Length
  double precision::FF(5,5,5,ninter-1),FF1(4,4,4),A(nband+1,N),B(nband+1,N),G(N),Gtest(9)
  double precision::evec(ms+1,(Ls+1)*N,Z),kinM(Ls+1,N,N),occn(ms+1,Z)

 integer::Ng,Na,indicator
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision::cxc

  double precision::M(nband+1,N),hk(ninter),rk(ninter)

 t0=0.d0
 !kinetic energy
 Ekin=0.d0
  do  mn=0,ms    !!will be changed later
      do k=1,Z       !!will be changed later
        if(occn(mn+1,k).ne. 0) then
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+evec(mn+1,l*N+i,k)*kinM(l+1,i,j)*evec(mn+1,l*N+j,k)
               enddo
             enddo
           enddo
           Ekin=Ekin+occn(mn+1,k)*temp
        endif
      enddo
   enddo
 Ekin=0.5d0*Ekin



  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho0,Ql0,rhoF0,Robin,Length,hk,M,rk)
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho,Ql1,rhoF0,Robin,Length,hk,M,rk) 
   Ql1=Ql1-Ql0
   rho1=rho-rho0

  Gtest=0.d0
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         Gtest(4)=Gtest(4)+rho1(1,i,j)*CC(i,j) !!Gtest(4)=\int -z/r(\rho-\rho_0)
         Gtest(6)=Gtest(6)+rho1(1,i,j)*(CC(i,j)-KKC(i,j))
         if(cef.ne.0)Gtest(7)=Gtest(7)+rho1(2,i,j)*EE(i,j)
      enddo
    enddo
   Gtest(4)=-2.d0*sqrt(pi)*dble(Zp)*Gtest(4)
   Gtest(6)=2.d0*sqrt(pi)*dble(Z)*Gtest(6)
   Gtest(7)=-cef*Gtest(7)

    do l=0,lrho
       do i=1,N
          do j=max(i-nband,1),min(i+nband,N)
            Gtest(5)=Gtest(5)+Ql0(l+1,i)*KinM(l+1,i,j)*Ql1(l+1,j)
            Gtest(3)=Gtest(3)+Ql1(l+1,i)*KinM(l+1,i,j)*Ql1(l+1,j)
         enddo
       enddo
    enddo
  Gtest(5)=Gtest(5)/(4.d0*pi)!!Gtest(5)+Gtest(6)=D(\rho_0,\rho-\rho_0)
  Gtest(3)=Gtest(3)/(4.d0*pi)!!Gtest(3)=D(\rho_1,\rho_1)

  Gtest(1)=Gtest(4)+Gtest(5)+Gtest(6)+Gtest(7)
  Gtest(2)=Gtest(3)

 if(wVxc==0)then
  t0=-(Ekin-Ekin0+Gtest(1))/Gtest(2)
  if(t0<0 .and. indicator==0)then
       print*,(Ekin-Ekin0+Gtest(1)),Gtest(2),'the step legth is negative!! the program is stopped','t0=',t0
       write(3,*)(Ekin-Ekin0+Gtest(1)),Gtest(2),'the step legth is negative!! the program is stopped','t0=',t0
       stop
  elseif(t0<0 .and. indicator==1)then
       t0=0.d0
  endif
  if(indicator==0)then!call from the program
     print*,'step length without exchange potential',t0!,&
            !'D(\rho-\rho_0,\rho-\rho_0)=',Gtest(2),Ekin-Ekin0+Gtest(1)!,Ekin,Ekin0,Gtest(1)
     write(3,*) 'the step legth is', t0 
  endif
 elseif(wVxc==1)then
    Gtest(3)=0.d0
    Gtest(4)=1.d0
    call MNBRAK(Gtest(3),Gtest(4),Gtest(5),Gtest(6),Gtest(7),Gtest(8),(Ekin-Ekin0+Gtest(1)),Gtest(2),&
              rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
              FormMat,tFormMat,LegMat,thr13,cxc)
    CALL golden(Gtest(3),Gtest(4),Gtest(5),1.d-5,t0,(Ekin-Ekin0+Gtest(1)),Gtest(2),&
                rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,Gtest(9))
   if(indicator==0)then!call from the program
      print*,'step length using golden',t0
      write(3,*)'step length using golden',t0
      if( t0<0)then!
         t0=0.7d0!before 20oct2016 I used t0=0.7 and did not use the bisection
         call bisection(t0,(Ekin-Ekin0+Gtest(1)),Gtest(2),rho,rho0,N,nband,ninter,&
                        Ng,lrho,Na,Ls,wg,xg,wa,hr1,FormMat,tFormMat,LegMat,thr13,cxc,&
                        indicator,Z,wVxc)
         print*,'step length using bisection',t0
         write(3,*)'step length using bisection',t0
      endif
   endif!indicator==0
 endif!wVxc==0
    if(t0<0)t0=0.d0
    if(t0>1)t0=1.d0
    if(indicator==0)write(3,*)'The step length is',t0
   rho=(1.d0-t0)*rho0+t0*rho
 end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine xcfunc(t,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,f)
 Implicit none
 integer::N,nband,ninter,Ng,lrho,Na,Ls,i,j,l,wVxc
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng),AA(N,N),AB(nband+1,N)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision, dimension(:,:,:)::Vxc(2*Ls+1,nband+1,N)
 double precision::f,t,a,b,cxc,temp
 double precision, dimension(lrho+1,N,N)::rho,rho0,rho1

 wVxc=2
 f=t*a+t**2/2*b
 rho1=(1.d0-t)*rho0+t*rho
 call exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho1,&
                       lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,wVxc)
  Vxc=-cxc*Vxc
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
 if(temp>0)temp=-temp
 f=f+3.d0/4.d0*temp
END subroutine

