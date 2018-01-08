 Subroutine bisection(x,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,&
                      wa,hr1,FormMat,tFormMat,LegMat,thr13,cxc,indicator,Z,&
                      wVxc)

 Implicit none
 integer::N,nband,ninter,Ng,lrho,Na,Ls,i,j,l,wVxc,k,indicator,Z
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng),AA(N,N),AB(nband+1,N)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision, dimension(:,:,:)::Vxc(2*Ls+1,nband+1,N)
 double precision::a,b,cxc,temp
 double precision, dimension(lrho+1,N,N)::rho,rho0,rho1
 double precision::x0,x1,x,y,y0


!!Bisection method to find the root of the  derivative of the energy
 !Subroutine Bisect(e,m,x,x0,x1)
 ! integer m  
 ! real*8 e,x,x0,x1,y0,yy,Y
 x0=0.d0
 x1=1.d0

 if(indicator==1)then
  if(wVxc==1)then
!    x1=2.d0/(dble(Z)-18.d0)!!!!!1/(N_p/2) !!!!!!!!!!!!!!!pay attention to this probably should be changed
  else
    print*,'?????????????????,',Z
  endif
 endif
  k=0
13 continue
 rho1=rho-rho0

 call exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho0+x0*rho1,&
                       lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,2)

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

  y0=a+x0*b+temp  !!y0=Y(x0)
  x=(x0+x1)/2.d0 !!x=(x0+x1)/2.d0

  call exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho0+x*rho1,&
                       lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,2)
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

  y=a+x*b+temp!!yy=Y(x) 
  k=k+1
  if(k.ge. 80)return
 ! print*,k
 ! print*,x0,y0,x,y
 ! if (dabs(y*y0).le. 1E-20) goto 6!!(dabs(yy*y0).eq.0) return
 ! the above condition should be commented out since in the case of radial 
 !\rho_1=0 and then will exit from the first iteration
  if ((y*y0)<0.d0) x1=x!!x1=x
  if ((y*y0)>0.d0) x0=x!!x0=x
  if (dabs(x1-x0)>1E-10) goto 13
 end Subroutine
