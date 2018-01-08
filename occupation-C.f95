 Subroutine occupationC(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,rho1,evec,occn,&
                FF,FF1,A,B,G,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,Robin,Length,constantM,hk,M,rk)
 implicit none
 integer::lrho,N,Ls,Z,ms,nband,ninter,Ng,Na,wVxc,Robin
 double precision::constantM(Ls+1,Ls+1,2*Ls+1,ms+1),evec(ms+1,(Ls+1)*N,Z),Length,Zp
 double precision::occn(ms+1,Z),kinM(Ls+1,N,N),FF(5,5,5,ninter-1),FF1(4,4,4),pi
 double precision::A(nband+1,N),B(nband+1,N),KKC(N,N),EE(N,N),CC(N,N),cxc,cef
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter),G(N)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)

 integer::i,j,k,l,l1,l2,mn
 double precision::rho0(lrho+1,N,N),rho1(lrho+1,N,N)
 double precision::t1,temp,Ekin0_C,occn1(ms+1,Z),occn0(ms+1,Z)
  double precision::M(nband+1,N),hk(ninter),rk(ninter)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! occupation number for the atom 6 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            occn0=0.d0
            occn0(1,1)=2.d0 
            occn0(1,2)=2.d0
            occn0(1,3)=(dble(z)-4.d0)!2.d0
            occn1=0.d0 
            occn1(1,1)=2.d0 
            occn1(1,2)=2.d0
            occn1(2,1)=(dble(z)-4.d0)!2.d0

 rho1=0.d0
 rho0=0.d0
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,ms    !!will be changed later
       do k=1,Z
           if(occn1(mn+1,k).ne.0 .or. occn0(mn+1,k).ne.0) then
          do l1=mn,Ls
            do l2=max(mn,abs(l-l1)),min(Ls,l+l1)!mn,Ls!
            rho1(l+1,i,j)=rho1(l+1,i,j)+occn1(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                 evec(mn+1,l1*N+i,k)*evec(mn+1,l2*N+j,k) 
            rho0(l+1,i,j)=rho0(l+1,i,j)+occn0(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                 evec(mn+1,l1*N+i,k)*evec(mn+1,l2*N+j,k)   
          enddo
        enddo 
            endif
      enddo
    enddo
   enddo
  enddo
 enddo

 Ekin0_C=0.d0
  do  mn=0,ms    
      do k=1,Z       
        if(occn0(mn+1,k).ne. 0) then
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+evec(mn+1,l*N+i,k)*kinM(l+1,i,j)*evec(mn+1,l*N+j,k)
               enddo
             enddo
           enddo
           Ekin0_C=Ekin0_C+occn0(mn+1,k)*temp
        endif
      enddo
   enddo
 Ekin0_C=0.5d0*Ekin0_C

 call  oda(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,rho1,rho0,evec,occn1,Ekin0_C,&
                FF,FF1,A,B,G,t1,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,1,Robin,Length,hk,M,rk)
  
  occn=(1.d0-t1)*occn0+t1*occn1
  if((sum(occn)-dble(Z)).ge.1.d-5)then
     print*, 'occupation-C error'
     stop
  endif

 write(3,*)'on m=0 we have',occn(1,3),'electrons'
 write(3,*)'on m=1 we have',occn(2,1),'electrons'     
 print*,'on m=0 we have',occn(1,3),'electrons' 
 print*,'on m=1 we have',occn(2,1),'electrons' 
end subroutine
