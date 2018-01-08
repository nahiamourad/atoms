 Subroutine accidentalplus(lrho,mmax,Ls,ninter,nband,N,Z,pi,cef,evec,rho,occn,&
                FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,constantM,Robin,Length)
 
  Implicit none
  integer::ninter,nband,N,Z,mn,k,l,i,j,lrho,mmax,Ls,wVxc,normal,config,mi,Robin
  double precision, dimension(lrho+1,N,N)::rho,rho0,rho1
  double precision, dimension(N,N)::KKC,EE,CC
  double precision::pi,t0,cef,temp,Ekin0,Zp,N_p,Length
  double precision::FF(5,5,5,ninter-1),FF1(4,4,4),A(nband+1,N),B(nband+1,N),G(N)
  double precision::evec(mmax+1,(Ls+1)*N,Z),kinM(Ls+1,N,N)
  double precision, dimension(mmax+1,Z)::occn,occn1,occn0

 integer::Ng,Na,l1,l2
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter) 
 double precision::cxc,constantM(Ls+1,Ls+1,2*Ls+1,mmax+1)

 normal=0
 config=0
 N_p=0.d0
 occn1=occn
 occn0=occn
 mi=0
 if(wVxc==0 .and. (Z==41 .or. Z==42))mi=1!! mi=1 change the degenracy bt 5s and 4d into 6s and 4d
 if(wVxc==1)then
        N_p=dble(Z)-36.d0
        occn1(0+1,9+mi)=2.d0*1.d0!5s
        occn1(0+1,10+mi)=(N_p-2.d0)/5.d0!4d_0
        occn1(1+1,5)=2.d0*occn1(0+1,10+mi)!4d_1
        occn1(2+1,2)=2.d0*occn1(0+1,10+mi)!4d_2

       if(N_p<=10.d0)then
        occn0(0+1,9+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,9+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,9+mi)!4d_2
        occn0(0+1,10+mi)=0.d0
       else
        occn0(0+1,9+mi)=2.d0!4d_0
        occn0(1+1,5)=4.d0!4d_1
        occn0(2+1,2)=4.d0!4d_2
        occn0(0+1,10+mi)=N_p-10.d0
       endif

 !print*,'max of difference',maxval(abs(occn-occn0)), maxval(abs(occn-occn1))
 if(maxval(abs(occn-occn0)).ne.0 .and. maxval(abs(occn-occn1)).ne.0)then
        normal=1
        goto 11
 endif  
      if(maxval(abs(occn-occn0))==0)config=1
        occn0(0+1,9+mi)=0.d0
        occn0(1+1,5)=0.d0
        occn0(2+1,2)=0.d0
        occn0(0+1,10+mi)=0.d0
        occn1(0+1,9+mi)=0.d0
        occn1(1+1,5)=0.d0
        occn1(2+1,2)=0.d0
        occn1(0+1,10+mi)=0.d0


       if(occn(0+1,10+mi)==0)then!'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 4d'
        occn1(0+1,10+mi)=N_p!5s
        occn0(0+1,9+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,9+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,9+mi)!4d_2
       else                   !'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 5s 4d'
        occn1(0+1,9+mi)=N_p!5s
        occn0(0+1,10+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,10+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,10+mi)!4d_2
       endif

 elseif(wVxc==0)then
       N_p=dble(Z)-36.d0
       if(mi==1)N_p=dble(Z)-38.d0
        occn1(0+1,9+mi)=2.d0*1.d0!5s
        occn1(0+1,10+mi)=(N_p-2.d0)/5.d0!4d_0
        occn1(1+1,5)=2.d0*occn1(0+1,10+mi)!4d_1
        occn1(2+1,2)=2.d0*occn1(0+1,10+mi)!4d_2

       if(Z .le.46)then
        occn0(0+1,9+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,9+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,9+mi)!4d_2
        occn0(0+1,10+mi)=0.d0!5s
       elseif(Z==47)then
        occn0(0+1,9+mi)=2.d0!4d_0
        occn0(1+1,5)=4.d0!4d_1
        occn0(2+1,2)=4.d0!4d_2
        occn0(0+1,10+mi)=N_p-10.d0!5s
       endif

 !print*,'max of difference',maxval(abs(occn-occn0)), maxval(abs(occn-occn1))
 if(maxval(abs(occn-occn0)).ne.0 .and. maxval(abs(occn-occn1)).ne.0)then
        normal=1
        goto 11
 endif 
      if(maxval(abs(occn-occn0))==0)config=1
        occn0(0+1,9+mi)=0.d0
        occn0(1+1,5)=0.d0
        occn0(2+1,2)=0.d0
        occn0(0+1,10+mi)=0.d0
        occn1(0+1,9+mi)=0.d0
        occn1(1+1,5)=0.d0
        occn1(2+1,2)=0.d0
        occn1(0+1,10+mi)=0.d0


       if((occn(0+1,10+mi)==0 .and. z.le.46))then!'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 4d'
        occn1(0+1,10+mi)=N_p!5s
        occn0(0+1,9+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,9+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,9+mi)!4d_2
       elseif(occn(0+1,10+mi).ne.0 .and. z.le.46)then!'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 5s 4d'
        occn1(0+1,9+mi)=N_p!5s
        occn0(0+1,10+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,10+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,10+mi)!4d_2
       endif

       if((occn(0+1,10+mi)==1.d0 .and. z==47))then!'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s'
        occn1(0+1,10+mi)=N_p!5s
        occn0(0+1,9+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,9+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,9+mi)!4d_2
       elseif(occn(0+1,10+mi).ne.1.d0 .and. z==47)then!'the configuration was 1s 2s 2p 3s 3p 3d 4s 4p 5s 4d'
        occn1(0+1,9+mi)=N_p!5s
        occn0(0+1,10+mi)=N_p/5.d0!4d_0
        occn0(1+1,5)=2.d0*occn0(0+1,10+mi)!4d_1
        occn0(2+1,2)=2.d0*occn0(0+1,10+mi)!4d_2
       endif
 endif!wVxc
  !print*,'sum of occupation numbers',sum(occn1),sum(occn0)
 rho1=0.d0
 rho0=0.d0
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,mmax    !!will be changed later
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

 Ekin0=0.d0
  do  mn=0,mmax    
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
           Ekin0=Ekin0+occn0(mn+1,k)*temp
        endif
      enddo
   enddo
 Ekin0=0.5d0*Ekin0
 
 call oda(lrho,mmax,Ls,ninter,nband,N,Z,pi,cef,rho1,rho0,evec,occn1,Ekin0,&
                FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,1,Robin,Length)!rho=(1.d0-t0)*rho0+t0*rho
  !if(t0==0)t0=0.d0!max(0,1.d0-5.d0/4.d0)
  if(N_p*t0.ge. 2.d0)t0=2.d0/N_p

       rho=rho1
       if(wVxc==1)then  
            if(config==1)then
             occn(0+1,9+mi)=N_p/5.d0*(1.d0-t0)!4d_0
             occn(1+1,5)=2.d0*occn(0+1,9+mi)!4d_1
             occn(2+1,2)=2.d0*occn(0+1,9+mi)!4d_2
             occn(0+1,10+mi)=N_p*t0!5s       
            else
             occn(0+1,9+mi)=N_p*t0!5s
             occn(0+1,10+mi)=N_p/5.d0*(1.d0-t0)!4d_0
             occn(1+1,5)=2.d0*occn(0+1,10+mi)!4d_1
             occn(2+1,2)=2.d0*occn(0+1,10+mi)!4d_2
            endif
        elseif(wVxc==0)then
           if(config==1)then
            occn(0+1,9+mi)=N_p/5.d0*(1.d0-t0)!4d_0
            occn(1+1,5)=2.d0*occn(0+1,9+mi)!4d_1
            occn(2+1,2)=2.d0*occn(0+1,9+mi)!4d_2
            occn(0+1,10+mi)=N_p*t0!5s
           else
            occn(0+1,9+mi)=N_p*t0!5s
            occn(0+1,10+mi)=N_p/5.d0*(1.d0-t0)!4d_0
            occn(1+1,5)=2.d0*occn(0+1,10+mi)!4d_1
            occn(2+1,2)=2.d0*occn(0+1,10+mi)!4d_2
           endif!config=1
        endif
   print*,'---------------------------------------accidental degeneracy on s-shell',N_p*t0,'Z=',sum(occn)
   write(3,*)'------------------------------------accidental degeneracy on s-shell',N_p*t0
 
11 continue
 if(normal==1)then
 rho=0.d0
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,mmax
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
 endif
 End Subroutine
