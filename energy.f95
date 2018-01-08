!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine energy(nband,mmax,Z,evec,N,Ls,nev,kinM,occn,rhoF,Ql,&
                   lrho,KKC,pi,rho,AA,Vxc,Ekin,E,CC,EE,cef,Ekin0,t0,Zp,vs)
 implicit none
 integer:: i,j,k,l,mn,mmax,Z,N,Ls,nev,lrho,nband
 double precision::Ekin,E,temp,pi,cef,Ekin0,t0,Zp,vs
 double precision::evec(mmax+1,(Ls+1)*N,nev),kinM(Ls+1,N,N),AA(N,N),occn(mmax+1,Z),AB(nband+1,N)
 double precision::Ql(lrho+1,N),rhoF(lrho+1,N),KKC(N,N),rho(lrho+1,N,N),Vxc(2*Ls+1,nband+1,N)
 double precision::CC(N,N),EE(N,N)
 E=0.d0
 Ekin=0.d0
  do 7 mn=0,mmax    !!will be changed later
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
7    enddo
 Ekin=0.5d0*Ekin
 Ekin=(1.d0-t0)*Ekin0+t0*Ekin
 write(3,*)'Kinetic Energy=', Ekin
 E=Ekin
   temp=0.d0
    do l=0,lrho
       do k=1,N
          temp=temp+rhoF(l+1,k)*Ql(l+1,k)
       enddo
    enddo
 E=E+0.5d0*temp
   temp=0.d0
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         temp=temp+rho(1,i,j)*(KKC(i,j)+CC(i,j))
      enddo
    enddo
   temp=-sqrt(pi)*dble(Z)*temp
 E=E+temp
  if(Z.ne.Zp)then
   temp=0.d0
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         temp=temp+rho(1,i,j)*CC(i,j)
      enddo
    enddo
   temp=2.d0*sqrt(pi)*(dble(Z)-Zp)*temp
 E=E+temp
 endif
 write(3,*)'Sum of Coulumb Energy=', E-Ekin

 temp=0.d0
 do l=0,lrho
 AB=Vxc(l+1,:,:)
 call construct(nband,N,AB,AA)                                                    
 do i=1,N
    do j=max(1,i-nband),min(i+nband,N)
        temp=temp+rho(l+1,i,j)*AA(i,j)
    enddo
 enddo
 enddo
 write(3,*)'Exchange Energy=',3.d0/4.d0*temp
 E=E+3.d0/4.d0*temp

 temp=0.d0
 if(cef.ne.0)then
    do i=1,N
      do j=max(1,i-nband),min(i+nband,N)
         temp=temp+rho(2,i,j)*EE(i,j)
      enddo
    enddo
 temp=cef*temp
 endif
 
 vs=2*Ekin+(E-Ekin)+temp
 write(3,*) 'Electric field Energy=',-temp
 write(3,*) 'Virial sum',  vs
 E=E-temp
 write(3,*)'Total Energy= ', E

 end subroutine
