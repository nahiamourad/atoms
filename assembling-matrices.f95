!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! elementary assembling matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Subroutine assemblingmatrices(N,ninter,nband,hk,rk,Ng,xg,wg,eta,A,B,C,M,E,KC,G,&
                               CC,EE,MM,Ls,KinM,KKC,FormMat,tFormMat,BB,AA,hr1,hr2,&
                               Robin,Length)
 Implicit none
 integer::i,j,k,Lambda,N,ninter,nband,Ng,Ls,l,Robin
 double precision:: integrate,eta,xg(Ng),wg(Ng),G(N),Length,temp
 double precision, dimension(:)::rk(ninter),hk(ninter),hr1(ninter),hr2(ninter)
 double precision, dimension(5,5)::alpha,mu,beta
 double precision, dimension(nband+1,N)::A,B,C,M,E,KC
 double precision, dimension(N,N)::AA,BB,CC,MM,EE,KKC
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:,:)::kinM(Ls+1,N,N)


!!  Mass and radial kinetic energy elementary assembling matrices
!tol=1.d-10;
 alpha=0.d0
 alpha(1,2)=-1.d0
 alpha(3,4)=128.d0/45.d0
 alpha(3,5)=128.d0/45.d0
 alpha(4,5)=5888.d0/945.d0
 alpha=alpha+transpose(alpha)
 alpha(1,1)=1.d0
 alpha(2,2)=1.d0
 alpha(3,3)=16.d0/3.d0
 alpha(4,4)=3328.d0/189.d0
 alpha(5,5)=3328.d0/189.d0

 mu=0.d0
 mu(1,2)=1.d0/6.d0
 mu(1,3)=1.d0/3.d0
 mu(1,4)=4.d0/15.d0
 mu(1,5)=4.d0/45.d0
 mu(2,3)=1.d0/3.d0
 mu(2,4)=4.d0/45.d0
 mu(2,5)=4.d0/15.d0
 mu(3,4)=64.d0/315.d0
 mu(3,5)=64.d0/315.d0
 mu(4,5)=128.d0/2835.d0
 mu=mu+transpose(mu)
 mu(1,1)=1.d0/3.d0
 mu(2,2)=1.d0/3.d0
 mu(3,3)=8.d0/15.d0
 mu(4,4)=128.d0/405.d0
 mu(5,5)=128.d0/405.d0

 beta=0.d0
 beta(1,2)=1.d0/12.d0
 beta(1,3)=2.d0/15.d0
 beta(1,4)=16.d0/315.d0
 beta(1,5)=16.d0/315.d0
 beta(2,3)=1.d0/5.d0                                 
 beta(2,4)=4.d0/105.d0
 beta(2,5)=68.d0/315.d0
 beta(3,4)=16.d0/315.d0
 beta(3,5)=16.d0/105.d0
 beta(4,5)=64.d0/2835.d0
 beta=beta+transpose(beta)
 beta(1,1)=1.d0/12.d0
 beta(2,2)=1.d0/4.d0
 beta(3,3)=4.d0/15.d0
 beta(4,4)=64.d0/945.d0
 beta(5,5)=704.d0/2835.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Kinetic energy (radial and orbital), mass, and 
!!  external potential matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Radial kinetic energy matrix
 AA=0.d0
 do i=1,3
  do j=0,i-1
     AA(4-i,4-j)=alpha(Lambda(i),Lambda(j))/hk(1)
  enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      AA(4*k-i,4*k-j)=alpha(Lambda(i),Lambda(j))/hk(k)
   enddo
  enddo
 enddo 

 do i=2,4
  do j=1,i-1
     AA(4*ninter-i,4*ninter-j)=alpha(Lambda(i),Lambda(j))/hk(ninter)
  enddo
 enddo

 if(Robin==1 .or. Robin==2)then
   do i=1,4
       AA(4*ninter-i,4*ninter)=alpha(Lambda(i),2)/hk(ninter)
   enddo
 endif

 AA=AA+transpose(AA)
 do k=2,ninter
    AA(4*k-4,4*k-4)=alpha(2,2)/hk(k-1)+alpha(1,1)/hk(k)
 enddo

 if(Robin==1 .or. Robin==2)AA(4*ninter,4*ninter)=alpha(2,2)/hk(ninter)!!k=ninter+1
 if(Robin==2)AA(4*ninter,4*ninter)=AA(4*ninter,4*ninter)+1.d0/(3.d0*Length)

  do k=1,ninter
   do i=1,3
    AA(4*k-i,4*k-i)=alpha(Lambda(i),Lambda(i))/hk(k)
   enddo
  enddo

 A=0.d0
   do   j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) A(nband+1+i-j,j)=AA(i,j)
    enddo
   enddo 


!! Orbital kinetic energy matrix
 BB=0.d0 
 do i=1,3
   do j=0,i-1
      BB(4-i,4-j)=integrate('B',Ng,xg,wg,hk(1),hr1,hr2,1,Lambda(i),Lambda(j),0,&
                                    FormMat,tFormMat,ninter)
   enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      BB(4*k-i,4*k-j)=integrate('B',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(j),0,&
                                    FormMat,tFormMat,ninter)
   enddo
  enddo
 enddo
 do i=2,4
  do j=1,i-1
     BB(4*ninter-i,4*ninter-j)=integrate('B',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                        Lambda(i),Lambda(j),0,FormMat,tFormMat,ninter)
  enddo
 enddo

 if(Robin==1 .or. Robin==2)then
  do i=1,4
     BB(4*ninter-i,4*ninter)=integrate('B',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                        Lambda(i),2,0,FormMat,tFormMat,ninter)
  enddo
 endif

 BB=BB+transpose(BB)
 do k=2,ninter
    BB(4*k-4,4*k-4)=integrate('B',Ng,xg,wg,hk(k-1),hr1,hr2,k-1,&
                                2,2,0,FormMat,tFormMat,ninter)+&
                   integrate('B',Ng,xg,wg,hk(k),hr1,hr2,k,1,1,0,&
                                        FormMat,tFormMat,ninter)
 enddo
 if(Robin==1 .or. Robin==2)BB(4*ninter,4*ninter)=integrate('B',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                2,2,0,FormMat,tFormMat,ninter)
 if(Robin==2)BB(4*ninter,4*ninter)=BB(4*ninter,4*ninter)+1.d0/(3.d0*Length)

 do k=1,ninter
  do i=1,3
    BB(4*k-i,4*k-i)=integrate('B',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(i),0,&
                                FormMat,tFormMat,ninter)
  enddo
 enddo

 B=0.d0
   do j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) B(nband+1+i-j,j)=BB(i,j)
    enddo
   enddo 

 kinM=0.d0
 do l=0,Ls
 kinM(l+1,:,:)=AA+dble(l*(l+1))*BB
 enddo

!! Mass matrix
 MM=0.d0
 do i=1,3
  do j=0,i-1
     MM(4-i,4-j)=mu(Lambda(i),Lambda(j))*hk(1)
  enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      MM(4*k-i,4*k-j)=mu(Lambda(i),Lambda(j))*hk(k)
   enddo
  enddo
 enddo 
 do i=2,4
  do j=1,i-1
      MM(4*ninter-i,4*ninter-j)=mu(Lambda(i),Lambda(j))*hk(ninter)
  enddo
 enddo

 if(Robin==1 .or.Robin==2)then
   do i=1,4
        MM(4*ninter-i,4*ninter)=mu(Lambda(i),2)*hk(ninter)
    enddo
 endif

 MM=MM+transpose(MM)
 do  k=2,ninter
    MM(4*k-4,4*k-4)=mu(2,2)*hk(k-1)+mu(1,1)*hk(k)
 enddo
if(Robin==1 .or. Robin==2) MM(4*ninter,4*ninter)=mu(2,2)*hk(ninter)
if(Robin==2) MM(4*ninter,4*ninter)=MM(4*ninter,4*ninter)+Length
 
 do k=1,ninter
  do i=1,3 
    MM(4*k-i,4*k-i)=mu(Lambda(i),Lambda(i))*hk(k)
  enddo
 enddo
 

 M=0.d0
   do j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) M(nband+1+i-j,j)=MM(i,j)
    enddo
   enddo 

! External potential matrix
 CC=0.d0
 do i=1,3
   do j=0,i-1
      CC(4-i,4-j)=integrate('C',Ng,xg,wg,hk(1),hr1,hr2,1,Lambda(i),Lambda(j),0,&
                                FormMat,tFormMat,ninter)
   enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      CC(4*k-i,4*k-j)=integrate('C',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(j),0,&
                                FormMat,tFormMat,ninter)
   enddo
  enddo
 enddo
 do i=2,4
  do j=1,i-1
     CC(4*ninter-i,4*ninter-j)=integrate('C',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                Lambda(i),Lambda(j),0,FormMat,tFormMat,ninter)
  enddo
 enddo
 if(Robin==1 .or. Robin==2)then
 do i=1,4
     CC(4*ninter-i,4*ninter)=integrate('C',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                Lambda(i),2,0,FormMat,tFormMat,ninter)
 enddo
 endif

 CC=CC+transpose(CC)
 do k=2,ninter
    CC(4*k-4,4*k-4)=integrate('C',Ng,xg,wg,hk(k-1),hr1,hr2,k-1,2,2,0,FormMat,tFormMat,ninter)+&
                   integrate('C',Ng,xg,wg,hk(k),hr1,hr2,k,1,1,0,FormMat,tFormMat,ninter)
 enddo
 if(Robin==1 .or. Robin==2)CC(4*ninter,4*ninter)=integrate('C',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                                           2,2,0,FormMat,tFormMat,ninter)
 if(Robin==2)CC(4*ninter,4*ninter)=CC(4*ninter,4*ninter)+0.5d0

 do k=1,ninter
  do i=1,3
     CC(4*k-i,4*k-i)=integrate('C',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(i),0,&
                               FormMat,tFormMat,ninter)
  enddo
 enddo

 C=0.d0
   do j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) C(nband+1+i-j,j)=CC(i,j)
    enddo
   enddo 


!! Radial stark effect energy matrix
 EE=0.d0
 do i=1,3
  do j=0,i-1
     EE(4-i,4-j)=beta(Lambda(i),Lambda(j))*hk(1)**2+&
                hk(1)*rk(1)*mu(Lambda(i),Lambda(j))
  enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      EE(4*k-i,4*k-j)=beta(Lambda(i),Lambda(j))*hk(k)**2+&
                     hk(k)*rk(k)*mu(Lambda(i),Lambda(j))
   enddo
  enddo
 enddo 
 do i=2,4
  do j=1,i-1
     EE(4*ninter-i,4*ninter-j)=beta(Lambda(i),Lambda(j))*hk(ninter)**2+&
                              hk(ninter)*rk(ninter)*mu(Lambda(i),Lambda(j))
  enddo
 enddo
 if(Robin==1)then
 do i=1,4
     EE(4*ninter-i,4*ninter)=beta(Lambda(i),2)*hk(ninter)**2+&
                              hk(ninter)*rk(ninter)*mu(Lambda(i),2)
 enddo
 endif

 EE=EE+transpose(EE)
 do k=2,ninter
    EE(4*k-4,4*k-4)=beta(2,2)*hk(k-1)**2+hk(k-1)*rk(k-1)*mu(2,2)+&
                   beta(1,1)*hk(k)**2+hk(k)*rk(k)*mu(1,1)
 enddo
 if(Robin==1)EE(4*ninter,4*ninter)=beta(2,2)*hk(ninter)**2+hk(ninter)*rk(ninter)*mu(2,2)

 do k=1,ninter
  do i=1,3
    EE(4*k-i,4*k-i)=beta(Lambda(i),Lambda(i))*hk(k)**2+&
                   hk(k)*rk(k)*mu(Lambda(i),Lambda(i))
  enddo
 enddo

 E=0.d0
   do  j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) E(nband+1+i-j,j)=EE(i,j)
    enddo
   enddo 

!! Additional potential energy matrix
 KKC=0.d0
 do i=1,3
  do j=0,i-1
     KKC(4-i,4-j)=integrate('K',Ng,xg,wg,hk(1),hr1,hr2,1,Lambda(i),Lambda(j),eta,&
                            FormMat,tFormMat,ninter)
  enddo
 enddo
 do k=2,ninter-1
  do i=1,4
   do j=0,i-1
      KKC(4*k-i,4*k-j)=integrate('K',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(j),eta,&
                                FormMat,tFormMat,ninter)*exp(-eta*rk(k))
   enddo
  enddo
 enddo
 do i=2,4
  do j=1,i-1
     KKC(4*ninter-i,4*ninter-j)=integrate('K',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                Lambda(i),Lambda(j),eta,FormMat,tFormMat,ninter)*&
                                 exp(-eta*rk(ninter))!!there was a mistake in i and j it is corrected on 24 Aug2015
  enddo
 enddo
 if(Robin==1 .or. Robin==2)then
 do i=1,4
     KKC(4*ninter-i,4*ninter)=integrate('K',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,&
                                Lambda(i),2,eta,FormMat,tFormMat,ninter)*&
                                 exp(-eta*rk(ninter))
 enddo
 endif

 KKC=KKC+transpose(KKC)
 do k=2,ninter
    KKC(4*k-4,4*k-4)=integrate('K',Ng,xg,wg,hk(k-1),hr1,hr2,k-1,2,2,eta,&
                             FormMat,tFormMat,ninter)*exp(-eta*rk(k-1))+&
                     integrate('K',Ng,xg,wg,hk(k),hr1,hr2,k,1,1,eta,&
                              FormMat,tFormMat,ninter)*exp(-eta*rk(k))
 enddo
 if(Robin==1 .or. Robin==2)KKC(4*ninter,4*ninter)=integrate('K',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,2,2,eta,&
                             FormMat,tFormMat,ninter)*exp(-eta*rk(ninter))
 if(Robin==2)then
     temp=0.d0
     do i=1,Ng
        temp=temp+wg(i)*exp(-eta*Length/xg(i))*xg(i)
     enddo
   KKC(4*ninter,4*ninter)=KKC(4*ninter,4*ninter)-temp
 endif

 do k=1,ninter
  do i=1,3
    KKC(4*k-i,4*k-i)=integrate('K',Ng,xg,wg,hk(k),hr1,hr2,k,Lambda(i),Lambda(i),eta,&
                               FormMat,tFormMat,ninter)*exp(-eta*rk(k))
  enddo
 enddo

 KC=0.d0
   do  j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) KC(nband+1+i-j,j)=KKC(i,j)
    enddo
   enddo 

!!!!!The addition vactor G
 G=0.d0
   do k=1,ninter
     if(k>=2) then
      G(4*(k-1))=integrate('G',Ng,xg,wg,hk(k-1),hr1,hr2,k-1,2,0,eta,&
                         FormMat,tFormMat,ninter)*exp(-eta*rk(k-1))+&
                 integrate('G',Ng,xg,wg,hk(k),hr1,hr2,k,1,0,eta,&
                             FormMat,tFormMat,ninter)*exp(-eta*rk(k))
     endif
      do i=1,3
      G(4*(k-1)+i)=integrate('G',Ng,xg,wg,hk(k),hr1,hr2,k,i+2,0,eta,&
                            FormMat,tFormMat,ninter)*exp(-eta*rk(k))
      enddo
   enddo    
      
 if(Robin==1 .or. Robin==2)G(4*ninter)=integrate('G',Ng,xg,wg,hk(ninter),hr1,hr2,ninter,2,0,eta,&
                         FormMat,tFormMat,ninter)*exp(-eta*rk(ninter))
 if(Robin==2)then
      temp=0.d0
    do i=1,Ng
       temp=temp+wg(i)*exp(-eta*Length/xg(i))/xg(i)
    enddo
       temp=temp*Length
    G(4*ninter)=G(4*ninter)+temp
 endif
 !open(unit=1,file='Matrices.dat')
 !  do  j=1,N
 !   do i=1,j
 !      if(max(1,j-nband)<=i) write(1,*)A(nband+1+i-j,j),B(nband+1+i-j,j),C(nband+1+i-j,j),M(nband+1+i-j,j)
 !   enddo
 !  enddo 
 
end subroutine

