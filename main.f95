 SUBROUTINE  main(Zi,Zf,wVxc,beta,Ng,Na,Ls,ms,lrho,Robin,steplength,logic_oda,grid_div,scal0,Length,Length0,ninter,ninter0,&
           epsE,eps,vc,continuation,norm,Rayleih,my_method,pert,diff_NI,step_NI,final_NI,&
           diff_R,step_R,final_R,diff_beta,step_beta,final_beta,threshold,final_eps,final_epsE,exceed,max_iter)
 implicit none
 integer::N,i,j,k,l,INFO,ILAENV,LWORK,mn,l1,l2,nn,lambda,ii,Z,diff
 integer::counter1,acc,counter,acc_pd,acc41
 real:: time1,time2
 double precision::Legendre,formfunction,tformfunction,Dlamch,cef,pi,constant,temp,F
 double precision::Ekin,Ekin0,Etot,Etot0, density,intG(5),cxc,x0,t0,density1,eta,Zp                           
 double precision, allocatable, dimension(:)::rk,hk,xg,wg,xa,wa,hr1,hr2,CGcoef
 double precision, allocatable, dimension(:,:)::A,B,C,M,E,KC,LegMat,FormMat,tFormMat,AA,MM,BB
 double precision, allocatable, dimension(:,:)::rhoF,KKC,CC,thr13,H1_norm
 double precision, allocatable, dimension(:,:,:)::Vxc,rho,rho0,evec,FF1,kinM,evec0
 double precision, allocatable, dimension(:,:,:,:)::FF,constantM
 integer, allocatable, dimension(:,:)::CGcoefint
 character(len=50)::filename,filename1,filename2,filename3,filename4,filename5
 INTEGER,allocatable,dimension(:)::IFAIL, IWORK
 double precision, allocatable, dimension(:):: W, WORK,G,Gtest
 double precision,allocatable,dimension(:,:):: AB, BB1, P,MM1,EE,HH,ev,ev0,occn,occn0,Ql,Hhat,Vxhat
 
 integer,parameter::nband=4,nev=50,nwigner=588 

 integer::Ng,Na,ninter,ninter0,step_NI,final_NI,Zi,Zf,exceed,max_iter
 double precision::scal0,steplength,epsE,eps,vc,Length,Length0,beta
 double precision::step_R,final_R,step_beta,final_beta,final_eps,final_epsE
 integer::logic_oda,wVxc,diff_NI,diff_R,Ls,ms,lrho,diff_beta,continuation
 integer::norm,pert,Robin,threshold,grid_div,Rayleih,my_method

               
 double precision::vs,Res,initial_eps
!!when I turn on the electric field I should pay attention to: Length,ninter,Ls,ms

 call cpu_time(time1)
 allocate(CGcoef(nwigner),CGcoefint(nwigner,4))
  CGcoefint=0
  CGcoef=0.d0
 open(unit=13, file='particularCG.dat', status='old')
  do i=1,nwigner
     read(13,*) CGcoefint(i,1),CGcoefint(i,2),CGcoefint(i,3),CGcoefint(i,4),CGcoef(i)
  enddo
 close(13)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! check the parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 diff=diff_NI+diff_R+diff_beta+threshold
 if(diff.ge.2)then
  print*,'multiple variables program stopped'
  write(3,*) 'multiple variables program stopped'
  stop
 endif
 if(Robin.ge.1.and.(wVxc.ne.0.or.beta.ne.0))then
  print*,'basis should be zero'
  write(3,*)'basis should be zero'
  stop
 endif  
 if(Ls.gt.6)then
  print*,'No data for Ls>6'
  write(3,*)'No data for Ls>6' 
  stop
 endif
 counter=0
 initial_eps=eps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(threshold.eq.1)then 
      open(unit=17,file='threshold.txt')
      write(17,*) '        eps         ','        epsE          ','              n_{0,3}         ',&
                  '            E           ','           e_{0,3)        ','            e_{1,1)   '
 endif

 if(diff_NI.eq.1)then
   if(norm.eq.1)open(unit=15,file='NI-E-norm.txt')
   open(unit=16,file='NI-E-EV.txt')
 endif
 if(diff_R.eq.1)then
   if(norm.eq.1)open(unit=15,file='R-E-norm.txt')
   open(unit=16,file='R-E-EV.txt')
 endif

 do 13 
  counter=counter+1
 if(diff_NI.eq.1)then
    if(abs(ninter-final_NI).le.1.d-10)stop
    if(counter.ge.2)ninter=ninter+step_NI
 endif
 if(diff_R.eq.1)then
    if(abs(Length-final_R).le.1.d-10)stop
    if(counter.ge.2)Length=Length+step_R
 endif
 if(diff.eq.0)then
    if(counter.eq.2)stop
 endif

 if(threshold.eq.1)then
    eps=eps*1.d-1
   if(abs(eps-final_eps).le. 1.d-15)then
      if(abs(epsE-final_epsE).le. 1.d-15)stop
      write(17,*) 
      eps=initial_eps
      epsE=epsE*1.d-1
    endif
    print*, 'eps=',eps, 'epsE=',epsE
    goto 21
 endif
21 continue


 if(Robin.eq.0)then
    N=4*ninter-1  !! number of dof
 elseif(Robin.eq.1.or.Robin.eq.2)then
    N=4*ninter
 endif

 allocate(rk(ninter),hk(ninter),hr1(ninter),hr2(ninter),xg(Ng),wg(Ng),xa(Na),wa(Na))
 allocate(FormMat(5,Ng),tFormMat(4,Ng),LegMat(2*Ls+1,Na))
 allocate(FF(5,5,5,ninter-1),FF1(4,4,4),constantM(Ls+1,Ls+1,2*Ls+1,ms+1))
 allocate(AB(nband+1,N), BB1(nband+1,N),G(N),Gtest(N))
 allocate(A(nband+1,N),B(nband+1,N),C(nband+1,N),M(nband+1,N),E(nband+1,N),KC(nband+1,N))
 allocate(BB(N,N),CC(N,N),EE(N,N),AA(N,N),MM(N,N),KKC(N,N))
 allocate(Vxc(2*Ls+1,nband+1,N),kinM(Ls+1,N,N),thr13(Ng,ninter))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! constructing the Grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(grid_div.eq.0)then
    if(scal0.eq.1.d0)then
      hk(ninter)=Length/dble(ninter)
    else 
      hk(ninter)=Length*(1.d0-scal0)/(1.d0-scal0**ninter)
    endif

    do j=1,ninter-1
       hk(ninter-j)=scal0*hk(ninter-j+1)
    enddo
 elseif(grid_div.eq.1)then
    hk(ninter0)=Length0*(1.d0-scal0)/(1.d0-scal0**ninter0)
    do j=1,ninter0-1
       hk(ninter0-j)=scal0*hk(ninter0-j+1)
    enddo
      hk(ninter)=(Length-Length0)/dble(ninter-ninter0)
    do j=ninter0+1,ninter
       hk(j)=hk(ninter)
    enddo
 endif!grid_div=0

    rk(1)=0.d0
    do k=1,ninter-1
       rk(k+1)=rk(k)+hk(k)
    enddo

     hr1(1)=0.d0
     hr2(1)=0.d0
    do k=2,ninter
       hr1(k)=hk(k)/rk(k)
       hr2(k)=hk(k)/rk(k)**2
    enddo
 
 open(unit=4,file='grid.txt')
 do i=1,ninter
 write(4,*) hk(i),rk(i),hr1(i),hr2(i)
 enddo
 write(4,*)  '    ',rk(ninter)+hk(ninter)
 close(4)

 call gauleg(0.d0,1.d0,xg,wg,Ng)
 call gauleg(-1.d0,1.d0,xa,wa,Na)

 pi=4.d0*atan(1.d0)
 cxc=2.d0*pi*(3.d0/pi)**(1.d0/3.d0)
    !2.d0*pi*4.d0/3.d0! this changed recently I looked at the note given by keiron

 do k=1,ninter
   do i=1,Ng
      thr13(i,k)= (xg(i)*hk(k)+rk(k))**(1.d0/3.d0) 
   enddo
 enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrices used for integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 LegMat=0.d0
 do l=0,2*Ls
   do i=1,Na
      LegMat(l+1,i)=sqrt(2.d0*dble(l)+1.d0)*Legendre(l,xa(i))
   enddo
 enddo
  LegMat=1.d0/sqrt(4.d0*pi)*LegMat

 FormMat=0.d0
  do i=1,5
   do k=1,Ng
      FormMat(i,k)=formfunction(i,xg(k))
   enddo
  enddo

 tFormMat=0.d0
  do i=2,5
   do k=1,Ng
      tFormMat(i-1,k)=tformfunction(i,xg(k))
   enddo
  enddo
   

  constantM=0.d0
  do l1=0,Ls
    do l2=0,Ls
      do l=0,2*Ls
        do mn=0,ms
           if(l1.ge.l2) then
              constantM(l1+1,l2+1,l+1,mn+1)=constant(l1,l2,l,mn,pi,CGcoef,CGcoefint,nwigner)
           else 
              constantM(l1+1,l2+1,l+1,mn+1)=constant(l2,l1,l,mn,pi,CGcoef,CGcoefint,nwigner)
            endif
              constantM(l1+1,l2+1,l+1,mn+1)=constantM(l1+1,l2+1,l+1,mn+1)*(-1)**(mod(mn,2))
         enddo
      enddo
    enddo
  enddo

    open(unit=5,file='test-for-matrices.txt')
      Gtest(3)=0.d0
      Gtest(4)=0.d0
        do l=0,Ls
           Gtest(3)=Gtest(3)+sqrt(dble(2*l+1)/(4.d0*pi))
         enddo
        do mn=0,ms!! changed on 18sep2016  it was 0,3
         temp=0.d0
         do l1=mn,Ls
           do l2=mn,Ls
              do l=abs(l1-l2),l1+l2
                  temp=temp+constantM(l1+1,l2+1,l+1,mn+1)*sqrt(dble(2*l+1)/(4.d0*pi))
              enddo
            enddo
          enddo
          if(mn.eq.0)then
           write(5,*)mn, temp,'it should be',Gtest(3)**2 ,'test for the Wigner 3j symbol'
          else
           write(5,*)mn, temp,'it should be', Gtest(4),'test for the Wigner 3j symbol'
          endif
          enddo

  FF=0.d0
  do i=0,4
    do j=0,4
      do nn=0,4
         do k=2,ninter
            FF(i+1,j+1,nn+1,k-1)=F(Lambda(i),Lambda(j),Lambda(nn),hk(k),rk(k),wg,xg,Ng,FormMat,tFormMat)
         enddo
      enddo
    enddo
  enddo
 
  FF1=0.d0
  do i=0,3
    do j=0,3
      do nn=0,3
          FF1(i+1,j+1,nn+1)=F(Lambda(i),Lambda(j),Lambda(nn),hk(1),rk(1),wg,xg,Ng,FormMat,tFormMat)
      enddo
    enddo
  enddo

 !!!I will do iteration for different Z at the same time 
 do 18 Z=Zi,Zf
 
 acc=0
 acc_pd=0
 acc41=0
 if(wVxc.eq.1 .and. ((Z.ge. 23 .and. Z.le. 28).or.(Z.ge. 41 .and. Z.le. 44)))acc=1
 if(wVxc.eq.0 .and. ((Z.ge. 23 .and. Z.le. 26).or.(Z.ge. 46 .and. Z.le. 47)))acc=1
 if(wVxc.eq.0 .and. (Z.eq.21 .or. Z.eq.22 .or. Z.eq.40))acc_pd=1
 if(wVxc.eq.0 .and. (Z.eq.41 .or. Z.eq.42))acc41=1!
  Zp=dble(Z)!+1.d0
 allocate(occn(ms+1,Z),evec(ms+1,(Ls+1)*N,Z),ev(Z,ms+1),ev0(Z,ms+1))
 allocate(H1_norm(ms+1,int(dble(Z)/2.d0)+1),occn0(ms+1,Z),evec0(1,1,1))
 evec0=0.d0
 if(norm.eq.1)then
  deallocate(evec0)
  allocate(evec0(ms+1,(Ls+1)*N,Z))
  evec0=0.d0
 endif
 occn0=0.d0

!!Find suitable eta
 eta=0.5d0
 do 22
  if(exp(-eta*Length).le.epsE)goto 23
    eta=eta+0.1d0
22 enddo 
23 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Assembling the Matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call assemblingmatrices(N,ninter,nband,hk,rk,Ng,xg,wg,eta,A,B,C,M,E,KC,G,&
                         CC,EE,MM,Ls,kinM,KKC,FormMat,tFormMat,BB,AA,hr1,hr2,&
                         Robin,Length)


   G=dble(Z)*eta**2*G 

  if(Robin.eq.1)then
    do l=0,Ls
       kinM(l+1,N,N)=kinM(l+1,N,N)+1.d0/Length
    enddo
  endif


 !call test(N,nband,A,Ng,xg,wg,ninter,hk,rk,Robin,Length,FormMat)
 ! test for the basis if they can capture a solution which behaves like 1/r at infinity

!  do i=1,N
!  do j=1,N
!     if(MM(i,j).ne.MM(j,i)) print*, 'M is not symetric see entry',i,j
!     if(EE(i,j).ne.EE(j,i)) print*, 'EE is not symetric see entry',i,j
!  enddo
! enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Test for the vector 4\pi G involving the matrices M,B,C,A and E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 intG(1)= 0.5d0*dble(Z)**2*eta**3*&
          (1.d0-exp(-2.d0*eta*Length))!is the integral (4\pi G)^2
 intG(2)=0.5d0*dble(Z)**2*eta**2*&
         (0.5d0-(Length*eta+0.5d0)*exp(-2.d0*eta*Length))!is the integral r(4\pi G)^2
 intG(3)=0.5d0*dble(Z)**2* &
          (0.5d0*eta-(Length**2*eta**3+Length*eta**2+0.5d0*eta)&
          *exp(-2.d0*eta*Length)) !is the integral  of r^2(4\pi G)^2
 intG(5)=0.25d0*dble(Z)**2/eta* &
          (3.d0-(2.d0*Length**4*eta**4+4.d0*Length**3*eta**3+&
                  6.d0*Length**2*eta**2+6.d0*Length*eta+3.d0)&
          *exp(-2.d0*eta*Length))!is the integral  of r^4(4\pi G)^2
 intG(4)=0.5d0*dble(Z)**2* &
          (0.75d0-(Length**3*eta**3+&
                  3.d0/2.d0*Length**2*eta**2+3.d0/2.d0*Length*eta+3.d0/4.d0)&
          *exp(-2.d0*eta*Length))!is the integral  of r^3(4\pi G)^2

 AB=C
 Gtest=G
 call DPBSV('U', N,nband, 1, AB, nband+1,Gtest, N, INFO ) 
 if(info.ne.0) write(5,*) info,'it should be zero'  
 !Test for the Matrix KC
 temp=0.d0
 do i=1,N
  do j=max(1,i-nband),min(N,i+nband)
     temp=temp+Gtest(i)*KKC(i,j)*Gtest(j)
  enddo
 enddo
 write(5,*)    temp, 'it should be' , dble(Z)**2/3.d0*&
           (eta**2/3.d0-(Length*eta**3+eta**2/3.d0)*exp(-3.d0*eta*Length)),&
           'is a test for the K matrix... test for G & C & K'
 !! end up with accuracy up to order 10^-11

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*MM(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be',intG(3),&
          'test for G & C & M'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*BB(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be',intG(1),&
          'test for G & C & B'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*CC(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be', intG(2),&
          'test for G & C'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*EE(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be', intG(4),&
          'test for G & E'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*AA(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp,'it should be', intG(1)-2.d0*eta*intG(2)+&
           eta**2*intG(3),'test for G & A' 
          !'is the integral of(1-eta*r)^2(4\pi G)^2'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*(AA(i,j)-BB(i,j)-MM(i,j)*eta**2+&
              2.d0*CC(i,j)*eta)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be', 'zero' ,'test for G & A & B & C & M'
 
 AB=B
 Gtest=G
 call DPBSV('U', N,nband, 1, AB, nband+1,Gtest, N, INFO ) 
 if(info.ne.0) write(5,*) info,'it should be zero'  
 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*BB(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be',intG(3),'test for G & B'
 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*MM(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be',intG(5),'test for G & B & M'

 temp=0.d0
 do i=1,N
    do j=max(1,i-nband),min(N,j+nband)
        temp=temp+Gtest(i)*CC(i,j)*Gtest(j)
    enddo
 enddo
 write(5,*)  temp, 'it should be',intG(4),'test for G & B & C'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The mass matrix is coercive 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      AB=M
      BB1=0.d0
      do i=1,N
         BB1(nband+1,i)=1.d0
      enddo
      ALLOCATE(W(N),Ql(1,N),P(1,N),WORK(7*N),IFAIL(nev), IWORK(5*N))
      W=0.d0
      call DSBGVX('N','I','U', N,nband, nband, AB,nband+1, BB1,&        
              nband+1, Ql,1,0.d0,0.d0,1,nev,2*DLAMCH('S'),i, W, P,&
               N, WORK, IWORK, IFAIL, INFO )
      !!DSBGVX: computes eigen values of A*x=(lambda)*B*x where A and B are symmetic banded matrices
      write(5,*)'The smallest eigenvalue of M is', W(1)
      deallocate(W,P,WORK,IFAIL,IWORK,Ql)
 close(5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Self-consistent iteration 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(pert.eq.1)goto 50

 allocate(rho(lrho+1,N,N),rho0(lrho+1,N,N))
 allocate(rhoF(lrho+1,N),Ql(lrho+1,N))

 counter1=0
 do 24 
   counter1=counter1+1
 if(diff_beta.eq.1)then
    if(abs(beta-final_beta).le.1.d-15)stop
    if(counter1.ge.2)beta=beta+step_beta
 else
    if(counter1.eq.2)goto 51
 endif


 if(beta.ne.0)then
  write(filename3,'(f11.9,a,i1,a,i1,a)')beta,"_results-oda=",logic_oda,"-xc=",wVxc,".txt"
  write(filename4,'(f11.9,a)')beta,"_E-and-n.txt" 
  write(filename5,'(f11.9,a)')beta,"_evec.txt"
 else
  write(filename3,'(a,i1,a,i1,a)')"results-oda=",logic_oda,"-xc=",wVxc,".txt"
  write(filename4,'(a)')"E-and-n.txt" 
  write(filename5,'(i2,a)')Z,"_evec.txt"
 endif
 open(unit=3,file=filename3)
 open(unit=10,file=filename4)

 write(3,*) 'Length of the interval is',Length,'number of grid points is',ninter
 write(3,*)'number of spherical harmonics used',Ls
 if(norm.eq.1)write(3,*) 'the H1 norm is calculated'
 if(Robin.eq.1) write(3,*) '###Robin boundary condition is considered'
 if(Robin.eq.2) write(3,*) '###infinite base is considered'
 if(grid_div.eq.1) write(3,'(a,f4.1,a1,f4.1,a,i2,a1,i2)') '###2 part grid Le=',&
                   Length0,'+',Length-Length0,'  NI=',ninter0,'+',ninter-ninter0

 print*,'beta=',beta
 write(3,*) 'electric field modulos is', beta
 print*,'eta=',eta,'Length=',Length,'Ninter=',ninter
 write(3,*)'eta=',eta
 write(3,*) 'exp(-eta*Length)=', exp(-eta*Length)

 cef=2.d0*sqrt(pi/3.d0)*beta

 if(continuation.eq.1.and.counter1.ne.1)then
 Etot0=Etot
 Ekin0=Ekin
 rho0=rho
 occn0=occn
  goto 25
 endif
 Ql=0.d0
 rho=0.d0
 rho0=0.d0
 Ekin0=0.d0
 Ekin=0.d0
 Etot=0.d0
 Etot0=0.d0
 t0=1.d0
 Vxc=0.d0
25 continue

 ii=0
4  continue
 if(exceed.eq.1.and. ii.gt.max_iter)then
  print*,'the iteration exceeds',max_iter,', the program stopped'
  write(3,*)'the iteration exceeds',max_iter,', the program stopped'
 endif
    ii=ii+1 
   ev=0.d0
   evec=0.d0
 do 3 mn=0,ms
  LWORK= ILAENV( 1, 'DSYTRD', 'U', (Ls-mn+1)*N, -1, -1, -1 )
  LWORK=(LWORK+3)*(Ls-mn+1)*N
  allocate(HH((Ls-mn+1)*N,(Ls-mn+1)*N),MM1((Ls-mn+1)*N,(Ls-mn+1)*N), WORK(LWORK))
  allocate(W((Ls-mn+1)*N), P((Ls-mn+1)*N,nev),IFAIL((Ls-mn+1)*N), IWORK(5*(Ls-mn+1)*N))
  allocate(Hhat((Ls-mn+1)*N,(Ls-mn+1)*N),Vxhat((Ls-mn+1)*N,(Ls-mn+1)*N))
  HH=0.d0
  MM1=0.d0
  Hhat=0.d0
  Vxhat=0.d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Computing the Hartree and the exchange matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  if(maxval(Ql).eq.0.and. minval(Ql).eq.0) goto 5
  call Hartreepotential(N,mn,lrho,ninter,Ql,Ls,ms,FF,FF1,constantM,Hhat,nband,&
                        Robin)
 if(wVxc.ne.0)then
 call exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho,&
                       lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,wVxc,cxc)
  Vxc=-cxc*Vxc

  do l2=mn,Ls
   do l1=mn,l2
     AB=0.d0
      do l=0,2*Ls 
       if(constantM(l1+1,l2+1,l+1,mn+1).ne.0) AB=AB+constantM(l1+1,l2+1,l+1,mn+1)*Vxc(l+1,:,:)
       enddo
         call construct(nband,N,AB,AA)
     do i=1,N                                
      do j=1,N 
         Vxhat((l1-mn)*N+i,(l2-mn)*N+j)=AA(i,j)
      enddo
     enddo
   enddo
  enddo
 endif
5 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Diagonalization of the Hamiltonian matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!I will just construct the upper part of \cH  and Hhat
  do 2 l=mn,Ls
     if(l.ne.Ls)temp=-cef*constantM(l+2,l+1,2,mn+1)
     if(maxval(abs(rho)).eq.0) then!ii.eq.1
      !temp=0.d0
      AA=0.5d0*kinM(l+1,:,:)-Zp*CC
     else
      AA=0.5d0*kinM(l+1,:,:)-dble(Z)*KKC+(dble(Z)-dble(Zp))*CC
     endif  
     do i=1,N                                
      do j=max(i-nband,1),min(N,i+nband) 
         HH((l-mn)*N+i,(l-mn)*N+j)=AA(i,j)  !! diagonal of \cH  
         MM1((l-mn)*N+i,(l-mn)*N+j)=MM(i,j)
         if(l.ne.Ls)HH((l-mn)*N+i,(l-mn+1)*N+j)=temp*EE(i,j)  !! off diagonal of \cH 
      enddo
     enddo
2 enddo

      HH=HH+Hhat+Vxhat 
      call DSYGVX( 1,'V','I', 'U', (Ls-mn+1)*N, HH, (Ls-mn+1)*N, MM1, (Ls-mn+1)*N,&
                   0.d0,0.d0,1,nev, 2*DLAMCH('S'),i,W,P,(Ls-mn+1)*N, WORK,&
                   LWORK, IWORK, IFAIL, INFO)
     !!DSYGVX: computes eigen values of A*x=(lambda)*B*x where A and B are symmetic matrices
     !print*,'for m= ',k,W(1),W(2),W(3),W(4),W(5),N,(Ls-mn+1)*N,Info!,W(6),i,Info,N!,WORK(1)!!contains the value of LWORK
   if(info.ne. 0) write(3,*) 'warning info is equal to ', info
   
   !!!Saving the eigenvalues and eigenvectors for each m
   do k=1,Z
   ev(k,mn+1)=W(k)   
   do j=1,(Ls-mn+1)*N
   evec(mn+1,mn*N+j,k)=P(j,k)
   enddo
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!arranging the eigenvectors in the first iteration 
!!because they are degenerate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 if(beta.eq.0 .and. Z.le.36)then
  if(ii.eq.1)then
   if(Z.ge. 1 .and. Z.le. 2)then
      goto 100 
   elseif(Z.ge.3.and.Z.le.10)then
      if(mn.eq.0) then
       evec(0+1,:,2)=P(:,3)
       evec(0+1,:,3)=P(:,2)
      endif
   elseif(Z.ge.11.and.Z.le.36)then
      if(mn.eq.0) then
       evec(0+1,:,2)=P(:,3)
       evec(0+1,:,3)=P(:,2)
        do k=4,6  !!n=3
          if(P(1,k).ne.0)then
             evec(0+1,:,4)=P(:,k)
          elseif(P(1,k).eq.0.and.P(N+1,k).ne.0)then
             evec(0+1,:,5)=P(:,k)
          elseif(P(1,k).eq.0.and.P(N+1,k).eq.0.and.P(2*N+1,k).ne.0)then
             evec(0+1,:,6)=P(:,k)
          else
           print*,'could not find the orbitals at first iteration'
          endif
        enddo
          if(Z.ge. 19)then
           do k=7,10  !!n=4
             if(P(1,k).ne.0)then
                evec(0+1,:,7)=P(:,k)
             elseif(P(1,k).eq.0.and.P(N+1,k).ne.0)then
                evec(0+1,:,8)=P(:,k)
             elseif(P(1,k).eq.0.and.P(N+1,k).eq.0.and.P(2*N+1,k).ne.0)then
                evec(0+1,:,9)=P(:,k)
             else
                print*,'could not find the orbitals at first iteration!!',k
             endif
            enddo
          endif!Z>=19
      endif!mn=0
        if(Z.ge.13.and. Z.le.30)then !!this is working automatically
          if(mn.eq.1) then
             do j=1,(Ls-mn+1)*N
                evec(1+1,mn*N+j,2)=P(j,3)
                evec(1+1,mn*N+j,3)=P(j,2)
             enddo       
          endif
         endif
        if(Z.ge. 31.and. mn.eq.1)then
             do j=1,(Ls-mn+1)*N
                evec(1+1,mn*N+j,2)=P(j,3)
                evec(1+1,mn*N+j,3)=P(j,2)
                evec(1+1,mn*N+j,4)=P(j,5)
                evec(1+1,mn*N+j,5)=P(j,4)
             enddo           
        endif
goto 100
   elseif(Z.eq.29) then 
    if(mn.eq.0) then
     evec(0+1,:,2)=P(:,3)
     evec(0+1,:,3)=P(:,2)
     evec(0+1,:,4)=P(:,6)
     evec(0+1,:,5)=P(:,5)
     evec(0+1,:,6)=P(:,4)
     evec(0+1,:,7)=P(:,9)         
    endif
    if(mn.eq.1) then
         do j=1,(Ls-mn+1)*N
            evec(1+1,mn*N+j,2)=P(j,3)
            evec(1+1,mn*N+j,3)=P(j,2)
         enddo       
    endif
   else
    print*,'there is something wrong in the main program'
   endif
100 continue
 endif
 endif! end of beta.eq.0 and Z<=36
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! writing the eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  deallocate(HH,MM1,WORK,W,P,IFAIL,IWORK,Hhat,Vxhat)
3 enddo
  write(3,*)' iteration number',ii
  print*,' iteration number',ii
 do mn=0,ms
   write(3,*) mn, (ev(i,mn+1),i=1,5)
   if(Z.ge.13)write(3,*) mn, (ev(i,mn+1),i=6,10)
   if(Z.ge.40 .and. mn.eq.0)write(3,*) mn, (ev(i,mn+1),i=11,13)
 enddo
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! occupation number 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 if(my_method.eq.1.and.(Z.eq.6.or. z.eq.5))then!beta.ne.0.and.
  if(my_method.eq.1 .and.maxval(abs(rho)).ne.0)then!.and.ii.ne.1 !.and.&(abs(ev(3,1)-ev(1,2)).le.1.d-3)
      if(lrho.ne.2*Ls)then
         print*, 'lrho should be', 2*Ls! for the minimization procedure to get the right quantity 
         stop
      endif
      call occupationC(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,rho,evec,occn,&
                FF,FF1,A,B,G,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,Robin,Length,constantM,hk,M,rk) 
      goto 20
   endif
     occn=0.d0
     occn(1,1)=2.d0
     occn(1,2)=2.d0
     if(ev(3,1).le.ev(1,2))then
       occn(1,3)=(dble(Z)-4.d0)
       occn(2,1)=0.d0
     else
       occn(1,3)=0.d0
       occn(2,1)=(dble(Z)-4.d0)
     endif

 elseif(beta.ne.0.and.Z.eq.4)then
    occn=0.d0
    occn(1,1)=2.d0
    if(ev(2,1)<ev(1,2))then
      occn(1,2)=2.d0
    elseif(ev(2,1)>ev(1,2))then
      occn(2,1)=2.d0
      print*, 'I entered here'
    else
      print*,'the problem is stopped'
      stop 
    endif 
 elseif(Z.le. 36)then
    call occupationnb(ii,ev,Z,ms,occn)
 elseif(Z.ge. 36)then!Z>36
    call occupationnbplus(ii,ev,Z,ms,occn)
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Computation of the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 rho=0.d0

 if(acc.eq.1.and. ((wVxc.eq.1 .and. dabs(ev(6,1)-ev(7,1)).le. vc).or. &!
                 (wVxc.eq.0 .and. dabs(ev(7,1)-ev(8,1)).le. vc .and. ii.ne.1)).and. Z.le. 36)then
      call accidental(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,rho,occn,&
                      FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                      tFormMat,LegMat,thr13,cxc,Zp,CC,constantM,Robin,Length)
 elseif(acc.eq.1.and. ((wVxc.eq.1 .and. dabs(ev(9,1)-ev(10,1)).le. vc).or. &!
                 (wVxc.eq.0 .and. dabs(ev(9,1)-ev(10,1)).le. vc .and. ii.ne.1)).and. Z>36)then
       call accidentalplus(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,rho,occn,&
                           FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                           tFormMat,LegMat,thr13,cxc,Zp,CC,constantM,Robin,Length)
 elseif(acc_pd.eq.1 .and. ((Z.le.22 .and.dabs(ev(7,1)-ev(8,1)).le. vc).or.&
                        (Z.ge.40 .and.dabs(ev(10,1)-ev(11,1)).le. vc) ).and. ii.ne.1)then
        call accidental_pd(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,rho,occn,&
                           FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                           tFormMat,LegMat,thr13,cxc,Zp,CC,constantM,Robin,Length)
 elseif(acc41.eq.1 .and.dabs(ev(10,1)-ev(11,1)).le. vc)then!!deg bt 6s and 4d
        call accidentalplus(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,rho,occn,&
                            FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                            tFormMat,LegMat,thr13,cxc,Zp,CC,constantM,Robin,Length)
 else
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,ms    !!will be changed later
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
20 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Computing the integral of the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 do l=0,lrho
  temp=0.d0
  do i=1,N
    do j=max(1,i-nband),min(i+nband,N)
       temp=temp+rho(l+1,i,j)*MM(i,j)
    enddo
  enddo
  write(3,*)'the integral of 2 sqrt(pi)r^2\rho_',l,' is',  2.d0*sqrt(pi)*temp
  if(l.ne.0.and. abs(temp)>eps .and. (beta.eq.0.and.my_method.eq.1)) then
   print*,beta, 'the ',l,' component of the density is not zero'
   !stop!!
  endif
 enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Mixing the densities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 if(logic_oda.eq.0)then !!no oda 
   t0=steplength
  if(minval(rho0).ne.0.or.maxval(rho0).ne.0) rho=(1.d0-t0)*rho+t0*rho0
 else         !!with oda
 if(ii.eq.1)then!maxval(abs(rho0)).eq.0
  t0=1.d0
 else
  call oda(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,rho,rho0,evec,occn,&
          Ekin0,FF,FF1,A,B,G,t0,KKC,EE,kinM,wVxc,Ng,Na,wg,xg,wa,hr1,FormMat,&
                tFormMat,LegMat,thr13,cxc,Zp,CC,0,Robin,Length,hk,M,rk)
  endif
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Solving the poisson equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho,Ql,rhoF,Robin,Length)!,hk,M,rk) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the energy 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call energy(nband,ms,Z,evec,N,Ls,nev,kinM,occn,rhoF,&
             Ql,lrho,KKC,pi,rho,AA,Vxc,Ekin,Etot,CC,EE,cef,Ekin0,t0,Zp,vs)
  print*, '                                                                                   total energy=',Etot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Stopping criteria
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 H1_norm=0.d0
 if(norm.eq.1)then
!!!!!!H_1 norm of the eigen functions
  do  mn=0,ms   
      do k=1,Z
          if(occn(mn+1,k).ne.0)then  
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+(evec(mn+1,l*N+i,k))*MM(i,j)*(evec0(mn+1,l*N+j,k))
               enddo
             enddo
           enddo

          if(temp.ge.0)then
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+(evec(mn+1,l*N+i,k)-evec0(mn+1,l*N+i,k))*(MM(i,j)+kinM(l+1,i,j))*&
                              (evec(mn+1,l*N+j,k)-evec0(mn+1,l*N+j,k))
               enddo
             enddo
           enddo
          elseif(temp<0)then
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+(evec(mn+1,l*N+i,k)+evec0(mn+1,l*N+i,k))*(MM(i,j)+kinM(l+1,i,j))*&
                              (evec(mn+1,l*N+j,k)+evec0(mn+1,l*N+j,k))
               enddo
             enddo
           enddo
          endif
          H1_norm(mn+1,k)=sqrt(abs(temp))
          endif
        enddo
    enddo
 endif!norm.eq.1 
   

   if(maxval(occn0).ne.0)occn=(1.d0-t0)*occn0+t0*occn

   if(my_method.eq.1.and.(z.eq.6 .or.z.eq.5))then!beta.ne.0 .and.
       print*, 'on m=0',occn(1,3),'on m=1',occn(2,1)
       write(3,*)'on m=0',occn(1,3),'on m=1',occn(2,1)
   endif

 if(abs(sum(occn)-Z)>eps) then
  print*,'error in the occupation number see main program'
  write(3,*)'error in the occupation number see main program'
  stop
 endif

 if(logic_oda.eq.0)then !!without oda
  if(maxval(H1_norm).le.eps .and. abs(Etot-Etot0).le. epsE)goto 12
 elseif(logic_oda.eq.1)then !!with oda
  if(Etot>Etot0)print*,'------------------------------------! energy increases!!'
  print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',abs(Etot-Etot0)
   if(abs(Etot-Etot0).le. epsE .and.maxval(H1_norm).le.eps.and.maxval(abs(ev-ev0)).le.vc)then
      print*,'energy converged now calculating the residual'
      call rayleighqoutient(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,occn,&
                FF,FF1,A,B,G,KKC,EE,kinM,constantM,ev,Robin,Length,Res,cxc,&
                Ng,xg,wg,Na,wa,LegMat,FormMat,tFormMat,hr1,thr13,wVxc)
      if(my_method.eq.1)then
        if(Res.le.1.d-5)goto 12
      else
        goto 12
      endif
   endif
 endif

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!! saving data for next interation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Etot0=Etot
 Ekin0=Ekin
 rho0=rho
 occn0=occn
 ev0=ev!!!this one just used in the stopping criterion
 if(norm.eq.1)evec0=evec
 write(3,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 goto 4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Saving converged Data in files  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
12 continue
 Gtest=0.d0
 if(beta.ne.0)then
   write(filename1,'(f11.9,a,i2.2,a,i1,a,i1,a)')beta,"_densityr2-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
   write(filename2,'(f11.9,a,i2.2,a,i1,a,i1,a)')beta,"_densityr2car-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
   write(filename,'(f11.9,a,i2.2,a,i1,a,i1,a)')beta,"_1Ddensity-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
 else
   write(filename1,'(a,i2.2,a,i1,a,i1,a)')"densityr2-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
   write(filename2,'(a,i2.2,a,i1,a,i1,a)')"densityr2car-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
   write(filename,'(a,i2.2,a,i1,a,i1,a)')"1Ddensity-Z=",Z,"-oda=",logic_oda,"-xc=",wVxc,".txt"
 endif

 open(unit=7,file=filename1)
 open(unit=9,file=filename2)
 open(unit=14,file=filename)!'density.txt')
       t0=1.d-4
 do k=1,ninter  !int(ninter/2) if I wanna do less disk
  do i=1,Ng
       x0=-1.d0
14     continue
       temp=density1(i,k,x0,rho,lrho,N,pi,FormMat,Ng,ninter) 
       write(7,*) xg(i)*hk(k)+rk(k),Acos(x0), temp
       write(9,*) (xg(i)*hk(k)+rk(k))*sin(Acos(x0)),(xg(i)*hk(k)+rk(k))*x0, temp
        x0=x0+t0
        if(x0.ge.xa(1))then
         goto 15
        else
          goto 14
        endif
15 continue
  
   do j=1,Na  
      temp=density(i,k,j,rho,lrho,N,&
                    LegMat,Ls,Na,FormMat,Ng,ninter)
       write(7,*) xg(i)*hk(k)+rk(k),Acos(xa(j)), temp
       write(9,*) (xg(i)*hk(k)+rk(k))*sin(Acos(xa(j))),(xg(i)*hk(k)+rk(k))*xa(j), temp
   enddo 
   !!!! saving the radial density for  fixed xa(Na) and any r
   if(k.eq.1)then
      Gtest(1)=0.d0
         do l=0,lrho
            do ii=0,3
               do j=0,3
                   Gtest(1)=Gtest(1)+rho(l+1,4-ii,4-j)*tFormMat(Lambda(ii)-1,i)*&
                                  tFormMat(Lambda(j)-1,i)*LegMat(l+1,Na)
               enddo
             enddo
         enddo
      Gtest(1)=Gtest(1)/hk(k)**2!!j---> Na
       else
        Gtest(1)=temp/(xg(i)*hk(k)+rk(k))**2
    endif

   !!!! saving the radial density for  fixed xa(1) and any r
   if(k.eq.1)then
      Gtest(2)=0.d0
         do l=0,lrho
            do ii=0,3
               do j=0,3
                   Gtest(2)=Gtest(2)+rho(l+1,4-ii,4-j)*tFormMat(Lambda(ii)-1,i)*&
                                  tFormMat(Lambda(j)-1,i)*LegMat(l+1,1)
               enddo
             enddo
         enddo
      Gtest(2)=Gtest(2)/hk(k)**2!!j---> Na
       else
        Gtest(2)=density(i,k,1,rho,lrho,N,LegMat,Ls,Na,FormMat,Ng,ninter)
        Gtest(2)=Gtest(2)/(xg(i)*hk(k)+rk(k))**2
    endif
     write(14,*) xg(i)*hk(k)+rk(k), Gtest(1), Gtest(2), temp

       x0=1.d0
16     continue
       temp=density1(i,k,x0,rho,lrho,N,pi,FormMat,Ng,ninter) 
       write(7,*) xg(i)*hk(k)+rk(k),Acos(x0), temp
       write(9,*) (xg(i)*hk(k)+rk(k))*sin(Acos(x0)),(xg(i)*hk(k)+rk(k))*x0, temp
        x0=x0-t0
        if(x0.le.xa(Na))then
          goto 17
        else
          goto 16
        endif
17 continue
    write(7,*) 
    write(9,*) 
  enddo
 enddo
 close(7)
 close(9)
 close(14)

 if(beta.eq.0 .and.((Z.ge. 5 .and.Z.le.10).or.(Z.eq.11 .or. Z.eq.12).or.&
                       (Z.ge. 23 .and. Z.le. 26).or.(Z.ge.41.and.Z.le.42)))  then
 call eigenfunctions(evec,ms,Ls,N,Z,Ng,xg,hk,rk,ninter,FormMat,Robin)
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Calculating the Rayleigh qoutient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(Zp.eq.Z)then
   write(10,*) 'eigenvalues for atom', Z
 else
   write(10,*) 'eigenvalues for ion with nucleus charge', Zp, 'and', Z, 'electons' 
 endif
 write(10,*)  '   m=0                     ','   m=1                     ','   m=2                       '
 do j=1,min(Z,10)
 write(10,*) (ev(j,mn+1),mn=0,ms)
 enddo
 
 if(Rayleih.eq.1 .and. Zp.eq.Z)then
      call rayleighqoutient(lrho,ms,Ls,ninter,nband,N,Z,pi,cef,evec,occn,&
                FF,FF1,A,B,G,KKC,EE,kinM,constantM,ev,Robin,Length,Res,cxc,&
                Ng,xg,wg,Na,wa,LegMat,FormMat,tFormMat,hr1,thr13,wVxc)
 write(10,*) 'eigenvalues calculated by Rayleigh qoutient' 
 write(10,*)  '   m=0                     ','   m=1                     ','   m=2                       '
 do j=1,min(Z,10)
 write(10,*) (ev(j,mn+1),mn=0,ms)
 enddo
 endif

 write(10,*) 'occupation numbers'
 write(10,*)  '   m=0                     ','   m=1                     ','   m=2                       '
 do j=1,min(Z,10)
 write(10,*)  (occn(mn+1,j), mn=0,ms)
 enddo

 write(10,*) 'Energy=',Etot

 if(diff_R.eq.1 .or. diff_NI.eq.1)then
!!!!!!H_1 norm of the eigen functions 
 H1_norm=0.d0
 if(norm.eq.1)then
  do  mn=0,ms   
      do k=1,Z
         if(occn(mn+1,k).ne.0)then     
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+evec(mn+1,l*N+i,k)*(MM(i,j)+kinM(l+1,i,j))*evec(mn+1,l*N+j,k)
               enddo
             enddo
           enddo
          H1_norm(mn+1,k)=sqrt(temp)
         endif
        enddo
    enddo
  endif

 if(Z<=54)then  
    i=0
   do k=1,Z
     if(occn(0+1,k).ne.0)then
        i=i+1
     else
        goto 19 
     endif
    enddo
19 continue
   if(diff_R.eq.1)then
      if(norm.eq.1)write(15,*)ninter,Length,Etot, (H1_norm(0+1,k),k=1,i)
      if(acc.eq.1.or.acc_pd.eq.1.or.acc41.eq.1)then
            if(Z<36)nn=1
            if(Z>36)nn=2
            write(16,*)Length,Etot, (ev(j,0+1),j=1,min(i+3,Z)),'elec. on d-shell:',&
                                    occn(2+1,nn)*5.d0/2.d0,'virial sum=',vs         
      else
           write(16,*)Length,Etot, (ev(j,0+1),j=1,min(i+3,Z)),'virial sum=',vs
      endif
   endif
   if(diff_NI.eq.1)then
      if(norm.eq.1)write(15,*)ninter,Etot, (H1_norm(0+1,k),k=1,i)
      if(beta.ne.0)then
            write(16,*)ninter,Etot, (ev(j,0+1),j=1,4),'ev on m=1',ev(1,2),'elec. on m=0:',&
                                     occn(0+1,3)         
      elseif(acc.eq.1.or.acc_pd.eq.1.or.acc41.eq.1)then
            if(Z<36)nn=1
            if(Z>36)nn=2
            write(16,*)ninter,Etot, (ev(j,0+1),j=1,min(i+3,Z)),'elec. on d-shell:',&
                                     occn(2+1,nn)*5.d0/2.d0,'virial sum=',vs
      else
            write(16,*)ninter,Etot, (ev(j,0+1),j=1,min(i+3,Z)),'virial sum=',vs
      endif
   endif
  else
    print*,'need to be defined see see program'
    stop
 endif
endif!!diff_NI.eq.1 .or. diff_R.eq.1


 if(Z.eq.6 .or. (Z.ge. 23 .and. Z.le. 26) .or. (Z.ge.41 .and. Z.le.42))then
 open(unit=11,file=filename5)
 if(Z.eq.6)then
   do i=1,(Ls+1)*N
    write(11,*) evec(0+1,i,1),evec(0+1,i,2),evec(0+1,i,3), evec(1+1,i,1)
   enddo
 endif
 if(wVxc.eq.0 .and. (Z.ge. 23 .and. Z.le. 26))then
   do i=1,(Ls+1)*N
    write(11,*) evec(0+1,i,6),evec(0+1,i,7), evec(0+1,i,8)!4s 5s 3d
   enddo
 endif
 if(wVxc.eq.0 .and. (Z.ge. 41 .and. Z.le. 42))then
   do i=1,(Ls+1)*N
    write(11,*) evec(0+1,i,9),evec(0+1,i,10), evec(0+1,i,11)!4s 5s 3d
   enddo
 endif
 close(11) 
 endif

 if(threshold.eq.1)then
   write(17,*) eps, epsE, occn(1,3),Etot, ev(3,1),ev(1,2)
 endif

24 enddo
51 continue
 deallocate(H1_norm)
 deallocate(occn,occn0,evec,ev,ev0)
 deallocate(rho,rho0,rhoF,Ql)
 if(norm.eq.1)deallocate(evec0)
18 enddo
50 continue

 if(pert.eq.1) then 
 call perturbation(nband,N,ms,Z,Ls,lrho,ninter,evec,occn,constantM,FF,FF1,pi,G,A,B,KinM,KKC,MM,EE,ev,beta,&
                   wVxc,Robin,Length)
 endif

 deallocate(rk,hk,hr1,hr2,xg,wg,xa,wa)
 deallocate(FormMat,tFormMat,LegMat)
 deallocate(FF,FF1,constantM)
 deallocate(AB, BB1,G,Gtest)
 deallocate(A,B,C,M,E,KC)
 deallocate(BB,CC,EE,AA,MM,KKC)
 deallocate(Vxc,kinM,thr13)

 call cpu_time(time2)
 print*, 'Elapsed CPU time = ', (time2-time1)/60.d0,'min'
 write(3,*)'Elapsed CPU time = ', (time2-time1)/60.d0,'min'
13 enddo
 
 close(3)
 close(10)!results,E-and-n
 if((diff_NI.eq.1.or.diff_R.eq.1).and.norm.eq.1)close(15)
 if(diff_NI.eq.1.or.diff_R.eq.1)close(16)
 if(threshold.eq.1)close(17) 

End SUBROUTINE
