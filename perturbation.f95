 subroutine perturbation(nband,N,mmax,Z,Ls,lrho,ninter,evec,occn,constantM,FF,&
                          FF1,pi,G,A,B,KinM,KKC,MM,EE,ev,beta,wVxc,Robin,Length)
 Implicit none
 integer:: nband,N,lrho,ninter,mmax,Z,Ls,wVxc,Robin,Length
 integer::l,i,j,mn,k,l1,l2
 double precision:: pi,G(N),cef,KKC(N,N),AA(N,N),MM(N,N),temp,EE(N,N),beta,G1(N),norm(4)
 double precision::A(nband+1,N),B(nband+1,N),ev(Z,mmax+1),ev0(Z,mmax+1),ev1(Z,mmax+1),tempv(4)
 double precision::occn(mmax+1,Z),occn0(mmax+1,Z),occn1(mmax+1,Z)
 double precision, dimension(lrho+1,N)::rhoF,Ql,Ql001,Ql100,Qlpert
 double precision, dimension(lrho+1,N,N)::rho,rho001,rho100,rhopert
 double precision::evec(mmax+1,(Ls+1)*N,Z),evec0(mmax+1,(Ls+1)*N,Z),FF1(4,4,4),kinM(Ls+1,N,N),evec1(mmax+1,(Ls+1)*N,Z)
 double precision, allocatable, dimension(:,:)::Hhat,Hhat001,Hhat100,Hhatpert
 double precision::constantM(Ls+1,Ls+1,2*Ls+1,mmax+1),FF(5,5,5,ninter-1)
 double precision::test(4),E(3)
 character(len=60)::energy,path
 character(len=100)::filename1,filename2,filename3,filename4

 beta=0.d0
 write(filename1,'(a,f11.9,a)')'electric-field-Z=06/xc=0/my_method/oda=0/',beta,"_E-and-n.txt" 
 write(filename2,'(a,f11.9,a)')'electric-field-Z=06/xc=0/my_method/oda=0/',beta,"_evec.txt"

 print*, 'entered the perturbation file'

 occn0=0.d0
 ev0=0.d0
 open(unit=9,file=filename1,status='old') !!!!!this should be removed later
 read (9,*)
 read (9,*)
 read(9,*)ev0(1,1),ev0(1,2)
 read(9,*)ev0(2,1)
 read(9,*)ev0(3,1)
 if(z.eq.6)then
 read(9,*)
 read(9,*)
 endif
 read(9,*)
 read(9,*)
 read(9,*)
 read(9,*)occn0(1,1),occn0(2,1)
 read(9,*)occn0(1,2)
 read(9,*)occn0(1,3)
 read(9,*)
 if(z.eq.6)then
 read(9,*)
 read(9,*)
 endif
 read(9,*)energy,E(1)
 close(9)

 evec0=0.d0
 open(unit=12, file=filename2, status='old')
      do i=1,(Ls+1)*N
          if(Z==6)then
            read(12,*) evec0(1,i,1), evec0(1,i,2), evec0(1,i,3), evec0(2,i,1)  !!evec0(mn+1,mn*N+j,k)
          elseif(Z==4)then
            read(12,*) evec0(1,i,1), evec0(1,i,2)
          endif
      enddo
 close(12)

  open(unit=15, file='first-order-perturbation/beta-E-pert.txt')
  open(unit=16, file='first-order-perturbation/beta-ev-pert.txt')
  open(unit=17, file='first-order-perturbation/beta-pert.txt')
  open(unit=13, file='first-order-perturbation/pert.txt')

 cef=2.d0*sqrt(pi/3.d0)!*beta/beta !! this cef is without the beta

 beta=1.d-6
 do 7
  if(beta.ge.1.d-4)then
   if(beta.ge.1.d-4)stop
   beta=beta+1.d-4
   goto 8
  endif
  if(beta.ge.1.d-5)then
   beta=beta+1.d-5
   goto 8
  endif
   beta=beta+1.d-6
8 continue
   print*,'--------------------------------------------------------------------------',beta

    if((beta.ge. 3.d-4 .and.beta.le. 8.d-3).or.(beta.ge. 1.4d-2 .and.beta.le. 2.7d-2))then
      write(path,'(a)')'/normal_oda/'
    else
      write(path,'(a)')'/my_method/oda=0/'
    endif

 write(filename3,'(a,i1,a,f11.9,a)')'electric-field-Z=06/xc=',wVxc,trim(path),beta,"_E-and-n.txt" 
 write(filename4,'(a,i1,a,f11.9,a)')'electric-field-Z=06/xc=',wVxc,trim(path),beta,"_evec.txt"
 occn=0.d0
 ev=0.d0
 open(unit=10,file=filename3,status='old',action='read') 
      read(10,*)
      read(10,*)
      read(10,*) ev(1,1),ev(1,2)
      read(10,*) ev(2,1)
      read(10,*) ev(3,1)
      if(Z==6)then
        read(10,*)
        read(10,*)
      endif
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) occn(1,1),occn(2,1)
      read(10,*) occn(1,2)
      read(10,*) occn(1,3)
      if(Z==6)then
        read(10,*)
        read(10,*)
      endif
      read(10,*)
      read(10,*)energy,E(2)
 close(10)

 evec=0.d0
 open(unit=11,file=filename4,status='old',action='read') 
     do i=1,(Ls+1)*N
        if(Z==6)then
          read(11,*) evec(1,i,1),evec(1,i,2),evec(1,i,3), evec(2,i,1)
        elseif(Z==4)then
          read(11,*) evec(1,i,1),evec(1,i,2)
        endif
     enddo
 close(11)


 occn1=(occn-occn0)/beta
 ev1=(ev-ev0)/beta
 E(3)=(E(2)-E(1))/beta


 if(Z==6)then
 print*, 'the occupation number wo pert', occn0(1,1),occn0(1,2),occn0(1,3),occn0(2,1)
 print*, 'the occupation number w  pert', occn(1,1),occn(1,2),occn(1,3),occn(2,1)
 print*, 'the occupation number deriv  ', occn1(1,1),occn1(1,2),occn1(1,3),occn1(2,1)
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 print*, 'the eigenvalues wo pert',ev0(1,1),ev0(2,1),ev0(3,1),ev0(1,2)
 print*, 'the eigenvalues w  pert',ev(1,1),ev(2,1),ev(3,1),ev(1,2)
 print*, 'the eigenvalues  deriv ',ev1(1,1),ev1(2,1),ev1(3,1),ev1(1,2)
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 print*, 'sum of occupation number', sum(occn0), sum(occn), sum(occn1)
 elseif(Z==4)then
 print*, 'the occupation number wo pert', occn0(1,1),occn0(1,2)
 print*, 'the occupation number w  pert', occn(1,1),occn(1,2)
 print*, 'the occupation number deriv  ', occn1(1,1),occn1(1,2)
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 print*, 'the eigenvalues wo pert',ev0(1,1),ev0(2,1)
 print*, 'the eigenvalues w  pert',ev(1,1),ev(2,1)
 print*, 'the eigenvalues  deriv ',ev1(1,1),ev1(2,1)
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 print*, 'sum of occupation number', sum(occn0), sum(occn), sum(occn1)
 endif
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


 do mn=0,mmax
    do k=1,Z
       if(occn0(mn+1,k).ne.0)then
          temp=0.d0
           do l=mn,Ls
             do i=1,N
               do j=max(i-nband,1),min(i+nband,N)
                    temp=temp+(evec(mn+1,l*N+i,k))*MM(i,j)*(evec0(mn+1,l*N+j,k))
               enddo
             enddo
           enddo
           print*,mn,k,temp
          if(temp<0)evec(mn+1,:,k)=-evec(mn+1,:,k)
       endif
    enddo
 enddo
 evec1=(evec-evec0)/beta

!!!!!! test for the obtained eigenvectors
 do mn=0,mmax
  do k=1,Z
    if(occn0(mn+1,k).ne.0) then
   temp=0.d0
   tempv=0.d0
   do l=mn,Ls
      do i=1,N
        do j=max(i-nband,1),min(N,i+nband)
           tempv(1)=tempv(1)+evec1(mn+1,l*N+i,k)*MM(i,j)*evec1(mn+1,l*N+j,k)
           tempv(2)=tempv(2)+evec(mn+1,l*N+i,k)*MM(i,j)*evec0(mn+1,l*N+j,k)
           tempv(3)=tempv(3)+evec0(mn+1,l*N+i,k)*MM(i,j)*(evec(mn+1,l*N+j,k)-evec0(mn+1,l*N+j,k))
           tempv(4)=tempv(4)+evec0(mn+1,l*N+i,k)*MM(i,j)*evec1(mn+1,l*N+j,k)
        enddo
      enddo
   enddo
   !print*,'mn=',mn,'the',k, 'eigenvectors norm', sqrt(tempv(1)),sqrt(temp), sqrt(tempv(2)) !
   !print*,'mn=',mn,'the',k, 'scalar product', tempv(3)!, tempv(4),(temp+tempv(1)-2.d0*tempv(4))/beta**2
    print*,'<phi_1,phi_1>,<phi_0,phi>,<phi_0,phi-phi_0>,<phi_0,phi_1> ',(tempv(l1),l1=1,4)
    if(mn==0)norm(k)=tempv(4)
    if(mn==1)norm(4)=tempv(4)
   endif
 enddo
 enddo
 print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Computation of the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 rho001=0.d0
 rho100=0.d0
 rho=0.d0
 rhopert=0.d0
 do l=0,lrho
  do i=1,N
   do j=1,N
     do mn=0,mmax    !!will be changed later
       do k=1,Z
         if(occn0(mn+1,k).ne.0) then
          do l1=mn,Ls
            do l2=max(mn,abs(l-l1)),min(Ls,l+l1)!mn,Ls!
            rho001(l+1,i,j)=rho001(l+1,i,j)+occn0(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                            evec0(mn+1,l1*N+i,k)*evec1(mn+1,l2*N+j,k)  
            rho(l+1,i,j)=rho(l+1,i,j)+occn0(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                            evec0(mn+1,l1*N+i,k)*evec0(mn+1,l2*N+j,k)  
            rhopert(l+1,i,j)=rhopert(l+1,i,j)+occn(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                            evec(mn+1,l1*N+i,k)*evec(mn+1,l2*N+j,k)
            enddo
          enddo 
          endif
         if(occn1(mn+1,k).ne.0) then
          do l1=mn,Ls
            do l2=max(mn,abs(l-l1)),min(Ls,l+l1)!mn,Ls!
            rho100(l+1,i,j)=rho100(l+1,i,j)+occn1(mn+1,k)*constantM(l1+1,l2+1,l+1,mn+1)*& 
                            evec0(mn+1,l1*N+i,k)*evec0(mn+1,l2*N+j,k)  
          enddo
        enddo 
          endif
      enddo
    enddo
   enddo
  enddo
 enddo 
!!!! test for the rhos
  do l=0,lrho
  tempv=0.d0
  do i=1,N
    do j=max(1,i-nband),min(i+nband,N)
       tempv(1)=tempv(1)+rho001(l+1,i,j)*MM(i,j)
       tempv(2)=tempv(2)+rho(l+1,i,j)*MM(i,j)
       tempv(3)=tempv(3)+rho100(l+1,i,j)*MM(i,j)
       tempv(4)=tempv(4)+rhopert(l+1,i,j)*MM(i,j)
    enddo
  enddo
   if(Z==6)then
     print*,'the integral of 2 sqrt(pi)r^2\rho_',l,' is',  (2.d0*sqrt(pi)*tempv(i), i=1,4)
   elseif(Z==4)then
     print*,'the integral of 2 sqrt(pi)r^2\rho_',l,' is',  (2.d0*sqrt(pi)*tempv(i), i=1,2),2.d0*sqrt(pi)*tempv(4)
   endif
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Solving the poisson equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rhopert,Qlpert,rhoF,Robin,Length)
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho,Ql,rhoF,Robin,Length)
  G1=0.d0
  call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G1,A,B,rho001,Ql001,rhoF,Robin,Length)
  if(Z==6)call  poisson(lrho,ninter,nband,N,FF,FF1,pi,G1,A,B,rho100,Ql100,rhoF,Robin,Length)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! To be labeled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tempv(4)=0.d0
  test=0.d0
  mn=0
1 continue
  if(mn>mmax) goto 6
  allocate(Hhat((Ls-mn+1)*N,(Ls-mn+1)*N),Hhat001((Ls-mn+1)*N,(Ls-mn+1)*N),Hhat100((Ls-mn+1)*N,(Ls-mn+1)*N))
  allocate(Hhatpert((Ls-mn+1)*N,(Ls-mn+1)*N))
  do 5 k=1,Z
  if(occn0(mn+1,k).ne.0)then

  call Hartreepotential(N,mn,lrho,ninter,Qlpert,Ls,mmax,FF,FF1,constantM,Hhatpert,nband)
  call Hartreepotential(N,mn,lrho,ninter,Ql,Ls,mmax,FF,FF1,constantM,Hhat,nband)
  call Hartreepotential(N,mn,lrho,ninter,Ql001,Ls,mmax,FF,FF1,constantM,Hhat001,nband)
  call Hartreepotential(N,mn,lrho,ninter,Ql100,Ls,mmax,FF,FF1,constantM,Hhat100,nband)
  
  do 3 l=mn,Ls
      AA=0.5d0*kinM(l+1,:,:)-dble(Z)*KKC  
    do 4 i=1,N
      
    tempv(1)=0.d0
    test=0.d0
    do l1=mn,Ls     
      if(l1.le.l)then                           
      do j=max(1,i-nband),min(N,i+nband) 
         tempv(1)=tempv(1)+Hhat((l1-mn)*N+i,(l-mn)*N+j)*evec1(mn+1,(l1)*N+j,k)+&
                  2.d0*Hhat001((l1-mn)*N+i,(l-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)+&
                  Hhat100((l1-mn)*N+i,(l-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)
         test(1)=test(1)+Hhat((l1-mn)*N+i,(l-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)
         test(2)=test(2)+Hhatpert((l1-mn)*N+i,(l-mn)*N+j)*evec(mn+1,(l1)*N+j,k)        
      enddo
       elseif(l1>l)then
      do j=max(1,i-nband),min(N,i+nband) 
         tempv(1)=tempv(1)+Hhat((l-mn)*N+i,(l1-mn)*N+j)*evec1(mn+1,(l1)*N+j,k)+&
                  2.d0*Hhat001((l-mn)*N+i,(l1-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)+&
                  Hhat100((l-mn)*N+i,(l1-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)
         test(1)=test(1)+Hhat((l-mn)*N+i,(l1-mn)*N+j)*evec0(mn+1,(l1)*N+j,k)
         test(2)=test(2)+Hhatpert((l-mn)*N+i,(l1-mn)*N+j)*evec(mn+1,(l1)*N+j,k)        
      enddo
       endif
    enddo

  do l1=1,N
  do j=1,N
     if(AA(l1,j).ne.AA(j,l1)) print*, 'AA is not symetric see entry',l1,j
  enddo
 enddo    
      
      tempv(2)=0.d0                              
      do j=max(i-nband,1),min(N,i+nband) 
          tempv(2)=tempv(2)+AA(i,j)*evec1(mn+1,(l)*N+j,k)-&
                   ev1(k,mn+1)*MM(i,j)*evec0(mn+1,(l)*N+j,k)-&!!!! pay attaention to this!!!!!!!!!!!!!!!!!!
                   ev0(k,mn+1)*MM(i,j)*evec1(mn+1,(l)*N+j,k)
          test(1)=test(1)+AA(i,j)*evec0(mn+1,(l)*N+j,k)-ev0(k,mn+1)*MM(i,j)*evec0(mn+1,(l)*N+j,k)
          test(2)=test(2)+AA(i,j)*evec(mn+1,(l)*N+j,k)-ev(k,mn+1)*MM(i,j)*evec(mn+1,(l)*N+j,k)
      enddo
          
      tempv(3)=0.d0
      do j=max(i-nband,1),min(N,i+nband) 
         if(l.ne.Ls)tempv(3)=tempv(3)-cef*constantM(l+2,l+1,2,mn+1)*EE(i,j)*evec0(mn+1,(l+1)*N+j,k)
         if(l.ne.mn)tempv(3)=tempv(3)-cef*constantM(l+1,l,2,mn+1)*EE(i,j)*evec0(mn+1,(l-1)*N+j,k)
         if(l.ne.Ls)test(2)=test(2)-cef*constantM(l+2,l+1,2,mn+1)*EE(i,j)*evec(mn+1,(l+1)*N+j,k)*beta
         if(l.ne.mn)test(2)=test(2)-cef*constantM(l+1,l,2,mn+1)*EE(i,j)*evec(mn+1,(l-1)*N+j,k)*beta
      enddo

      write(13,*)tempv(3),(tempv(1)+tempv(2)+tempv(3))! mn,k,l,i,tempv(1),tempv(2),tempv(3)
      tempv(4)=max(tempv(4),abs(tempv(1)+tempv(2)+tempv(3)))
       test(3)=max(test(3),abs(test(1)))
       test(4)=max(test(4),abs(test(2))) 
4 enddo
3 enddo
   endif
5 enddo
  mn=mn+1
  deallocate(Hhat,Hhat001,Hhat100,Hhatpert)
  goto 1
6 continue
   print*,'norm of the vector generated from the eqn satisfied by phi_1 ,phi_0 ,phi respectively ',tempv(4),(test(i),i=3,4)
        write(15,*)beta,E(1),E(2),E(2)-E(1),E(3),E(3)/beta,&
                   (occn0(1,k),occn(1,k),occn(1,k)-occn0(1,k),occn1(1,k),occn1(1,k)/beta,k=1,3)
        write(16,*)beta,(ev0(k,1),ev(k,1),ev(k,1)-ev0(k,1),ev1(k,1),ev1(k,1)/beta,k=1,3),&
                         ev0(1,2),ev(1,2),ev(1,2)-ev0(1,2),ev1(1,2),ev1(1,2)/beta
        write(17,*)beta,tempv(4)
7 enddo
end subroutine
