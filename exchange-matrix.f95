 Subroutine exchangematrix(N,nband,ninter,Ng,xg,wg,Na,wa,rho,&
                           lrho,Ls,LegMat,FormMat,tFormMat,Vxc,hr1,thr13,wVxc)
 Implicit none
 integer::l,i,j,k,N,nband,ninter,Ng,lrho,Na,Lambda,Ls,p,q,wVxc
 double precision:: integrateVxc,density
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::Matrix(N,N),FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision, dimension(:,:,:)::rho(lrho+1,N,N),Vxc(2*Ls+1,nband+1,N),DensityM(Ng,ninter,Na)

 Vxc=0.d0

!! Density Matrix
if(wVxc==0 .or. wVxc==1)then
 DensityM=0.d0
  do p=1,Ng
    do k=1,ninter
      do q=1,Na
        DensityM(p,k,q)=density(p,k,q,rho,lrho,N,LegMat,Ls,Na,FormMat,Ng,ninter)
        if(DensityM(p,k,q)<0) then
           if(abs(DensityM(p,k,q))>10E-11)then
            print*, 'the density at point (k,p,q)',k,p,q,'is not positive !!',DensityM(p,k,q)
            stop
            endif
           DensityM(p,k,q)=0.d0
         else
         DensityM(p,k,q)=(DensityM(p,k,q))**(1.d0/3.d0)
         endif
      enddo
     enddo
   enddo
 elseif(wVxc==2)then !!this case is done just for the oda
 DensityM=0.d0
  do p=1,Ng
    do k=1,ninter
      do q=1,Na
        DensityM(p,k,q)=density(p,k,q,rho,lrho,N,LegMat,Ls,Na,FormMat,Ng,ninter)
        DensityM(p,k,q)=abs(DensityM(p,k,q))**(1.d0/3.d0)
      enddo
     enddo
   enddo
 endif
 if(wVxc==0)return

 !do k=1,ninter
 !print*,k, DensityM(1,k,1),DensityM(Ng,k,1)
 !enddo
 !stop
!! Exchange correlation energy Matrix
  do 3 l=0,2*Ls
   Matrix=0.d0
   do i=1,3
    do j=0,i-1
       Matrix(4-i,4-j)=integrateVxc(Na,wa,Ng,xg,wg,1,&
                      Lambda(i),Lambda(j),l,LegMat,Ls,&
                      FormMat,tFormMat,ninter,thr13,hr1,DensityM)
    enddo
   enddo
   do k=2,ninter-1
    do i=1,4
     do j=0,i-1
        Matrix(4*k-i,4*k-j)=integrateVxc(Na,wa,Ng,xg,wg,k,&
                           Lambda(i),Lambda(j),l,LegMat,Ls,&
                           FormMat,tFormMat,ninter,thr13,hr1,DensityM)
     enddo
    enddo
   enddo
   do i=2,4
    do j=1,i-1
       Matrix(4*ninter-i,4*ninter-j)=integrateVxc(Na,wa,Ng,xg,wg,&
                                    ninter,Lambda(i),Lambda(j),l,LegMat,Ls,&
                                    FormMat,tFormMat,ninter,thr13,hr1,DensityM)
    enddo
   enddo

   do k=2,ninter
      Matrix(4*k-4,4*k-4)=integrateVxc(Na,wa,Ng,xg,wg,&
                          k-1,2,2,l,LegMat,Ls,&
                          FormMat,tFormMat,ninter,thr13,hr1,DensityM)+&
                          integrateVxc(Na,wa,Ng,xg,wg,&
                          k,1,1,l,LegMat,Ls,&
                          FormMat,tFormMat,ninter,thr13,hr1,DensityM)
   enddo
   do k=1,ninter
    do i=1,3
      Matrix(4*k-i,4*k-i)=integrateVxc(Na,wa,Ng,xg,wg,k,&
                          Lambda(i),Lambda(i),l,LegMat,Ls,&
                          FormMat,tFormMat,ninter,thr13,hr1,DensityM)
    enddo
   enddo


   do  j=1,N
    do i=1,j
       if(max(1,j-nband)<=i) Vxc(l+1,nband+1+i-j,j)=Matrix(i,j)
    enddo
   enddo
3 enddo
 end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function integrateVxc(Na,wa,Ng,xg,wg,k,i,j,l,LegMat,Ls,&
                       FormMat,tFormMat,ninter,thr13,hr1,DensityM)
 Implicit none
 integer::i,j,l,p,q,Na,Ng,Ls,k,ninter
 double precision::integrateVxc
 double precision, dimension(:)::wa(Na),wg(Ng),xg(Ng),hr1(ninter)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),FormMat(5,Ng),tFormMat(4,Ng),thr13(Ng,ninter)
 double precision, dimension(:,:,:)::DensityM(Ng,ninter,Na)

 integrateVxc=0.d0
 if(k==1) then
 do p=1,Ng
  do q=1,Na
     integrateVxc=integrateVxc+wg(p)*wa(q)*DensityM(p,k,q)*&
                  LegMat(l+1,q)*tFormMat(i-1,p)*&
                  FormMat(j,p)*thr13(p,1)
  enddo
 enddo
 else
 do p=1,Ng
  do q=1,Na
     integrateVxc=integrateVxc+wg(p)*wa(q)*DensityM(p,k,q)*&
                  LegMat(l+1,q)*FormMat(i,p)*&
                  FormMat(j,p)*thr13(p,k)/(xg(p)*hr1(k)+1.d0)
  enddo
 enddo
 integrateVxc=hr1(k)*integrateVxc
 endif
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function density(p,k,q,rho,lrho,N,LegMat,Ls,Na,FormMat,Ng,ninter)
 Implicit none
 integer::l,i,j,k,Lambda,lrho,Ls,Na,q,p,Ng,ninter,N
 double precision:: density
 double precision,dimension(:,:)::LegMat(2*Ls+1,Na),FormMat(5,Ng)
 double precision,dimension(:,:,:)::rho(lrho+1,N,N) 
 density=0.d0
 if(k==1) then
   do l=0,lrho
    do i=0,3
     do j=0,3
        density=density+rho(l+1,4-i,4-j)*FormMat(Lambda(i),p)*&
                 FormMat(Lambda(j),p)*LegMat(l+1,q)
     enddo
    enddo
   enddo
  elseif(k==ninter)then
   do l=0,lrho
    do i=1,4
     do j=1,4
        density=density+rho(l+1,4*ninter-i,4*ninter-j)*FormMat(Lambda(i),p)*&
                 FormMat(Lambda(j),p)*LegMat(l+1,q)
     enddo
    enddo
   enddo
  else
   !if(p.ne. 1 .and. p.ne. Ng) then
     do l=0,lrho
      do i=0,4
       do j=0,4
          density=density+rho(l+1,4*k-i,4*k-j)*FormMat(Lambda(i),p)*&
                   FormMat(Lambda(j),p)*LegMat(l+1,q)
       enddo
      enddo
     enddo 
    !elseif(p==1) then
    ! do l=0,lrho
    !      density=density+rho(l+1,4*k-4,4*k-4)*4.d0*LegMat(l+1,q)
    ! enddo
    !elseif(p==Ng) then
    ! do l=0,lrho
    !      density=density+rho(l+1,4*k,4*k)*4.d0*LegMat(l+1,q)
    ! enddo
    !endif   
  endif  
 end function   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function density1(p,k,q,rho,lrho,N,pi,FormMat,Ng,ninter) 
!! this function is buit just for plotting and the number q is passed now as a double precision poit not as an index
 Implicit none
 integer::l,i,j,k,Lambda,lrho,p,Ng,ninter,N
 double precision:: density,pi,q,density1,legendre
 double precision,dimension(:,:)::FormMat(5,Ng)
 double precision,dimension(:,:,:)::rho(lrho+1,N,N) 

 density=0.d0
 if(k==1) then
   do l=0,lrho
    do i=0,3
     do j=0,3
        density=density+rho(l+1,4-i,4-j)*FormMat(Lambda(i),p)*&
                 FormMat(Lambda(j),p)*sqrt(2.d0*dble(l)+1.d0)*Legendre(l,q)
     enddo
    enddo
   enddo
  elseif(k==ninter)then
   do l=0,lrho
    do i=1,4
     do j=1,4
        density=density+rho(l+1,4*ninter-i,4*ninter-j)*FormMat(Lambda(i),p)*&
                 FormMat(Lambda(j),p)*sqrt(2.d0*dble(l)+1.d0)*Legendre(l,q)
     enddo
    enddo
   enddo
  else
     do l=0,lrho
      do i=0,4
       do j=0,4
          density=density+rho(l+1,4*k-i,4*k-j)*FormMat(Lambda(i),p)*&
                   FormMat(Lambda(j),p)*sqrt(2.d0*dble(l)+1.d0)*Legendre(l,q)
       enddo
      enddo
     enddo  
  endif  
 density1=density/sqrt(4.d0*pi)
 end function 
