  Subroutine poisson(lrho,ninter,nband,N,FF,FF1,pi,G,A,B,rho,Ql,rhoF,Robin,Length)!,hk,M,rk)!,F0)
  Implicit none
  integer::lrho,N,nn,k,rn,i,j,ninter,nband,info,l,Robin,kappa
  double precision:: pi,rHS(N),Ql(lrho+1,N),FF1(4,4,4),G(N),rho(lrho+1,N,N),Length
!  double precision::F0(N),M(nband+1,N),hk(ninter),rk(ninter),temp
  double precision, dimension(:,:)::A(nband+1,N),B(nband+1,N),AB(nband+1,N),rhoF(lrho+1,N)
  double precision, dimension(:,:,:,:)::FF(5,5,5,ninter-1)

 kappa=1
 if(Robin==1 .or. Robin==2)kappa=0

 rhoF=0.d0
 Ql=0.d0
 do 1 l=0,lrho                       
  do nn=1,N
     rHS(nn)=0.d0
     rn=4-mod(nn,4)
     k=1+(nn-mod(nn,4))/4
     !print*, rn,k
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if(k==1) then
       do i=0,3
         do j=0,3
             rHS(nn)=rHS(nn)+FF1(i+1,j+1,rn+1)*rho(l+1,4-i,4-j)
         enddo
       enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     elseif(k==ninter)then
      if(rn.ne. 4)then
       do i=kappa,4
         do j=kappa,4
             rHS(nn)=rHS(nn)+FF(i+1,j+1,rn+1,ninter-1)*&
                     rho(l+1,4*ninter-i,4*ninter-j)
         enddo
       enddo  
      elseif(rn==4) then    
       do i=0,4
         do j=0,4
             rHS(nn)=rHS(nn)+FF(1,i+1,j+1,ninter-2)*&
                     rho(l+1,4*(ninter-1)-i,4*(ninter-1)-j)
         enddo
        enddo
        do i=kappa,4
           do j=kappa,4
              rHS(nn)=rHS(nn)+FF(5,i+1,j+1,ninter-1)*&
                      rho(l+1,4*ninter-i,4*ninter-j)
           enddo
        enddo
     endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     elseif(2.le.k .and. k.le.ninter-1)then
      if(rn.ne. 4)then
       do i=0,4
         do j=0,4
             rHS(nn)=rHS(nn)+FF(i+1,j+1,rn+1,k-1)*rho(l+1,4*k-i,4*k-j)
         enddo
       enddo       
      elseif(rn==4 .and. k.ne. 2) then

       do i=0,4
         do j=0,4
             rHS(nn)=rHS(nn)+FF(1,i+1,j+1,k-2)*rho(l+1,4*k-4-i,4*k-4-j)+&
                     FF(5,i+1,j+1,k-1)*rho(l+1,4*k-i,4*k-j)
         enddo
       enddo
      elseif(rn==4 .and. k==2) then      !! for k=2 I have a special situation with k-1 
       do i=0,3
         do j=0,3
             rHS(nn)=rHS(nn)+FF1(1,i+1,j+1)*rho(l+1,4-i,4-j)
         enddo
        enddo
        do i=0,4
           do j=0,4
              rHS(nn)=rHS(nn)+FF(5,i+1,j+1,1)*rho(l+1,8-i,8-j)
           enddo
        enddo
      endif
     elseif((Robin==1.or.Robin==2) .and. k==ninter+1 .and. rn==4)then
       do i=kappa,4
         do j=kappa,4
             rHS(nn)=rHS(nn)+FF(1,i+1,j+1,ninter-1)*&
                     rho(l+1,4*ninter-i,4*ninter-j)
         enddo
       enddo        
       if(Robin==2)rHS(nn)=rHS(nn)+rho(l+1,4*ninter,4*ninter)/3.d0
     else
       print*,'!!!!!!!!see poisson equation file, there is a missing case, k=',k,'q=',rn,'Robin',Robin
     endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  enddo
!!!!! Test of the right hand side in the poisson equation where [F_0]=\int r\Xi_i dr
! F0=0.d0
! F0(1)=(hk(1)**2)*1.d0/3.d0
! F0(2)=(hk(1)**2)*4.d0/45.d0
! F0(3)=(hk(1)**2)*4.d0/15.d0
!   do k=2,ninter
!      F0(4*(k-1))=(hk(k-1)**2)/3.d0+hk(k-1)*rk(k-1)*0.5d0+&
!                 (hk(k)**2)/6.d0+hk(k)*rk(k)*0.5d0
!      
!      F0(4*(k-1)+1)=(hk(k)**2)*1.d0/3.d0+hk(k)*rk(k)*2.d0/3.d0
!      F0(4*(k-1)+2)=(hk(k)**2)*4.d0/45.d0+hk(k)*rk(k)*16.d0/45.d0
!      F0(4*(k-1)+3)=(hk(k)**2)*4.d0/15.d0+hk(k)*rk(k)*16.d0/45.d0
!   enddo 
! !if(Robin==1)F0(4*ninter)=(hk(ninter)**2)/3.d0+hk(ninter)*rk(ninter)*0.5d0
! AB=M
! call DPBSV('U', N,nband,1, AB, nband+1,F0, N, INFO )
! !! after this step [F_0]=M^{-1}[F_0]
!
! temp=0.d0
! do i=1,N
!    temp=temp+F0(i)*rHS(i)
! enddo
! if(l==0) then
! print*,'----------------------------------------------------------------------------'
! print*, 'the integral of 2.d0*sqrt(pi)*rF_',l, 'is',2.d0*sqrt(pi)*temp, 'it should be',&
!                                                       2.d0*sqrt(pi)*6.d0/(sqrt(4.d0*pi))
! else
! print*, 'the integral of 2.d0*sqrt(pi)*rF_',l ,'is',2.d0*sqrt(pi)*temp
! endif
! !stop
! ! the above test insures that if I got correct rho then the F is correct


  rhoF(l+1,:)=rHS !!this vector is needed to calculate the energy

  rHS=4.d0*pi*rHS
  if(l==0)rHS=rHS-(2.d0*sqrt(pi))*G    !!G is already multiplied by 4\pi !!??????

  AB=A+dble(l*(l+1))*B
  if(Robin==1)AB(nband+1,N)=AB(nband+1,N)+1.d0/Length
 call DPBSV('U', N,nband, 1, AB, nband+1,rHS, N, INFO )
 Ql(l+1,:)=rHS
 if(info.ne.0)then
   print*,'in solving the poisson equation...info=',info
   stop
 endif
1 enddo

end subroutine
