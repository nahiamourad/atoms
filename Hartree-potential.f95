  Subroutine Hartreepotential(N,mn,lrho,ninter,Q,Ls,mmax,FF,FF1,constantM,Hhat,nband,&
                             Robin)
  Implicit none
  integer::N,i,j,k,mn,lrho,l1,l2,l,nn,ninter,Ls,mmax,nband,Robin,kappa
  double precision::temp
  double precision,dimension(:,:)::Matrix(N,N),Q(lrho+1,N),Hhat((Ls-mn+1)*N,(Ls-mn+1)*N)
  double precision,dimension(:,:,:)::FF1(4,4,4)
  double precision,dimension(:,:,:,:)::FF(5,5,5,ninter-1),constantM(Ls+1,Ls+1,2*Ls+1,mmax+1)
  
 kappa=1
 if(Robin==1 .or. Robin==2)kappa=0

  Hhat=0.d0
  do 10 l2=mn,Ls
    do 11 l1=mn,l2        !!! construct just the upper triangular  of Hhat l2>=l1
         Matrix=0.d0
         do i=1,3
          do j=0,i-1
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 1
               temp=0.d0
               do nn=0,3
                 temp=temp+FF1(i+1,j+1,nn+1)*Q(l+1,4-nn)                  
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4-i,4-j)=Matrix(4-i,4-j)+temp
                     !alpha(Lambda(i),Lambda(j))/hk(1)
1             continue
             enddo
          enddo
         enddo
         do k=2,ninter-1
          do i=1,4
           do j=0,i-1
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 2
              temp=0.d0
              do nn=0,4
               temp=temp+FF(i+1,j+1,nn+1,k-1)*Q(l+1,4*k-nn)                      
                                  !alpha(Lambda(i),Lambda(j))/hk(k)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*k-i,4*k-j)=Matrix(4*k-i,4*k-j)+temp
2             continue
             enddo
           enddo
          enddo
         enddo 

         do i=2,4
          do j=1,i-1
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 3
              temp=0.D0
              do nn=kappa,4
                temp=temp+FF(i+1,j+1,nn+1,ninter-1)*Q(l+1,4*ninter-nn)  
                                 !alpha(Lambda(i),Lambda(j))/hk(ninter)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                Matrix(4*ninter-i,4*ninter-j)=Matrix(4*ninter-i,4*ninter-j)+temp
3             continue
             enddo
          enddo
         enddo
	
     if(Robin==1 .or. Robin==2)then
         do i=1,4
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 13
              temp=0.D0
              do nn=kappa,4
                temp=temp+FF(i+1,1,nn+1,ninter-1)*Q(l+1,4*ninter-nn)  
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                Matrix(4*ninter-i,4*ninter)=Matrix(4*ninter-i,4*ninter)+temp
13            continue
             enddo
          enddo
     endif

         Matrix=Matrix+transpose(Matrix)
         
          do k=3,ninter-1
            do l=0,lrho
              if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 4
              temp=0.d0
              do nn=0,4
                 temp=temp+(FF(1,1,nn+1,k-2)*Q(l+1,4*k-4-nn)+FF(5,5,nn+1,k-1)*Q(l+1,4*k-nn))
                                  !!Lambda(0)=2; Lambda(4)=1 !alpha(2,2)/hk(k-1)+alpha(1,1)/hk(k)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*k-4,4*k-4)=Matrix(4*k-4,4*k-4)+temp
4             continue
             enddo
         enddo
           !for k=2 I have a special situation with k-1
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 5
               temp=0.d0
              do nn=0,3
                 temp=temp+(FF1(1,nn+1,1)*Q(l+1,4-nn)+FF(5,5,nn+1,1)*Q(l+1,8-nn))
                                  !alpha(2,2)/hk(k-1)+alpha(1,1)/hk(k)
              enddo
                 temp=temp+FF(5,5,5,1)*Q(l+1,4)
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp 
                 Matrix(4,4)=Matrix(4,4)+temp
5             continue
             enddo
        
           !for k=ninter
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 6
               temp=0.d0
              do nn=kappa,4
                 temp=temp+(FF(1,1,nn+1,ninter-2)*Q(l+1,4*ninter-4-nn)+&
                            FF(5,5,nn+1,ninter-1)*Q(l+1,4*ninter-nn))
                            !alpha(2,2)/hk(k-1)+alpha(1,1)/hk(k)
              enddo
                if(kappa==1)temp=temp+FF(1,1,1,ninter-2)*Q(l+1,4*ninter-4)
                temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*ninter-4,4*ninter-4)=Matrix(4*ninter-4,4*ninter-4)+temp
6             continue
             enddo

         if(Robin==1 .or. Robin==2)then
               do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 16
               temp=0.d0
              do nn=kappa,4
                 temp=temp+FF(1,1,nn+1,ninter-1)*Q(l+1,4*ninter-nn)
                            !alpha(2,2)/hk(k-1)+alpha(1,1)/hk(k)
              enddo
                temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*ninter,4*ninter)=Matrix(4*ninter,4*ninter)+temp
16             continue
             enddo            
         endif

         if(Robin==2)then
               do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 116
                temp=constantM(l1+1,l2+1,l+1,mn+1)*Q(l+1,4*ninter)/3.d0
                 Matrix(4*ninter,4*ninter)=Matrix(4*ninter,4*ninter)+temp
116             continue
             enddo            
         endif


          do k=2,ninter-1
           do i=1,3
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 7
                temp=0.d0
              do nn=0,4
                 temp=temp+FF(i+1,i+1,nn+1,k-1)*Q(l+1,4*k-nn)
                                     !alpha(Lambda(i),Lambda(i))/hk(k)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*k-i,4*k-i)=Matrix(4*k-i,4*k-i)+temp
7             continue
             enddo
           enddo
          enddo
           !for k=1
           do i=1,3
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 8
               temp=0.d0
              do nn=0,3
                 temp=temp+FF1(i+1,i+1,nn+1)*Q(l+1,4-nn)
                                 !alpha(Lambda(i),Lambda(i))/hk(k)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4-i,4-i)=Matrix(4-i,4-i)+temp
8             continue
             enddo
           enddo
           !for k=ninter
           do i=1,3
            do l=0,lrho
               if(constantM(l1+1,l2+1,l+1,mn+1)==0) goto 9
                 temp=0.d0
              do nn=kappa,4
                 temp=temp+FF(i+1,i+1,nn+1,ninter-1)*Q(l+1,4*ninter-nn)
                                               !alpha(Lambda(i),Lambda(i))/hk(k)
              enddo
                 temp=constantM(l1+1,l2+1,l+1,mn+1)*temp
                 Matrix(4*ninter-i,4*ninter-i)=Matrix(4*ninter-i,4*ninter-i)+temp
9             continue
             enddo
           enddo
    
     do i=1,N                                
      do j=max(1,i-nband),min(N,i+nband) 
         Hhat((l1-mn)*N+i,(l2-mn)*N+j)=Matrix(i,j)
      enddo
     enddo
         
11     enddo
10      enddo
 End Subroutine

