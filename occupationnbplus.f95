 Subroutine occupationnbplus(ii,ev,Z,mmax,occn)

 implicit none
 integer::Z,mmax,ii
 double precision::ev(Z,mmax+1),occn(mmax+1,Z)

 !internal symbols
 integer::i,k1m_1,k2m_1,k3m_1,k4m_1,k5m_1,k6m_1,k1m_2,k2m_2,k3m_2,k1m_3,k2m_3
 double precision::vector3(3),vector4(4),vector6(6),vector8(8),vector7(7)

 occn=0.d0
 if(Z.ge.37.and.Z.le.54) then
  if(ii==1)then
        occn(0+1,1)=1.d0!1s  
        occn(0+1,2)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,3)=1.d0!2s
        occn(0+1,4)=1.d0!3d_0
        occn(1+1,2)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2  
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,3)=2.d0!3p_1 
        occn(0+1,6)=1.d0!3s 
        if(Z.ge.37 .and. Z.le.42)then   
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!' 
        occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,4)=2.d0*occn(0+1,7)!4f_1
        occn(2+1,2)=2.d0*occn(0+1,7)!4f_2
        occn(3+1,1)=2.d0*occn(0+1,7)!4f_3
        elseif(Z.ge.43 .and. Z.le.52)then 
        occn(0+1,7)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
        if((dble(Z)/2.d0-sum(occn))>5.d0)print*,'error!!' 
        occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,8)!4d_1
        occn(2+1,3)=2.d0*occn(0+1,8)!4d_2  
        else
        occn(0+1,7)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
        occn(0+1,8)=1.d0!4d_0
        occn(1+1,5)=2.d0!4d_1
        occn(2+1,3)=2.d0!4d_2        
        if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!' 
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
        occn(1+1,6)=2.d0*occn(0+1,9)!4p_1         
        endif
  else
             !!where is the 2p
            do i=2,4
               vector3(i-1)=abs(ev(i,0+1)-ev(1,1+1))
            enddo
            do k1m_1=3,5
               if(abs(ev(k1m_1,1)-ev(1,2))==minval(vector3)) goto 11
            enddo
11     continue
             !!where is the 3p
            do i=5,7
               vector3(i-4)=abs(ev(i,0+1)-ev(2,1+1))
            enddo
            do k2m_1=5,7
               if(abs(ev(k2m_1,0+1)-ev(2,1+1))==minval(vector3)) goto 12
            enddo
12     continue
             !!index for 4p or 3d
            do i=6,8
               vector3(i-5)=abs(ev(i,0+1)-ev(3,1+1))
            enddo
            do k3m_1=6,8
               if(abs(ev(k3m_1,0+1)-ev(3,1+1))==minval(vector3)) goto 13
            enddo
13     continue
            !!index for 4p and 5s
            do i=7,10
               vector4(i-6)=abs(ev(i,0+1)-ev(4,1+1))
            enddo
            do k4m_1=7,10
               if(abs(ev(k4m_1,0+1)-ev(4,1+1))==minval(vector4)) goto 14
            enddo
14     continue   
     
            !!where is the 3d
            do i=5,10
               vector6(i-4)=abs(ev(i,0+1)-ev(1,2+1))
            enddo
            do k1m_2=5,10
               if(abs(ev(k1m_2,1)-ev(1,2+1))==minval(vector6)) goto 15
            enddo
15     continue

            !!where is the 4d
            do i=5,12
               vector8(i-4)=abs(ev(i,0+1)-ev(2,2+1))
            enddo
            do k2m_2=5,12
               if(abs(ev(k2m_2,1)-ev(2,2+1))==minval(vector8)) goto 16
            enddo
16     continue

            !!where is the 4f
            do i=7,13
               vector7(i-6)=abs(ev(i,0+1)-ev(1,3+1))
            enddo
            do k1m_3=7,13
               if(abs(ev(k1m_3,1)-ev(1,3+1))==minval(vector7)) goto 17
            enddo
17     continue

            !!where is the ??
            do i=8,13
               vector6(i-7)=abs(ev(i,0+1)-ev(5,1+1))
            enddo
            do k5m_1=8,13
               if(abs(ev(k5m_1,1)-ev(5,1+1))==minval(vector6)) goto 18
            enddo
18     continue

            !!where is the ??
            do i=9,14
               vector6(i-8)=abs(ev(i,0+1)-ev(6,1+1))
            enddo
            do k6m_1=9,14
               if(abs(ev(k6m_1,1)-ev(6,1+1))==minval(vector6)) goto 19
            enddo
19     continue

            !!where is the ??
          if(Z.ge.45) then
            do i=9,14
               vector6(i-8)=abs(ev(i,0+1)-ev(3,2+1))
            enddo
            do k3m_2=9,14
               if(abs(ev(k3m_2,1)-ev(3,2+1))==minval(vector6)) goto 20
            enddo
20         continue
           endif
          
            if(Z.ge. 51)then
            do i=9,14
               vector6(i-8)=abs(ev(i,0+1)-ev(2,3+1))
            enddo
            do k2m_3=9,14
               if(abs(ev(k2m_3,1)-ev(2,3+1))==minval(vector6)) goto 21
            enddo
21         continue
           endif
     
      if(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and. k1m_2==8 .and. Z.le.38)then
        print*,   'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s'
        write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!4s              
        occn(0+1,7)=1.d0!4p_0 
        occn(1+1,3)=2.d0!4p_1  
        occn(0+1,8)=1.d0!3d_0
        occn(1+1,4)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))!5s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k1m_2==6 .and. k2m_2==9 .and. Z.le.46)then
        print*,   'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d'
        write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,9)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,9)!4d_2       
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and. k1m_2==8 .and. k2m_2==10 .and. &
              k1m_3==11 .and. Z.ge.39 .and. Z.le. 48)then
        print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 4s 4p 3d 5s 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!4s
        occn(0+1,7)=1.d0!4p_0 
        occn(1+1,3)=2.d0!4p_1  
        occn(0+1,8)=1.d0!3d_0      
        occn(1+1,4)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,10)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,10)!4d_2  
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k1m_2==6 .and. k2m_2>=10 .and. Z.le.38)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))!5s   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and. k1m_2==8 .and. k2m_2==10 .and. &
              k1m_3==10 .and. Z.ge.39 .and. Z.le. 52)then
        print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 4s 4p 3d 5s 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_
        occn(0+1,6)=1.d0!4s
        occn(0+1,7)=1.d0!4p_0 
        occn(1+1,3)=2.d0!4p_1  
        occn(0+1,8)=1.d0!3d_0      
        occn(1+1,4)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,5)=2.d0*occn(0+1,10)!4f_1
        occn(2+1,3)=2.d0*occn(0+1,10)!4f_2 
        occn(3+1,1)=2.d0*occn(0+1,10)!4f_2   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and.k1m_2==8 .and. k2m_2==11 .and. &
              k1m_3==11 .and. k5m_1==10 .and. Z.ge.39 .and. Z.le. 44) then
        print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s 5p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 4s 4p 3d 5s 5p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!4s
        occn(0+1,7)=1.d0!4p_0 
        occn(1+1,3)=2.d0!4p_1  
        occn(0+1,8)=1.d0!3d_0      
        occn(1+1,4)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/3.d0!5p_0
        occn(1+1,5)=2.d0*occn(0+1,10)!5p_1   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and.k1m_2==8 .and. k2m_2==11 .and. &
              k1m_3==11 .and. k5m_1==13 .and. Z.ge.39 .and. Z.le. 40) then
        print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d 5s 6s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 4s 4p 3d 5s 6s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!4s
        occn(0+1,7)=1.d0!4p_0 
        occn(1+1,3)=2.d0!4p_1  
        occn(0+1,8)=1.d0!3d_0      
        occn(1+1,4)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))!6s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. &
              (k5m_1==10.or.(k5m_1==11 .and. Z==42)).and. k1m_2==6 .and. &!!!!(k5m_1==11 .and. Z==42) this one is wrong but I face it with atom 42 because of the
                k2m_2==10 .and. k1m_3>=11 .and. Z.ge.39 .and. Z.le. 48)then!!! accidental degeneracy it was confused
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0      
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,10)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,10)!4d_2   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k1m_2==6 .and. &
               k2m_2==8 .and. k1m_3==8 .and. Z.ge.31 .and. Z.le. 44)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0      
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,4)=2.d0*occn(0+1,8)!4f_1
        occn(2+1,2)=2.d0*occn(0+1,8)!4f_2   
        occn(3+1,1)=2.d0*occn(0+1,8)!4f_3  
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==10 .and.&
              k1m_2==6 .and. k2m_2>=11 .and. Z.ge. 39 .and. Z.le. 44)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        occn(0+1,9)=1.d0!5s   
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/3.d0!5p_0
        occn(1+1,5)=2.d0*occn(0+1,10)!5p_1
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==11 .and.&
              k1m_2==6 .and. k2m_2>=12 .and. Z.ge. 41 .and. Z.le. 46)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0   
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/3.d0!5p_0
        occn(1+1,5)=2.d0*occn(0+1,11)!5p_1
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==9 .and. k5m_1==10 .and.&
              k1m_2==6 .and. k2m_2==10 .and. k1m_3.ge.11 .and. Z.ge. 39 .and.Z.le. 48)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 5s 4p 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 5s 4p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!5s
        occn(0+1,9)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,10)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,10)!4d_2
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and.&
              k1m_2==6 .and. k2m_2==8 .and. k3m_2>=10 .and. k1m_3==8 .and. Z.ge. 45 .and.Z.le. 50)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4f 4p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4f 4p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2 
        occn(3+1,1)=2.d0!4f_3
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
        occn(1+1,5)=2.d0*occn(0+1,9)!4p_1   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==10 .and.&
              k1m_2==6 .and. k2m_2==11 .and. k1m_3>=11 .and. Z.ge. 44 .and.Z.le. 54)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0!5p_0
        occn(1+1,5)=2.d0!5p_1
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,6)=2.d0*occn(0+1,11)!4d_1  
        occn(2+1,2)=2.d0*occn(0+1,11)!4d_2   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==10 .and. k6m_1==12 .and.&
              k1m_2==6 .and. k2m_2==12 .and. k1m_3>=12 .and. Z.ge. 45 .and.Z.le. 46)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p 6s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 5p 6s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0!5p_0
        occn(1+1,5)=2.d0!5p_1
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))!6s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k6m_1==10 .and.&
              k1m_2==6 .and. k2m_2==9 .and. k1m_3==10 .and. Z.ge. 46 .and.Z.le. 60)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!4d_0
        occn(1+1,5)=2.d0!4d_1
        occn(2+1,2)=2.d0!4d_2
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,6)=2.d0*occn(0+1,10)!4f_1
        occn(2+1,3)=2.d0*occn(0+1,10)!4f_2
        occn(3+1,1)=2.d0*occn(0+1,10)!4f_3
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k6m_1>=11 .and.&
              k1m_2==6 .and. k2m_2==9 .and. k1m_3>=11 .and. Z.ge. 47 .and.Z.le. 48)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!4d_0
        occn(1+1,5)=2.d0!4d_1
        occn(2+1,2)=2.d0!4d_2
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))!5s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==11 .and. k6m_1==12 .and.&
              k1m_2==6 .and. k2m_2==12 .and. k1m_3>=13 .and. Z.ge. 47 .and.Z.le. 56)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0!6s
        occn(0+1,11)=1.d0!5p_0
        occn(1+1,5)=2.d0!5p_1
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,12)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,6)=2.d0*occn(0+1,12)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,12)!4d_2
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==11 .and. k6m_1>12 .and.&
              k1m_2==6 .and. k2m_2>=12 .and. Z.ge. 47 .and.Z.le. 48)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p 7s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 5p 7s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0!6s
        occn(0+1,11)=1.d0!5p_0
        occn(1+1,5)=2.d0!5p_1
        if((dble(Z)/2.d0-sum(occn))>1.d0)print*,'error!!'
        occn(0+1,12)=(dble(Z)/2.d0-sum(occn))!7s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k6m_1==11 .and.&
              k1m_2==6 .and. k2m_2==9 .and. k3m_2==11 .and. k1m_3==11 .and. Z.ge. 49 .and.Z.le. 62)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!4d_0
        occn(1+1,5)=2.d0!4d_1
        occn(2+1,2)=2.d0!4d_2
        occn(0+1,10)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,6)=2.d0*occn(0+1,11)!4f_1
        occn(2+1,3)=2.d0*occn(0+1,11)!4f_2
        occn(3+1,1)=2.d0*occn(0+1,11)!4f_3
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k6m_1==11 .and.&
              k1m_2==6 .and. k2m_2==9 .and. k3m_2>=12 .and. k1m_3>=12 .and. Z.ge. 49 .and.Z.le. 54)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!4d_0
        occn(1+1,5)=2.d0!4d_1
        occn(2+1,2)=2.d0!4d_2
        occn(0+1,10)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/3.d0!5p_0
        occn(1+1,6)=2.d0*occn(0+1,11)!5p_1
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==9 .and. k6m_1==10 .and.&
              k1m_2==6 .and. k2m_2==8 .and. k3m_2==10 .and. k1m_3==8 .and. Z.ge. 51 .and.Z.le. 60)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4f 4p 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4f 4p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2 
        occn(3+1,1)=2.d0!4f_3
        occn(0+1,9)=1.d0!4p_0
        occn(1+1,5)=2.d0!4p_1 
        if((dble(Z)/2.d0-sum(occn))>5.d0)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,6)=2.d0*occn(0+1,10)!4d_1
        occn(2+1,3)=2.d0*occn(0+1,10)!4d_2
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1>=11 .and. k6m_1>=11 .and.&
              k1m_2==6 .and. k2m_2==12 .and. k1m_3>=11 .and. Z.ge. 39 .and.Z.le. 40)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))!6s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1==11 .and. k6m_1>=12 .and.&
              k1m_2==6 .and. k2m_2==11 .and. k1m_3>=12 .and. Z.ge. 41 .and.Z.le. 50)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0
        occn(1+1,4)=2.d0!4p_1
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0!6s
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,5)=2.d0*occn(0+1,11)!4d_1
        occn(2+1,2)=2.d0*occn(0+1,11)!4d_2
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==7 .and. k4m_1==8 .and. k5m_1>=9 .and. k6m_1>=9 .and.&
              k1m_2==7 .and. k2m_2==8 .and. k1m_3==8 .and. Z.ge. 31 .and.Z.le. 44)then
        print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 4s 3d 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!4s
        occn(0+1,7)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,4)=2.d0*occn(0+1,8)!4f_1
        occn(2+1,2)=2.d0*occn(0+1,8)!4f_2
        occn(3+1,1)=2.d0*occn(0+1,8)!4f_3
       elseif(k1m_1==3 .and. k2m_1==6 .and. k3m_1==7 .and. k4m_1==8 .and. k5m_1>=10 .and. k6m_1>=10 .and.&
              k1m_2==7 .and. k2m_2==8 .and. k1m_3==8 .and. Z.ge. 45 .and.Z.le. 46)then
        print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 4f 5s'
        write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 4f 5s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!4s
        occn(0+1,6)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,7)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,8)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))!5s
       elseif(k1m_1==3 .and. k2m_1==6 .and. k3m_1==7 .and. k4m_1==9 .and. k5m_1>=10 .and. k6m_1>=10 .and.&
              k1m_2==7 .and. k2m_2==9 .and. k1m_3==9 .and. Z==46)then
        print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 5s 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 5s 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!4s
        occn(0+1,6)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,7)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,8)=1.d0!5s
        occn(0+1,9)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
       elseif(k1m_1==3 .and. k2m_1==6 .and. k3m_1==7 .and. k4m_1==9 .and. k5m_1==10 .and. k6m_1>=11 .and.&
              k1m_2==7 .and. k2m_2==9 .and. k1m_3==9 .and. Z.ge. 47 .and. Z.le. 52)then
        print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 5s 4f 4p'
        write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 5s 4f 4p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!4s
        occn(0+1,6)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,7)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,8)=1.d0!5s
        occn(0+1,9)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
        occn(1+1,5)=2.d0*occn(0+1,10)!4p_1
       elseif(k1m_1==3 .and. k2m_1==6 .and. k3m_1==7 .and. k4m_1==9 .and. k5m_1==10 .and. k6m_1==11 .and.&
              k1m_2==7 .and. k2m_2==9 .and. k3m_2==11.and. k1m_3==9 .and. k2m_3>=12 .and.&
               Z.ge. 53 .and. Z.le. 62)then
        print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 5s 4f 4p 4d'
        write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 5s 4f 4p 4d'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!4s
        occn(0+1,6)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,7)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,8)=1.d0!5s
        occn(0+1,9)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2
        occn(3+1,1)=2.d0!4f_3
        occn(0+1,10)=1.d0!4p_0
        occn(1+1,5)=2.d0!4p_1
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0
        occn(1+1,6)=2.d0*occn(0+1,11)!4d_1
        occn(2+1,3)=2.d0*occn(0+1,11)!4d_2
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==7 .and.&
              k1m_2==6 .and. k2m_2==7 .and. k1m_3==7 .and. Z.ge. 37 .and. Z.le. 42)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4f'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4f'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        if((dble(Z)/2.d0-sum(occn))>7)print*,'error!!'
        occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/7.d0!4f_0
        occn(1+1,4)=2.d0*occn(0+1,7)!4f_1
        occn(2+1,2)=2.d0*occn(0+1,7)!4f_2 
        occn(3+1,1)=2.d0*occn(0+1,7)!4f_3     
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==7 .and. k5m_1==9 .and.&
              k1m_2==6 .and. k2m_2==7 .and. k3m_2>=10 .and. k1m_3==7 .and. Z.ge. 45 .and. Z.le. 50)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4f 4s 4p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4f 4s 4p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2 
        occn(3+1,1)=2.d0!4f_3
        occn(0+1,8)=1.d0!4s
        if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
        occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0        
        occn(1+1,5)=2.d0*occn(0+1,9)!4p_1   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==7 .and. k5m_1==9 .and.&
              k1m_2==6 .and. k2m_2==7 .and. k1m_3==7 .and. Z.ge. 43 .and. Z.le. 44)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4f 4s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4f 4s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2 
        occn(3+1,1)=2.d0!4f_3
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,8)=(dble(Z)/2.d0-sum(occn))!4s
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==7 .and. k5m_1==9 .and.&
              k1m_2==6 .and. k2m_2==7 .and. k3m_2==10 .and. k1m_3==7 .and. k2m_3>=11 .and.&
              Z.ge. 51 .and. Z.le. 60)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4f 4s 4p'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4f 4s 4p'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4f_0
        occn(1+1,4)=2.d0!4f_1
        occn(2+1,2)=2.d0!4f_2 
        occn(3+1,1)=2.d0!4f_3
        occn(0+1,8)=1.d0!4s
        occn(0+1,9)=1.d0!4p_0        
        occn(1+1,5)=2.d0!4p_1 
        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
        occn(0+1,10)=(dble(Z)/2.d0-sum(occn))/5.d0!4d_0        
        occn(1+1,6)=2.d0*occn(0+1,10)!4d_1
        occn(2+1,3)=2.d0*occn(0+1,10)!4d_2   
       elseif(k1m_1==3 .and. k2m_1==5 .and. k3m_1==6 .and. k4m_1==8 .and. k5m_1>11 .and.&
              k1m_2==6 .and. k2m_2>11 .and. Z.ge. 41 .and. Z.le. 42)then
        print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 7s'
        write(3,*)'the configuration 1s 2s 2p 3s 3p 3d 4s 4p 5s 6s 7s'
        occn(0+1,1)=1.d0!1s
        occn(0+1,2)=1.d0!2s 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1
        occn(0+1,6)=1.d0!3d_0
        occn(1+1,3)=2.d0!3d_1
        occn(2+1,1)=2.d0!3d_2
        occn(0+1,7)=1.d0!4s
        occn(0+1,8)=1.d0!4p_0 
        occn(1+1,4)=2.d0!4p_1  
        occn(0+1,9)=1.d0!5s
        occn(0+1,10)=1.d0   
        if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
        occn(0+1,11)=(dble(Z)/2.d0-sum(occn))!7s
       else      
         print*,'m=1',k1m_1,k2m_1,k3m_1,k4m_1,k5m_1,k6m_1,'m=2',k1m_2,k2m_2,k3m_2,'m=3',k1m_3,k2m_3
         stop
      endif!!cases
  endif!i==1
 else
  print*,'did not reach this Z yet'
  stop
 endif!36<Z<55

  if(maxval(occn)>2)then
  print*,maxval(occn), 'there is something wrong with the occupation number plus'
  stop
  endif

 occn=2.d0*occn
  print*,ii,'sum of occupation number=', sum(occn)

 if(abs(sum(occn)-Z)>0.1) then
  print*,'error in the occupationnbplus'
  stop
 endif
end Subroutine
