 Subroutine occupationnb(ii,ev,Z,mmax,occn)

 implicit none
 integer::ii,i,k1m_1,k2m_1,k3m_1,k4m_1,k1m_2,Z,mmax,lpd
 double precision::ev(Z,mmax+1),occn(mmax+1,Z),B(2),B3(3),B6(6),B4(4),Sm!test(4),test1(8),
  
 Sm=sum(occn)
  occn=0.d0 
  k4m_1=0

  if(Z==1) then!!!! Z=1 Hydrogen
      occn(1,1)=0.5d0

  elseif(Z==2) then!!!! Z=2 Helium 
      occn(1,1)=1.d0

  elseif(Z==3) then!!!! Z=3 Lithium  
      occn(0+1,1)=1.d0
      occn(0+1,2)=0.5d0

  elseif(Z==4) then!!!! Z=4 Beryllium 
      occn(0+1,1)=1.d0
      occn(0+1,2)=1.d0

  elseif(Z.ge.5.and.Z.le.10)then
      occn(0+1,1)=1.d0!1s
      occn(0+1,2)=1.d0!2s
     if(ii.ne.1) then
            do k1m_2=3,4
               B(k1m_2-2)=abs(ev(k1m_2,1)-ev(1,2))
            enddo
            do i=3,4
               if(abs(ev(i,1)-ev(1,2))==minval(B)) goto 12
            enddo
12     continue
        if(i==3) then
           occn(0+1,3)=(Z/2.d0-sum(occn))/3.d0!2p_0
           occn(1+1,1)=2.d0*occn(0+1,3) !2p_1
        else
           print*,' the 3s goes below the 2p'
           occn(0+1,3)=1.d0!3s
           occn(0+1,4)=(Z/2.d0-sum(occn))/3.d0!2p_0
           occn(1+1,1)=2.d0*occn(0+1,4)!2p_1
        endif
      else !!ii==1
      occn(0+1,3)=(Z/2.d0-sum(occn))/3.d0!2p_0
      occn(1+1,1)=2.d0*occn(0+1,3)!2p_1
     endif
   
  elseif(Z.ge.11.and.Z.le.12) then
      occn(0+1,1)=1.d0!1s
      occn(0+1,2)=1.d0!2s  
      occn(0+1,3)=1.d0!2p_0
      occn(1+1,1)=2.d0!2p_1
      if(ii.ne.1) then
             !!where is the 2p
            do k1m_1=3,4
               B(k1m_1-2)=abs(ev(k1m_1,0+1)-ev(1,1+1))
            enddo
            do i=3,4
               if(abs(ev(i,1)-ev(1,2))==minval(B)) goto 22
            enddo
22     continue
            !!where is the 3d
            do k1m_1=4,6
               B(k1m_1-3)=abs(ev(k1m_1,0+1)-ev(1,2+1))
            enddo
            do k1m_1=4,6
               if(abs(ev(k1m_1,1)-ev(1,2+1))==minval(B)) goto 21
            enddo
21     continue

         if(k1m_1==4 .and. i==3) then
           print*,'the configuration is 1s 2s 2p 3d'
           write(3,*)'the configuration is 1s 2s 2p 3d'
           occn(0+1,4)=(Z/2.d0-sum(occn))/5.d0!3d_0
           occn(1+1,2)=2.d0*occn(0+1,4)!3d_1
           occn(2+1,1)=2.d0*occn(0+1,4)!3d_2
          elseif(i==3 .and.k1m_1>4)then
           print*,'the configuration is 1s 2s 2p 3s'
           write(3,*)'the configuration is 1s 2s 2p 3s'
           if((Z/2.d0-sum(occn))>1.d0)print*,'error!!!!!!!'
           occn(0+1,4)=(Z/2.d0-sum(occn))!3s 
          elseif(i==4)then
           print*,'the configuration is 1s 2s 3s 2p'
           write(3,*)'the configuration is 1s 2s 3s 2p'
           occn(0+1,3)=0.d0! the 3s goes below the 2p
           occn(1+1,1)=0.d0! the 3s goes below the 2p
           occn(0+1,3)=1.d0!3s
           if((Z/2.d0-sum(occn))>3.d0)print*,'error!!!!!!!'
           occn(0+1,4)=(Z/2.d0-sum(occn))/3.d0!2p_0
           occn(1+1,1)=2.d0*occn(0+1,4)!2p_1
           else
            print*,'something missing'
         endif
      else
       occn(0+1,4)=(Z/2.d0-sum(occn))!3s      
      endif

  elseif(Z.ge.13.and.Z.le.18) then
       occn(0+1,1)=1.d0!1s
       occn(0+1,2)=1.d0!2s 
      if(ii.ne.1)then
             !!where is the 2p
            do k1m_1=3,4
               B(k1m_1-2)=abs(ev(k1m_1,0+1)-ev(1,1+1))
            enddo
            do i=3,4
               if(abs(ev(i,1)-ev(1,2))==minval(B)) goto 13
            enddo
13     continue
            !!where is the 3d
            do k1m_2=5,6
               B(k1m_2-4)=abs(ev(k1m_2,0+1)-ev(1,2+1))
            enddo
            do k1m_2=5,6
               if(abs(ev(k1m_2,1)-ev(1,2+1))==minval(B)) goto 14
            enddo
14     continue
             !!where is the 3p
            do k2m_1=5,7
               B3(k2m_1-4)=abs(ev(k2m_1,0+1)-ev(2,1+1))
            enddo
            do k2m_1=5,7
               if(abs(ev(k2m_1,0+1)-ev(2,1+1))==minval(B3)) goto 15
            enddo
15     continue
         if(i==3.and.k2m_1==5) then
           occn(0+1,3)=1.d0!2p_0
           occn(1+1,1)=2.d0!2p_1
           occn(0+1,4)=1.d0!3s 
           occn(0+1,5)=(Z/2.d0-sum(occn))/3.d0!3p_0
           occn(1+1,2)=2.d0*occn(0+1,5)!3p_1  
         elseif(i==3.and.k2m_1==6) then
            print*,'the 4s is below the 3p'
           occn(0+1,3)=1.d0!2p_0
           occn(1+1,1)=2.d0!2p_1
           occn(0+1,4)=1.d0!3s 
            if(Z.le.14)then
              occn(0+1,5)=(Z/2.d0-sum(occn))!4s  
            else
              occn(0+1,5)=1.d0!4s               
              occn(0+1,6)=(Z/2.d0-sum(occn))/3.d0!3p
              occn(1+1,2)=2.d0*occn(0+1,6)!3p_1  
            endif
         elseif(i==4.and. k1m_2==5)then
            print*,'the 3s is below the 2p and the 3d below the 3p'
           occn(0+1,3)=1.d0!3s 
           occn(0+1,4)=1.d0!2p_0
           occn(1+1,1)=2.d0!2p_1
           occn(0+1,5)=(Z/2.d0-sum(occn))/5.d0!3d_0
           occn(1+1,2)=2.d0*occn(0+1,5)!3d_1 
           occn(2+1,1)=2.d0*occn(0+1,5)!3d_2
          else
           print*, 'could not find this case'
         endif    
      else !!ii==1
       occn(0+1,3)=1.d0!2p_0
       occn(1+1,1)=2.d0!2p_1
       occn(0+1,4)=1.d0!3s 
       occn(0+1,5)=(Z/2.d0-sum(occn))/3.d0!3p_0
       occn(1+1,2)=2.d0*occn(0+1,5)!3p_1   
      endif 

  elseif(Z.ge.19.and.Z.le.36) then
      occn(0+1,1)=1.d0!1s
      occn(0+1,2)=1.d0!2s  
      if(ii.ne.1)then
             !!where is the 2p
            do i=3,4
               B(i-2)=abs(ev(i,0+1)-ev(1,1+1))
            enddo
            do i=3,4
               if(abs(ev(i,1)-ev(1,2))==minval(B)) goto 16
            enddo
16     continue
             !!where is the 3p
            do k2m_1=5,7
               B3(k2m_1-4)=abs(ev(k2m_1,0+1)-ev(2,1+1))
            enddo
            do k2m_1=5,7
               if(abs(ev(k2m_1,0+1)-ev(2,1+1))==minval(B3)) goto 18
            enddo
18     continue
             !!index for 4p or 3d
            do k3m_1=6,8
               B3(k3m_1-5)=abs(ev(k3m_1,0+1)-ev(3,1+1))
            enddo
            do k3m_1=6,8
               if(abs(ev(k3m_1,0+1)-ev(3,1+1))==minval(B3)) goto 19
            enddo
19     continue

        if(Z.ge.21)then
            !!index for 4p and 5s
            do k4m_1=7,10
               B4(k4m_1-6)=abs(ev(k4m_1,0+1)-ev(4,1+1))
            enddo
            do k4m_1=7,10
               if(abs(ev(k4m_1,0+1)-ev(4,1+1))==minval(B4)) goto 20
            enddo
20     continue   
        endif          
            !!where is the 3d k1m
            do k1m_2=5,10
               B6(k1m_2-4)=abs(ev(k1m_2,0+1)-ev(1,2+1))
            enddo
            do k1m_2=5,10
               if(abs(ev(k1m_2,1)-ev(1,2+1))==minval(B6)) goto 17
            enddo
17     continue

         if(i==4 .and. k2m_1==5 .and. k1m_2==5) then
            occn(0+1,3)=1.d0!3s
            occn(0+1,4)=1.d0!2p_0
            occn(1+1,1)=2.d0!2p_1
                if(Z.le.22)then
                  print*,'the configuration is 1s 2s 3s 2p 3d'
                  write(3,*)'the configuration is 1s 2s 3s 2p 3d'
                  if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
                  occn(0+1,5)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                  occn(1+1,2)=2.d0*occn(0+1,5)!3d_1
                  occn(2+1,1)=2.d0*occn(0+1,5)!3d_2
                elseif(Z.ge.23 .and. Z.le.24)then
                  print*,'the configuration is 1s 2s 3s 2p 3d 4s'
                  write(3,*)'the configuration is 1s 2s 3s 2p 3d 4s'
                  occn(0+1,5)=1.d0!3d_0
                  occn(1+1,2)=2.d0!3d_1
                  occn(2+1,1)=2.d0!3d_2
                  if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!!!'
                  occn(0+1,6)=(dble(Z)/2.d0-sum(occn))!4s   
                 else
                  print*, 'need more shells'
                  stop               
                endif
          elseif(i==3 .and. k2m_1==5 .and. k1m_2==6)then
            occn(0+1,3)=1.d0!2p_0
            occn(1+1,1)=2.d0!2p_1
            occn(0+1,4)=1.d0!3s
            occn(0+1,5)=1.d0!3p_0
            occn(1+1,2)=2.d0!3p_1
               if(Z.le. 28) then
                  print*,'the configuration is 1s 2s 2p 3s 3p 3d'
                  write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d'
                  if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
                  occn(0+1,6)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                  occn(1+1,3)=2.d0*occn(0+1,6)!3d_1
                  occn(2+1,1)=2.d0*occn(0+1,6)!3d_2
               elseif(Z.ge. 28 .and.Z.le.30)then
                  print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s'! I am not sure that the last one is really 4s it works as far as the density is radial'
                  write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s'
                  occn(0+1,6)=1.d0!3d_0
                  occn(1+1,3)=2.d0!3d_1
                  occn(2+1,1)=2.d0!3d_2
                  if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!!!'
                  occn(0+1,7)=(dble(Z)/2.d0-sum(occn))!4s
               elseif(Z.ge.31.and. Z.le. 38)then 
                      occn(0+1,6)=1.d0!3d_0
                      occn(1+1,3)=2.d0!3d_1
                      occn(2+1,1)=2.d0!3d_2
                      occn(0+1,7)=1.d0!4s
                   if(k4m_1==8)then
                      print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p'
                      write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p'
                      if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!!!'
                      occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p
                      occn(1+1,4)=2.d0*occn(0+1,8)
                    elseif(k4m_1==9)then
                      if(Z==31)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 5s'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s 5s'
                          if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!!!'
                          occn(0+1,8)=(dble(Z)/2.d0-sum(occn))!5s
                       else!!Z>32 and Z<38
                             print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 5s 4p'
                             write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s 5s 4p'
                             occn(0+1,8)=1.d0!5s
                             if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!!!'
                             occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                             occn(1+1,4)=2.d0*occn(0+1,9)
                       endif
                    else
                      print*,'precise this'
                    endif
                  else
                    print*,'need more shells defined till 38'
                    stop
               endif 

          elseif(i==3 .and. k2m_1==5 .and. k1m_2.ge. 7 .and. k3m_1.ge. 7)then    
            occn(0+1,3)=1.d0!2p_0
            occn(1+1,1)=2.d0!2p_1
            occn(0+1,4)=1.d0!3s
            occn(0+1,5)=1.d0!3p_0
            occn(1+1,2)=2.d0!3p_1
            if(Z.le.20)then
               print*,'the configuration is 1s 2s 2p 3s 3p 4s'
               if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'
               occn(0+1,6)=(dble(Z)/2.d0-sum(occn))!4s
            else!!Z>20
               occn(0+1,6)=1.d0!4s
                if(k1m_2==7.and. k3m_1==7)then
                    if(Z.le.30)then
                         !!!!!!special case faced only when there is degeneracy for which it is not able to distinguish d from p 
                         if(Z==21 .and. k4m_1==7)then
                                      do lpd=3,5
                                           B3(lpd-2)=abs(ev(lpd,1+1)-ev(1,2+1))
                                      enddo
                                       do lpd=3,5
                                         if(abs(ev(lpd,1+1)-ev(1,2+1))==minval(B3)) goto 24
                                      enddo
24     continue
                         if(lpd==3)then
                           goto 25
                         elseif(lpd==4)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p(degenracy)'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p(degeneracy)'
                          if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'  
                          occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                          occn(1+1,3)=2.d0*occn(0+1,7)!4p_1
                          goto 27
                         else
                         print*,'something to be corrscted here'
                         write(3,*)'something to be corrscted here'
                         stop
                         endif
                        endif
25     continue 
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d'
                        write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 3d'
                        if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
                        occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                        occn(1+1,3)=2.d0*occn(0+1,7)!3d_1
                        occn(2+1,1)=2.d0*occn(0+1,7)!3d_2 
                     elseif(Z.ge.31 .and.Z.le. 36) then
                        if(k4m_1==8)then
                           print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d 4p'
                           write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 3d 4p'
                           occn(0+1,7)=1.d0!3d_0
                           occn(1+1,3)=2.d0!3d_1
                           occn(2+1,1)=2.d0!3d_2
                           if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!!!'
                           occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                           occn(1+1,4)=2.d0*occn(0+1,8)!4p_1
                        else
                           print*,'5s may be here'
                        endif
                      else
                        print*,' need more shells defined for 36 only'                
                     endif                  
                elseif(k3m_1==7 .and. k1m_2>7)then 
                      if(Z.le.26)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p'
                          if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'  
                          occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                          occn(1+1,3)=2.d0*occn(0+1,7)!4p_1
                      elseif(Z.ge. 27 .and. Z.le. 38) then
                          if(k1m_2==8)then
                           print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d'
                           write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 3d'
                           occn(0+1,7)=1.d0!4p_0 
                           occn(1+1,3)=2.d0!4p_1
                           if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'  
                           occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                           occn(1+1,4)=2.d0*occn(0+1,8)!3d_1
                           occn(2+1,1)=2.d0*occn(0+1,8)!3d_2
                          elseif(k1m_2==9)then
                                if(Z.le. 28)then
                                   print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s'
                                   write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s'
                                   occn(0+1,7)=1.d0!4p_0 
                                   occn(1+1,3)=2.d0!4p_1
                                   if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'  
                                   occn(0+1,8)=(dble(Z)/2.d0-sum(occn))!5s
                                elseif(Z.ge. 29 .and.Z.le. 38)then
                                   print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s 3d'
                                   write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s 3d'
                                   occn(0+1,7)=1.d0!4p_0 
                                   occn(1+1,3)=2.d0!4p_1 
                                   occn(0+1,8)=1.d0!5s
                                   if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!' 
                                   occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                                   occn(1+1,4)=2.d0*occn(0+1,9)!3d_1
                                   occn(2+1,1)=2.d0*occn(0+1,9)!3d_2
                                 else
                                   print*,'need more shells defined for at most 28'
                                   stop
                                endif
                          elseif(k1m_2==10)then
                                 if(Z.ge. 29.and. Z.le. 34)then !!this is seen for Z.ge.29
                                   print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s 5p'
                                   write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s 5p'
                                   occn(0+1,7)=1.d0!4p_0 
                                   occn(1+1,3)=2.d0!4p_1 
                                   occn(0+1,8)=1.d0!5s
                                   if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!' 
                                   occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!3p_0
                                   occn(1+1,4)=2.d0*occn(0+1,9)!3p_1
                                  elseif(Z.le. 29)then!!Z<29
                                   print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s'
                                   write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p 5s'
                                   occn(0+1,7)=1.d0!4p_0 
                                   occn(1+1,3)=2.d0!4p_1
                                   if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'  
                                   occn(0+1,8)=(dble(Z)/2.d0-sum(occn))!5s   
                                  else
                                    print*, ' need to define more shells !!'
                                    stop                                
                                  endif
                          else
                           print*,'should define this'
                               print*,k1m_2,k2m_1
                               stop
                          endif
                      endif
                elseif(k3m_1==8 .and. k1m_2>7)then 
                    if(Z.le.22)then
                       if(k3m_1.ne.k4m_1)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s'
                          if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'  
                          occn(0+1,7)=(dble(Z)/2.d0-sum(occn))!5s
                         !!!!!!special case faced only when there is degeneracy for which it is not able to distinguish d from p 
                       elseif(k3m_1==k4m_1)then
                                      do lpd=3,5
                                           B3(lpd-2)=abs(ev(lpd,1+1)-ev(1,2+1))
                                      enddo
                                       do lpd=3,5
                                         if(abs(ev(lpd,1+1)-ev(1,2+1))==minval(B3)) goto 26
                                      enddo
26     continue
                         if(lpd==4)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 4s 4p (degeneracy)'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 4p(degeneracy)'
                          if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'  
                          occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                          occn(1+1,3)=2.d0*occn(0+1,7)!4p_1    
                         elseif(lpd==3)then
                          print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d(degeneracy)'
                          write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 3d(degeneracy)'
                          if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
                           occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                           occn(1+1,3)=2.d0*occn(0+1,7)!3d_1
                           occn(2+1,1)=2.d0*occn(0+1,7)!3d_2 
                         else
                         print*,'something to be corrscted here'
                         write(3,*)'something to be corrscted here'
                         stop
                         endif
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
                       endif
                    else!!Z>22
                         if(k3m_1==8 .and. k1m_2>8)then
                             if(Z.le. 28)then
                                 print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p'
                                 write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p'
                                 occn(0+1,7)=1.d0!5s
                                 if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'  
                                 occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p                   
                                 occn(1+1,3)=2.d0*occn(0+1,8)!4p_1 
                              elseif(k1m_2==9 .and. Z>28)then
                                 print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 3d'
                                 write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 3d'
                                 occn(0+1,7)=1.d0!5s
                                 occn(0+1,8)=1.d0!4p_0
                                 occn(1+1,3)=2.d0!4p_1
                                 if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'  
                                 occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0                   
                                 occn(1+1,4)=2.d0*occn(0+1,9)!3d_1
                                 occn(2+1,1)=2.d0*occn(0+1,9)!3d_2
                              elseif(k1m_2>9 .and. k4m_1==9 .and. Z>28 .and. Z.le. 34)then
                                 print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 5p'
                                 write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 5p'
                                 occn(0+1,7)=1.d0!5s
                                 occn(0+1,8)=1.d0!4p_0
                                 occn(1+1,3)=2.d0!4p_1
                                 if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'  
                                 occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!5p_0                   
                                 occn(1+1,4)=2.d0*occn(0+1,9)!5p_1
                              elseif(k1m_2>9 .and. k4m_1>9 .and. Z>28 .and. Z.le. 30)then
                                 print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 6s'
                                 write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s 4p 6s'
                                 occn(0+1,7)=1.d0!5s
                                 occn(0+1,8)=1.d0!4p_0
                                 occn(1+1,3)=2.d0!4p_1
                                 if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!'  
                                 occn(0+1,9)=(dble(Z)/2.d0-sum(occn))!6p_0  
                               else
                                print*, 'missing cases??' ,k3m_1,k4m_1,k1m_2
                              endif
                         elseif(k3m_1==8 .and. k1m_2==8)then
                           print*,'the configuration is 1s 2s 2p 3s 3p 4s 5s 3d'
                           write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 5s 3d'
                           occn(0+1,7)=1.d0!5s
                           if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'  
                           occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0                   
                           occn(1+1,3)=2.d0*occn(0+1,8)!3d_1
                           occn(2+1,1)=2.d0*occn(0+1,8)!3d_2
                         endif  
                    endif
                elseif(k3m_1==8 .and. k1m_2==7)then
                          if(Z.le.30)then
                           print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d'
                           write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 3d'
                           if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'  
                           occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0                   
                           occn(1+1,3)=2.d0*occn(0+1,7)!3d_1
                           occn(2+1,1)=2.d0*occn(0+1,7)!3d_2
                          else!Z>30 and Z.le. 36
                           print*,'the configuration is 1s 2s 2p 3s 3p 4s 3d 4p'
                           write(3,*)'the configuration is 1s 2s 2p 3s 3p 4s 3d 4p'  
                           occn(0+1,7)=1.d0!3d_0                   
                           occn(1+1,3)=2.d0!3d_1
                           occn(2+1,1)=2.d0!3d_2  
                           if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!'
                           occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                           occn(1+1,4)=2.d0*occn(0+1,8)!4p_1                         
                          endif             
                else
                  print*,'dont know',i,k2m_1,k3m_1,k4m_1,k1m_2
                endif            
            endif
          elseif(i==4 .and. k2m_1==5 .and. k1m_2==6 .and. Z.ge. 24) then
               print*,'the configuration is 1s 2s 3s 2p 3p 3d'
               write(3,*)'the configuration 1s 2s 3s 2p 3p 3d'
               occn(0+1,3)=1.d0!3s
               occn(0+1,4)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,5)=1.d0!3p_0
               occn(1+1,2)=2.d0!3p_1
               if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
               occn(0+1,6)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
               occn(1+1,3)=2.d0*occn(0+1,6)!3d_1
               occn(2+1,1)=2.d0*occn(0+1,6)!3d_2       
          elseif(i==3 .and. k2m_1==6 .and. k1m_2==7 .and. Z.ge. 25 .and. Z.le. 30) then
               print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d'
               write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d'
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!4s
               occn(0+1,6)=1.d0!3p_0
               occn(1+1,2)=2.d0!3p_1
               if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
               occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
               occn(1+1,3)=2.d0*occn(0+1,7)!3d_1
               occn(2+1,1)=2.d0*occn(0+1,7)!3d_2    
          elseif(i==3 .and. k2m_1==6 .and. k1m_2==7 .and. Z.ge. 31 .and. Z.le. 36) then
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!4s
               occn(0+1,6)=1.d0!3p_0
               occn(1+1,2)=2.d0!3p_1
               occn(0+1,7)=1.d0!3d_0
               occn(1+1,3)=2.d0!3d_1
               occn(2+1,1)=2.d0!3d_2 
              if(k4m_1==8)then
               print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 4p'
               write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 4p'
               if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!!!'
               occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
               occn(1+1,4)=2.d0*occn(0+1,8)!3p_1
              elseif(k4m_1==9)then
                   if(Z.le. 32)then
                      print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 5s'
                      write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 5s'
                      if((dble(Z)/2.d0-sum(occn))>1.d0)print*,'error!!!!'
                      occn(0+1,8)=(dble(Z)/2.d0-sum(occn))!5s
                   else!!Z>32 and Z<36
                      print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d 5s 4p'
                      write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d 5s 4p'
                      occn(0+1,8)=1.d0!5s
                      if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!!!'
                      occn(0+1,9)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
                      occn(1+1,4)=2.d0*occn(0+1,9)!4p_1
                   endif 
                else !!k4m_1>9
                   print*, 'define for k4m_1=', k4m_1
                   stop                   
              endif
          elseif(i==3 .and.k1m_2==8 .and. k2m_1==6 .and. k3m_1==8 .and. Z.ge.23.and. Z.le. 32)then
              print*,'the configuration is 1s 2s 2p 3s 4s 3p 5s 3d'
              write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 5s 3d'
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!4s
               occn(0+1,6)=1.d0!3p_0
               occn(1+1,2)=2.d0!3p_1
               occn(0+1,7)=1.d0!5s
               if((dble(Z)/2.d0-sum(occn))>5.d0)print*,'error!!!!'
               occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
               occn(1+1,3)=2.d0*occn(0+1,8)!3d_1
               occn(2+1,1)=2.d0*occn(0+1,8)!3d_2 
          elseif(i==3 .and.k1m_2==10 .and. k2m_1==7 .and. k3m_1==8 .and. Z.ge.23.and. Z.le. 28)then
              print*,'the configuration is 1s 2s 2p 3s 4s 5s 3p 4p'
              write(3,*)'the configuration 1s 2s 2p 3s 4s 5s 3p 4p'
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!4s
               occn(0+1,6)=1.d0!5s
               occn(0+1,7)=1.d0!3p_0
               occn(1+1,2)=2.d0!3p_1
               if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!!!'
               occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!3p_0
               occn(1+1,3)=2.d0*occn(0+1,8)!3p_1
          elseif(i==3 .and.k1m_2==5 .and. k2m_1==5 .and. k3m_1==6 .and. Z.ge. 23 .and. Z.le. 28)then
              print*,'the configuration is 1s 2s 2p 3s 3d 3p'
              write(3,*)'the configuration 1s 2s 2p 3s 3d 3p'
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!3d_0
               occn(1+1,2)=2.d0!3d_1
               occn(2+1,1)=2.d0!3d_2
               if((dble(Z)/2.d0-sum(occn))>3.d0)print*,'error!!!!'
               occn(0+1,6)=(dble(Z)/2.d0-sum(occn))/3.d0!3p_0
               occn(1+1,3)=2.d0*occn(0+1,6)!3p_1
          elseif(i==3 .and.k1m_2==7 .and. k2m_1==6 .and. k3m_1==7 .and. Z.ge. 21 .and. Z.le. 30)then
              print*,'the configuration is 1s 2s 2p 3s 4s 3p 3d'
              write(3,*)'the configuration 1s 2s 2p 3s 4s 3p 3d'
               occn(0+1,3)=1.d0!2p_0
               occn(1+1,1)=2.d0!2p_1  
               occn(0+1,4)=1.d0!3s
               occn(0+1,5)=1.d0!4s
               occn(0+1,6)=1.d0!3p_0               
               occn(1+1,2)=2.d0!3p_1
               if((dble(Z)/2.d0-sum(occn))>5.d0)print*,'error!!!!'
               occn(0+1,7)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
               occn(1+1,3)=2.d0*occn(0+1,7)!3d_1
               occn(2+1,1)=2.d0*occn(0+1,7)!3d_2
           elseif(i==3 .and. k2m_1==5 .and. k3m_1==6 .and. k1m_2==7 .and. k4m_1>7 .and.&
                 (Z.ge.23 .and. Z.le.28))then 
               !!!this is a wrong case bust faced in the degeneracy case when vxc=1   
               !! I fitted to be 1s 2s 2p 3s 3p 4s 3d but may be it is 1s 2s 2p 3s 3p 3d
                occn(0+1,3)=1.d0!2p_0
                occn(1+1,1)=2.d0!2p_1
                occn(0+1,4)=1.d0!3s
                occn(0+1,5)=1.d0!3p_0
                occn(1+1,2)=2.d0!3p_1
                  print*,'the configuration is 1s 2s 2p 3s 3p 3d (degeneracy 4s3d)',sum(occn),9
                  write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d (degeneracy 4s3d)'
                  if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!!!'
                  occn(0+1,6)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
                  occn(1+1,3)=2.d0*occn(0+1,6)!3d_1
                  occn(2+1,1)=2.d0*occn(0+1,6)!3d_2      
          else
            print*,i,k2m_1,k3m_1,k4m_1,k1m_2,'the configuration is not defined case'
          endif
      elseif(ii==1)then 
        occn(0+1,3)=1.d0!2p_0
        occn(1+1,1)=2.d0!2p_1
        occn(0+1,4)=1.d0!3s 
        occn(0+1,5)=1.d0!3p_0
        occn(1+1,2)=2.d0!3p_1 
         if(Z.le. 28)then
             print*,'the configuration is 1s 2s 2p 3s 3p 3d'
             write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d'
            if((dble(Z)/2.d0-sum(occn))>5)print*,'error!!'
              occn(0+1,6)=(dble(Z)/2.d0-sum(occn))/5.d0!3d_0
              occn(1+1,3)=2.d0*occn(0+1,6)!3d_1
              occn(2+1,1)=2.d0*occn(0+1,6)!3d_2
         elseif(Z>28 .and. Z.le. 30)then
             print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s'
             write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s'
              occn(0+1,6)=1.d0!3d_0
              occn(1+1,3)=2.d0!3d_1
              occn(2+1,1)=2.d0!3d_2    
            if((dble(Z)/2.d0-sum(occn))>1)print*,'error!!' 
              occn(0+1,7)=(dble(Z)/2.d0-sum(occn))!4s
          else!!Z>30
             print*,'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p'
             write(3,*)'the configuration is 1s 2s 2p 3s 3p 3d 4s 4p'
              occn(0+1,6)=1.d0!3d_0
              occn(1+1,3)=2.d0!3d_1
              occn(2+1,1)=2.d0!3d_2    
              occn(0+1,7)=1.d0!4s
              if((dble(Z)/2.d0-sum(occn))>3)print*,'error!!' 
              occn(0+1,8)=(dble(Z)/2.d0-sum(occn))/3.d0!4p_0
              occn(1+1,4)=2.d0*occn(0+1,8)!4p_1
         endif       
      endif

   else
   print*,'could not find the Z in the occupationnb file'
  endif
  if(maxval(occn)>2)then
  print*,maxval(occn), 'there is something wrong with the occupation number'
  stop
  endif
27   continue
  occn=2.d0*occn     
  print*,ii,'sum of occupation number=', sum(occn)

  if(abs(sum(occn)-Sm)>0.1 .and. ii.ne.1)then
    print*,'different occupation number!!',Sm,sum(occn)
    stop
  endif
end subroutine
