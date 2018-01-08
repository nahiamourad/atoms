 function Legendre(l,t)
 integer:: l,j
 double precision:: Legendre,t,p1,p2,p3
     p1=1.d0
     p2=0.d0
 if(l>=1) then
    do j=1,l
       p3=p2
       p2=p1
       p1=((2.d0*dble(j)-1.d0)*t*p2-(dble(j)-1.d0)*p3)/dble(j)
    enddo
 endif
  Legendre=p1
 end function 
