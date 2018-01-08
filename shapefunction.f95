!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Shape functions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function formfunction(i,x)
  integer::i
  double precision::formfunction
  double precision::x
  formfunction=0.d0
  if(i>5) print*,'there is an error in calling formfunction'
  if(i==1)formfunction=1.d0-x
  if(i==2)formfunction=x
  if(i==3)formfunction=-4.d0*x**2+4.d0*x
  if(i==4)formfunction=-128.d0/3.d0*(x**4-9.d0/4.d0*x**3+13.d0/8.d0*x**2-3.d0/8.d0*x)
  if(i==5)formfunction=-128.d0/3.d0*(x**4-7.d0/4.d0*x**3+7.d0/8.d0*x**2-1.d0/8.d0*x)
 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function Lambda(i)
  integer:: Lambda,i
  if(i>4 .or. i<0) print*,'there is an error in calling the Lambda function'
  if(i==4)Lambda=1
  if(i==3)Lambda=3
  if(i==2)Lambda=4
  if(i==1)Lambda=5
  if(i==0)Lambda=2
 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function tformfunction(i,x)
  integer::i
  double precision::tformfunction,x
  tformfunction=0.d0
  if(i==1 .or. i>5) print*, 'there is something wrong see file shapefunction'
  if(i==2)tformfunction=1.d0
  if(i==3)tformfunction=-4.d0*x+4.d0
  if(i==4)tformfunction=-128.d0/3.d0*(x**3-9.d0/4.d0*x**2+13.d0/8.d0*x-3.d0/8.d0)
  if(i==5)tformfunction=-128.d0/3.d0*(x**3-7.d0/4.d0*x**2+7.d0/8.d0*x-1.d0/8.d0)
 end function


