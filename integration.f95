!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Itegrations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function integrate(cha,Ng,xg,wg,h,hr1,hr2,k,i,j,eta,FormMat,tFormMat,ninter)
  implicit none
  character:: cha(1)
  integer:: Ng,i,j,nn,k,ninter
  logical:: LSAME,matrixB,matrixC,matrixK,vectorG
  double precision::h,FormMat(5,Ng),tFormMat(4,Ng),integrate,eta
  double precision, dimension(Ng):: xg,wg
  double precision, dimension(ninter)::hr1,hr2

  integrate=0.d0
  matrixB=LSAME(cha,'B')
  matrixC=LSAME(cha,'C')   
  matrixK=LSAME(cha,'K')
  vectorG=LSAME(cha,'G')
  

  if(matrixB) then 
  if(k==1) then
   do nn=1,Ng
   integrate=integrate+wg(nn)*tFormMat(i-1,nn)*tFormMat(j-1,nn)
   enddo
   integrate=integrate/h
  else
   do nn=1,Ng
      integrate=integrate+wg(nn)*FormMat(i,nn)*FormMat(j,nn)/((1.d0+xg(nn)*hr1(k))**2)
   enddo
      integrate=integrate*hr2(k)
   endif
  endif

 if(matrixC) then  
  if(k==1) then
   do nn=1,Ng
      integrate=integrate+wg(nn)*tFormMat(i-1,nn)*FormMat(j,nn)
   enddo
  else
   do nn=1,Ng
      integrate=integrate+wg(nn)*FormMat(i,nn)*FormMat(j,nn)/(1.d0+xg(nn)*hr1(k))
   enddo
      integrate=integrate*hr1(k)
   endif
 endif

  if(matrixK) then
    if(k.ne. 1) then
     do nn=1,Ng
        integrate=integrate+wg(nn)*FormMat(i,nn)*FormMat(j,nn)*&
                  exp(-eta*h*xg(nn))/(1.d0+xg(nn)*hr1(k))
     enddo
        integrate=hr1(k)*integrate
    else
     do nn=1,Ng
        integrate=integrate+wg(nn)*tFormMat(i-1,nn)*FormMat(j,nn)*&
                  exp(-eta*h*xg(nn))
     enddo
    endif
  endif



  if(vectorG) then
    do nn=1,Ng
       integrate=integrate+wg(nn)*FormMat(i,nn)*exp(-eta*h*xg(nn)) 
    enddo
       integrate=integrate*h 
  endif
 end function

