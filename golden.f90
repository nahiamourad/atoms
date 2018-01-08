SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc)
!Given a function FUNC(X), and given distinct initial points AX and
!BX, this routine searches in the downhill direction (defined by the
!function as evaluated at the initial points) and returns new points
!AX, BX, CX which bracket a minimum of the function. Also returned
!are the function values at the three points, FA, FB and FC.
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 Implicit none
 integer::N,nband,ninter,Ng,lrho,Na,Ls
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision::A,B,cxc
 double precision, dimension(lrho+1,N,N)::rho,rho0

 double precision::AX,BX,CX,DUM,FA,FB,FC,FU,R,Q,U,ULIM
 DOUBLE PRECISION, PARAMETER::GOLD=1.618034,GLIMIT=100.,TINY=1.D-20
!The first parameter is the default ratio by which successive intervals
!are magnified; the second is the maximum magnification allowed for
!a parabolic-fit step.
 call xcFUNC(AX,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FA)
 call xcFUNC(BX,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FB)
IF(FB.GT.FA) THEN
  DUM=AX
  AX=BX
  BX=DUM
  DUM=FB
  FB=FA
  FA=DUM
ENDIF
 CX=BX+GOLD*(BX-AX)
 CALL xcFUNC(CX,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FC)
1 IF(FB.GE.FC) THEN
  R=(BX-AX)*(FB-FC)
  Q=(BX-CX)*(FB-FA)
  U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
  ULIM=BX+GLIMIT*(CX-BX)
  IF((BX-U)*(U-CX).GT.0) THEN
    CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
    IF(FU.LT.FC) THEN
          AX=BX
          FA=FB
          BX=U
          FB=FU
          GOTO 1
    ELSE IF(FU.GT.FB) THEN
          CX=U
          FC=FU
          GOTO 1
    ENDIF
        U=CX+GOLD*(CX-BX)
        CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
  ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
    CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
        IF(FU.LT.FC) THEN
          BX=CX
          CX=U
          U=CX+GOLD*(CX-BX)
          FB=FC
          FC=FU
          CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
    ENDIF
  ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
    U=ULIM
          CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
  ELSE
    U=CX+GOLD*(CX-BX)
        CALL xcFUNC(U,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,FU)
  ENDIF
  AX=BX
  BX=CX
  CX=U
  FA=FB
  FB=FC
  FC=FU
  GOTO 1
ENDIF
RETURN
END

SUBROUTINE golden(ax,bx,cx,tol,xmin,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,GN) 
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!Implicit Integer (l-n)
 implicit none
 integer::N,nband,ninter,Ng,lrho,Na,Ls
 double precision, dimension(:)::wg(Ng),xg(Ng),wa(Na),hr1(ninter)
 double precision, dimension(:,:)::FormMat(5,Ng),tFormMat(4,Ng)
 double precision, dimension(:,:)::LegMat(2*Ls+1,Na),thr13(Ng,ninter)
 double precision::a,b,cxc
 double precision, dimension(lrho+1,N,N)::rho,rho0

!DOUBLE PRECISION golden,ax,bx,cx,tol,xmin,xcfunc,R,C,a,b 
!EXTERNAL xcfunc 
double precision, PARAMETER::R=.61803399,C=1.-R 
double precision::ax,bx,cx,tol,xmin,GN
DOUBLE PRECISION:: f1,f2,x0,x1,x2,x3 
x0=ax 
x3=cx 
if(abs(cx-bx)>abs(bx-ax)) then 
  x1=bx 
  x2=bx+C*(cx-bx) 
else 
  x2=bx 
  x1=bx-C*(bx-ax) 
endif 
 call xcfunc(x1,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,f1)!f(x1) 
 call xcfunc(x2,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,f2)!f(x2) 
do 
  if(abs(x3-x0)<=tol*(abs(x1)+abs(x2))) exit 
  if(f2<f1) then 
    x0=x1 
    x1=x2 
    x2=R*x1+C*x3 
    f1=f2 
    call xcfunc(x2,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,f2)!f(x2) 
  else 
    x3=x2 
    x2=x1 
    x1=R*x2+C*x0 
    f2=f1 
    call xcfunc(x1,a,b,rho,rho0,N,nband,ninter,Ng,lrho,Na,Ls,wg,xg,wa,hr1,&
                FormMat,tFormMat,LegMat,thr13,cxc,f1)!f(x1) 
  endif 
end do 
if(f1<f2) then 
  GN=f1 
  xmin=x1 
else 
  GN=f2 
  xmin=x2 
endif 
END SUBROUTINE golden
