MODULE TPA
IMPLICIT NONE
REAL(kind=8)::Rp !母核半径
REAL(kind=8)::VN,VC,Vl  !核势。库伦势，离心势
REAL(kind=8)::V0  !V0=192.42 V0=194.37 !势井深度
REAL(kind=8),PARAMETER::a0=0.5958  !核势扩散度 a0=0.5958 a0=0.5966
REAL(kind=8)::Ap,Zp,Np !母核质量数，质子数，中子数
REAL(kind=8)::Ad,Zd,Nd !子核质量数，质子数，中子数
REAL(kind=8)::Ma=4d0,Za=2d0  !alpha粒子的质量数和质子数
REAL(kind=8)::L !轨道角动量
REAL(kind=8)::u  !约化质量
REAL(kind=8)::Qa  !衰变能
REAL(kind=8)::GAma
REAL(kind=8),PARAMETER::ee=1.4399652                    !参数
REAL(kind=8),PARAMETER::h=4.135667E-21                  !普朗克常量Mev    C2_hpi2=38991.855
REAL(kind=8),PARAMETER::hbar=6.58211889E-22             !约化普朗克常量Mev
REAL(kind=8),PARAMETER::c=3.0E23                        !光速 fm/s
REAL(kind=8),PARAMETER::Pa=1.0 !预形成因子
REAL(kind=8)::R1,R2,R3									!!!!经典转折点R1、R2、R3
REAL(kind=8)::V1,V2,V3	
CONTAINS


REAL FUNCTION V(r) !!! considering the isospin  effect
IMPLICIT NONE
REAL(kind=8)::r
VN=-V0*(1+cosh(Rp/a0))/(cosh(r/a0)+cosh(Rp/a0))!Rp=1.28*Ap**(1/3)-0.76+0.8*Ap**(-1/3)
IF (r<Rp)THEN
    VC=(Zd*Za*ee*(3d0-(r/Rp)**2))/(2d0*Rp)
ELSE
    VC=(Zd*Za*ee)/r
ENDIF
Vl=((c*hbar)**2)*((L+0.5d0)**2)/(2d0*u*r*r)
!Vl=(L*(L+1)*(c*hbar)**2)/(2d0*u*r**2)
V=VN+VC+Vl
RETURN 
END FUNCTION V


REAL FUNCTION K(r)
IMPLICIT NONE 
REAL(kind=8)::r
K=SQRT(2d0*u*abs(Qa-V(r))/((c*hbar)**2))
RETURN
END FUNCTION K


!!!!!!!!!!!!!!!!!!!!!简约的辛普森法求穿透概率P!!!!!!!!!!!!!!!!!!!
REAL FUNCTION P_(X2,X3)
IMPLICIT NONE
REAL(kind=8)::X2,X3
REAL(kind=8)::sum,i,dr
INTEGER::N=1E5
dr=(X3-X2)/dfloat(N)
sum=0.0
DO i=X2,X3,dr
    sum=sum+K(i+0.5*dr)*dr*2d0/3d0+(K(i)+K(i+dr))*dr/6d0
ENDDO
P_=EXP(-2d0*sum)
RETURN
END FUNCTION P_


!!!!!!!!!!!!!!!!!!!!!简约的辛普森法求归一化因子F!!!!!!!!!!!!!!!!!!!
REAL FUNCTION F_(X1,X2)
IMPLICIT NONE
REAL(kind=8)::X1,X2
REAL(kind=8)::sum,i,dr
INTEGER::N=1E5
dr=(X2-X1)/dfloat(N)
sum=0.0
DO i=X1,X2,dr
    if(K(i)<1E-20.or.K(i+0.5*dr)<1E-20.or.K(i+dr)<1E-20)THEN
        sum=sum+0 
    else
        sum=sum+2d0*dr/(3d0*K(i+0.5*dr))+(dr/K(i)+dr/K(i+dr))/6d0
    endif
ENDDO
F_=2d0/sum
RETURN
END FUNCTION F_


!!!!!!!!!!!!!!!!!!!!!!!半衰期T1/2!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION Thalf_(Pa,F,P)
REAL(kind=8)::Pa,F,P
GAma=((c*hbar)*(c*hbar))*Pa*P*F/(4d0*u)
Thalf_=((hbar)*log(2d0))/GAma
!Thalf_=log10(((hbar)*log(2d0))/GAma)
RETURN
END FUNCTION Thalf_
END MODULE TPA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!主程序!!!!!!!!!!!!!!!!!!!!!!
PROGRAM main
USE TPA
IMPLICIT NONE
integer,allocatable::ApS(:),ZpS(:)
real(kind=8),allocatable::LS(:),QaS(:)
real(kind=8)::X1,X2,iii,VV
real(kind=8)::P,F,Thalf
integer::ii,line=535
allocate(ApS(line),ZpS(line),LS(line),QaS(line))
open(10,file='input535.txt')
do ii=1,line
    read(10,*)ApS(ii),ZpS(ii),LS(ii),QaS(ii)
enddo
!write(10,"('Ap=',F13.7/'Zp=',13.8/'L=',F13.9/)")ApS(ii),ZpS(ii),LS(ii),QaS(ii))
close(10)
open(11,file='output535.txt')
!open(12,file='Thalf.txt')
write(11,"(3X,'Ap',3X,'Zp',5X,'Qa',10X,'Thalf',12X)")
do ii=1,line
    Qa=QaS(ii)
    L=LS(ii)
    Ap=dfloat(ApS(ii))
    Zp=dfloat(ZpS(ii))
    Ad=Ap-4d0
    Zd=Zp-2d0
    Nd=Ad-Zd
    Np=Ap-Zp
    V0=192.42+31.059*(Nd-Zd)/Ad     !V0=194.37+31.059*(Nd-Zd)/Ad
    Rp=1.28*Ap**(1d0/3d0)-0.76+0.8*Ap**(-1d0/3d0)
    u=931.5d0*Ma*Ad/Ap
    write(*,"('母核质量数：Ap=',F13.7/'母核质子数：Zp=',F13.8/'势井深度：V0=',F13.7/'母核半径：Rp=',F13.9/'折合质量：u=',F18.9/,'衰变能：Qa=',F13.9/,'角动量： L=',F13.7)")Ap,Zp,V0,Rp,u,Qa,L
!!!!!!!!!!!!!!!!求转折点和极点!!!!!!!!!!!!!!!!!!!!!!!
   X1=0d0
   X2=5d0
   R1=(X1+X2)/2d0
   V1=V(R1)
   do while(abs(V1-Qa)>1E-6)
       if(V1>Qa)then
           X1=R1
       else
          X2=R1
       endif
         R1=(X1+X2)/2d0
         V1=V(R1)
   enddo
   write(*,"('第一个转折点：'/'R1=',F13.8)")R1
   write(*,"('第一个转折处的势能：'/'V1=',F13.8)")V1
   X1=5d0
   X2=15d0
   R2=(X1+X2)/2d0
   V2=V(R2)
   do while(abs(V2-Qa)>1E-6)
      if(V2>Qa)then
          X2=R2
      else
          X1=R2
      endif
          R2=(X1+X2)/2d0
          V2=V(R2)
   enddo
   write(*,"('第二个转折点：'/'R2=',F13.8)")R2
   write(*,"('第二个转折处的势能：'/'V2=',F13.8)")V2
   X1=15d0
   X2=120d0
   R3=(X1+X2)/2d0
   V3=V(R3)
   do while(abs(V3-Qa)>1E-6)
       if(V3>Qa)then
          X1=R3
       else
          X2=R3
       endif
          R3=(X1+X2)/2d0
          V3=V(R3)
   enddo
   write(*,"('第三个转折点：'/'R3=',F13.8)")R3
   write(*,"('第三个转折处的势能：'/'V3=',F13.8)")V3

   P=P_(R2,R3)
   write(*,*)'穿透概率:P=',P
   F=F_(R1,R2)
   write(*,*)'归一化因子:F=',F
   Thalf=Thalf_(Pa,F,P)
   write(*,*)'半衰期：T1/2=',Thalf
   
   write(11,"(2X,I3,2X,I3,2X,F13.8,2X,E15.6,2X)")&
   &IDINT(Ap),IDINT(Zp),Qa,Thalf   
enddo
 
close(11) 
END