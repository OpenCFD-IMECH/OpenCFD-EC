 ! ------------------------------------------------------------------------------------------
 ! 数值格式（重构）  (二维、三维通用) ; 用界面左右两侧多个网格点的值，构造界面上的左、右值
 !  UL, UR 界面（半点）的值; u(:) 界面左、右两侧各LAP个点的值, LAP=3时共6个点；界面位于1/2点；
 ! 2014-3-14: 增加了 MUSCL with k=-1， 去掉了WENO3 
 !---正通量重构 (左值)-------------------------------------------------------------------------
    subroutine scheme_fP(uL,u,Iflag_Scheme)
     use   const_var
     implicit none
     real(PRE_EC):: uL,u(1-LAP:LAP),minmod
     real(PRE_EC):: IS1,IS2,a1,a2,w1,w2,r1,r2,f1,f
     real(PRE_EC),parameter:: k2c=1.d0,k2u=-1.d0,b2=2.d0
     real(PRE_EC),parameter::k3=1.d0/3.d0,ep=1.d-6
	 real(PRE_EC),parameter:: eta=0.05    ! OMUSCL2中的耗散系数 >0, 值越大，耗散越大
     real(PRE_EC):: up,um,s
     integer:: Iflag_Scheme

!    UL为i=1/2点	 
   Select case (Iflag_Scheme)
     case (Scheme_CD2)     ! 2阶中心 （第1个边界网格使用）
      UL=0.5d0*(u(0)+u(1))
	 case (Scheme_UD1)          ! 1阶迎风 
        UL=u(0)
     case ( Scheme_UD3 )        ! 3rd-order-Upwind
        UL=(-u(-1)+5.d0*u(0)+2.d0*u(1)  )/6.d0        
     case(Scheme_NND2 )         ! 2nd-order-NND
        UL=u(0)+0.5d0*minmod(u(0)-u(-1),u(1)-u(0))            
     case (Scheme_MUSCL2U)        ! MUSCL2 with k=-1 
       UL=u(0)+0.25d0*((1.d0-k2u)*minmod(u(0)-u(-1),b2*(u(1)-u(0)))+(1.d0+k2u)*minmod(u(1)-u(0),b2*(u(0)-u(-1)))) ! MUSCL2, k=-1
     case (Scheme_MUSCL2C ) 
       UL=u(0)+0.25d0*((1.d0-k2c)*minmod(u(0)-u(-1),b2*(u(1)-u(0)))+(1.d0+k2c)*minmod(u(1)-u(0),b2*(u(0)-u(-1)))) ! MUSCL2, k=1
     case (Scheme_MUSCL3 )      ! 3rd-order MUSCL (with Van Albada limiter)
       up=u(1)-u(0); um=u(0)-u(-1)                           ! 1阶 前差、后差
       s=(2.d0*up*um+ep)/(up*up+um*um+ep)                    ! Van Albada限制器 （光滑区，前差与后差接近，该值接近1）
       UL=u(0)+0.25d0*s*((1.d0-k3*s)*um+(1.d0+k3*s)*up)      ! 3阶MUSCL (光滑区逼近3阶迎风)

     case (Scheme_OMUSCL2 )          ! 2阶优化的MUSCL 格式 (Developed by Leng Yan)
         r1=(u(0)-u(-1)+ep)/(u(1)-u(0)+ep); r2=(u(1)-u(0)+ep)/(u(2)-u(1)+ep)
!		 f1=0.8d0-0.175d0/r2+0.375d0*r1
!        f1=0.95-0.25d0/r2+0.3d0*r1        !eita=0.05
         f1=(1.d0-eta)+(eta-0.55d0)*0.5d0/r2 + (eta+0.55d0)*0.5d0*r1
		 f=max(0.d0,min(2.d0,f1,2.d0*r1))
         UL=u(0)+0.5d0*f*(u(1)-u(0))                               

     case (Scheme_WENO5)
        call weno5_L(u,UL)
     case (Scheme_WENO7)
        call weno7_L(u,UL)
	
	 case (Scheme_UD5)        ! 5th-order upwind
        UL=(2.d0*u(-2)-13.d0*u(-1)+47.d0*u(0)+27.d0*u(1)-3.d0*u(2))/60.d0
     end Select 

    end
!---负通量重构----------------------------------------------------
  subroutine scheme_fm(uR,u,Iflag_Scheme)
     use   const_var
     implicit none
     real(PRE_EC):: uR,u(1-LAP:LAP),minmod
     real(PRE_EC):: IS1,IS2,a1,a2,w1,w2,r1,r2,f1,f
     real(PRE_EC),parameter:: k2c=1.d0,k2u=-1.d0,b2=2.d0
     real(PRE_EC),parameter::k3=1.d0/3.d0,ep=1.d-6
	 real(PRE_EC),parameter:: eta=0.05    ! OMUSCL2中的耗散系数 >0, 值越大，耗散越大
     real(PRE_EC):: up,um,s
     integer:: Iflag_Scheme
   Select case (Iflag_Scheme)
     case (Scheme_CD2)     ! 2阶中心 （第1个边界网格使用）
      UR=0.5d0*(u(0)+u(1))
	 case (Scheme_UD1) 
	   UR=u(1)                                                   ! 1阶迎风
     case (Scheme_UD3) 
       UR=(2.d0*u(0)+5.d0*u(1)-u(2)  )/6.d0                          !UD 3nd order
     case (Scheme_NND2) 
       UR=u(1)-0.5d0*minmod(u(1)-u(0),u(2)-u(1))                         ! NND 2nd order
     case (Scheme_MUSCL2U)
       uR=u(1)-0.25d0*((1.d0-k2u)*minmod(u(2)-u(1),b2*(u(1)-u(0)))+(1.d0+k2u)*minmod(u(1)-u(0),b2*(u(2)-u(1))))   ! MUSCL
     case (Scheme_MUSCL2C) 
       uR=u(1)-0.25d0*((1.d0-k2c)*minmod(u(2)-u(1),b2*(u(1)-u(0)))+(1.d0+k2c)*minmod(u(1)-u(0),b2*(u(2)-u(1))))   ! MUSCL
     case ( Scheme_MUSCL3 )             ! 3阶MUSCL (Van Albada限制器)
       up=u(2)-u(1) ; um=u(1)-u(0)                                     !前差、后差
       s=(2.d0*up*um+ep)/(up*up+um*um+ep)
       UR=u(1)-0.25d0*s*((1.d0-k3*s)*up+(1.d0+k3*s)*um)
    case (Scheme_OMUSCL2) 
        r1=(u(2)-u(1)+ep)/(u(1)-u(0)+ep); r2=(u(1)-u(0)+ep)/(u(0)-u(-1)+ep)
!	 	 f1=0.8d0-0.175d0/r2+0.375d0*r1
         f1=(1.d0-eta)+(eta-0.55d0)*0.5d0/r2 + (eta+0.55d0)*0.5d0*r1
	 	 f=max(0.d0,min(2.d0,f1,2.d0*r1))
         UR=u(1)-0.5d0*f*(u(1)-u(0))                               ! 2阶优化的MUSCL: OMUSCL2

    case (Scheme_WENO5)
       call weno5_R(u,UR)
    case (Scheme_WENO7)
       call weno7_R(u,UR)
    case (Scheme_UD5) 
        UR=(-3.d0*u(-1)+27.d0*u(0)+47.d0*u(1)-13.d0*u(2)+2.d0*u(3))/60.d0
   end Select
  end


!-------------------------------------------------------------------------------------------------------   
! 5阶WENO格式 (左值) 
     subroutine weno5_L(u,UL)
     use  const_var
     implicit none
     real(PRE_EC)::  u(1-LAP:LAP), UL
     real(PRE_EC)::  S0,S1,S2,a0,a1,a2,am,q03,q13,q23
     real(PRE_EC),parameter::  ep=1.d-6, C03=3.d0/10.d0, C13=3.d0/5.d0, C23=1.d0/10.d0
         
		 S0=13.d0/12.d0*(u(0)-2.d0*u(1)+u(2))**2+  1.d0/4.d0*(3.d0*u(0)-4.d0*u(1)+u(2))**2
         S1=13.d0/12.d0*(u(-1)-2.d0*u(0)+u(1))**2+  1.d0/4.d0*(u(-1)-u(1))**2
         S2=13.d0/12.d0*(u(-2)-2.d0*u(-1)+u(0))**2+  1.d0/4.d0*(u(-2)-4.d0*u(-1)+3.d0*u(0))**2
         a0=C03/((ep+S0)**2)
         a1=C13/((ep+S1)**2)
         a2=C23/((ep+S2)**2)
         am=a0+a1+a2
         q03=1.d0/3.d0*u(0)+5.d0/6.d0*u(1)-1.d0/6.d0*u(2)
         q13=-1.d0/6.d0*u(-1)+5.d0/6.d0*u(0)+1.d0/3.d0*u(1)
         q23=1.d0/3.d0*u(-2)-7.d0/6.d0*u(-1)+11.d0/6.d0*u(0)
         UL=(a0*q03+a1*q13+a2*q23)/am
     end
!---------------------------------------------
! 5阶WENO格式 (右值) 
     subroutine weno5_R(u,UR)
     use  const_var
     implicit none
     real(PRE_EC)::  u(1-LAP:LAP), UR
     real(PRE_EC):: S0,S1,S2,a0,a1,a2,am,q03,q13,q23
     real(PRE_EC),parameter::  ep=1.d-6, C03=3.d0/10.d0, C13=3.d0/5.d0, C23=1.d0/10.d0

       S0=13.d0/12.d0*(u(1)-2.d0*u(0)+u(-1))**2+  1.d0/4.d0*(3.d0*u(1)-4.d0*u(0)+u(-1))**2
       S1=13.d0/12.d0*(u(2)-2.d0*u(1)+u(0))**2+  1.d0/4.d0*(u(2)-u(0))**2
       S2=13.d0/12.d0*(u(3)-2.d0*u(2)+u(1))**2+  1.d0/4.d0*(u(3)-4.d0*u(2)+3.d0*u(1))**2
       a0=C03/((ep+S0)**2)
       a1=C13/((ep+S1)**2)
       a2=C23/((ep+S2)**2)
       am=a0+a1+a2
       q03=1.d0/3.d0*u(1)+5.d0/6.d0*u(0)-1.d0/6.d0*u(-1)
       q13=-1.d0/6.d0*u(2)+5.d0/6.d0*u(1)+1.d0/3.d0*u(0)
       q23=1.d0/3.d0*u(3)-7.d0/6.d0*u(2)+11.d0/6.d0*u(1)
       UR=(a0*q03+a1*q13+a2*q23)/am
     end
!----------------------------------------------------------------------

! 7阶WENO格式 (左值) UL=u(1/2) 
     subroutine weno7_L(v,UL)
     use  const_var
     implicit none
     real(PRE_EC)::  v(1-LAP:LAP), UL
      
      real(PRE_EC)::  S0,S1,S2,S3, s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
             a0,a1,a2,a3,am,q0,q1,q2,q3
      real(PRE_EC),parameter::                             &
         C0=1.d0/35.d0, C1=12.d0/35.d0, C2=18.d0/35.d0,  C3=4.d0/35.d0, &
         a11=-2.d0/6.d0,a12=9.d0/6.d0,a13=-18.d0/6.d0,a14=11.d0/6.d0, &
         a21=1.d0/6.d0,               a23=3.d0/6.d0,a24=2.d0/6.d0, &
         a31=-2.d0/6.d0,a32=-3.d0/6.d0,            a34=-1.d0/6.d0,   &
         a41=-11.d0/6.d0,a42=18.d0/6.d0,a43=-9.d0/6.d0,a44=2.d0/6.d0,  &
         b12=4.d0,b13=-5.d0,b14=2.d0,  b22= -2.d0,      &   
         b41=2.d0,b42=-5.d0,b43=4.d0,    c12=3.d0, &
         d12=13.d0/12.d0,d13=1043.d0/960.d0,d14=1.d0/12.d0
	  real(PRE_EC),parameter::                                    &
         e11=-3.d0/12.d0, e12=13.d0/12.d0, e13=-23.d0/12.d0, e14=25.d0/12.d0,  &
         e21=1.d0/12.d0, e22=-5.d0/12.d0, e23=13.d0/12.d0,  e24=3.d0/12.d0,    &
         e31=-1.d0/12.d0, e32=7.d0/12.d0, e33=7.d0/12.d0,    e34=-1.d0/12.d0,  &  
         e41=3.d0/12.d0,  e42=13.d0/12.d0, e43=-5.d0/12.d0,  e44=1.d0/12.d0

      real(PRE_EC),parameter:: ep=1.d-8    !! WENO-JS

! 1  阶导数  
         S10=a11*v(-3)+a12*v(-2)+a13*v(-1) +a14*v(0)
         S11=a21*v(-2) -   v(-1)+a23*v(0)   +a24*v(1)
         S12=a31*v(-1)+a32*v(0)  +    v(1) +a34*v(2)
         S13=a41*v(0)  +a42*v(1)+a43*v(2) +a44*v(3)
 ! 2 阶导数
         S20=-v(-3)+b12*v(-2)+b13*v(-1)+b14*v(0)             
         S21=             v(-1)+b22*v(0)  +v(1)         
         S22=             v(0)  +b22*v(1)+v(2)         
         S23=b41*v(0)+b42*v(1)+b43*v(2)-v(3)         
! 3 阶导数
         S30=-v(-3)+c12*(v(-2)-v(-1)) +v(0)                                   
         S31=-v(-2)+c12*(v(-1)-v(0))   +v(1)                 
         S32=-v(-1)+c12*(v(0)-v(1))   +v(2)                 
         S33=-v(0)  +c12*(v(1)-v(2)) +v(3)                 

       S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30
       S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31
       S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32
       S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33

!-------WENO J-S----------------------
       a0=C0/((ep+S0)**2)
       a1=C1/((ep+S1)**2)
       a2=C2/((ep+S2)**2)
       a3=C3/((ep+S3)**2)

!-----------------------------------------------
     am=a0+a1+a2+a3

!  4阶差分格式的通量
     q0=e11*v(-3)+e12*v(-2)+e13*v(-1) +e14*v(0)
     q1=e21*v(-2)+e22*v(-1)+e23*v(0)   +e24*v(1)
     q2=e31*v(-1)+e32*v(0)  +e33*v(1) +e34*v(2)
     q3=e41*v(0)  +e42*v(1)+e43*v(2) +e44*v(3)

!  由4个4阶差分格式组合成1个7阶差分格式
     UL=(a0*q0+a1*q1+a2*q2+a3*q3)/am
   end   




! 7阶WENO格式 (右值)  UR=U(1/2)
     subroutine weno7_R(v,UR)
     use  const_var
     implicit none
     real(PRE_EC)::  v(1-LAP:LAP), UR

      real(PRE_EC)::  S0,S1,S2,S3, s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
             a0,a1,a2,a3,am,q0,q1,q2,q3
      real(PRE_EC),parameter::                             &
         C0=1.d0/35.d0, C1=12.d0/35.d0, C2=18.d0/35.d0,  C3=4.d0/35.d0, &
         a11=-2.d0/6.d0,a12=9.d0/6.d0,a13=-18.d0/6.d0,a14=11.d0/6.d0, &
         a21=1.d0/6.d0,               a23=3.d0/6.d0,a24=2.d0/6.d0, &
         a31=-2.d0/6.d0,a32=-3.d0/6.d0,            a34=-1.d0/6.d0,   &
         a41=-11.d0/6.d0,a42=18.d0/6.d0,a43=-9.d0/6.d0,a44=2.d0/6.d0,  &
         b12=4.d0,b13=-5.d0,b14=2.d0,  b22= -2.d0,      &   
         b41=2.d0,b42=-5.d0,b43=4.d0,    c12=3.d0, &
         d12=13.d0/12.d0,d13=1043.d0/960.d0,d14=1.d0/12.d0
	  real(PRE_EC),parameter::                                    &
         e11=-3.d0/12.d0, e12=13.d0/12.d0, e13=-23.d0/12.d0, e14=25.d0/12.d0,  &
         e21=1.d0/12.d0, e22=-5.d0/12.d0, e23=13.d0/12.d0,  e24=3.d0/12.d0,    &
         e31=-1.d0/12.d0, e32=7.d0/12.d0, e33=7.d0/12.d0,    e34=-1.d0/12.d0,  &  
         e41=3.d0/12.d0,  e42=13.d0/12.d0, e43=-5.d0/12.d0,  e44=1.d0/12.d0
      real(PRE_EC),parameter:: ep=1.d-8 

 
    

!      7th order WENO scheme
! 1  阶导数
         S10=a11*v(4)+a12*v(3)+a13*v(2)  +a14*v(1)
         S11=a21*v(3)-    v(2) +a23*v(1)    +a24*v(0)
         S12=a31*v(2)+a32*v(1)   +    v(0)  +a34*v(-1)
         S13=a41*v(1)  +a42*v(0)+a43*v(-1)  +a44*v(-2)
! 2 阶导数
         S20=-v(4)+b12*v(3)+b13*v(2)+b14*v(1)              
         S21=             v(2) +b22*v(1)  +v(0)         
         S22=             v(1)   +b22*v(0)+v(-1)         
         S23=b41*v(1)+b42*v(0)+b43*v(-1)-v(-2)         
! 3 阶导数 
         S30=-v(4)+c12*(v(3)-v(2))+v(1)                                  
         S31=-v(3)+c12*(v(2)-v(1))+  v(0)                 
         S32=-v(2)+c12*(v(1)  -v(0))+v(-1)                 
         S33=-v(1)+  c12*(v(0)-v(-1))+v(-2)                 

       S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30
       S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31
       S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32
       S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33

        a0=C0/((ep+S0)**2)
        a1=C1/((ep+S1)**2)
        a2=C2/((ep+S2)**2)
        a3=C3/((ep+S3)**2)

!-----------------------------------------------

     am=a0+a1+a2+a3

!  4阶差分格式的通量
     q0=e11*v(4)+e12*v(3)+e13*v(2)+e14*v(1)
     q1=e21*v(3)+e22*v(2)+e23*v(1)  +e24*v(0)
     q2=e31*v(2)+e32*v(1)  +e33*v(0)+e34*v(-1)
     q3=e41*v(1)+  e42*v(0)+e43*v(-1)+e44*v(-2)

!  由4个4阶差分格式组合成1个7阶差分格式
     UR=(a0*q0+a1*q1+a2*q2+a3*q3)/am
    end


!------------minmod 函数------------------------------------------
      function minmod(a,b)
	  use precision_EC
      implicit none
      real(PRE_EC) a,b,minmod
      if(a*b .le. 0.d0) then
       minmod=0.d0
      else if (abs(a) .le. abs(b)) then
       minmod=a
      else
       minmod=b
      endif
      end function minmod
!--------------------------------------------------------------------------------
