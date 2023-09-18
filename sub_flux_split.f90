!--OpenCFD EC Ver 0.6d----------------------------------------------------------------
! 通量分裂模块，包含HLL/HLLC, Roe, Steger-Warming, Van Leer, AUSM+ 等多种方法
! Developed by Li Xinliang and Leng Yan
!--------------------------------------------------------------------------------------
!  Flux by using  HLL or  HLLC approximate Riemann slover
!  Copyright by Li Xinliang, Institute of Mechanics, CAS; lixl@imech.ac.cn
!  Revised by Leng Yan (In HLLC, the second kind method in Toro's book is used; the first method is unstable)
!  Ref. Toro : Riemann Solvers and Numerical Methods for Fluid Dynamics,  p330       
! 2011-3-21: A bug in Flux_Roe_1D() is removed
! 2011-10-8: A bug in AUSM-PW() is removed
!--------------------------------------------------

    subroutine Flux_HLL_HLLC_1D(QL,QR,Flux,gamma,Iflag_flux)   
     use const_var
     implicit none
     integer:: Iflag_flux,IFlag_Reconstruction
     real(PRE_EC):: QL(5),QR(5),UL(5),UR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar    ! uu,vv,ww velocity
     real(PRE_EC):: p_pvrs,p_star,qql,qqr,Sl,Sr,Fl(5),Fr(5),S_star,tmpl,tmpr,tmp,F_starl(5),F_starr(5)
      
     dl=QL(1); uul=QL(2); vvl=QL(3);wwl=QL(4); pl=QL(5)
     dr=QR(1); uur=QR(2); vvr=QR(3);wwr=QR(4); pr=QR(5)
     UL(1)=dl; UL(2)=dl*uul; UL(3)=dl*vvl; UL(4)=dl*wwl; UL(5)=pl/(gamma-1.d0)+dl*(uul*uul+vvl*vvl+wwl*wwl)*0.5d0 
     UR(1)=dr; UR(2)=dr*uur; UR(3)=dr*vvr; UR(4)=dr*wwr; UR(5)=pr/(gamma-1.d0)+dr*(uur*uur+vvr*vvr+wwr*wwr)*0.5d0 
     al=sqrt(gamma*pl/dl); ar=sqrt(gamma*pr/dr)
     p_pvrs=0.5d0*(pl+pr)-(uur-uul)*(dl+dr)*(al+ar)*0.125d0
     p_star=max(0.d0,p_pvrs)  
    
     if(p_star .le. pl) then
       qql=1
     else
       qql=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pl-1.d0) )
     endif
     if(p_star .le. pr) then
       qqr=1
     else
       qqr=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pr-1.d0) )
     endif
     Sl=uul-al*qql; Sr=uur+ar*qqr    ! speed of the lift and the right shockwaves
       
     Fl(1)=UL(2); Fl(2)=UL(2)*uul+pl; Fl(3)=UL(3)*uul; Fl(4)=UL(4)*uul; Fl(5)=uul*(UL(5)+pl)
     Fr(1)=UR(2); Fr(2)=UR(2)*uur+pr; Fr(3)=UR(3)*uur; Fr(4)=UR(4)*uur; Fr(5)=uur*(UR(5)+pr)
!----HLL---------------------------------
     
     if(Iflag_flux .eq. Flux_HLL) then  ! HLL Flux
       if( Sl .ge. 0 ) then
         Flux=Fl
       else if (Sr .le. 0) then
         Flux=Fr
       else
         Flux=(Sr*Fl-Sl*Fr+Sl*Sr*(UR-UL))/(Sr-Sl)
       endif
     else  
!------------HLLC Flux------------------------------------------------------------- 
       S_star=(pr-pl+dl*uul*(Sl-uul)-dr*uur*(Sr-uur))/(dl*(Sl-uul)-dr*(Sr-uur))
       tmpl=dl*(Sl-uul)/(Sl-S_star) ; tmpr=dr*(Sr-uur)/(Sr-S_star)
       F_starl(1)=Fl(1)+Sl*(tmpl-UL(1))
       F_starl(2)=Fl(2)+Sl*(tmpl*S_star-UL(2))
       F_starl(3)=Fl(3)+Sl*(tmpl*vvl-UL(3))
       F_starl(4)=Fl(4)+Sl*(tmpl*wwl-UL(4))
       tmp=UL(5)/dl+(S_star-uul)*(S_star+pl/(dl*(Sl-uul)))
       F_starl(5)=Fl(5)+Sl*(tmpl*tmp-UL(5))

       F_starr(1)=Fr(1)+Sr*(tmpr-UR(1))
       F_starr(2)=Fr(2)+Sr*(tmpr*S_star-UR(2))
       F_starr(3)=Fr(3)+Sr*(tmpr*vvr-UR(3))
       F_starr(4)=Fr(4)+Sr*(tmpr*vvr-UR(4))
       tmp=UR(5)/dr+(S_star-uur)*(S_star+pr/(dr*(Sr-uur)))
       F_starr(5)=Fr(5)+Sr*(tmpr*tmp-UR(5))

! Revised by Leng Yan  
!     S_star=(pr-pl+dl*uul*(Sl-uul)-dr*uur*(Sr-uur))/(dl*(Sl-uul)-dr*(Sr-uur))
!     F_starl(1)=S_star*(Sl*UL(1)-Fl(1))/(Sl-S_star)
!     F_starl(2)=(S_star*(Sl*UL(2)-Fl(2))+Sl*(pl+dl*(Sl-uul)*(S_star-uul)))/(Sl-S_star)
!     F_starl(3)=S_star*(Sl*UL(3)-Fl(3))/(Sl-S_star)
!     F_starl(4)=S_star*(Sl*UL(4)-Fl(4))/(Sl-S_star)
!     F_starl(5)=(S_star*(Sl*UL(5)-Fl(5))+Sl*(pl+dl*(Sl-uul)*(S_star-uul))*S_star)/(Sl-S_star)
    
!     F_starr(1)=S_star*(Sr*UR(1)-Fr(1))/(Sr-S_star)
!     F_starr(2)=(S_star*(Sr*UR(2)-Fr(2))+Sr*(pr+dr*(Sr-uur)*(S_star-uur)))/(Sr-S_star)
!     F_starr(3)=S_star*(Sr*UR(3)-Fr(3))/(Sr-S_star)
!     F_starr(4)=S_star*(Sr*UR(4)-Fr(4))/(Sr-S_star)
!     F_starr(5)=(S_star*(Sr*UR(5)-Fr(5))+Sr*(pr+dr*(Sr-uur)*(S_star-uur))*S_star)/(Sr-S_star)
   
       if( Sl .ge. 0 ) then
         Flux=Fl
       else if (Sr .le. 0) then
         Flux=Fr
       else if (S_star .ge. 0) then
         Flux=F_starl
       else
         Flux=F_starr
       endif
     endif
    end

!c========================================================
! Developed By Li Xinliang
!-------Extend 1D steger_warming FLux for finite volume method
! F=F+ (uL) + F- (uR)	

    subroutine Flux_steger_warming_1Da(QL,QR,Flux,gamma)   
     use const_var
     implicit none
     integer:: IFlag_Reconstruction
     real(PRE_EC):: QL(5),QR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar ,pr1   ! uu velocity
     real(PRE_EC):: tmp0,tmp1,tmp2,tmp3,E1P,E2P,E3P,E1M,E2M,E3M,fp(5),fm(5)

 
     dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5) 
     dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)	
	 al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
     ar=sqrt(gamma*pr/dr)  ! find a bug, removed
      
	 tmp1=2.d0*(gamma-1.d0)
     tmp3=(3.d0-gamma)/(2.d0*(gamma-1.d0)) 
! eigenvalues---------	      
     E1P=(uul+abs(uul))*0.5d0
     E2P=(uul-al+abs(uul-al))*0.5d0
     E3P=(uul+al+abs(uul+al))*0.5d0
     tmp0=dl/(2.d0*gamma) 
     fp(1)=tmp0*(tmp1*E1P+E2P+E3P)
     fp(2)=tmp0*(tmp1*E1P*uul+E2P*(uul-al)+E3P*(uul+al))
     fp(3)=tmp0*(tmp1*E1P*vvl+E2P*vvl+E3P*vvl)
     fp(4)=tmp0*(tmp1*E1P*wwl+E2P*wwl+E3P*wwl)
     fp(5)=tmp0*(E1P*(gamma-1.d0)*(uul*uul+vvl*vvl+wwl*wwl)+E2P*((uul-al)**2+vvl*vvl+wwl*wwl)*0.5d0    &
	       +E3P*((uul+al)**2+vvl*vvl+wwl*wwl)*0.5d0+tmp3*al*al*(E2P+E3P))
     
     E1M=(uur-abs(uur))*0.5d0
     E2M=(uur-ar-abs(uur-ar))*0.5d0
     E3M=(uur+ar-abs(uur+ar))*0.5d0
     tmp0=dr/(2.d0*gamma) 
     fm(1)=tmp0*(tmp1*E1M+E2M+E3M)
     fm(2)=tmp0*(tmp1*E1M*uur+E2M*(uur-ar)+E3M*(uur+ar))
     fm(3)=tmp0*(tmp1*E1M*vvr+E2M*vvr+E3M*vvr)
     fm(4)=tmp0*(tmp1*E1M*wwr+E2M*wwr+E3M*wwr)
     fm(5)=tmp0*(E1M*(gamma-1.d0)*(uur*uur+vvr*vvr+wwr*wwr)+E2M*((uur-ar)**2+vvr*vvr+wwr*wwr)*0.5d0    &
           +E3M*((uur+ar)**2+vvr*vvr+wwr*wwr)*0.5d0+tmp3*ar*ar*(E2M+E3M))
     Flux=fp+fm
    end
!-----------------------------------------------------------
!-----------------------------------------------------------
! Flux by using Roe approximate Riemann slover
! 代码由冷岩开发
! A bug is removed, 2011-5-6
! A bug is removed, 2015-8-31
    subroutine Flux_Roe_1D(QL,QR,Flux,gamma)  
     use const_var 
     implicit none
     real(PRE_EC):: QL(5),QR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar    ! uu velocity
     real(PRE_EC):: avd,avu,avv,avw,avh,ava
     real(PRE_EC):: lamda(3),UL(5),UR(5),Fl(5),Fr(5)
     real(PRE_EC):: arr(5)
     real(PRE_EC):: hl,hr,tmpp,delt=0.1d0
     dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5)
     dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
     al=sqrt(gamma*pl/dl); ar=sqrt(gamma*pr/dr)
     hl=gamma*pl/((gamma-1.d0)*dl)+0.5*(uul*uul+vvl*vvl+wwl*wwl)
     hr=gamma*pr/((gamma-1.d0)*dr)+0.5*(uur*uur+vvr*vvr+wwr*wwr)   

     UL(1)=dl; UL(2)=dl*uul; UL(3)=dl*vvl; UL(4)=dl*wwl; UL(5)=pl/(gamma-1.d0)+dl*(uul*uul+vvl*vvl+wwl*wwl)*0.5d0 
     UR(1)=dr; UR(2)=dr*uur; UR(3)=dr*vvr; UR(4)=dr*wwr; UR(5)=pr/(gamma-1.d0)+dr*(uur*uur+vvr*vvr+wwr*wwr)*0.5d0 
     Fl(1)=UL(2); Fl(2)=UL(2)*uul+pl; Fl(3)=UL(3)*uul; Fl(4)=UL(4)*uul; Fl(5)=uul*(UL(5)+pl)
     Fr(1)=UR(2); Fr(2)=UR(2)*uur+pr; Fr(3)=UR(3)*uur; Fr(4)=UR(4)*uur; Fr(5)=uur*(UR(5)+pr)
!----Roe平均---------------------
     tmpp=sqrt(dr/dl)
     avd=0.25d0*dl*(1.d0+tmpp)*(1.d0+tmpp)
     avu=(uul+tmpp*uur)/(1.d0+tmpp)
     avv=(vvl+tmpp*vvr)/(1.d0+tmpp)
     avw=(wwl+tmpp*wwr)/(1.d0+tmpp)
     avh=(hl+tmpp*hr)/(1.d0+tmpp)
     ava=sqrt((gamma-1.)*(avh-0.5d0*(avu*avu+avv*avv+avw*avw)))

!  A bug is removed ------- 2011-5-6
! modified by Li Xinliang 2011-5-6
     lamda(1)=abs(avu-ava)
     lamda(2)=abs(avu)
     lamda(3)=abs(avu+ava)

! Harten型熵修正 
     if(lamda(1) < delt) then
       lamda(1)=( lamda(1)**2 +delt*delt)/(2.d0*delt)
     end if
     if(lamda(2) < delt) then
       lamda(2)=( lamda(2)**2 +delt*delt)/(2.d0*delt)
     end if
     if(lamda(3) < delt) then
       lamda(3)=( lamda(3)**2 +delt*delt)/(2.d0*delt)
     end if


!----------------------------------------------
     arr(2)=(UR(3)-UL(3))/(avv+1.e-8)-(dr-dl)
!     arr(3)=(UR(4)-UL(4))/(avw+1.e-8)-(UR(2)-UL(2))   ! Bug Bug Bug ???
      arr(3)=(UR(4)-UL(4))/(avw+1.e-8)-(dr-dl)   ! Bug Bug Bug ???

     arr(4)=(gamma-1.d0)*((dr-dl)*(avh-avu*avu-avv*avv-avw*avw)+avu*(UR(2)-UL(2))+avv*(UR(3)-UL(3)) & 
	               +avw*(UR(4)-UL(4))-(UR(5)-UL(5)))/(ava*ava)
     arr(1)=((avu+ava)*(dr-dl)-UR(2)+UL(2)-ava*arr(4))/(2.d0*ava)
     arr(5)=dr-dl-(arr(1)+arr(4))
     flux(1)=0.5d0*(Fl(1)+Fr(1))-0.5d0*(lamda(1)*arr(1)+lamda(2)*arr(4)+lamda(3)*arr(5))
     flux(2)=0.5d0*(Fl(2)+Fr(2))-0.5d0*(lamda(1)*arr(1)*(avu-ava)+lamda(2)*arr(4)*avu+lamda(3)*arr(5)*(avu+ava))
     flux(3)=0.5d0*(Fl(3)+Fr(3))-0.5d0*(lamda(1)*arr(1)*avv+lamda(2)*arr(2)*avv+lamda(2)*arr(4)*avv+lamda(3)*arr(5)*avv)
     flux(4)=0.5d0*(Fl(4)+Fr(4))-0.5d0*(lamda(1)*arr(1)*avw+lamda(2)*arr(3)*avw+lamda(2)*arr(4)*avw+lamda(3)*arr(5)*avw)
     flux(5)=0.5d0*(Fl(5)+Fr(5))-0.5d0*(lamda(1)*arr(1)*(avh-avu*ava)+lamda(2)*arr(2)*avv*avv+lamda(2)*arr(3)*avw*avw  & 
	               +lamda(2)*arr(4)*0.5d0*(avu*avu+avv*avv+avw*avw)+lamda(3)*arr(5)*(avh+avu*ava))
    
	end

!c========================================================
! Van Leer 流通矢量分裂 （代码由冷岩开发）
!-------Extend 1D Van Leer FLux for finite volume method
! F=F+ (uL) + F- (uR)	

    subroutine Flux_Van_Leer_1Da(QL,QR,Flux,gamma)   
     use const_var
     implicit none
     real(PRE_EC):: QL(5),QR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar ,Ml,Mr,Mp,Mm  ! uu velocity
     real(PRE_EC):: tmp0,fp(5),fm(5)

     dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5) 
     dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
     al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
     ar=sqrt(gamma*pr/dr)  
     Ml=uul/(al); Mr=uur/(ar)
     if(Ml>=1.d0) then
       fp(1)=dl*uul
       fp(2)=dl*uul*uul+pl
       fp(3)=dl*uul*vvl
       fp(4)=dl*uul*wwl
       fp(5)=uul*(gamma*pl/(gamma-1.d0)+0.5d0*dl*(uul*uul+vvl*vvl+wwl*wwl))
     else if(abs(Ml)<1.d0) then 
       Mp=0.25d0*(1.d0+Ml)*(1.d0+Ml)
       tmp0=dl*al*Mp
       fp(1)=tmp0
       fp(2)=tmp0*((gamma-1.d0)*uul+2.d0*al)/gamma
       fp(3)=tmp0*vvl
       fp(4)=tmp0*wwl
       fp(5)=tmp0*(((gamma-1.d0)*uul+2.d0*al)*((gamma-1.d0)*uul+2.d0*al)*0.5d0/(gamma*gamma-1.d0)+0.5d0*(vvl*vvl+wwl*wwl))
      else if(Ml<=-1.d0)   then
        fp(1)=0.d0
        fp(2)=0.d0
        fp(3)=0.d0
        fp(4)=0.d0
        fp(5)=0.d0
      end if	
		
      if(Mr>= 1.d0) then
        fm(1)=0.d0
        fm(2)=0.d0
        fm(3)=0.d0
        fm(4)=0.d0
        fm(5)=0.d0
      else if(abs(Mr) < 1.d0) then 
        Mm=-0.25d0*(Mr-1.d0)*(Mr-1.d0)
        tmp0=dr*ar*Mm
        fm(1)=tmp0
        fm(2)=tmp0*((gamma-1.d0)*uur-2.d0*ar)/gamma
        fm(3)=tmp0*vvr
        fm(4)=tmp0*wwr
        fm(5)=tmp0*(((gamma-1.d0)*uur-2.d0*ar)*((gamma-1.d0)*uur-2.d0*ar)*0.5d0/(gamma*gamma-1.d0)+0.5d0*(vvr*vvr+wwr*wwr))   
      else if(Mr<=-1.d0)    then
        fm(1)=dr*uur
        fm(2)=dr*uur*uur+pr
        fm(3)=dr*uur*vvr
        fm(4)=dr*uur*wwr
        fm(5)=uur*(gamma*pr/(gamma-1.d0)+0.5d0*dr*(uur*uur+vvr*vvr+wwr*wwr))
      end if
	  
	  Flux=fp+fm
      
	  end

!------------------------------------------------------------------------------------------------
! AUSM+ 通量分裂 （代码由冷岩开发）
!-------Extend 1D Ausm FLux for finite volume method
! F=F+ (uL) + F- (uR)	

    subroutine Flux_Ausm_1Da(QL,QR,Flux,gamma)   
     use const_var
     implicit none
     real(PRE_EC):: QL(5),QR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, hl,dr,uur,vvr,wwr,pr,ar,hr  ! uu velocity
     real(PRE_EC):: a,Ml,Mr,Mp4,Mm4,Pp5,Pm5,M,mp,mm,p,dm,xm
     real(PRE_EC):: fp(5),fm(5)

     dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5)
     dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
     hl=gamma*pl/((gamma-1.d0)*dl)+0.5*(uul*uul+vvl*vvl+wwl*wwl)
     hr=gamma*pr/((gamma-1.d0)*dr)+0.5*(uur*uur+vvr*vvr+wwr*wwr)
	    
     al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
     ar=sqrt(gamma*pr/dr) 
     a=0.5d0*(al+ar)
     Ml=uul/a; Mr=uur/a
 
     if(abs(Ml)>=1.d0) then
       Mp4=0.5d0*(Ml+abs(Ml))
     else
       Mp4=0.25*(Ml+1.d0)*(Ml+1.d0)+0.125d0*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)
     end if
     if(abs(Mr)>=1.d0) then
       Mm4=0.5d0*(Mr-abs(Mr))
     else
       Mm4=-0.25*(Mr-1.d0)*(Mr-1.d0)-0.125d0*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)
     end if
        
     M=Mp4+Mm4
     mp=dl*a*max(0.d0,M)
     mm=dr*a*min(0.d0,M)
     xm=mp+mm
     if(M>=0.d0) then
       dm=a*abs(M)*dl
     else
       dm=a*abs(M)*dr
     end if

     if(abs(Ml)>=1.d0) then
       Pp5=0.5d0*(Ml+abs(Ml))/Ml
     else
       Pp5=0.25*(Ml+1.d0)*(Ml+1.d0)*(2.d0-Ml)+3.d0*Ml*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)/16.d0
     end if
     
	 if(abs(Mr)>=1.d0) then
       Pm5=0.5d0*(Mr-abs(Mr))/Mr
     else
       Pm5=0.25*(Mr-1.d0)*(Mr-1.d0)*(2.d0+Mr)-3.d0*Mr*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)/16.d0
     end if
     p=pl*Pp5+pr*Pm5
    
     fp(1)=1.d0
     fp(2)=uul
     fp(3)=vvl
     fp(4)=wwl
     fp(5)=hl
        
     fm(1)=1.d0
     fm(2)=uur
     fm(3)=vvr
     fm(4)=wwr
     fm(5)=hr
        
     Flux(1)=0.5d0*xm*(fp(1)+fm(1))-0.5d0*dm*(fm(1)-fp(1))
     Flux(2)=0.5d0*xm*(fp(2)+fm(2))-0.5d0*dm*(fm(2)-fp(2))+p
     Flux(3)=0.5d0*xm*(fp(3)+fm(3))-0.5d0*dm*(fm(3)-fp(3))
     Flux(4)=0.5d0*xm*(fp(4)+fm(4))-0.5d0*dm*(fm(4)-fp(4))
     Flux(5)=0.5d0*xm*(fp(5)+fm(5))-0.5d0*dm*(fm(5)-fp(5))
    
	end
!-----------------------------------------------------------
!AUSMPW-----------------------------------------------------
! Code by Leng Yan -----------------------------------------
    subroutine Flux_Ausmpw_1Da(QL,QR,Flux,gamma)   
     use const_var
     implicit none
     real(PRE_EC):: QL(5),QR(5),Flux(5),gamma
     real(PRE_EC):: dl,uul,vvl,wwl,pl,al, hl,dr,uur,vvr,wwr,pr,ar,hr  ! uu velocity
     real(PRE_EC):: a,Ml,Mr,Mp4,Mm4,Pp5,Pm5,M,wpl,plpl,plpr,MLp,MRm,FL,FR,p
     real(PRE_EC):: fp(5),fm(5)
 
     dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5)
     dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
     hl=gamma*pl/((gamma-1.d0)*dl)+0.5*(uul*uul+vvl*vvl+wwl*wwl)
     hr=gamma*pr/((gamma-1.d0)*dr)+0.5*(uur*uur+vvr*vvr+wwr*wwr)
	    
     al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
     ar=sqrt(gamma*pr/dr) 
     a=0.5d0*(al+ar)
     Ml=uul/a; Mr=uur/a
 
     wpl=1.d0-min(pl/pr,pr/pl)**3.d0
     if((0.75d0 .le. min(pl/pr,pr/pl)) .and.(1.d0 .gt. min(pl/pr,pr/pl))) then
	   plpl=4.d0*min(pl/pr,pr/pl)-3.d0
	 else
	   plpl=0.d0
	 endif
     if((0.75d0 .le. min(pr/pl,pl/pr)) .and. (1.d0 .gt. min(pr/pl,pl/pr))) then
	   plpr=4.d0*min(pr/pl,pl/pr)-3.d0
	 else
	   plpr=0.d0
	 endif
     if(abs(Ml)>=1.d0) then
       Pp5=0.5d0*(Ml+abs(Ml))/Ml
     else
       Pp5=0.25*(Ml+1.d0)*(Ml+1.d0)*(2.d0-Ml)+3.d0*Ml*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)/16.d0
     end if
     if(abs(Mr)>=1.d0) then
       Pm5=0.5d0*(Mr-abs(Mr))/Mr
     else
       Pm5=0.25*(Mr-1.d0)*(Mr-1.d0)*(2.d0+Mr)-3.d0*Mr*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)/16.d0
     end if
     p=pl*Pp5+pr*Pm5

     if(abs(Ml)>=1.d0) then
       Mp4=0.5d0*(Ml+abs(Ml))
	   fL=0.d0
     else
       Mp4=0.25*(Ml+1.d0)*(Ml+1.d0)+0.125d0*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)
	   fL=(pl/p-1.d0)*Plpl*abs(0.25*(Ml+1.d0)*(Ml+1.d0))*min(1.d0,(sqrt(uul*uul+vvl*vvl+wwl*wwl)/a)**0.25d0)   !!! 
     end if
     if(abs(Mr)>=1.d0) then
       Mm4=0.5d0*(Mr-abs(Mr))
	   fR=0.d0
     else
       Mm4=-0.25*(Mr-1.d0)*(Mr-1.d0)-0.125d0*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)
	   fR=(pr/p-1.d0)*Plpr*abs(-0.25*(Mr-1.d0)*(Mr-1.d0))*min(1.d0,(sqrt(uur*uur+vvr*vvr+wwr*wwr)/a)**0.25d0)  ! 
     end if
        
     M=Mp4+Mm4
	 if(M>=0.d0) then
	   MLp=Mp4+Mm4-Mm4*wpl*(1.d0+fR)+(fL*Mp4+fR*Mm4)
     else
       MLp=Mp4*wpl*(1.d0+fL)
	 end if
     if(M<0.d0) then
	   MRm=Mp4+Mm4-Mp4*wpl*(1.d0+fL)+(fL*Mp4+fR*Mm4)
	 else
	   MRm=Mm4*wpl*(1.d0+fR)
	 end if
    
     fp(1)=dl
     fp(2)=dl*uul
     fp(3)=dl*vvl
     fp(4)=dl*wwl
     fp(5)=dl*hl

        
     fm(1)=dr
     fm(2)=dr*uur
     fm(3)=dr*vvr
     fm(4)=dr*wwr
     fm(5)=dr*hr
        
     Flux(1)=MLp*a*fp(1)+MRm*a*fm(1)
     Flux(2)=MLp*a*fp(2)+MRm*a*fm(2)+p
     Flux(3)=MLp*a*fp(3)+MRm*a*fm(3)
     Flux(4)=MLp*a*fp(4)+MRm*a*fm(4)
     Flux(5)=MLp*a*fp(5)+MRm*a*fm(5)
    end

