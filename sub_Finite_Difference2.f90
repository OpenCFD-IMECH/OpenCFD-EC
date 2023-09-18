


! 用差分法计算通量 （相当于补丁程序）

   subroutine Residual_FDM(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   Use FDM_data
   implicit none
   integer:: nMesh,mBlock,nx,ny,nz,NVAR1,i,j,k,m

   Type (Block_TYPE),pointer:: B
   TYPE (Mesh_TYPE),pointer:: MP
   Type (FDM_Block_Type),pointer:: Bm   
     MP=> Mesh(nMesh)
     B => MP%Block(mBlock)
	 Bm=>FDM_Mesh(nMesh)%Block(mBlock)
     nx=B%nx-1 ; ny=B%ny-1 ; nz=B%nz-1   ! 网格中心点的数目   
     NVAR1=MP%NVAR
!	print*, "------------------------------------------------"
!   主干程序与差分法子程序（补丁程序）之间的接口
   call Residual_FDM_local(NVAR1,nx,ny,nz,B%Res(1,-1,-1,-1),    &
       d(1-LAP,1-LAP,1-LAP), uu(1-LAP,1-LAP,1-LAP), v(1-LAP,1-LAP,1-LAP),   & 
	   w(1-LAP,1-LAP,1-LAP),p(1-LAP,1-LAP,1-LAP), T(1-LAP,1-LAP,1-LAP),  &
	   B%mu(-1,-1,-1), B%mu_t(-1,-1,-1),        &
       Bm%ix(1,1,1),Bm%iy(1,1,1),Bm%iz(1,1,1),Bm%jx(1,1,1),Bm%jy(1,1,1),Bm%jz(1,1,1), &
	   Bm%kx(1,1,1),Bm%ky(1,1,1),Bm%kz(1,1,1),Bm%Jac(1,1,1),Cp,PrL,Prt,  &
	   FD_Flux,FD_scheme)
   
     do k=4,nz-3
     do j=4,ny-3
     do i=4,nx-3
     do m=1,5
       B%Res(m,i,j,k)=B%Res(m,i,j,k)*B%Vol(i,j,k)*Bm%Jac(i,j,k)                       
     enddo
     enddo
     enddo
     enddo

   end subroutine Residual_FDM


!----- i- direction --------------------------------------------------------------------------


! 利用差分方法计算通量
! 接口简单，仅需传入几何量（坐标与Jocabian变换系数）及物理量(包括层流及湍流粘性系数)，返回流通量



   subroutine Residual_FDM_local(NVAR,nx,ny,nz,Res,d,u,v,w,p,T,mu,mut,ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac,Cp,Pr,Prt,FD_Flux,FD_scheme)
   use   const_var
   implicit none
   integer:: NVAR,nx,ny,nz   ! 储存变量点的数目（如变量储存在网格中心，则为网格中心点的数目）
!                           注，这里的nx,ny,nz与主干程序中的B%nx, B%ny, B%nz 不同。 主干程序为网格节点的数目，本子程序nx,ny,nz为网格中心点的数目
   
 !  real(PRE_EC),dimension(-1:nx+2,-1:ny+2,-1:nz+2):: d,u,v,w,p,T  ! 与主干程序 (OpenCFD-EC) 结构一致， 如与其他主干程序衔接，可修改
   real(PRE_EC),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP):: d,u,v,w,p,T  ! 与主干程序 (OpenCFD-EC) 结构一致， 如与其他主干程序衔接，可修改
   real(PRE_EC),dimension(-1:nx+2,-1:ny+2,-1:nz+2):: mu,mut         ! 层流及湍流粘性系数, 
   real(PRE_EC),dimension(NVAR,-1:nx+2,-1:ny+2,-1:nz+2 ):: Res     ! 残差（右端项）, 数据结构与主干程序(opencfd-ec)一致

   real(PRE_EC):: fluxi(nx,5),fluxj(ny,5),fluxk(nz,5)
   real(PRE_EC),dimension(nx,ny,nz):: ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac
   real(PRE_EC),allocatable,dimension(:,:,:,:)::  EV1,EV2,EV3
   real(PRE_EC),dimension(nx,5):: fx1,fx2,hi1,hi2
   real(PRE_EC),dimension(ny,5):: fy1,fy2,hj1,hj2
   real(PRE_EC),dimension(nz,5):: fz1,fz2,hk1,hk2
   integer:: i,j,k,m,FD_Flux,FD_scheme
   real(PRE_EC):: Cp,Pr,Prt,Jac1,Amu,Amk

   real(PRE_EC):: ui,vi,wi,Ti,uj,vj,wj,Tj,uk,vk,wk,Tk,ux,vx,wx,Tx,uy,vy,wy,Ty,uz,vz,wz,Tz,t11,t12,t13,t22,t23,t33,E1,E2,E3
   real(PRE_EC):: A1,A2,A3

! 无粘通量的计算
!--------i- -----------------------------------
      
       do k=4,nz-3
       do j=4,ny-3
         do i=1,nx
           Jac1=1.d0/Jac(i,j,k)
		   if(FD_Flux .eq. Flux_Steger_Warming ) then
             A1=ix(i,j,k)*Jac1 ; A2=iy(i,j,k)*Jac1 ; A3=iz(i,j,k)*Jac1
             call split_Steger_Warming(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fx1(i,:),fx2(i,:),A1,A2,A3)
           else if(FD_Flux .eq. Flux_Van_Leer ) then
             A1=ix(i,j,k) ; A2=iy(i,j,k) ; A3=iz(i,j,k)
		     call split_Van_Leer(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fx1(i,:),fx2(i,:),A1,A2,A3,Jac1)
		   else
		     print*, "Flux Splitting Method is not supported"
			 print*, "In this version, only Steger-Warming or Van Leer Splitting method is supported !!!"
			 stop
		   endif
  	     
		 enddo
        do m=1,5     
	     call fp(nx,fx1(1,m),hi1(1,m),FD_scheme)
	     call fm(nx,fx2(1,m),hi2(1,m),FD_scheme)
  	      do i=3,nx-3
	      Fluxi(i,m)=-(hi1(i,m)+hi2(i,m))
	      enddo
          do i=4,nx-3
		  Res(m,i,j,k)=Fluxi(i,m)-Fluxi(i-1,m)              ! fx=f(i+1/2)-f(i-1/2)
		  enddo
	    enddo
      


	    enddo
        enddo
      

!-------j- -------------------------------------
       do k=4,nz-3
       do i=4,nx-3
	     do j=1,ny
	     Jac1=1.d0/Jac(i,j,k)
		 if(FD_Flux .eq. Flux_Steger_Warming ) then
            A1=jx(i,j,k)*Jac1 ; A2=jy(i,j,k)*Jac1 ; A3=jz(i,j,k)*Jac1
            call split_Steger_Warming(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fy1(j,:),fy2(j,:),A1,A2,A3)
		 else if(FD_Flux .eq. Flux_Van_Leer ) then
           A1=jx(i,j,k) ; A2=jy(i,j,k) ; A3=jz(i,j,k)
		   call split_Van_Leer(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fy1(j,:),fy2(j,:),A1,A2,A3,Jac1)
		 endif
  	     enddo

        do m=1,5     
	     call fp(ny,fy1(1,m),hj1(1,m),FD_scheme)
	     call fm(ny,fy2(1,m),hj2(1,m),FD_scheme)
  	      do j=3,ny-3
	      Fluxj(j,m)=-(hj1(j,m)+hj2(j,m))
	      enddo
          do j=4,ny-3
          Res(m,i,j,k)=Res(m,i,j,k)+Fluxj(j,m)-Fluxj(j-1,m)
          enddo

        enddo

	    enddo
        enddo

 
 !-------k- -------------------------------------
       do j=4,ny-3
       do i=4,nx-3
	   do k=1,nz
	   Jac1=1.d0/Jac(i,j,k)
	     if(FD_Flux .eq. Flux_Steger_Warming ) then
         A1=kx(i,j,k)*Jac1 ; A2=ky(i,j,k)*Jac1 ; A3=kz(i,j,k)*Jac1
         call split_Steger_Warming(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fz1(k,:),fz2(k,:),A1,A2,A3)
		 else if(FD_Flux .eq. Flux_Van_Leer ) then
           A1=kx(i,j,k) ; A2=ky(i,j,k) ; A3=kz(i,j,k)
		   call split_Van_Leer(d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),fz1(k,:),fz2(k,:),A1,A2,A3,Jac1)
		 endif
  	   enddo

        do m=1,5     
	     call fp(nz,fz1(1,m),hk1(1,m),FD_scheme)
	     call fm(nz,fz2(1,m),hk2(1,m),FD_scheme)
  	      do k=3,nz-3
	      Fluxk(k,m)=-(hk1(k,m)+hk2(k,m))
	      enddo
		  do k=4,nz-3
		  Res(m,i,j,k)=Res(m,i,j,k)+Fluxk(k,m)-Fluxk(k-1,m)
		  enddo
        enddo
	    enddo
        enddo
 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!   粘性通量的计算
    allocate(EV1(nx,ny,nz,4),Ev2(nx,ny,nz,4),EV3(nx,ny,nz,4))

!  计算应力张量与热流项  (2阶中心差分)
   do k=2,nz-1
   do j=2,ny-1
   do i=2,nx-1

     ui=0.5d0*(u(i+1,j,k)-u(i-1,j,k))
     vi=0.5d0*(v(i+1,j,k)-v(i-1,j,k))
     wi=0.5d0*(w(i+1,j,k)-w(i-1,j,k))
 	 Ti=0.5d0*(T(i+1,j,k)-T(i-1,j,k))
     
     uj=0.5d0*(u(i,j+1,k)-u(i,j-1,k))
     vj=0.5d0*(v(i,j+1,k)-v(i,j-1,k))
     wj=0.5d0*(w(i,j+1,k)-w(i,j-1,k))
     Tj=0.5d0*(T(i,j+1,k)-T(i,j-1,k))
     
     uk=0.5d0*(u(i,j,k+1)-u(i,j,k-1))
     vk=0.5d0*(v(i,j,k+1)-v(i,j,k-1))
     wk=0.5d0*(w(i,j,k+1)-w(i,j,k-1))
     Tk=0.5d0*(T(i,j,k+1)-T(i,j,k-1))

     ux=ui*ix(i,j,k)+uj*jx(i,j,k)+uk*kx(i,j,k)
     vx=vi*ix(i,j,k)+vj*jx(i,j,k)+vk*kx(i,j,k)
     wx=wi*ix(i,j,k)+wj*jx(i,j,k)+wk*kx(i,j,k)
     Tx=Ti*ix(i,j,k)+Tj*jx(i,j,k)+Tk*kx(i,j,k)

     uy=ui*iy(i,j,k)+uj*jy(i,j,k)+uk*ky(i,j,k)
     vy=vi*iy(i,j,k)+vj*jy(i,j,k)+vk*ky(i,j,k)
     wy=wi*iy(i,j,k)+wj*jy(i,j,k)+wk*ky(i,j,k)
     Ty=Ti*iy(i,j,k)+Tj*jy(i,j,k)+Tk*ky(i,j,k)

     uz=ui*iz(i,j,k)+uj*jz(i,j,k)+uk*kz(i,j,k)
     vz=vi*iz(i,j,k)+vj*jz(i,j,k)+vk*kz(i,j,k)
     wz=wi*iz(i,j,k)+wj*jz(i,j,k)+wk*kz(i,j,k)
     Tz=Ti*iz(i,j,k)+Tj*jz(i,j,k)+Tk*kz(i,j,k)

     Amu=mu(i,j,k)+mut(i,j,k)                      ! 层流+湍流粘性系数
	 Amk=Cp*(mu(i,j,k)/Pr + mut(i,j,k)/Prt)        ! 层流+湍流热传导系数

     t11=(4.d0/3.d0*ux-2.d0/3.d0*(vy+wz))*Amu
     t22=(4.d0/3.d0*vy-2.d0/3.d0*(ux+wz))*Amu
     t33=(4.d0/3.d0*wz-2.d0/3.d0*(ux+vy))*Amu
     t12=(uy+vx)*Amu
     t13=(uz+wx)*Amu
     t23=(vz+wy)*Amu

     E1=u(i,j,k)*t11+v(i,j,k)*t12+w(i,j,k)*t13+Amk*Tx
     E2=u(i,j,k)*t12+v(i,j,k)*t22+w(i,j,k)*t23+Amk*Ty
     E3=u(i,j,k)*t13+v(i,j,k)*t23+w(i,j,k)*t33+Amk*Tz

     Jac1=1.d0/Jac(i,j,k)
 	 Ev1(i,j,k,1)=(ix(i,j,k)*t11+iy(i,j,k)*t12+iz(i,j,k)*t13)*Jac1
 	 Ev1(i,j,k,2)=(ix(i,j,k)*t12+iy(i,j,k)*t22+iz(i,j,k)*t23)*Jac1
 	 Ev1(i,j,k,3)=(ix(i,j,k)*t13+iy(i,j,k)*t23+iz(i,j,k)*t33)*Jac1
 	 Ev1(i,j,k,4)=(ix(i,j,k)*E1 +iy(i,j,k)*E2 +iz(i,j,k)*E3)*Jac1
  
 	 Ev2(i,j,k,1)=(jx(i,j,k)*t11+jy(i,j,k)*t12+jz(i,j,k)*t13)*Jac1
 	 Ev2(i,j,k,2)=(jx(i,j,k)*t12+jy(i,j,k)*t22+jz(i,j,k)*t23)*Jac1
 	 Ev2(i,j,k,3)=(jx(i,j,k)*t13+jy(i,j,k)*t23+jz(i,j,k)*t33)*Jac1
 	 Ev2(i,j,k,4)=(jx(i,j,k)*E1+ jy(i,j,k)*E2+ jz(i,j,k)*E3)*Jac1
 	 
     Ev3(i,j,k,1)=(kx(i,j,k)*t11+ky(i,j,k)*t12+kz(i,j,k)*t13)*Jac1
 	 Ev3(i,j,k,2)=(kx(i,j,k)*t12+ky(i,j,k)*t22+kz(i,j,k)*t23)*Jac1
 	 Ev3(i,j,k,3)=(kx(i,j,k)*t13+ky(i,j,k)*t23+kz(i,j,k)*t33)*Jac1
 	 Ev3(i,j,k,4)=(kx(i,j,k)*E1 +ky(i,j,k)*E2+ kz(i,j,k)*E3)*Jac1

   enddo
   enddo
   enddo

!  计算粘性通量
   do k=4,nz-3
   do j=4,ny-3
   do i=4,nx-3
   do m=2,5
    Res(m,i,j,k)=Res(m,i,j,k) +  0.5*(Ev1(i+1,j,k,m-1)-Ev1(i-1,j,k,m-1)   &
	    +Ev2(i,j+1,k,m-1)-Ev2(i,j-1,k,m-1)+Ev3(i,j,k+1,m-1)-Ev3(i,j,k-1,m-1)  )                     ! 2阶中心差分   
   enddo
   enddo
   enddo
   enddo

    deallocate(EV1,EV2,EV3)


  end  

!---------------------------------------------------------------------------------------


    subroutine split_Steger_Warming(d,u,v,w,p,fP,fm, A1,A2,A3)
    use precision_EC
    implicit none
    real(PRE_EC),parameter:: gamma=1.4d0
    integer nx,ny,nz,i,j,k
    real(PRE_EC):: d,u,v,w,p,cc,A1,A2,A3,fp(5),fm(5)
    real(PRE_EC):: tmp0,tmp1,tmp2,tmp3,ss,ak1,ak2,ak3,    &
             E1,E2,E3,E1P,E2P,E3P,E1M,E2M,E3M,   &
          vs,uc1,uc2,vc1,vc2,wc1,wc2,vvc1,vvc2,vv,W2,P2 
    real(PRE_EC),parameter:: epsl=1.d-10      ! 
       
!c El 为特征值,其中El(:,1)为x方向的特征值（5个， u, u, u, u+c, u-c)
        
      tmp1=2.d0*(gamma-1.d0)
      tmp2=1.d0/(2.d0*gamma)
      tmp3=(3.d0-gamma)/(2.d0*(gamma-1.d0)) 

        ss=sqrt(A1*A1+A2*A2 +A3*A3)
        ak1=A1/ss
		ak2=A2/ss
		ak3=A3/ss

        cc=sqrt(gamma*p/d)
	    vs=A1*u+A2*v+A3*w  
!c E1 is lamda1, lamda2 and lamda3; E2 is lamda4; E3 is lamda5 
        E1=vs
        E2=vs-cc*ss
        E3=vs+cc*ss
!----------------------------------------
!      E1P=(E1+sqrt(E1*E1+epsl*epsl))/2.d0
!      E2P=(E2+sqrt(E2*E2+epsl*epsl))/2.d0
!      E3P=(E3+sqrt(E3*E3+epsl*epsl))/2.d0

      E1P=(E1+abs(E1))/2.d0
      E2P=(E2+abs(E2))/2.d0
      E3P=(E3+abs(E3))/2.d0


      E1M=E1-E1P
      E2M=E2-E2P
      E3M=E3-E3P
!----------------------------------------
      tmp0=d/(2.d0*gamma) 
      uc1=u-cc*ak1
      uc2=u+cc*ak1
      vc1=v-cc*ak2
      vc2=v+cc*ak2
      wc1=w-cc*ak3
      wc2=w+cc*ak3
      vvc1=(uc1*uc1+vc1*vc1+wc1*wc1)/2.d0
      vvc2=(uc2*uc2+vc2*vc2+wc2*wc2)/2.d0
      vv=(gamma-1.d0)*(u*u+v*v +w*w)   
      W2=tmp3*cc*cc
! --------The equation seems wrong ! P2 should be zero !
!      P2=tmp1*d(i,j,k)*ak1(i,j)*(ak2(i,j)*w(i,j,k)-ak3(i,j)*v(i,j,k))    !???????
!--------------------------------------------------------
      fp(1)=tmp0*(tmp1*E1P+E2P+E3P)
      fp(2)=tmp0*(tmp1*E1P*u+E2P*uc1+E3P*uc2)
      fp(3)=tmp0*(tmp1*E1P*v+E2P*vc1+E3P*vc2)
      fp(4)=tmp0*(tmp1*E1P*w+E2P*wc1+E3P*wc2)
      fp(5)=tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P))
        
      fm(1)=tmp0*(tmp1*E1M+E2M+E3M)
      fm(2)=tmp0*(tmp1*E1M*u+E2M*uc1+E3M*uc2)
      fm(3)=tmp0*(tmp1*E1M*v+E2M*vc1+E3M*vc2)
      fm(4)=tmp0*(tmp1*E1M*w+E2M*wc1+E3M*wc2)
      fm(5)=tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M))

    return 
   end


!c---------------------------------------------------------------

    subroutine Split_Van_Leer(d,u,v,w,p,fP,fm,A1,A2,A3,Jac1)
    use precision_EC
    implicit none
    real(PRE_EC),parameter:: gamma=1.4d0
    real(PRE_EC):: d,u,v,w,p,cc,Ms,A1,A2,A3,Jac1,fp(5),fm(5)
    real(PRE_EC):: tmp0,ss,ak1,ak2,ak3,vs,vvc1,vvc2,W1,W2

        ss=sqrt(A1*A1+A2*A2+A3*A3)
        ak1=A1/ss
		ak2=A2/ss
		ak3=A3/ss

        cc=sqrt(gamma*p/d)
	    vs=(A1*u+A2*v+A3*w)/ss
		Ms=vs/cc
		
		

		if(Ms .ge. 1.d0) then
		  fp(1)=Jac1*d*vs*ss
          fp(2)=Jac1*(A1*p+d*u*vs*ss)
          fp(3)=Jac1*(A2*p+d*v*vs*ss)
          fp(4)=Jac1*(A3*p+d*w*vs*ss)
          fp(5)=Jac1*vs*ss*(gamma*p/(gamma-1.d0)+0.5d0*d*(u*u+v*v+w*w))
        
          fm(1)=0.d0
          fm(2)=0.d0
          fm(3)=0.d0
          fm(4)=0.d0
          fm(5)=0.d0
		else if(Ms .le. -1.d0) then
		  fp(1)=0.d0
          fp(2)=0.d0
          fp(3)=0.d0
          fp(4)=0.d0
          fp(5)=0.d0
        
          fm(1)=Jac1*d*vs*ss
          fm(2)=Jac1*(A1*p+d*u*vs*ss)
          fm(3)=Jac1*(A2*p+d*v*vs*ss)
          fm(4)=Jac1*(A3*p+d*w*vs*ss)
          fm(5)=Jac1*vs*ss*(gamma*p/(gamma-1.d0)+0.5d0*d*(u*u+v*v+w*w))
		else if(abs(Ms) .lt. 1.d0) then
		  tmp0=ss*Jac1
		  vvc1=d*cc*(Ms+1)*(Ms+1)/4.d0
		  vvc2=-d*cc*(Ms-1)*(Ms-1)/4.d0
!		  W1=vvc1*(0.5d0*((1.d0-gamma)*vs*vs+2.d0*(gamma-1.d0)*vs*cc+2.d0*cc*cc)/(gamma*gamma-1.d0)+0.5d0*(u*u+v*v+w*w))
!		  W2=vvc2*(0.5d0*((1.d0-gamma)*vs*vs-2.d0*(gamma-1.d0)*vs*cc+2.d0*cc*cc)/(gamma*gamma-1.d0)+0.5d0*(u*u+v*v+w*w))
		  W1=vvc1*(((1.d0-gamma)*vs*vs+2.d0*(gamma-1.d0)*vs*cc+2.d0*cc*cc)/(gamma*gamma-1.d0)+0.5d0*(u*u+v*v+w*w))
		  W2=vvc2*(((1.d0-gamma)*vs*vs-2.d0*(gamma-1.d0)*vs*cc+2.d0*cc*cc)/(gamma*gamma-1.d0)+0.5d0*(u*u+v*v+w*w))
!		  W1=vvc1*(cc*cc/(gamma-1.d0)+0.5d0*(u*u+v*v+w*w))  !可使定常流中总焓守恒
!		  W2=vvc2*(cc*cc/(gamma-1.d0)+0.5d0*(u*u+v*v+w*w))
		  fp(1)=tmp0*vvc1
          fp(2)=tmp0*vvc1*(u+ak1*(-vs+2.d0*cc)/gamma)
          fp(3)=tmp0*vvc1*(v+ak2*(-vs+2.d0*cc)/gamma)
          fp(4)=tmp0*vvc1*(w+ak3*(-vs+2.d0*cc)/gamma)
          fp(5)=tmp0*W1
        
          fm(1)=tmp0*vvc2
          fm(2)=tmp0*vvc2*(u+ak1*(-vs-2.d0*cc)/gamma)
          fm(3)=tmp0*vvc2*(v+ak2*(-vs-2.d0*cc)/gamma)
          fm(4)=tmp0*vvc2*(w+ak3*(-vs-2.d0*cc)/gamma)
          fm(5)=tmp0*W2
		endif
    return 
   end


!c---------------------------------------------------------------
!    Positive flux     
	 subroutine fp(nx,v,hj,FD_scheme)
       use const_var
	   implicit none
       integer:: nx,FD_scheme
	   real(PRE_EC):: v(nx),hj(nx)
     if(FD_scheme .eq. FD_WENO5) then
       call fp_weno5(nx,v,hj)           ! 采用5阶WENO
	 else if(FD_scheme .eq. FD_WENO7) then
	   call fp_weno7(nx,v,hj)           ! 内点采用WENO 7
       call fp_weno5_onepoint(nx,v,hj,3)  ! 近左边界点(k=3)仍采用WENO 5
       call fp_weno5_onepoint(nx,v,hj,nx-2)  ! 近右边界点(k=nx-2)人采用WENO 5
     else if (FD_scheme .eq. FD_OMP6) then   ! 内点采用OMP6
	   call fp_OMP6(nx,v,hj)
       call fp_weno5_onepoint(nx,v,hj,3)  ! 近左边界点(k=3)仍采用WENO 5
       call fp_weno5_onepoint(nx,v,hj,nx-2)  ! 近右边界点(k=nx-2)人采用WENO 5
       call fp_weno5_onepoint(nx,v,hj,nx-3)  ! 近右边界点(k=nx-3)人采用WENO 5
     endif
	 end


!    Negative flux     
	 subroutine fm(nx,v,hj,FD_scheme)
       use const_var
	   implicit none
       integer:: nx,FD_scheme
	   real(PRE_EC):: v(nx),hj(nx)
     if(FD_scheme .eq. FD_WENO5) then
       call fm_weno5(nx,v,hj)           ! 采用5阶WENO
	 else if(FD_scheme .eq. FD_WENO7) then
	   call fm_weno7(nx,v,hj)           ! 内点采用WENO 7
       call fm_weno5_onepoint(nx,v,hj,3)  ! 近左边界点(k=3)仍采用WENO 5
       call fm_weno5_onepoint(nx,v,hj,nx-2)  ! 近右边界点(k=nx-2)人采用WENO 5
     else if (FD_scheme .eq. FD_OMP6) then   ! 内点采用OMP6
	   call fm_OMP6(nx,v,hj)
       call fm_weno5_onepoint(nx,v,hj,3)  ! 近左边界点(k=3)仍采用WENO 5
       call fm_weno5_onepoint(nx,v,hj,4)  ! 近左边界点(k=4)仍采用WENO 5
       call fm_weno5_onepoint(nx,v,hj,nx-2)  ! 近右边界点(k=nx-2)人采用WENO 5
     endif
	 end



!----------------------------------------------------------------

 ! 5阶WENO格式计算正通量  
       subroutine fp_weno5(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,k
	   real(PRE_EC):: v(nx),hj(nx)
 
      do k=3,nx-2
      call fp_weno5_onepoint(nx,v,hj,k)
	  enddo

	 return
	 end
!-------------------------------------------------
 ! 5阶WENO格式计算负通量  
       subroutine fm_weno5(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,k
	   real(PRE_EC):: v(nx),hj(nx)
       do k=3,nx-2
       call fm_weno5_onepoint(nx,v,hj,k)
	   enddo
       return
	 end	  
    


!---------------------------------------------------------------
 
 
 ! 5阶WENO格式计算正通量 (单个点) 

       subroutine fp_weno5_onepoint(nx,v,hj,k)
       use precision_EC
	   implicit none
       integer:: nx,k
	   real(PRE_EC):: v(nx),hj(nx)
       real(PRE_EC):: S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
   	   real(PRE_EC),parameter:: ep=1.d-6,C03=3.d0/10.d0,    C13=3.d0/5.d0,    C23=1.d0/10.d0

         S0=13.d0/12.d0*(v(k)-2.d0*v(k+1)+v(k+2))**2+  1.d0/4.d0*(3.d0*v(k)-4.d0*v(k+1)+v(k+2))**2
         S1=13.d0/12.d0*(v(k-1)-2.d0*v(k)+v(k+1))**2+  1.d0/4.d0*(v(k-1)-v(k+1))**2
         S2=13.d0/12.d0*(v(k-2)-2.d0*v(k-1)+v(k))**2+  1.d0/4.d0*(v(k-2)-4.d0*v(k-1)+3.*v(k))**2

      a0=C03/((ep+S0)**2)
	  a1=C13/((ep+S1)**2)
	  a2=C23/((ep+S2)**2)

	  W0=a0/(a0+a1+a2)
      W1=a1/(a0+a1+a2)
	  W2=a2/(a0+a1+a2)

	  q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k+1)-1.d0/6.d0*v(k+2)
	  q13=-1.d0/6.d0*v(k-1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k+1)
	  q23=1.d0/3.d0*v(k-2)-7.d0/6.d0*v(k-1)+11.d0/6.d0*v(k)
	  hj(k)=W0*q03+W1*q13+W2*q23

	 return
	 end
!-------------------------------------------------
 ! 5阶WENO格式计算负通量 (单个点) 
       subroutine fm_weno5_onepoint(nx,v,hj,k)
       use precision_EC
	   implicit none
       integer:: nx,k
	   real(PRE_EC):: v(nx),hj(nx)
       real(PRE_EC):: S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
   	   real(PRE_EC),parameter:: ep=1.d-6,   C03=3.d0/10.d0 ,    C13=3.d0/5.d0,    C23=1.d0/10.d0

       S0=13.d0/12.d0*(v(k)-2.d0*v(k-1)+v(k-2))**2+ 1.d0/4.d0*(3.d0*v(k)-4.d0*v(k-1)+v(k-2))**2
       S1=13.d0/12.d0*(v(k+1)-2.d0*v(k)+v(k-1))**2+ 1.d0/4.d0*(v(k+1)-v(k-1))**2
       S2=13.d0/12.d0*(v(k+2)-2.d0*v(k+1)+v(k))**2+ 1.d0/4.d0*(v(k+2)-4.d0*v(k+1)+3.d0*v(k))**2

      a0=C03/((ep+S0)**2)
	  a1=C13/((ep+S1)**2)
	  a2=C23/((ep+S2)**2)

	  W0=a0/(a0+a1+a2)
      W1=a1/(a0+a1+a2)
	  W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k-1)-1.d0/6.d0*v(k-2)
	 q13=-1.d0/6.d0*v(k+1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k-1)
	 q23=1.d0/3.d0*v(k+2)-7.d0/6.d0*v(k+1)+11.d0/6.d0*v(k)

	 hj(k-1)=W0*q03+W1*q13+W2*q23

     return

	 end	  
    
!c---------------------------------------------------------------
! 7阶WENO格式 WENO-Z (正通量)
       subroutine fp_weno7(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,i
	   real(PRE_EC):: v(nx),hj(nx)

      real(PRE_EC)::  S0,S1,S2,S3, s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
             a0,a1,a2,a3,am,q0,q1,q2,q3,tao7
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
         e41=3.d0/12.d0,  e42=13.d0/12.d0, e43=-5.d0/12.d0,  e44=1.d0/12.d0,   &
		 ep=1.d-20   !! WENO-Z

       do i=4,nx-3
! 7th order WENO scheme
! 1  阶导数  
         S10=a11*v(i-3)+a12*v(i-2)+a13*v(i-1) +a14*v(i)
         S11=a21*v(i-2) -   v(i-1)+a23*v(i)   +a24*v(i+1)
         S12=a31*v(i-1)+a32*v(i)  +    v(i+1) +a34*v(i+2)
         S13=a41*v(i)  +a42*v(i+1)+a43*v(i+2) +a44*v(i+3)
 ! 2 阶导数
         S20=-v(i-3)+b12*v(i-2)+b13*v(i-1)+b14*v(i)             
         S21=             v(i-1)+b22*v(i)  +v(i+1)         
         S22=             v(i)  +b22*v(i+1)+v(i+2)         
         S23=b41*v(i)+b42*v(i+1)+b43*v(i+2)-v(i+3)         
! 3 阶导数
         S30=-v(i-3)+c12*(v(i-2)-v(i-1)) +v(i)                                   
         S31=-v(i-2)+c12*(v(i-1)-v(i))   +v(i+1)                 
         S32=-v(i-1)+c12*(v(i)-v(i+1))   +v(i+2)                 
         S33=-v(i)  +c12*(v(i+1)-v(i+2)) +v(i+3)                 

        S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30
        S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31
        S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32
        S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33

!----------WENO-Z-------------------------
      tao7=abs(S0-S3)
      a0=C0*(1.d0+(tao7/(S0+ep))**2)    ! WENO-Z 
      a1=C1*(1.d0+(tao7/(S1+ep))**2)
      a2=C2*(1.d0+(tao7/(S2+ep))**2)
      a3=C3*(1.d0+(tao7/(S3+ep))**2)

!-----------------------------------------------
     am=a0+a1+a2+a3

!  4阶差分格式的通量
     q0=e11*v(i-3)+e12*v(i-2)+e13*v(i-1) +e14*v(i)
     q1=e21*v(i-2)+e22*v(i-1)+e23*v(i)   +e24*v(i+1)
     q2=e31*v(i-1)+e32*v(i)  +e33*v(i+1) +e34*v(i+2)
     q3=e41*v(i)  +e42*v(i+1)+e43*v(i+2) +e44*v(i+3)

!  由4个4阶差分格式组合成1个7阶差分格式
     hj(i)=(a0*q0+a1*q1+a2*q2+a3*q3)/am
    enddo
   end

!=========================================================================
! 7阶WENO (WENO-Z)     
      subroutine fm_weno7(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,i
	   real(PRE_EC):: v(nx),hj(nx)

      real(PRE_EC)::  S0,S1,S2,S3, s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
             a0,a1,a2,a3,am,q0,q1,q2,q3,tao7
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
         e41=3.d0/12.d0,  e42=13.d0/12.d0, e43=-5.d0/12.d0,  e44=1.d0/12.d0, &
		 ep=1.d-20   !! WENO-Z
   
       do i=4,nx-3

!     7th order WENO scheme
! 1  阶导数
         S10=a11*v(i+3)+a12*v(i+2)+a13*v(i+1)  +a14*v(i)
         S11=a21*v(i+2)-    v(i+1) +a23*v(i)    +a24*v(i-1)
         S12=a31*v(i+1)+a32*v(i)   +    v(i-1)  +a34*v(i-2)
         S13=a41*v(i)  +a42*v(i-1)+a43*v(i-2)  +a44*v(i-3)
! 2 阶导数
         S20=-v(i+3)+b12*v(i+2)+b13*v(i+1)+b14*v(i)              
         S21=             v(i+1) +b22*v(i)  +v(i-1)         
         S22=             v(i)   +b22*v(i-1)+v(i-2)         
         S23=b41*v(i)+b42*v(i-1)+b43*v(i-2)-v(i-3)         
! 3 阶导数 
         S30=-v(i+3)+c12*(v(i+2)-v(i+1))+v(i)                                  
         S31=-v(i+2)+c12*(v(i+1)-v(i))+  v(i-1)                 
         S32=-v(i+1)+c12*(v(i)  -v(i-1))+v(i-2)                 
         S33=-v(i)+  c12*(v(i-1)-v(i-2))+v(i-3)                 

       S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30
       S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31
       S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32
       S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33

       tao7=abs(S0-S3)
       a0=C0*(1.d0+(tao7/(S0+ep))**2)    ! WENO-Z  
       a1=C1*(1.d0+(tao7/(S1+ep))**2)
       a2=C2*(1.d0+(tao7/(S2+ep))**2)
       a3=C3*(1.d0+(tao7/(S3+ep))**2)

!-----------------------------------------------

     am=a0+a1+a2+a3

!  4阶差分格式的通量
     q0=e11*v(i+3)+e12*v(i+2)+e13*v(i+1)+e14*v(i)
     q1=e21*v(i+2)+e22*v(i+1)+e23*v(i)  +e24*v(i-1)
     q2=e31*v(i+1)+e32*v(i)  +e33*v(i-1)+e34*v(i-2)
     q3=e41*v(i)+  e42*v(i-1)+e43*v(i-2)+e44*v(i-3)

!  由4个4阶差分格式组合成1个7阶差分格式
     hj(i-1)=(a0*q0+a1*q1+a2*q2+a3*q3)/am
    enddo

    end
!--------------------------------------------

! Optimized 6th order Monotonicity-Preserving Schemes (Li XL et al.)
! 网格基属于中心格式

!================================================================================
       subroutine fp_OMP6(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,i
	   real(PRE_EC):: v(nx),hj(nx)
  	   real(PRE_EC)::mid_nf
   	   real(PRE_EC)::minmod2,minmod4
       real(PRE_EC)::d1,d2,d3,ul,md,lc,mp,fmax,fmin
       real(PRE_EC),parameter ::kappa=4.d0, ep=1.e-10 
       real(PRE_EC),parameter:: n=-1.3d0/140.d0,  m=0.001d0   ! n: diffusition coefficient; m: dissipation coefficient (larger m means larger dissipation)
 !     real(PRE_EC),parameter:: n=0.d0,  m=0.015d0   ! n: diffusition coefficient; m: dissipation coefficient (larger m means larger dissipation)
       real(PRE_EC),parameter:: a=0.5d0*(m+n),  b=0.5d0*(m-n)	
	   real(PRE_EC),parameter:: a1=a, a2=1.d0/60.d0-b-6.d0*a, a3= -2.d0/15.d0+6.d0*b+15.d0*a, a4=37.d0/60.d0-15.d0*b-20.d0*a, &
					a5=37.d0/60.d0+20.d0*b+15.d0*a,a6=-2.d0/15.d0-6.d0*a-15.d0*b,a7=1.d0/60.d0+a+6.d0*b,a8=-b
           
        do i=4,nx-4 
		    mid_nf=a1*v(i+4)+a2*v(i+3)+a3*v(i+2)+a4*v(i+1)+a5*v(i)+a6*v(i-1)+a7*v(i-2)+a8*v(i-3)
		    mp=v(i)+minmod2((v(i+1)-v(i)),kappa*(v(i)-v(i-1)))
			if((mid_nf-v(i))*(mid_nf-mp) .ge. ep) then
			 d1=v(i-2)+v(i)-2.d0*v(i-1)
			 d2=v(i-1)+v(i+1)-2.d0*v(i)
			 d3=v(i)+v(i+2)-2.d0*v(i+1)
			 ul=v(i)+kappa*(v(i)-v(i-1))
			 md=0.5d0*(v(i)+v(i+1))-0.5d0*minmod4(4.0*d2-d3,4.0*d3-d2,d2,d3)
			 lc=v(i)+0.5d0*(v(i)-v(i-1))+kappa*minmod4(4.0*d1-d2,4.0*d2-d1,d2,d1) /3.d0
			 fmin=max(min(v(i),v(i+1),md),min(v(i),ul,lc))
			 fmax=min(max(v(i),v(i+1),md),max(v(i),ul,lc))
			 mid_nf=mid_nf+minmod2(fmax-mid_nf,fmin-mid_nf)
		    endif
		  hj(i)=mid_nf
        enddo
	  return
	  endsubroutine

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Optimized 6th order Monotonicity-Preserving Schemes (Li XL et al.)
! 网格基属于中心格式

       subroutine fm_OMP6(nx,v,hj)
       use precision_EC
	   implicit none
       integer:: nx,i
	   real(PRE_EC):: v(nx),hj(nx)
  	   real(PRE_EC)::mid_nf
  	   real(PRE_EC)::minmod2,minmod4
       real(PRE_EC)::d1,d2,d3,ul,md,lc,mp,fmax,fmin
       real(PRE_EC),parameter::kappa=4.d0, ep=1.e-10   !(Kappa=1 is more robust)
!       real(PRE_EC),parameter:: n=0.d0,  m=0.015d0   ! n: diffusition coefficient; m: dissipation coefficient (larger m means larger dissipation)
       real(PRE_EC),parameter:: n=-1.3d0/140.d0,  m=0.001d0   ! n: diffusition coefficient; m: dissipation coefficient (larger m means larger dissipation)
       real(PRE_EC),parameter:: a=0.5d0*(m+n),  b=0.5d0*(m-n)	
       real(PRE_EC),parameter:: a1=a, a2=1.d0/60.d0-b-6.d0*a, a3= -2.d0/15.d0+6.d0*b+15.d0*a, a4=37.d0/60.d0-15.d0*b-20.d0*a, &
             a5=37.d0/60.d0+20.d0*b+15.d0*a,a6=-2.d0/15.d0-6.d0*a-15.d0*b,a7=1.d0/60.d0+a+6.d0*b,a8=-b

         do i=4,nx-4
			mid_nf=a1*v(i-3)+a2*v(i-2)+a3*v(i-1)+a4*v(i) +a5*v(i+1)+a6*v(i+2)+a7*v(i+3)+a8*v(i+4)
		    mp=v(i+1)+minmod2((v(i)-v(i+1)),kappa*(v(i+1)-v(i+2)))
			if((mid_nf-v(i+1))*(mid_nf-mp) .ge. ep) then  
			  d1=v(i+3)+v(i+1)-2.d0*v(i+2) 
			  d2=v(i+2)+v(i)-2.d0*v(i+1) 
			  d3=v(i+1)+v(i-1)-2.d0*v(i)
			  !----
		   	  ul=v(i+1)+kappa*(v(i+1)-v(i+2))
			  md=0.5d0*(v(i+1)+v(i))-0.5d0*minmod4(4.0*d2-d3,4.0*d3-d2,d2,d3)
			  lc=v(i+1)+0.5d0*(v(i+1)-v(i))+kappa*minmod4(4.0*d1-d2,4.0*d2-d1,d2,d1) /3.d0
			  fmin=max(min(v(i+1),v(i),md),min(v(i+1),ul,lc))
			  fmax=min(max(v(i+1),v(i),md),max(v(i+1),ul,lc))
			  mid_nf=mid_nf+minmod2(fmax-mid_nf,fmin-mid_nf)
			endif
			hj(i)=mid_nf
           enddo
	  	 return
		endsubroutine


! ==========================================
    function minmod2(x1,x2)
    use precision_EC
    implicit none
	real(PRE_EC)::x1,x2,minmod2
	minmod2=0.5d0*(sign(1.d0,x1)+sign(1.d0,x2))*min(abs(x1),abs(x2))
	return
    endfunction
!=========================================================
    function minmod4(x1,x2,x3,x4)
    use precision_EC
    implicit none
	real(PRE_EC)::x1,x2,x3,x4,minmod4

	minmod4=0.5d0*(sign(1.d0,x1)+sign(1.d0,x2))
	minmod4=minmod4*abs(0.5d0*(sign(1.d0,x1)+sign(1.d0,x3)))
	minmod4=minmod4*abs(0.5d0*(sign(1.d0,x1)+sign(1.d0,x4)))
	minmod4=minmod4*min(abs(x1),abs(x2),abs(x3),abs(x4))

	return
    endfunction

