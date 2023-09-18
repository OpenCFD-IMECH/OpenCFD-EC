!------------------------------------------------------------------------------
! SST model  (See: J. Blazek's book, section 7.2.3) 
! Do not consider transition (full turbulence)
! See:  CFL3D 5.0 manual
! 2018-9-29:  A bug is removed;

  subroutine  Turbulence_model_SST(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,nMesh
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz
   real(PRE_EC):: s1x,s1y,s1z
   real(PRE_EC):: ui,vi,wi,Kti,Wti,uj,vj,wj,Ktj,Wtj,uk,vk,wk,Ktk,Wtk,  &
                  ux,vx,wx,Ktx,Wtx,uy,vy,wy,Kty,Wty,uz,vz,wz,Ktz,Wtz
   real(PRE_EC):: omega,arg1,arg2,arg3,f2,Kws,CD_kw,Pk,Pk1,Pk0,Qk,Qw
   real(PRE_EC):: t11,t22,t33,t12,t13,t23,muk,muw,vn1,vn2,kfi,wfi
   real(PRE_EC):: Cw_SST,beta_SST,sigma_k_SST,sigma_w_SST
   
   real*8,parameter:: sigma_k1_SST=0.85d0,sigma_w1_SST=0.5d0,beta1_SST=0.075d0,Cw1_SST=0.533d0, &
                        sigma_k2_SST=1.d0,  sigma_w2_SST=0.856d0,beta2_SST=0.0828d0,Cw2_SST=0.440d0
   real*8,parameter::a1_SST=0.31d0, betas_SST=0.09d0
   real(PRE_EC),Pointer,dimension(:,:,:):: Kt,Wt,Fluxk,Fluxw  ! 湍能、湍能比耗散率；源项；通量项
   real(PRE_EC),Pointer,dimension(:,:,:):: f1

   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
 
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz
 ! 计算湍流粘性系数
   allocate(Kt(0:nx,0:ny,0:nz),Wt(0:nx,0:ny,0:nz),Fluxk(nx,ny,nz),Fluxw(nx,ny,nz))
   allocate(f1(nx,ny,nz))

! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Kt,Wt,Fluxk,Fluxw,d,uu,v,w,f1)
    
!$OMP DO   
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)*Re             ! 量纲转换
   enddo
   enddo
   enddo
!$OMP END DO   

!$OMP DO   
    do k=0,nz
	do j=0,ny
	do i=0,nx
	 Kt(i,j,k)=B%U(6,i,j,k)/B%U(1,i,j,k)
	 Wt(i,j,k)=B%U(7,i,j,k)/B%U(1,i,j,k)
    enddo
	enddo
	enddo
!$OMP END DO   
   

!$OMP DO   
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
! 计算涡量, 交叉对流项
    ui=uu(i+1,j,k)-uu(i-1,j,k)            
    vi=v(i+1,j,k)-v(i-1,j,k)  
    wi=w(i+1,j,k)-w(i-1,j,k)  
    Kti=Kt(i+1,j,k)-Kt(i-1,j,k)  
    wti=wt(i+1,j,k)-wt(i-1,j,k)  
 
    uj=uu(i,j+1,k)-uu(i,j-1,k)   
    vj=v(i,j+1,k)-v(i,j-1,k)
    wj=w(i,j+1,k)-w(i,j-1,k) 
    Ktj=Kt(i,j+1,k)-Kt(i,j-1,k) 
    wtj=wt(i,j+1,k)-wt(i,j-1,k) 
    
	uk=uu(i,j,k+1)-uu(i,j,k-1)  
    vk=v(i,j,k+1)-v(i,j,k-1)
    wk=w(i,j,k+1)-w(i,j,k-1)  
    Ktk=Kt(i,j,k+1)-Kt(i,j,k-1)  
    wtk=wt(i,j,k+1)-wt(i,j,k-1) 
	 
    ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
    jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
    kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)
   
    ux=ui*ix+uj*jx+uk*kx
    vx=vi*ix+vj*jx+vk*kx
    wx=wi*ix+wj*jx+wk*kx
    Ktx=kti*ix+ktj*jx+ktk*kx
    wtx=wti*ix+wtj*jx+wtk*kx

    uy=ui*iy+uj*jy+uk*ky
    vy=vi*iy+vj*jy+vk*ky
    wy=wi*iy+wj*jy+wk*ky
    Kty=kti*iy+ktj*jy+ktk*ky
    wty=wti*iy+wtj*jy+wtk*ky

    uz=ui*iz+uj*jz+uk*kz
    vz=vi*iz+vj*jz+vk*kz
    wz=wi*iz+wj*jz+wk*kz
    Ktz=kti*iz+ktj*jz+ktk*kz
    wtz=wti*iz+wtj*jz+wtk*kz

! 涡量
     omega=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
	 arg2=max( 2.d0* sqrt(abs(Kt(i,j,k)))/(0.09*Wt(i,j,k)*B%dw(i,j,k)*Re) , &
	          500.d0*B%mu(i,j,k)/(d(i,j,k)*Wt(i,j,k)*B%dw(i,j,k)**2 *Re*Re) )
     f2=tanh(arg2**2)
     B%mu_t(i,j,k)=a1_SST*d(i,j,k)*Kt(i,j,k)/max(a1_SST*Wt(i,j,k),f2*abs(omega)/Re)


 ! 计算f1 (识别是否为近壁区，近壁区趋近于1）      
     
 !    Kws=2.d0*(ktx*wtx+kty*wty+ktz*ktz)*d(i,j,k)*sigma_w2_SST/(Wt(i,j,k)+1.d-20)      ! Bug
      Kws=2.d0*(ktx*wtx+kty*wty+ktz*wtz)*d(i,j,k)*sigma_w2_SST/(Wt(i,j,k)+1.d-20)      ! 交叉输运项
    
     CD_kw=max(Kws,1.d-20)
     arg3=max(sqrt(abs(Kt(i,j,k)))/(0.09*Wt(i,j,k)*B%dw(i,j,k) *Re)  , &
	          500.d0*B%mu(i,j,k)/(d(i,j,k)*Wt(i,j,k)*B%dw(i,j,k)**2 *Re*Re) )
	 arg1=min(arg3,4.d0*d(i,j,k)*sigma_w2_SST*Kt(i,j,k)/(CD_kw*B%dw(i,j,k)**2 ))
     f1(i,j,k)=tanh(arg1**4)             ! 开关函数，近壁区趋近于1，远壁区趋近于0  （用来切换k-w及k-epsl方程)
     
     
!    湍应力 （使用了涡粘模型）     ! Blazek's Book, Eq. (7.25)
        
!	     t11=(4.d0/3.d0)*ux-(2.d0/3.d0)*(vy+wz)  
!         t22=(4.d0/3.d0)*vy-(2.d0/3.d0)*(ux+wz) 
!         t33=(4.d0/3.d0)*wz-(2.d0/3.d0)*(ux+vy) 
!         t12=uy+vx
!         t13=uz+wx
!         t23=vz+wy


!    湍能方程的源项（生成-耗散)     
!     Pk1=t11*ux+t22*vy+t33*wz+t12*(uy+vx)+t13*(uz+wx)+t23*(vz+wy)     
!	  Pk=B%mu_t(i,j,k)*Pk1                                                   ! 湍能生成项 （湍应力乘以应变率）       

      Pk=B%mu_t(i,j,k)*omega*omega
!     Pk0=min(Pk,20.d0*betas_SST*Kt(i,j,k)*Wt(i,j,k)*Re*Re)                        ! 对湍能生成项进行限制，防止湍能过大
   
     Pk0=Pk        ! 不进行限制

	 Qk=Pk0/Re-Re*betas_SST*d(i,j,k)*Wt(i,j,k)*Kt(i,j,k)    ! k方程的源项  （生成项-耗散项）

     Cw_SST=f1(i,j,k)*Cw1_SST+(1.d0-f1(i,j,k))*Cw2_SST    ! 模型系数，利用f1函数进行切换
     beta_SST=f1(i,j,k)*beta1_SST+(1.d0-f1(i,j,k))*beta2_SST    ! 模型系数，利用f1函数进行切换
!    Qw= Cw_SST*d(i,j,k)*Pk1          &
!	           -beta_SST*d(i,j,k)*Wt(i,j,k)**2+(1.d0-f1(i,j,k))*Kws     ! W方程的源项    
     Qw= Cw_SST*d(i,j,k)*omega*omega/Re          &
	           -Re*beta_SST*d(i,j,k)*Wt(i,j,k)**2+(1.d0-f1(i,j,k))*Kws/Re     ! W方程的源项    

!-------------------------------------------

	B%Res(6,i,j,k)=QK*B%vol(i,j,k)
	B%Res(7,i,j,k)=Qw*B%vol(i,j,k)
 
	enddo
	enddo
    enddo
!$OMP END DO   
!$OMP END PARALLEL
	

! 设定湍流粘性系数虚网格的值
! mut in Ghost Cell of the boundary
! 采用临近点的值

   B%mu_t(0,:,:)=B%mu_t(1,:,:)
   B%mu_t(nx,:,:)=B%mu_t(nx-1,:,:)
   B%mu_t(:,0,:)=B%mu_t(:,1,:)
   B%mu_t(:,ny,:)=B%mu_t(:,ny-1,:)
   B%mu_t(:,:,0)=B%mu_t(:,:,1)
   B%mu_t(:,:,nz)=B%mu_t(:,:,nz-1)
!  固壁上
   call  Amut_boundary(nMesh,mBlock)
 
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Kt,Wt,Fluxk,Fluxw,d,uu,v,w,f1)
!$OMP DO   
!------i- direcion ---------------------------
     do k=1,nz-1 
     do j=1,ny-1
     do i=1,nx
! 扩散项
! 扩散系数，界面上的值=两侧值的平均, 边界上的扩散系数=内侧的值   
       sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
       sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
       muk=(B%mu(i-1,j,k)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i-1,j,k)+B%mu_t(i,j,k)) )*0.5d0 /Re        ! 扩散系数 (k方程), 界面上的值=两侧值的平均
       muw=(B%mu(i-1,j,k)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i-1,j,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)
       s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
          Kti=Kt(i,j,k)-Kt(i-1,j,k)                               ! K
          Wti=wt(i,j,k)-wt(i-1,j,k)                               ! W
          Ktj=0.25d0*(Kt(i,j+1,k)-Kt(i,j-1,k)+Kt(i-1,j+1,k)-Kt(i-1,j-1,k))
          Wtj=0.25d0*(Wt(i,j+1,k)-Wt(i,j-1,k)+Wt(i-1,j+1,k)-Wt(i-1,j-1,k))
          Ktk=0.25d0*(Kt(i,j,k+1)-Kt(i,j,k-1)+Kt(i-1,j,k+1)-Kt(i-1,j,k-1))
          Wtk=0.25d0*(Wt(i,j,k+1)-Wt(i,j,k-1)+Wt(i-1,j,k+1)-wt(i-1,j,k-1))

          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)
          ktx=kti*ix+ktj*jx+ktk*kx                                ! Kt对坐标的导数
          wtx=wti*ix+wtj*jx+wtk*kx
          kty=kti*iy+ktj*jy+ktk*ky
          wty=wti*iy+wtj*jy+wtk*ky
          ktz=kti*iz+ktj*jz+ktk*kz
          wtz=wti*iz+wtj*jz+wtk*kz
!          对流项
          vn1=uu(i-1,j,k)*s1x+v(i-1,j,k)*s1y+w(i-1,j,k)*s1z   ! 法向速度
          vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
          kfi=0.5d0*((vn1+abs(vn1))*kt(i-1,j,k)+(vn2-abs(vn2))*kt(i,j,k))  ! 对流项，一阶 L-F格式
          wfi=0.5d0*((vn1+abs(vn1))*wt(i-1,j,k)+(vn2-abs(vn2))*wt(i,j,k))  ! 对流项，一阶 L-F格式
		  Fluxk(i,j,k)= (-kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))* B%Si(i,j,k)
		  Fluxw(i,j,k)= (-wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))* B%Si(i,j,k)
      enddo
     enddo
    enddo
!$OMP END DO   

!$OMP DO   
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
         B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i+1,j,k)-Fluxk(i,j,k)        
         B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i+1,j,k)-Fluxw(i,j,k)        
	  enddo
      enddo
      enddo
!$OMP END DO   

!-----j- direction ------------------------

!$OMP  DO   
    do k=1,nz-1 
    do j=1,ny
    do i=1,nx-1
       sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
       sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
       muk=(B%mu(i,j-1,k)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i,j-1,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (k方程), 界面上的值=两侧值的平均
       muw=(B%mu(i,j-1,k)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i,j-1,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)

      s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
      kti=0.25d0*(kt(i+1,j,k)-kt(i-1,j,k)+kt(i+1,j-1,k)-kt(i-1,j-1,k))
      ktj=kt(i,j,k)-kt(i,j-1,k)
      ktk=0.25d0*(kt(i,j,k+1)-kt(i,j,k-1)+kt(i,j-1,k+1)-kt(i,j-1,k-1))
      wti=0.25d0*(wt(i+1,j,k)-wt(i-1,j,k)+wt(i+1,j-1,k)-wt(i-1,j-1,k))
      wtj=wt(i,j,k)-wt(i,j-1,k)
      wtk=0.25d0*(wt(i,j,k+1)-wt(i,j,k-1)+wt(i,j-1,k+1)-wt(i,j-1,k-1))
     
	  ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)
      ktx=kti*ix+ktj*jx+ktk*kx
      kty=kti*iy+ktj*jy+ktk*ky
      ktz=kti*iz+ktj*jz+ktk*kz
      wtx=wti*ix+wtj*jx+wtk*kx
      wty=wti*iy+wtj*jy+wtk*ky
      wtz=wti*iz+wtj*jz+wtk*kz
      vn1=uu(i,j-1,k)*s1x+v(i,j-1,k)*s1y+w(i,j-1,k)*s1z   ! 法向速度
      vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
      kfi=0.5d0*((vn1+abs(vn1))*kt(i,j-1,k)+(vn2-abs(vn2))*kt(i,j,k))  ! 一阶 L-F格式
      wfi=0.5d0*((vn1+abs(vn1))*wt(i,j-1,k)+(vn2-abs(vn2))*wt(i,j,k))  ! 一阶 L-F格式
      
	  Fluxk(i,j,k)= (-kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))* B%Sj(i,j,k)
	  Fluxw(i,j,k)= (-wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))* B%Sj(i,j,k)

    enddo
    enddo
    enddo
!$OMP END DO 
!$OMP  DO 
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
      B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i,j+1,k)-Fluxk(i,j,k)           
      B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i,j+1,k)-Fluxw(i,j,k)           
	 enddo
     enddo
     enddo
!$OMP END DO 
!-----k- direction ------------------------
!$OMP  DO 
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
      sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
      sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
      muk=(B%mu(i,j,k-1)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i,j,k-1)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (k方程), 界面上的值=两侧值的平均
      muw=(B%mu(i,j,k-1)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i,j,k-1)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)

      s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
      kti=0.25d0*(kt(i+1,j,k)-kt(i-1,j,k)+kt(i+1,j,k-1)-kt(i-1,j,k-1))
      ktj=0.25d0*(kt(i,j+1,k)-kt(i,j-1,k)+kt(i,j+1,k-1)-kt(i,j-1,k-1))
      ktk=kt(i,j,k)-kt(i,j,k-1)
      wti=0.25d0*(wt(i+1,j,k)-wt(i-1,j,k)+wt(i+1,j,k-1)-wt(i-1,j,k-1))
      wtj=0.25d0*(wt(i,j+1,k)-wt(i,j-1,k)+wt(i,j+1,k-1)-wt(i,j-1,k-1))
      wtk=wt(i,j,k)-wt(i,j,k-1)
      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)
      ktx=kti*ix+ktj*jx+ktk*kx
      kty=kti*iy+ktj*jy+ktk*ky
      ktz=kti*iz+ktj*jz+ktk*kz
      wtx=wti*ix+wtj*jx+wtk*kx
      wty=wti*iy+wtj*jy+wtk*ky
      wtz=wti*iz+wtj*jz+wtk*kz

     vn1=uu(i,j,k-1)*s1x+v(i,j,k-1)*s1y+w(i,j,k-1)*s1z   ! 法向速度
     vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
     kfi=0.5d0*((vn1+abs(vn1))*kt(i,j,k-1)+(vn2-abs(vn2))*kt(i,j,k))  ! 一阶 L-F格式
     wfi=0.5d0*((vn1+abs(vn1))*wt(i,j,k-1)+(vn2-abs(vn2))*wt(i,j,k))  ! 一阶 L-F格式
    
	 Fluxk(i,j,k)= ( -kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
 	 Fluxw(i,j,k)= ( -wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
   
	enddo
    enddo
    enddo
!$OMP END DO
 
!$OMP  DO 
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i,j,k+1)-Fluxk(i,j,k)           
       B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i,j,k+1)-Fluxw(i,j,k)           
	 enddo
     enddo
     enddo
!$OMP END DO 

 !$OMP DO   
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)/Re
    B%mu_t(i,j,k)=B%mu_t(i,j,k)/Re
   enddo
   enddo
   enddo
!$OMP END DO   


!$OMP END PARALLEL 
!--------------源项---------------------------------------------------------
  deallocate(f1,Kt,Wt,Fluxk,Fluxw)

end  subroutine  Turbulence_model_SST

