!------------------------------------------------------------------------------
! Source term of SA model  (See: J. Blazek's book, P240-243) 
! Do not consider transition (full turbulence)
! Ver 0.81a, using Eq. (7.4.2) (J. Blazek's book)
! Ver 0.81b, limit for sorce term  ( not less than 0.3*Omega)
! Ver 0.98c, Nondimensional (See CFL3D manual)

  subroutine  Turbulence_model_NewSA(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,nMesh
   real(PRE_EC):: ui,vi,wi,uj,vj,wj,uk,vk,wk,ux,vx,wx,uy,vy,wy,uz,vz,wz
   real(PRE_EC):: s1x,s1y,s1z,v0,vti,vtj,vtk,vtx,vty,vtz
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz
   real(PRE_EC):: S, S1,X,fv1,fv2,fv3,ft1,ft2,r,g,fw,Q_SA,vn1,vn2,vfi,Q3
   real(PRE_EC):: ppi,pj,pk,px,py,pz,p_plus,Cb1,Cw1
  
   real(PRE_EC),parameter:: SA_sigma=2.d0/3.d0,Cv1=7.1d0,Cv2=5.d0,Cb2=0.622d0,SA_k=0.41d0
   real(PRE_EC),parameter:: Cw2=0.3d0, Cw3=2.d0, Ct1=1.d0, Ct2=2.d0, Ct3=1.2d0, Ct4=0.5d0  !Ct3=1.3d0
!   real(PRE_EC),parameter:: Cb1=0.1355d0, Cw1=Cb1/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma
   real(PRE_EC),parameter:: ep=1.d-6
   real(PRE_EC),Pointer,dimension(:,:,:):: vt,Fluxv,fluxv2
   Type (Block_TYPE),pointer:: B
 
  
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz



 ! 计算湍流粘性系数
   allocate(vt(0:nx,0:ny,0:nz),Fluxv(nx,ny,nz),fluxv2(nx,ny,nz))
	
! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)*Re             ! 量纲转换
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,vt,d)
 
    do k=0,nz
    do j=0,ny
    do i=0,nx
     vt(i,j,k)=B%U(6,i,j,k)
	 X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
     fv1=X**3/(X**3+Cv1**3)
     B%mu_t(i,j,k)=fv1*d(i,j,k)*vt(i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
   
   ! 限定湍流粘性系数
   call limit_mut(nMesh,mBlock)

! 设定湍流粘性系数虚网格的值
 

! 计算vt方程的残差 B%Res(6,:,:,:)

!$OMP PARALLEL  DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,Re,vt,d,uu,v,w,Fluxv,fluxv2,CP1_NSA,CP2_NSA)
!------i- direcion ---------------------------
!$OMP DO
    do k=1,nz-1 
     do j=1,ny-1
      do i=1,nx
          s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
          vti=vt(i,j,k)-vt(i-1,j,k)                               ! SA模型中的vt
          vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i-1,j+1,k)-vt(i-1,j-1,k))
          vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i-1,j,k+1)-vt(i-1,j,k-1))
          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)
          vtx=vti*ix+vtj*jx+vtk*kx
          vty=vti*iy+vtj*jy+vtk*ky
          vtz=vti*iz+vtj*jz+vtk*kz
          v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i-1,j,k)+B%mu(i-1,j,k)/d(i-1,j,k))  ! (I-1/2,J,K) 点的动力学粘性系数
 
          vn1=uu(i-1,j,k)*s1x+v(i-1,j,k)*s1y+w(i-1,j,k)*s1z   ! 法向速度
          vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
          vfi=0.5d0*((vn1+abs(vn1))*vt(i-1,j,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
          Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Si(i,j,k)    !!! Re
          Fluxv2(i,j,k)=(vtx*s1x+vty*s1y+vtz*s1z)*B%Si(i,j,k)
	  enddo
     enddo
    enddo
!$OMP END DO
!$OMP DO
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
	     v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
         B%Res(6,i,j,k)=Fluxv(i+1,j,k)-Fluxv(i,j,k) -v0/Re*(Fluxv2(i+1,j,k)-Fluxv2(i,j,k))    !!! Re     
      enddo
      enddo
      enddo
!$OMP END DO
!-----j- direction ------------------------
!$OMP DO
    do k=1,nz-1 
    do j=1,ny
    do i=1,nx-1
      s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
      vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j-1,k)-vt(i-1,j-1,k))
      vtj=vt(i,j,k)-vt(i,j-1,k)
      vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i,j-1,k+1)-vt(i,j-1,k-1))
      ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j-1,k)+B%mu(i,j-1,k)/d(i,j-1,k))  ! (I-1/2,J,K) 点的动力学粘性系数   s11=Amu1*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
      vn1=uu(i,j-1,k)*s1x+v(i,j-1,k)*s1y+w(i,j-1,k)*s1z   ! 法向速度
      vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
      vfi=0.5d0*((vn1+abs(vn1))*vt(i,j-1,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
      Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Sj(i,j,k)   !!! Re
      Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)* B%Sj(i,j,k)
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
	  v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
      B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j+1,k)-Fluxv(i,j,k)- v0/Re*(Fluxv2(i,j+1,k)-Fluxv2(i,j,k))   !!! Re        
     enddo
     enddo
     enddo
!$OMP END DO
!-----k- direction ------------------------
!$OMP DO
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
    s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
    vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j,k-1)-vt(i-1,j,k-1))
    vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i,j+1,k-1)-vt(i,j-1,k-1))
    vtk=vt(i,j,k)-vt(i,j,k-1)

      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j,k-1)+B%mu(i,j,k-1)/d(i,j,k-1))  ! (I,J,K-1/2) 点的动力学粘性系数
     vn1=uu(i,j,k-1)*s1x+v(i,j,k-1)*s1y+w(i,j,k-1)*s1z   ! 法向速度
     vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
     vfi=0.5d0*((vn1+abs(vn1))*vt(i,j,k-1)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
     Fluxv(i,j,k)= ( -vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
     Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)*B%Sk(i,j,k)   
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
 	   v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
       B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j,k+1)-Fluxv(i,j,k)-v0/Re*(Fluxv2(i,j,k+1)-Fluxv2(i,j,k))           
	 enddo
     enddo
     enddo
!$OMP END DO
!--------------源项---------------------------------------------------------
!$OMP DO
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
!----- get S (normal of vorticity)  S=sqrt(0.5*Omiga_ij*Omiga_ij) at the cell's center ------
!  计算涡量；计算湍流粘性系数的梯度
      

   ui=uu(i+1,j,k)-uu(i-1,j,k)            
   vi=v(i+1,j,k)-v(i-1,j,k)  
   wi=w(i+1,j,k)-w(i-1,j,k)  

   uj=uu(i,j+1,k)-uu(i,j-1,k)   
   vj=v(i,j+1,k)-v(i,j-1,k)
   wj=w(i,j+1,k)-w(i,j-1,k) 

   uk=uu(i,j,k+1)-uu(i,j,k-1)  
   vk=v(i,j,k+1)-v(i,j,k-1)
   wk=w(i,j,k+1)-w(i,j,k-1)  
  

  
   
   ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
   jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
   kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)

   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
 
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz


   ppi=p(i+1,j,k)-p(i-1,j,k)  
   pj=p(i,j+1,k)-p(i,j-1,k) 
   pk=w(i,j,k+1)-p(i,j,k-1)  
   px=ppi*ix+pj*jx+pk*kx
   py=ppi*iy+pj*jy+pk*ky
   pz=ppi*iz+pj*jz+pk*kz


! 涡量
   S=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
   X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
!--------------------------------------------------------------------------   

! 源项 采用 Blazek's book,  p241 (7.38), (7.39)


!   Blazek's book 的公式稳定性不好 (容易算出负的湍流粘性系数)
!   Source term, Blazek's Book section (7.2.1),  modified from original form
!   fv1=X**3/(X**3+Cv1**3)
!   fv2=1.d0/(1.d0+X/Cv2)**3
!   fv3=(1.d0+X*fv1)*(1.d0-fv2)/max(X,0.001d0)
!   S1=fv3*S+vt(i,j,k)*fv2/(SA_k*B%dw(i,j,k))**2           
!   r=vt(i,j,k)/(S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k))

!------------------------------------------------------------------------
! Source term, original form; See: http://turbmodels.larc.nasa.gov/spalart.html
   fv1=X**3/(X**3+Cv1**3)
   fv2=1.d0-X/(1.d0+X*fv1)
   ft2=Ct3*exp(-Ct4*X*X)            
   S1=max(S+fv2*vt(i,j,k)/(Re*(SA_k*B%dw(i,j,k))**2),0.d0)    !!! Re
   r=min(vt(i,j,k)/(Re*S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k)),10.0)    !!! Re
   g=r+Cw2*(r**6-r)
   fw=g*((1.d0+Cw3**6)/(g**6+Cw3**6))**(1.d0/6.d0)
!---------Source term, see: CFL3D manual ---------------------
   p_plus=B%mu(i,j,k)/(d(i,j,k)**3*(B%mu(i,j,k)/d(i,j,k)*S)**1.5d0+ep)*sqrt(px*px+py*py+pz*pz)
 
    Cb1=0.1355d0*(1.d0 +Cp1_NSA*(1.d0-exp(-p_plus/Cp2_NSA)) )
    Cw1=0.1355d0/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma
!   Cw1=Cb1/(SA_k*SA_k)*exp(-P_plus/Cp2)+(1.d0+Cb2)/SA_sigma



    Q3=Cb1*((1-ft2)*fv2+ft2)/(SA_k*SA_k)-Cw1*fw
	Q_SA=Cb1*(1.d0-ft2)*S1*vt(i,j,k)   &
	    +Q3*(vt(i,j,k)/B%dw(i,j,k))**2/Re   !!! Re


!---------------------------------------------------------------------------
   B%Res(6,i,j,k)=B%Res(6,i,j,k)+Q_SA*B%vol(i,j,k)     ! Bug removed

   enddo
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)/Re             ! 量纲转换
    B%mu_t(i,j,k)=B%mu_t(i,j,k)/Re
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



  deallocate(vt,Fluxv,fluxv2)

end  subroutine  Turbulence_model_NewSA

!-----------------------------------------------------------------------------------
!-----根据vt (v~) 计算出SA模型中的湍流粘性系数mut-----------------------------------
! Blazek's Book p241 (7.37)


