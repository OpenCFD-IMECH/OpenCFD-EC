!------------------------------------------------------------------------------
   subroutine turbulence_model_BL(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,ksub,i,j,k,kflag,i1,j1,k1,i2,j2,k2,i0,j0,k0,nMesh
   real(PRE_EC):: ui,vi,wi,uj,vj,wj,uk,vk,wk,ux,vx,wx,uy,vy,wy,uz,vz,wz
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz,x0,y0,z0
   real(PRE_EC),allocatable,dimension(:,:,:):: omiga
   real(PRE_EC),allocatable,dimension(:):: Amu1d,Amut1d,d1d,u1d,yy,omiga1d
   integer,allocatable:: flag1(:,:,:) 
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz
   B%mu_t(:,:,:)=0.d0

! test if the block contains wall    
   kflag=0
   do ksub=1, B%subface
     if(B%bc_msg(ksub)%bc .eq. BC_WALL) kflag=1
   enddo  
   if(Kflag .eq. 0) return    ! No wall in this block

!  This Block Contains Wall
   allocate(omiga(0:nx,0:ny,0:nz),flag1(0:nx,0:ny,0:nz))
   omiga=0.d0
   flag1=0

!----- get Omiga (vorticity)  omiga=sqrt(omigax**2+omigay**2+omigaz**2) at the cell's center ------
! 采用中心差分求解
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (nx,ny,nz,uu,v,w,omiga,B)
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1

! 物理量对于计算坐标（下标）的导数
 
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

!----对物理坐标的偏导数----------------------------------------------
   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky
   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz
!-------------------------------------------------
    omiga(i,j,k)=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
do ksub=1, B%subface
  Bc=> B%bc_msg(ksub)

 if(Bc%bc .eq. BC_WALL) then   ! Wall boundary
  ! 沿网格线一维处理（假设网格线垂直壁面）---------------------------------------
  if(Bc%face .eq. 1 .or. Bc%face .eq. 4) then   ! i+ or i- 
    allocate(yy(nx),Amu1d(nx),Amut1d(nx),d1d(nx),u1d(nx),omiga1d(nx))
    Amut1d=0.d0; Amu1d=0.d0  ! 初始化
   
    
    if(Bc%face .eq. 1) then
      i1=1; i2=0 
    else
      i1=nx-1; i2=nx 
    endif  
   
    do k=Bc%kb,Bc%ke-1
    do j=Bc%jb,Bc%je-1
       x0=(B%xc(i1,j,k)+B%xc(i2,j,k))*0.5d0 
       y0=(B%yc(i1,j,k)+B%yc(i2,j,k))*0.5d0
       z0=(B%zc(i1,j,k)+B%zc(i2,j,k))*0.5d0

     do i=1,nx-1
      if(Bc%face .eq. 1) then
       i0=i
      else
       i0=nx-i
      endif

       yy(i0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(i0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(i0)=d(i,j,k)  
       omiga1d(i0)=omiga(i,j,k)
       Amu1d(i0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(nx-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do i=1,nx-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 1) then
       i0=i
      else
       i0=nx-i
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(i0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(i0))
     endif
     enddo
   
    enddo
    enddo
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

! !!! To set Amu_t=0 in the wall      设置虚网格点上的mut(-1)=-mut(1), 以保证壁面上mut=0 (mut=0.5*(mut(-1)+mut(1))

!  设置壁面第1层网格上的湍流粘性系数为0
    if(Bc%face .eq. 1) then
     B%mu_t(0,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(1,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0   
    else
     B%mu_t(nx,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(nx-1,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
    endif


 else if(Bc%face .eq. 2 .or. Bc%face .eq. 4) then   ! face of j- or j+ 
    allocate(yy(ny),Amu1d(ny),Amut1d(ny),d1d(ny),u1d(ny),omiga1d(ny))
     Amut1d=0.d0; Amu1d=0.d0  ! 初始化

    
    if(Bc%face .eq. 2) then
      j1=1; j2=0 
    else
      j1=ny-1; j2=ny
    endif  
    do k=Bc%kb,Bc%ke-1
    do i=Bc%ib,Bc%ie-1
       x0=(B%xc(i,j1,k)+B%xc(i,j2,k))*0.5d0 
       y0=(B%yc(i,j1,k)+B%yc(i,j2,k))*0.5d0
       z0=(B%zc(i,j1,k)+B%zc(i,j2,k))*0.5d0

     do j=1,ny-1
      if(Bc%face .eq. 2) then
       j0=j
      else
       j0=ny-j
      endif

       yy(j0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(j0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(j0)=d(i,j,k)  
       omiga1d(j0)=omiga(i,j,k)
       Amu1d(j0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(ny-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do j=1,ny-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 2) then
       j0=j
      else
       j0=ny-j
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(j0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(j0))
     endif
     enddo
    enddo
    enddo
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

    if(Bc%face .eq. 2) then
     B%mu_t(Bc%ib:Bc%ie-1, 0,  Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, 1,  Bc%kb:Bc%ke-1)=0.d0   
    else
     B%mu_t(Bc%ib:Bc%ie-1, ny,  Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, ny-1, Bc%kb:Bc%ke-1)=0.d0
    endif




  else if(Bc%face .eq. 3 .or. Bc%face .eq. 6) then   ! face of k- or k+ 
    allocate(yy(nz),Amu1d(nz),Amut1d(nz),d1d(nz),u1d(nz),omiga1d(nz))
    Amut1d=0.d0; Amu1d=0.d0  ! 初始化

    
    if(Bc%face .eq. 3) then
      k1=1; k2=0 
    else
      k1=nz-1; j2=nz
    endif  
    do j=Bc%jb,Bc%je-1
    do i=Bc%ib,Bc%ie-1
       x0=(B%xc(i,j,k1)+B%xc(i,j,k2))*0.5d0 
       y0=(B%yc(i,j,k1)+B%yc(i,j,k2))*0.5d0
       z0=(B%zc(i,j,k1)+B%zc(i,j,k2))*0.5d0

     do k=1,nz-1
      if(Bc%face .eq. 3) then
        k0=k
      else
        k0=nz-k
      endif

       yy(k0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(k0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(k0)=d(i,j,k)  
       omiga1d(k0)=omiga(i,j,k)
       Amu1d(k0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(nz-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do k=1,nz-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 3) then
       k0=k
      else
       k0=nz-k
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(k0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(k0))
     endif
     enddo
   
    enddo
    enddo
  
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

    if(Bc%face .eq. 3) then
     B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, 0)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, 1)=0.d0   
    else
     B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, nz)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, nz-1)=0.d0
    endif
   
  endif
 
 endif
 enddo


   call Amut_boundary(nMesh,mBlock)
   deallocate(omiga,flag1)

end


!c------------------------------------------------------------------------
! B-L model of turbulence
! Ref:  Wilox DC. Turbulence Modeling for CFD (2nd Edition), p77
   subroutine BL_model_1d(ny,yy,Amu,Amu_t,d,u,omiga)
      use precision_EC
	  implicit none
      integer ny,j,Iflag
      real(PRE_EC),dimension(ny) :: yy, Amu,Amu_t,d,u,omiga
      real(PRE_EC),parameter::  AP=26.,Ccp=1.6,Ckleb=0.3,Cwk=0.25d0,AKT=0.4,AK=0.0168  ! Cwk= 1.d0
      real(PRE_EC):: Tw,Ret,Fmax,etamax,Udif,etap,FF,Fwak,bl,Fkleb,Visti,Visto
      
           TW=abs(Amu(1)*omiga(1))
            do j=1,Ny
             if(abs(Amu(j)*omiga(j)) .gt. TW) Tw=abs(Amu(j)*omiga(j))  
            enddo
             Ret=sqrt(d(1)*TW)/Amu(1)

           Fmax=0.d0 ;   etamax=0.d0 ;     Udif=0.d0
  

 !----------------------------------
           do j=1,Ny
            if(u(j) .gt. Udif) Udif=u(j)
             etap=yy(j)*Ret
             FF=yy(j)*abs(omiga(j))*(1.d0-exp(-etap/AP))

            if(FF.gt.Fmax) then
             Fmax=FF
             etamax=yy(j)   
            endif

!            if(FF.gt.Fmax) then
!             Fmax=FF
!             etamax=yy(j)    
!            else
!             goto 100    ! Find the first peak of F(y)   ! 只要第1个峰值         
!            endif

           enddo
100        continue

           IFlag=0
           Fwak=min(etamax*Fmax,Cwk*etamax*Udif*Udif/Fmax)
           do j=1,Ny
            etap=Ret*yy(j)
            bl=AKT*yy(j)*(1.d0-exp(-etap/AP))
            visti=d(j)*bl*bl*abs(omiga(j))
            Fkleb=1.d0/(1.d0+5.5d0*(Ckleb*yy(j)/etamax)**6)
            visto=AK*Ccp*d(j)*Fwak*Fkleb
            if(abs(visto).lt.abs(visti)) IFlag=1
            if(Iflag.eq.0) then
             Amu_t(j)=visti
            else
             Amu_t(j)=visto
            endif
           enddo 

  end



