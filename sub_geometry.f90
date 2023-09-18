!-----------------------------------------------------
!   计算几何量：控制体体积和中心点坐标，Jocabian系数  
!   2013-4-26:  可处理退化线 （面积为0的面） 
!   2013-5-3: 修改粘性项Jocabian系数计算方法，与物理量导数方法一致 

  subroutine Comput_Goemetric_var(nMesh)
   use   Global_Var
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nMesh
   real(PRE_EC):: t1x,t1y,t1z,t2x,t2y,t2z,s1x,s1y,s1z,xa,ya,za,ss
   real(PRE_EC):: xi,yi,zi,xj,yj,zj,xk,yk,zk,Jac,Jac1
   real(PRE_EC):: xi1,xi2,yi1,yi2,zi1,zi2,xj1,xj2,yj1,yj2,zj1,zj2,xk1,xk2,yk1,yk2,zk1,zk2
   real(PRE_EC),allocatable,dimension(:,:,:)::Vi,Vj,Vk
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
!  计算控制体的体积
!  计算控制体各表面的面积 （为了避免内存占用过多，表面的法方向、切方向在计算中求出，不进行存储）
   MP=>Mesh(nMesh)
   do m=1,MP%Num_Block   
     B => MP%Block(m)
     allocate(Vi(B%nx,B%ny,B%nz),Vj(B%nx,B%ny,B%nz),Vk(B%nx,B%ny,B%nz))
!  Area of surface i, j and k     
     do k=1,B%nz
       do j=1,B%ny
         do i=1,B%nx
!-------------------------------------------
           t1x=B%x(i,j+1,k)-B%x(i,j,k+1); t1y=B%y(i,j+1,k)-B%y(i,j,k+1); t1z=B%z(i,j+1,k)-B%z(i,j,k+1)   ! 对角线1
           t2x=B%x(i,j+1,k+1)-B%x(i,j,k); t2y=B%y(i,j+1,k+1)-B%y(i,j,k) ; t2z=B%z(i,j+1,k+1)-B%z(i,j,k)  ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)  ! 长度
           B%Si(i,j,k)=ss*0.5d0
           
		   if(ss .ge. Lim_Zero) then
		     B%ni1(i,j,k)=s1x/ss; B%ni2(i,j,k)=s1y/ss ; B%ni3(i,j,k)=s1z/ss
           else
		     B%ni1(i,j,k)=1.d0; B%ni2(i,j,k)=0.d0 ; B%ni3(i,j,k)=0.d0          ! 退化线；任意确定法方向
		   endif
		     

           xa=(B%x(i,j,k)+B%x(i,j+1,k)+B%x(i,j,k+1)+B%x(i,j+1,k+1))*0.25d0
           ya=(B%y(i,j,k)+B%y(i,j+1,k)+B%y(i,j,k+1)+B%y(i,j+1,k+1))*0.25d0
           za=(B%z(i,j,k)+B%z(i,j+1,k)+B%z(i,j,k+1)+B%z(i,j+1,k+1))*0.25d0
           Vi(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!----------------------------------------------------------------------------------
           t1x=B%x(i+1,j,k+1)-B%x(i,j,k); t1y=B%y(i+1,j,k+1)-B%y(i,j,k) ; t1z=B%z(i+1,j,k+1)-B%z(i,j,k)  ! 对角线1
           t2x=B%x(i+1,j,k)-B%x(i,j,k+1); t2y=B%y(i+1,j,k)-B%y(i,j,k+1);  t2z=B%z(i+1,j,k)-B%z(i,j,k+1)   ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)
           B%Sj(i,j,k)=ss*0.5d0
    
	   	   if(ss .ge. Lim_Zero) then
    	     B%nj1(i,j,k)=s1x/ss; B%nj2(i,j,k)=s1y/ss ; B%nj3(i,j,k)=s1z/ss
           else
    	     B%nj1(i,j,k)=1.d0; B%nj2(i,j,k)=0.d0 ; B%nj3(i,j,k)=0.d0
		   endif


           xa=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j,k+1)+B%x(i+1,j,k+1))*0.25d0
           ya=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j,k+1)+B%y(i+1,j,k+1))*0.25d0
           za=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j,k+1)+B%z(i+1,j,k+1))*0.25d0
           Vj(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!----------------------=----------------------------------------------------
           t1x=B%x(i+1,j+1,k)-B%x(i,j,k); t1y=B%y(i+1,j+1,k)-B%y(i,j,k) ; t1z=B%z(i+1,j+1,k)-B%z(i,j,k)  ! 对角线1
           t2x=B%x(i,j+1,k)-B%x(i+1,j,k); t2y=B%y(i,j+1,k)-B%y(i+1,j,k) ; t2z=B%z(i,j+1,k)-B%z(i+1,j,k)   ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)
           B%Sk(i,j,k)=ss*0.5d0  
           if(ss .ge. Lim_Zero) then
             B%nk1(i,j,k)=s1x/ss; B%nk2(i,j,k)=s1y/ss; B%nk3(i,j,k)=s1z/ss
           else
             B%nk1(i,j,k)=1.d0; B%nk2(i,j,k)=0.d0; B%nk3(i,j,k)=0.d0
		   endif

           xa=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j+1,k)+B%x(i+1,j+1,k))*0.25d0
           ya=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j+1,k)+B%y(i+1,j+1,k))*0.25d0
           za=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j+1,k)+B%z(i+1,j+1,k))*0.25d0
           Vk(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!-----------------------------------------------
         enddo
       enddo
     enddo

! 控制体 体积    
     do k=1,B%nz-1
       do j=1,B%ny-1
         do i=1,B%nx-1
           B%vol(i,j,k)=(Vi(i+1,j,k)-Vi(i,j,k)+Vj(i,j+1,k)-Vj(i,j,k)+Vk(i,j,k+1)-Vk(i,j,k))/3.d0
         enddo
       enddo
     enddo

!  网格中心点坐标
     do k=0,B%nz
       do j=0,B%ny
         do i=0,B%nx
           B%xc(i,j,k)=(B%x(i,j,k)+B%x(i,j+1,k)+B%x(i,j,k+1)+B%x(i,j+1,k+1)+ &
                       B%x(i+1,j,k)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1)+B%x(i+1,j+1,k+1))*0.125
           B%yc(i,j,k)=(B%y(i,j,k)+B%y(i,j+1,k)+B%y(i,j,k+1)+B%y(i,j+1,k+1)+ &
                       B%y(i+1,j,k)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1)+B%y(i+1,j+1,k+1))*0.125
           B%zc(i,j,k)=(B%z(i,j,k)+B%z(i,j+1,k)+B%z(i,j,k+1)+B%z(i,j+1,k+1)+ &
                       B%z(i+1,j,k)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1)+B%z(i+1,j+1,k+1))*0.125
         enddo
       enddo
     enddo              
     deallocate(Vi,Vj,Vk)
 ! -----------计算 Jocabian系数 （粘性项计算导数时使用）------------------------
 !  (I+1/2,J,K)点
 ! revised, 2013-5-3:  Jocabian系数与 物理量导数 计算方法一致，避免额外误差；
 ! revised, 2013-5-4:  避免使用角点坐标（计算域立方体的棱），以免出现不稳定性
       do k=1,B%nz-1 
       do j=1,B%ny-1
       do i=1,B%nx
        
		  xi=B%xc(i,j,k)-B%xc(i-1,j,k) 
          yi=B%yc(i,j,k)-B%yc(i-1,j,k)
          zi=B%zc(i,j,k)-B%zc(i-1,j,k)
		 

          if( (i==1 .or. i==B%nx) .and. (j==1 .or. j==B%ny-1) ) then
		   xj1=B%xc(i,j-1,k)
		   yj1=B%yc(i,j-1,k)
		   zj1=B%zc(i,j-1,k)
		   xj2=B%xc(i,j+1,k)
		   yj2=B%yc(i,j+1,k)
		   zj2=B%zc(i,j+1,k)
          else
		   xj1=0.5d0*(B%xc(i,j-1,k)+B%xc(i-1,j-1,k))
		   yj1=0.5d0*(B%yc(i,j-1,k)+B%yc(i-1,j-1,k))
		   zj1=0.5d0*(B%zc(i,j-1,k)+B%zc(i-1,j-1,k))
		   xj2=0.5d0*(B%xc(i,j+1,k)+B%xc(i-1,j+1,k))
		   yj2=0.5d0*(B%yc(i,j+1,k)+B%yc(i-1,j+1,k))
		   zj2=0.5d0*(B%zc(i,j+1,k)+B%zc(i-1,j+1,k))
          endif

          if( (i==1 .or. i==B%nx) .and. (k==1 .or. k==B%nz-1) ) then
            xk1=B%xc(i,j,k-1)
            yk1=B%yc(i,j,k-1)
            zk1=B%zc(i,j,k-1)
            xk2=B%xc(i,j,k+1)
            yk2=B%yc(i,j,k+1)
            zk2=B%zc(i,j,k+1)
		  else
            xk1=0.5d0*(B%xc(i,j,k-1)+B%xc(i-1,j,k-1))
            yk1=0.5d0*(B%yc(i,j,k-1)+B%yc(i-1,j,k-1))
            zk1=0.5d0*(B%zc(i,j,k-1)+B%zc(i-1,j,k-1))
            xk2=0.5d0*(B%xc(i,j,k+1)+B%xc(i-1,j,k+1))
            yk2=0.5d0*(B%yc(i,j,k+1)+B%yc(i-1,j,k+1))
            zk2=0.5d0*(B%zc(i,j,k+1)+B%zc(i-1,j,k+1))
          endif
	      
 		   xj=0.5d0*(xj2-xj1)
		   yj=0.5d0*(yj2-yj1)
		   zj=0.5d0*(zj2-zj1)
 		   xk=0.5d0*(xk2-xk1)
		   yk=0.5d0*(yk2-yk1)
		   zk=0.5d0*(zk2-zk1)

           Jac1=(xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk)
           if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
		   Jac=1.d0/Jac1

!   Jac=B%Jaci(i,j,k)
!   9个Jocabian变换系数    
          B%ix1(i,j,k)=Jac*(yj*zk-zj*yk)
          B%iy1(i,j,k)=Jac*(zj*xk-xj*zk)
          B%iz1(i,j,k)=Jac*(xj*yk-yj*xk)
          B%jx1(i,j,k)=Jac*(yk*zi-zk*yi)
          B%jy1(i,j,k)=Jac*(zk*xi-xk*zi)
          B%jz1(i,j,k)=Jac*(xk*yi-yk*xi)
          B%kx1(i,j,k)=Jac*(yi*zj-zi*yj)
          B%ky1(i,j,k)=Jac*(zi*xj-xi*zj)
          B%kz1(i,j,k)=Jac*(xi*yj-yi*xj)
      enddo
	  enddo
	  enddo
 
 ! (I,J-1/2,K) 点的值， 即 (i+1/2,j,k+1/2)点的值
 
! Revised, 2013-5-3, 坐标的导数与物理量的导数 计算方法相同
      do k=1,B%nz-1 
      do j=1,B%ny
      do i=1,B%nx-1
       xj=B%xc(i,j,k)-B%xc(i,j-1,k)
       yj=B%yc(i,j,k)-B%yc(i,j-1,k)
       zj=B%zc(i,j,k)-B%zc(i,j-1,k)
      
! Revised, 2013-5-4, 避免使用角点（棱）坐标
	   if( (j==1 .or. j==B%ny) .and. (i==1 .or. i==B%nx-1) ) then
		xi1=B%xc(i-1,j,k)
		yi1=B%yc(i-1,j,k)
		zi1=B%zc(i-1,j,k)
		xi2=B%xc(i+1,j,k)
		yi2=B%yc(i+1,j,k)
		zi2=B%zc(i+1,j,k)
	   else
		xi1=0.5d0*(B%xc(i-1,j,k)+B%xc(i-1,j-1,k))
		yi1=0.5d0*(B%yc(i-1,j,k)+B%yc(i-1,j-1,k))
		zi1=0.5d0*(B%zc(i-1,j,k)+B%zc(i-1,j-1,k))
		xi2=0.5d0*(B%xc(i+1,j,k)+B%xc(i+1,j-1,k))
		yi2=0.5d0*(B%yc(i+1,j,k)+B%yc(i+1,j-1,k))
		zi2=0.5d0*(B%zc(i+1,j,k)+B%zc(i+1,j-1,k))
       endif
	
	   if( (j==1 .or. j==B%ny) .and. (k==1 .or. k==B%nz-1) ) then
		 xk1=B%xc(i,j,k-1)
		 yk1=B%yc(i,j,k-1)
		 zk1=B%zc(i,j,k-1)
		 xk2=B%xc(i,j,k+1)
		 yk2=B%yc(i,j,k+1)
		 zk2=B%zc(i,j,k+1)
	   else	 
		 xk1=0.5d0*(B%xc(i,j,k-1)+B%xc(i,j-1,k-1))
		 yk1=0.5d0*(B%yc(i,j,k-1)+B%yc(i,j-1,k-1))
		 zk1=0.5d0*(B%zc(i,j,k-1)+B%zc(i,j-1,k-1))
		 xk2=0.5d0*(B%xc(i,j,k+1)+B%xc(i,j-1,k+1))
		 yk2=0.5d0*(B%yc(i,j,k+1)+B%yc(i,j-1,k+1))
		 zk2=0.5d0*(B%zc(i,j,k+1)+B%zc(i,j-1,k+1))
       endif
	    xi=0.5d0*(xi2-xi1)
	    yi=0.5d0*(yi2-yi1)
	    zi=0.5d0*(zi2-zi1)
	    xk=0.5d0*(xk2-xk1)
	    yk=0.5d0*(yk2-yk1)
	    zk=0.5d0*(zk2-zk1)
	
      
	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1
	 
	  B%ix2(i,j,k)=Jac*(yj*zk-zj*yk)
      B%iy2(i,j,k)=Jac*(zj*xk-xj*zk)
      B%iz2(i,j,k)=Jac*(xj*yk-yj*xk)
      B%jx2(i,j,k)=Jac*(yk*zi-zk*yi)
      B%jy2(i,j,k)=Jac*(zk*xi-xk*zi)
      B%jz2(i,j,k)=Jac*(xk*yi-yk*xi)
      B%kx2(i,j,k)=Jac*(yi*zj-zi*yj)
      B%ky2(i,j,k)=Jac*(zi*xj-xi*zj)
      B%kz2(i,j,k)=Jac*(xi*yj-yi*xj)
	 enddo
	 enddo
	 enddo
	
! (I,J,K-1/2) 点的值， 即 (i+1/2,j+1/2,k)点的值
! Revised, 2013-5-3 
 
     do k=1,B%nz 
     do j=1,B%ny-1
     do i=1,B%nx-1
      xk=B%xc(i,j,k)-B%xc(i,j,k-1)
      yk=B%yc(i,j,k)-B%yc(i,j,k-1)
      zk=B%zc(i,j,k)-B%zc(i,j,k-1)
 
 
 	 if( (k==1 .or. k==B%nz) .and. (i==1 .or. i==B%nx-1) ) then
	  xi1=B%xc(i-1,j,k)
	  yi1=B%yc(i-1,j,k)
	  zi1=B%zc(i-1,j,k)
	  xi2=B%xc(i+1,j,k)
	  yi2=B%yc(i+1,j,k)
	  zi2=B%zc(i+1,j,k)
     else
	  xi1=0.5d0*(B%xc(i-1,j,k)+B%xc(i-1,j,k-1))
	  yi1=0.5d0*(B%yc(i-1,j,k)+B%yc(i-1,j,k-1))
	  zi1=0.5d0*(B%zc(i-1,j,k)+B%zc(i-1,j,k-1))
	  xi2=0.5d0*(B%xc(i+1,j,k)+B%xc(i+1,j,k-1))
	  yi2=0.5d0*(B%yc(i+1,j,k)+B%yc(i+1,j,k-1))
	  zi2=0.5d0*(B%zc(i+1,j,k)+B%zc(i+1,j,k-1))
     endif

  	 if( (k==1 .or. k==B%nz) .and. (j==1 .or. j==B%ny-1) ) then
      xj1=B%xc(i,j-1,k)
      yj1=B%yc(i,j-1,k)
      zj1=B%zc(i,j-1,k)
      xj2=B%xc(i,j+1,k)
      yj2=B%yc(i,j+1,k)
      zj2=B%zc(i,j+1,k)
	 else
      xj1=0.5d0*(B%xc(i,j-1,k)+B%xc(i,j-1,k-1))
      yj1=0.5d0*(B%yc(i,j-1,k)+B%yc(i,j-1,k-1))
      zj1=0.5d0*(B%zc(i,j-1,k)+B%zc(i,j-1,k-1))
      xj2=0.5d0*(B%xc(i,j+1,k)+B%xc(i,j+1,k-1))
      yj2=0.5d0*(B%yc(i,j+1,k)+B%yc(i,j+1,k-1))
      zj2=0.5d0*(B%zc(i,j+1,k)+B%zc(i,j+1,k-1))
	 endif
	  xi=0.5d0*(xi2-xi1)
	  yi=0.5d0*(yi2-yi1)
	  zi=0.5d0*(zi2-zi1)
	  xj=0.5d0*(xj2-xj1)
	  yj=0.5d0*(yj2-yj1)
	  zj=0.5d0*(zj2-zj1)
    
	 

	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1
	 
	 B%ix3(i,j,k)=Jac*(yj*zk-zj*yk)
     B%iy3(i,j,k)=Jac*(zj*xk-xj*zk)
     B%iz3(i,j,k)=Jac*(xj*yk-yj*xk)
     B%jx3(i,j,k)=Jac*(yk*zi-zk*yi)
     B%jy3(i,j,k)=Jac*(zk*xi-xk*zi)
     B%jz3(i,j,k)=Jac*(xk*yi-yk*xi)
     B%kx3(i,j,k)=Jac*(yi*zj-zi*yj)
     B%ky3(i,j,k)=Jac*(zi*xj-xi*zj)
     B%kz3(i,j,k)=Jac*(xi*yj-yi*xj)
     enddo
	 enddo
	 enddo
!  (I,J,K)点的值, 标量方程的源项需要 (仅最密的网格使用)
    if(nMesh .eq. 1) then
    do k=1,B%nz-1
    do j=1,B%ny-1
    do i=1,B%nx-1
	 xi=B%xc(i+1,j,k)-B%xc(i-1,j,k) 
     yi=B%yc(i+1,j,k)-B%yc(i-1,j,k)
     zi=B%zc(i+1,j,k)-B%zc(i-1,j,k)
     xj=B%xc(i,j+1,k)-B%xc(i,j-1,k)   
     yj=B%yc(i,j+1,k)-B%yc(i,j-1,k)
     zj=B%zc(i,j+1,k)-B%zc(i,j-1,k)
     xk=B%xc(i,j,k+1)-B%xc(i,j,k-1) 
     yk=B%yc(i,j,k+1)-B%yc(i,j,k-1)
     zk=B%zc(i,j,k+1)-B%zc(i,j,k-1)
	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1

     B%ix0(i,j,k)=Jac*(yj*zk-zj*yk)
     B%iy0(i,j,k)=Jac*(zj*xk-xj*zk)
     B%iz0(i,j,k)=Jac*(xj*yk-yj*xk)
     B%jx0(i,j,k)=Jac*(yk*zi-zk*yi)
     B%jy0(i,j,k)=Jac*(zk*xi-xk*zi)
     B%jz0(i,j,k)=Jac*(xk*yi-yk*xi)
     B%kx0(i,j,k)=Jac*(yi*zj-zi*yj)
     B%ky0(i,j,k)=Jac*(zi*xj-xi*zj)
     B%kz0(i,j,k)=Jac*(xi*yj-yi*xj)
    enddo
    enddo
	enddo
	endif
	
	enddo

 
  end subroutine Comput_Goemetric_var

!-------------------------------------------------
  subroutine check_mesh_quality
   use   Global_Var
   implicit none
   integer::nMesh
   do nMesh=1,Num_Mesh
   call check_mesh_quality_onemesh(nMesh)
   enddo
   end


! 检查网格质量
! 检查方法： 网格的连续性 （体积的连续性、法方向的连续性）
  subroutine check_mesh_quality_onemesh(nMesh)
   use   Global_Var
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nMesh,i1,j1,k1,i2,j2,k2,im,jm,km,in,jn,kn
   real(PRE_EC):: fl,ft,flmax,ftmax,fl1,fl2,fl3,ft1,ft2,ft3,Af,Aa
   real(PRE_EC):: x1,x2,y1,y2,z1,z2,Vmax,Vmin
   character(len=50):: filename

   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   
   if(my_id .eq. 0)    print*, "Check mesh quality ......"

   MP=>Mesh(nMesh)
   do m=1,MP%Num_Block   
     B => MP%Block(m)
     nx=B%nx; ny=B%ny ; nz=B%nz
 
 ! ----统计最大、最小网格 ----------------------
      Vmax=B%vol(1,1,1)
	  Vmin=B%vol(1,1,1)
	  im=1; jm=1; km=1
	  in=1; jn=1; kn=1
      
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   if(B%vol(i,j,k) > Vmax) then
	    Vmax=B%vol(i,j,k)
		im=i ; jm=j ; km=k
	   endif
	   if(B%vol(i,j,k) < Vmin) then
	    Vmin=B%vol(i,j,k)
		in=i; jn=j; kn=k
	   endif
	   enddo
	   enddo
	   enddo



 !-----计算网格长度比、网格线的转角 ------------
 
 	 flmax=1.d0
	 ftmax=0.d0

!  长度比


     do k=1,nz-1
	 do j=1,ny-1
 	 do i=1,nx-1
      fl1=sqrt((B%xc(i,j,k)-B%xc(i-1,j,k))**2+(B%yc(i,j,k)-B%yc(i-1,j,k))**2+(B%zc(i,j,k)-B%zc(i-1,j,k))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i+1,j,k))**2+(B%yc(i,j,k)-B%yc(i+1,j,k))**2+(B%zc(i,j,k)-B%zc(i+1,j,k))**2)
      fl2=sqrt((B%xc(i,j,k)-B%xc(i,j-1,k))**2+(B%yc(i,j,k)-B%yc(i,j-1,k))**2+(B%zc(i,j,k)-B%zc(i,j-1,k))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i,j+1,k))**2+(B%yc(i,j,k)-B%yc(i,j+1,k))**2+(B%zc(i,j,k)-B%zc(i,j+1,k))**2)
      fl3=sqrt((B%xc(i,j,k)-B%xc(i,j,k-1))**2+(B%yc(i,j,k)-B%yc(i,j,k-1))**2+(B%zc(i,j,k)-B%zc(i,j,k-1))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i,j,k+1))**2+(B%yc(i,j,k)-B%yc(i,j,k+1))**2+(B%zc(i,j,k)-B%zc(i,j,k+1))**2)
	        
         if(fl1 .lt. 1.d0) fl1=1.d0/fl1
         if(fl2 .lt. 1.d0) fl2=1.d0/fl2
         if(fl3 .lt. 1.d0) fl3=1.d0/fl3
         fl=max(fl1,fl2,fl3)
     
	    x1=B%xc(i,j,k)-B%xc(i-1,j,k) ; x2= B%xc(i+1,j,k)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i-1,j,k) ; y2= B%yc(i+1,j,k)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i-1,j,k) ; z2= B%zc(i+1,j,k)-B%zc(i,j,k)
        ft1=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
    
	    x1=B%xc(i,j,k)-B%xc(i,j-1,k) ; x2= B%xc(i,j+1,k)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i,j-1,k) ; y2= B%yc(i,j+1,k)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i,j-1,k) ; z2= B%zc(i,j+1,k)-B%zc(i,j,k)
        ft2=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))

	    x1=B%xc(i,j,k)-B%xc(i,j,k-1) ; x2= B%xc(i,j,k+1)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i,j,k-1) ; y2= B%yc(i,j,k+1)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i,j,k-1) ; z2= B%zc(i,j,k+1)-B%zc(i,j,k)
        ft3=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
        ft=min(1.d0,1.d0*min(ft1,ft2,ft3))
        ft=acos(ft)           ! 网格线折角 (容易出现NaN)
        
        Af=0.9d0*exp(-4.d0*(fl-1.d0)**2)+0.1d0
        Aa=0.9d0*exp(-(4.d0/3.1415926535d0*ft)**2)+0.1d0
        B%dtime_mesh(i,j,k)=min(Af,Aa)
       
	    if(ISNAN(Af) .or. ISNAN(Aa)) then
		 print*, "--------Find bad grid ------------"
		 print*, "Af, Aa=", Af,Aa
		 print*, "Block_no, i,j,k=",B%block_no, i,j,k
		 write(*,"(7E30.20)") fl,ft,ft1,ft2,ft3,min(ft1,ft2,ft3),acos(min(ft1,ft2,ft3))
		endif

!---------找出质量最差的网格，输出-------------------	   
       if(fl .gt. flmax) then
	    flmax=fl
	    i1=i
	    j1=j
	    k1=k
	   endif
   
       if(ft .gt. ftmax) then
	   ftmax=ft
	   i2=i
	   j2=j
	   k2=k
	   endif


	  enddo
      enddo
	  enddo

       if(nMesh .eq. 1) then
        open(103,file="mesh-quality.dat",access="append")
        write(103,*) "------------------Block ", m, "---------------------------"
        write(103,*) "Cell Number=", (nx-1)*(ny-1)*(nz-1)
        write(103,*) "Max volume=", Vmax, im,jm,km
		write(103,*) "Min Volume=", Vmin, in,jn,kn
	    write(103,*) "============="
	    write(103,*) "max grid factor=",flmax, i1,j1,k1
		write(103,*) "dtime_mesh=", B%dtime_mesh(i1,j1,k1)
	    write(103,*) "max grid-line turn angle (degree)=",ftmax*180.d0/3.1415926535d0, i2,j2,k2
	    write(103,*) "dtime_mesh=", B%dtime_mesh(i2,j2,k2)
        close(103)
	   endif
     
    enddo

    if(nMesh .eq. 1) then
       if(IF_Debug .eq. 1) then
	    write(filename, "('mesh-quality-'I5.5'.dat')") my_id
		open(106, file=filename)
        write(106,*) "variables=x,y,z,dtfact"
	    do m=1,MP%Num_Block   
         B => MP%Block(m)
         nx=B%nx; ny=B%ny ; nz=B%nz
         write(106,*) "zone i=", nx-1, " j= ", ny-1, " k= ",nz-1
		do k=1,nz-1
		do j=1,ny-1
        do i=1,nx-1
		write(106,"(4f20.10)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%dtime_mesh(i,j,k)
	    enddo
		enddo
		enddo
	    enddo
	    close(106)
       endif

!	  close(103)
     if(my_id .eq. 0)    print*, "Check mesh quality OK"

    endif


	end  subroutine check_mesh_quality_onemesh


!------------------------------------------------------------------    
!  设定各重网格上的控制信息
  subroutine set_control_para
   use Global_var
   implicit none
   integer nMesh
   TYPE (Mesh_TYPE),pointer:: MP
   MP=>Mesh(1)            ! 最细的网格
!  最细网格上的控制参数与主控制参数相同
   MP%Iflag_turbulence_model=Iflag_turbulence_model
   MP%Iflag_Scheme=Iflag_Scheme
   MP%IFlag_flux=IFlag_flux
   MP%IFlag_Reconstruction=IFlag_Reconstruction
   MP%Bound_Scheme=Bound_scheme   !  边界格式

!  设定粗网格上的控制参数
   do nMesh=2,Num_Mesh
     MP=>Mesh(nMesh)
     MP%Iflag_turbulence_model=Turbulence_NONE    ! 粗网格不使用湍流模型
     MP%Iflag_Scheme=Scheme_UD1                   ! 粗网格使用1阶迎风格式
     MP%IFlag_flux=IFlag_flux                     ! 粗网格的通量分裂技术、时间推进近似及重构技术与细网格相同
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Bound_Scheme=Scheme_UD1                   ! 粗网格边界点使用1阶格式
   enddo
  
  end subroutine set_control_para




  subroutine check_mesh_multigrid 
   use Global_var
   implicit none
   integer,allocatable,dimension(:):: NI,NJ,NK
   integer:: NB,NST,m,k,NN,Km,Km_grid,N_Cell,Ntmp,Bsub,Ksub
   integer:: ib,ie,jb,je,kb,ke,bc,ist,iend,jst,jend,kst,kend
   print*, "Check if Multi-Grid can be used ..."

   if( Mesh_File_Format .eq. 1) then   ! 格式文件
     open(99,file="Mesh3d.dat")
     read(99,*) NB
     allocate(NI(NB),NJ(NB),NK(NB))
     read(99,*) (NI(m), NJ(m), NK(m), m=1,NB)
     close(99)
   else
     open(99,file="Mesh3d.dat",form="unformatted")
     read(99) NB
     allocate(NI(NB),NJ(NB),NK(NB))
     read(99) (NI(m), NJ(m), NK(m), m=1,NB)
     close(99)
  endif


   N_Cell=0
   Km_grid=NI(1)  ! 初始值    
   do m=1,NB 
	 N_Cell=N_Cell+(NI(m)-1)*(NJ(m)-1)*(NK(m)-1)  ! 统计总网格单元数 
!  判断可使用的网格重数      
 	 Km=1
	 NN=2
!  判断准则： 网格数-1 能被2**km 整除， 且最稀的网格单元数不小于2
     do while( mod((NI(m)-1),NN) .eq. 0 .and. (NI(m)-1)/NN .ge. 2     &
		     .and. mod((NJ(m)-1),NN) .eq. 0 .and. (NJ(m)-1)/NN .ge. 2    &
		     .and. mod((NK(m)-1),NN) .eq. 0 .and. (NK(m)-1)/NN .ge. 2) 
       Km=Km+1              ! 所允许的网格重数
	   NN=NN*2
     enddo
     Km_grid=min(Km_grid,Km)
   enddo
!   Print*, " Finished check Mesh3d.dat,  Most stage is ", Km_grid
!   print*,  "Check bc3d.inp ..." 
   open(88,file="bc3d.inp")
   read(88,*)
   read(88,*) 
   do m=1,NB
     read(88,*)
     read(88,*)
     read(88,*) Bsub    !number of the subface in the Block m
     do ksub=1, Bsub
       read(88,*)  ib,ie,jb,je,kb,ke,bc
	   if(bc .lt. 0) read(88,*)
	   ist=min(abs(ib),abs(ie)) ; iend=max(abs(ib),abs(ie))
       jst=min(abs(jb),abs(je));  jend=max(abs(jb),abs(je))
       kst=min(abs(kb),abs(ke));  kend=max(abs(kb),abs(ke))

	   NN=1
	   Km=1
       do while( mod((ist-1),NN) .eq. 0 .and. mod((iend-1),NN) .eq.0       &
		       .and. mod((jst-1),NN) .eq. 0 .and. mod((jend-1),NN) .eq. 0      &
		       .and. mod((kst-1),NN) .eq. 0 .and. mod((kend-1),NN) .eq. 0    ) 
         NN=NN*2
		 Km=Km+1
	   enddo
       Km_grid=min(Km_grid,Km)
     enddo
   enddo
   close(88)

!--------------------------------------------------------- 
   print*, "Total Block number is ", NB, "Total Cell number is " , N_Cell
!   print*, "Most stage number of multi-grid is ", Km_grid
!   print*, "-------------------------------------------------"
!-------------------------------------------------------------
   if(Num_Mesh .gt. Km_grid .or. Num_mesh .gt. 3) then
     print*, "Wrong !, Stage number of multi-grid error !!!"
	 stop 
   endif
   print*, "Check multigrid OK"

!-------------------------------------------------------
   deallocate(NI,NJ,NK)

  end subroutine check_mesh_multigrid

