!-----------------------------------------------------------
! Copyright by LiXinliang
! Comput Q
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
  end module Const_Variables

 
  module Global_Variables
  use Const_Variables
  implicit none

   TYPE Block_TYPE           !  variables for each block 
     integer :: Block_no,nx,ny,nz
	 real(PRE_EC),pointer,dimension(:,:,:)::  x,y,z
	 real(PRE_EC),pointer,dimension(:,:,:)::  xc,yc,zc
	 real(PRE_EC),pointer,dimension(:,:,:) :: dc,uc,vc,wc,Tc,Qc

   End TYPE Block_TYPE  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
   integer:: Num_Block, Mesh_File_Format
   integer,allocatable,dimension(:):: NI,NJ,NK

  end module Global_Variables
!------------------------------------------------------------------




!-------------------------------------------
     program PLOT_Q  
     use Global_Variables
     implicit none
     integer:: k,nx,ny,nz,m
	 Type (Block_TYPE),pointer:: B
 
     print*, "read Mesh3d.dat... (PLOT3D Format)"
     print*, " please input the format, 0 unformatted, 1 formatted"
      read(*,*) Mesh_File_Format
  
     if(Mesh_File_Format .eq. 0) then
       open(99,file="Mesh3d.dat",form="unformatted")
       read(99) Num_Block         
     else
      open(99,file="Mesh3d.dat")
      read(99,*) Num_Block        
     endif
     allocate(Block(Num_Block))             
     allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   
    if(Mesh_File_Format .eq. 0) then
     read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
    else
     read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
    endif

    open(100,file="flow3d.dat",form="unformatted")
    open(101,file="Q-tec.dat")
    write(101,*) "variables=x,y,z,Q"

!--------------------------------------------
    do m=1,Num_Block
      B => Block(m)
      B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   
	  nx=B%nx ; ny= B%ny ; nz=B%nz
       allocate(B%x(nx,ny,nz),B%y(nx,ny,nz),B%z(nx,ny,nz))
       allocate(B%xc(nx-1,ny-1,nz-1), B%yc(nx-1,ny-1,nz-1), B%zc(nx-1,ny-1,nz-1)) 
       allocate(B%dc(0:nx,0:ny,0:nz),B%uc(0:nx,0:ny,0:nz),B%vc(0:nx,0:ny,0:nz), &
             B%wc(0:nx,0:ny,0:nz),B%Tc(0:nx,0:ny,0:nz),B%Qc(0:nx,0:ny,0:nz))

  
	  call read_mesh_flow(m)
      call computQ (m)
      call PlotQ(m)
      deallocate(B%x,B%y,B%z,B%xc,B%yc,B%zc,B%dc,B%uc,B%vc,B%wc,B%Tc,B%Qc)
	enddo

    deallocate(NI,NJ,NK)
    close(99)
    close(100)
	close(101)

	 end

!=====================================================================================
   subroutine read_mesh_flow(m)
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B

! ---------node Coordinates----------------------------------------  
     B => Block(m)
     nx=B%nx ; ny= B%ny ; nz=B%nz

!------read mesh ---------------------------------------
   if(Mesh_File_Format .eq. 0) then
	read(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else
  	read(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                 (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                 (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif
 
!----------Comput the Cell-Center's coordinate-----------------------------------
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
     B%xc(i,j,k)=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j+1,k)+        &
	              B%x(i,j,k+1)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1) + &
                  B%x(i,j+1,k+1)+B%x(i+1,j+1,k+1))/8.d0
     B%yc(i,j,k)=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j+1,k)+        &
	              B%y(i,j,k+1)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1) + &
                  B%y(i,j+1,k+1)+B%y(i+1,j+1,k+1))/8.d0
     B%zc(i,j,k)=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j+1,k)+        &
	              B%z(i,j,k+1)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1) + &
                  B%z(i,j+1,k+1)+B%z(i+1,j+1,k+1))/8.d0
   
   enddo
   enddo
   enddo


!-----read flow (center value)---------------------
        read(100)   (((B%dc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                    (((B%uc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                    (((B%vc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                    (((B%wc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
		            (((B%Tc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 


 end 

!----------------------------------------------------

   subroutine computQ (m)
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: xi,yi,zi,xj,yj,zj,xk,yk,zk,Jac1,Jac, &
                  ix,iy,iz,jx,jy,jz,kx,ky,kz, &
				  ui,vi,wi,uj,vj,wj,uk,vk,wk, &
				  ux,vx,wx,uy,vy,wy,uz,vz,wz

   
   B => Block(m)
	  nx=B%nx ; ny= B%ny ; nz=B%nz


   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
!----- get S (normal of vorticity)  S=sqrt(0.5*Omiga_ij*Omiga_ij) at the cell's center ------
!  计算涡量；计算湍流粘性系数的梯度
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
     Jac=1.d0/Jac1

     ix=Jac*(yj*zk-zj*yk)
     iy=Jac*(zj*xk-xj*zk)
     iz=Jac*(xj*yk-yj*xk)
     jx=Jac*(yk*zi-zk*yi)
     jy=Jac*(zk*xi-xk*zi)
     jz=Jac*(xk*yi-yk*xi)
     kx=Jac*(yi*zj-zi*yj)
     ky=Jac*(zi*xj-xi*zj)
     kz=Jac*(xi*yj-yi*xj)



   ui=B%uc(i+1,j,k)-B%uc(i-1,j,k)            
   vi=B%vc(i+1,j,k)-B%vc(i-1,j,k)  
   wi=B%wc(i+1,j,k)-B%wc(i-1,j,k)  

   uj=B%uc(i,j+1,k)-B%uc(i,j-1,k)   
   vj=B%vc(i,j+1,k)-B%vc(i,j-1,k)
   wj=B%wc(i,j+1,k)-B%wc(i,j-1,k) 

   uk=B%uc(i,j,k+1)-B%uc(i,j,k-1)  
   vk=B%vc(i,j,k+1)-B%vc(i,j,k-1)
   wk=B%wc(i,j,k+1)-B%wc(i,j,k-1)  
   

   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
 
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz

! Q
   B%Qc(i,j,k)=ux*vy+ux*wz+vy*wz -uy*vx-uz*wx-vz*wy  !! Q=II(UX)

!--------------------------------------------------------------------------   
 enddo
 enddo
 enddo

end


   subroutine PlotQ (m)
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B

   
   B => Block(m)
    nx=B%nx ; ny= B%ny ; nz=B%nz
    write(101,*) "zone  i= ", nx-1, " j= ", ny-1, " k= ", nz-1
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
     write(101,'(4E18.9)') B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%Qc(i,j,k)
	enddo
	enddo
	enddo
 end
