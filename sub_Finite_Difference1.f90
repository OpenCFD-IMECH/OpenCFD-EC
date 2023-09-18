
!----差分-有限体积混合算法----------------
! Copyright by Li Xinliang (c) lixl@imech.ac.cn

 subroutine init_FDM
  use Global_var
  use FDM_data
  implicit none
   call read_para_FDM
   call comput_Jacobian_FDM
 end

 !---------------------------------------------------
   subroutine read_para_FDM
   use Global_var
   use FDM_data
   implicit none
   logical ext
   integer:: ntmp(10), nbk,k, m,ierr
   integer, allocatable :: BFDM(:)
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
  
    nbk=0
   if(my_id .eq. 0) then
    inquire(file="FDM.in",exist=ext)
	if(ext) then
	 print*, "Find FDM.in, read it ..."
	 open(99,file="FDM.in")
	 read(99,*)
	 read(99,*) 
	 read(99,*)   FD_Flux, FD_scheme   ! 内嵌差分法采用的通量方式、数值格式
     read(99,*)
	 read(99,*)   nbk    ! 内嵌差分法的块数
         allocate(BFDM(nbk))	 
	 read(99,*)
	 read(99,*)   (BFDM(k),k=1,nbk)
     close(99)
   	 else
	  print*, "No FDM.in, Donot use hybrid FVM/FDM"
	 endif 
    endif

      ntmp(1)=FD_Flux
      ntmp(2)=FD_Scheme
	  ntmp(3)=nbk
        call  MPI_bcast(ntmp(1),3,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
      FD_Flux=ntmp(1)
      FD_Scheme=ntmp(2)
	  nbk=ntmp(3)

   
    if(nbk .ne. 0) then
      if(my_id .ne. 0)  allocate(BFDM(nbk))
      call MPI_bcast(BFDM,nbk,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
	endif
  

   MP=>Mesh(1)         ! 仅最密的网格使用FDM
   do m=1,MP%Num_Block   
    B=> MP%Block(m)       ! 本块
    B%IFLAG_FVM_FDM=Method_FVM
Loop1:	do k=1,nbk
	      if( BFDM(k) == B%block_no ) then
	      B%IFLAG_FVM_FDM=Method_FDM
	      exit Loop1
	      endif
	    enddo Loop1
     print*, "my_id, block_no, FDM, FD_Flux, FD_Scheme=", my_id, B%block_no, B%IFLAG_FVM_FDM 
  enddo	   
  end

!-----计算差分法相关的Jacobian系数-------------------

  subroutine comput_Jacobian_FDM
  use Global_var
  use FDM_data
  implicit none
  Type (Block_TYPE),pointer:: B
  Type (FDM_Block_Type),pointer:: Bm

  real(PRE_EC),allocatable,dimension(:,:,:):: x,y,z
  integer:: i,j,k,m,nMesh,nB,nx,ny,nz
  character(len=50):: filename

!---------------------------------------------------------
! 依赖的全局数据：
! Num_Mesh  （整型变量） 网格的套数  （例如，采用3重网格，则该数为3）
! Mesh(k)%Num_Block  网格的块数；
! Flag_FDM(:) （整型数组）, 如果Flag_FDM(k)=0, 则第k块网格采用有限体积法， 如果为1，则采用有限差分法；
!----------------------------------------------------------

! 创建数据结构
   allocate(FDM_Mesh(Num_Mesh))
   do m=1,Num_Mesh
    Num_block=Mesh(m)%Num_Block
    allocate(FDM_Mesh(m)%Block(Num_block))
    do nB=1,Num_block
    Bm=>FDM_Mesh(m)%Block(nB)
	B=>Mesh(m)%Block(nB) 
 	  if(B%IFLAG_FVM_FDM .eq. Method_FDM)   then  ! 该块采用差分法计算
       nx=B%nx-1; ny=B%ny-1; nz=B%nz-1   
	  allocate(Bm%ix(nx,ny,nz), Bm%iy(nx,ny,nz), Bm%iz(nx,ny,nz), &
	          Bm%jx(nx,ny,nz), Bm%jy(nx,ny,nz), Bm%jz(nx,ny,nz), &
			  Bm%kx(nx,ny,nz), Bm%ky(nx,ny,nz), Bm%kz(nx,ny,nz), Bm%jac(nx,ny,nz))
      endif
	enddo
   enddo
!---------------------------------------------------------------
! Comput Jocabian coefficients
	do nMesh=1,Num_Mesh

!    write(filename,"('test-Jacobian-'I1.1'.dat')") nMesh
!    open(99,file=filename)
!    write(99,*) "variables=x,y,z,ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac" 

	do m=1,Mesh(nMesh)%Num_Block
     B  => Mesh(nMesh)%Block(m)
     Bm => FDM_Mesh(nMesh)%Block(m)
	 if(B%IFLAG_FVM_FDM .ne. Method_FDM ) cycle    ! only for FDM (Finite difference method)
 !    print*, " Comput Jocabian Coefficient, Block No.", m

     nx=B%nx-1  ! 网格点数（若物理量储存于网格中心，指的是网格中心的数目）
	 ny=B%ny-1 
	 nz=B%nz-1     


!   申请内存  储存坐标 
     allocate(x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz))

!  物理量所在点的坐标 （若物理量储存在网格中心，则为网格中心点的坐标）
     do k=1,nz
	 do j=1,ny
	 do i=1,nx
      x(i,j,k)=B%xc(i,j,k)
	  y(i,j,k)=B%yc(i,j,k)
	  z(i,j,k)=B%zc(i,j,k)
	 enddo
	 enddo
     enddo
!  计算Jacobian系数
    call grid_Jacobian(nx,ny,nz,   x(1,1,1),y(1,1,1),z(1,1,1), &
	                   Bm%ix(1,1,1),Bm%iy(1,1,1),Bm%iz(1,1,1), &
	                   Bm%jx(1,1,1),Bm%jy(1,1,1),Bm%jz(1,1,1), &
					   Bm%kx(1,1,1),Bm%ky(1,1,1),Bm%kz(1,1,1), &
					   Bm%Jac(1,1,1) )



!	 write(99,*) "zone i= ", nx, " j= ", ny , " k= ", nz
!	 do k=1,nz
!	 do j=1,ny
!	 do i=1,nx
!	    write(99,"(13f20.8)") x(i,j,k),y(i,j,k), z(i,j,k),  &
!		 Bm%ix(i,j,k),Bm%iy(i,j,k),Bm%iz(i,j,k),Bm%jx(i,j,k),Bm%jy(i,j,k),Bm%jz(i,j,k),  &
!		 Bm%kx(i,j,k),Bm%ky(i,j,k),Bm%kz(i,j,k),Bm%Jac(i,j,k)
!	 enddo
!	 enddo
!   enddo

     deallocate(x,y,z)
    enddo
!     close(99)
    enddo

	 if(my_id ==0) print*, "comput Jocabian OK ..."

	end   subroutine comput_Jacobian_FDM



!c================================================================
     subroutine grid_Jacobian(nx,ny,nz,x,y,z,ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac)
     use precision_EC
	 implicit none
     integer:: nx,ny,nz,i,j,k

     real(PRE_EC),dimension(nx,ny,nz):: x,y,z,ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac
	 real(PRE_EC),allocatable,dimension(:,:,:):: xi,xj,xk,yi,yj,yk,zi,zj,zk
 	 real(PRE_EC):: Jac1,xi1,xj1,xk1,yi1,yj1,yk1,zi1,zj1,zk1
	 allocate(xi(nx,ny,nz),xj(nx,ny,nz),xk(nx,ny,nz),  & 
	          yi(nx,ny,nz),yj(nx,ny,nz),yk(nx,ny,nz),  &
			  zi(nx,ny,nz),zj(nx,ny,nz),zk(nx,ny,nz))

	 call dx3d(nx,ny,nz,x,xi)
	 call dx3d(nx,ny,nz,y,yi)
	 call dx3d(nx,ny,nz,z,zi)
	 call dy3d(nx,ny,nz,x,xj)
	 call dy3d(nx,ny,nz,y,yj)
	 call dy3d(nx,ny,nz,z,zj)
	 call dz3d(nx,ny,nz,x,xk)
	 call dz3d(nx,ny,nz,y,yk)
	 call dz3d(nx,ny,nz,z,zk)
    
	 do k=1,nz
	 do j=1,ny
	 do i=1,nx
	  xi1=xi(i,j,k); xj1=xj(i,j,k); xk1=xk(i,j,k)
	  yi1=yi(i,j,k); yj1=yj(i,j,k); yk1=yk(i,j,k)
	  zi1=zi(i,j,k); zj1=zj(i,j,k); zk1=zk(i,j,k)
	  Jac1=1.d0/(xi1*yj1*zk1+yi1*zj1*xk1+zi1*xj1*yk1-zi1*yj1*xk1-yi1*xj1*zk1-xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k) 
      Jac(i,j,k)=Jac1
	  ix(i,j,k)=Jac1*(yj1*zk1-zj1*yk1)
	  iy(i,j,k)=Jac1*(zj1*xk1-xj1*zk1)
	  iz(i,j,k)=Jac1*(xj1*yk1-yj1*xk1)
	  jx(i,j,k)=Jac1*(yk1*zi1-zk1*yi1)
	  jy(i,j,k)=Jac1*(zk1*xi1-xk1*zi1)
	  jz(i,j,k)=Jac1*(xk1*yi1-yk1*xi1)
	  kx(i,j,k)=Jac1*(yi1*zj1-zi1*yj1)
	  ky(i,j,k)=Jac1*(zi1*xj1-xi1*zj1)
	  kz(i,j,k)=Jac1*(xi1*yj1-yi1*xj1)
     if(Jac1 .lt. 0) print*, " Jocabian < 0 !!!, i,j,k, Jac1"
  	enddo
	enddo
	enddo
    deallocate(xi,yi,zi,xj,yj,zj,xk,yk,zk)

    end subroutine grid_Jacobian

!c==========================================================================================
! 采用差分方法计算三个方向的导数 （内部6阶中心差分，边界处降阶）, 用于计算Jacobian系数 
       subroutine dx0(nx,f,fx)
        use precision_EC
		implicit none
		integer:: nx,i
		real(PRE_EC):: b1,b2,a1,a2,a3
        real(PRE_EC):: f(nx),fx(nx)
         b1=8.d0/(12.d0)
         b2=1.d0/(12.d0)
         a1=1.d0/(60.d0)
         a2=-3.d0/(20.d0)
         a3=3.d0/(4.d0)

         do i=4,nx-3
          fx(i)  =a1*(f(i+3)-f(i-3)) +a2*(f(i+2)-f(i-2)) +a3*(f(i+1)-f(i-1))  ! 6th centred
         enddo
 
          fx(1)=(-3.d0*f(1)+4.d0*f(2)-f(3))/2.d0            
          fx(2)=(-2.d0*f(1)-3.d0*f(2)+6.d0*f(3)-f(4)) /6.d0  
	      fx(3)=b1*(f(4)-f(2)) -b2*(f(5)-f(1))
          fx(nx-2)=b1*(f(nx-1)-f(nx-3)) -b2*(f(nx)-f(nx-4))
          fx(nx-1)=(f(nx-3)-6.d0*f(nx-2)+3.d0*f(nx-1) +2.d0*f(nx))/6.d0
          fx(nx)=(f(nx-2)-4.d0*f(nx-1)+3.d0*f(nx))/2.d0
       end
!-----------------------------------------------------------------------
    subroutine dx3d(nx,ny,nz,f,fx)
       use precision_EC
 	   implicit none
	   integer:: nx,ny,nz,i,j,k
	   real(PRE_EC):: f(nx,ny,nz),fx(nx,ny,nz),f1d(nx),fx1d(nx)
	   do k=1,nz
	   do j=1,ny
	   do i=1,nx
	    f1d(i)=f(i,j,k)
	   enddo
	   call dx0(nx,f1d,fx1d)
	   do i=1,nx
	   fx(i,j,k)=fx1d(i)
	   enddo
	   enddo
	   enddo
	end subroutine dx3d   

    subroutine dy3d(nx,ny,nz,f,fy)
       use precision_EC
	   implicit none
	   integer:: nx,ny,nz,i,j,k
	   real(PRE_EC):: f(nx,ny,nz),fy(nx,ny,nz),f1d(ny),fx1d(ny)
	   do k=1,nz
	   do i=1,nx
	   do j=1,ny
	    f1d(j)=f(i,j,k)
	   enddo
	   call dx0(ny,f1d,fx1d)
	   do j=1,ny
	   fy(i,j,k)=fx1d(j)
	   enddo
	   enddo
	   enddo
	end subroutine dy3d   

    subroutine dz3d(nx,ny,nz,f,fz)
       use precision_EC
	   implicit none
	   integer:: nx,ny,nz,i,j,k
	   real(PRE_EC):: f(nx,ny,nz),fz(nx,ny,nz),f1d(nz),fx1d(nz)
	   do j=1,ny
	   do i=1,nx
	   do k=1,nz
	    f1d(k)=f(i,j,k)
	   enddo
	   call dx0(nz,f1d,fx1d)
	   do k=1,nz
	   fz(i,j,k)=fx1d(k)
	   enddo
	   enddo
	   enddo
	end subroutine dz3d   