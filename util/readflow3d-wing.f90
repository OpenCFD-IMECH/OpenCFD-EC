
  module Global_Variables
!  use Const_Variables
  implicit none
  real*8,save:: Ma,Re,gamma,p00,AoA
  integer:: Num_Block
   TYPE Block_TYPE           !  variables for each block 
     integer:: nx,ny,nz
     real*8,pointer,dimension(:,:,:):: xc,yc,zc   
     real*8,pointer,dimension(:,:,:):: d,u,v,w,T,p
   End TYPE Block_TYPE 
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  end module Global_Variables
 
!-------------------------------------------
     program readflow   
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     integer::  i,j,k,m
     real*8:: xb,xe
     call init
     print*, "Plot 2d sections ..."
     print*, "please input k"
     read(*,*) k
     print*, "z=", Block(1)%zc(1,1,k)
     open(100,file="flow2d-z.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p"
     do m=1,Num_Block
     B=>Block(m)
     write(100,*) "zone i=", B%nx+1, " j= ",B%ny+1 
     do j=0,B%ny
     do i=0,B%nx
     write(100,"(9f16.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%d(i,j,k),B%u(i,j,k),B%v(i,j,k), &
                     B%w(i,j,k),B%T(i,j,k),B%p(i,j,k)
     enddo
     enddo
     enddo
     close(100)

     print*, "write Cp ..."

     open(100,file="Cp.dat")
     write(100,*) "variables=x,Cp"
     do m=2,3
     B=>Block(m)
     write(100,*) "zone i=", B%nx+1 
     if(m .eq. 2) then
       xb=B%xc(B%nx-1,1,k) ; xe=B%xc(1,1,k)
     else
       xb= B%xc(1,1,k) ; xe=B%xc(B%nx-1,1,k)
     endif

     do i=0,B%nx
     write(100,"(9f16.9)") (B%xc(i,1,k)-xb)/(xe-xb),-2.d0*(B%p(i,1,k)-p00)
     enddo
     enddo
     close(100)

     print*, "please input j"
     read(*,*) j
     open(100,file="flow2d-j.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p"
     do m=2,3
     B=>Block(m)
     write(100,*) "zone i=", B%nx+1, " j= ",B%nz+1 
     do k=0,B%nz
     do i=0,B%nx
     write(100,"(9f16.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%d(i,j,k),B%u(i,j,k),B%v(i,j,k), &
                     B%w(i,j,k),B%T(i,j,k),B%p(i,j,k)
     enddo
     enddo
     enddo
     close(100)


    
     end


!------------------------------------------------------------------------------     
! Read the message of the mesh and the initial flow;
! 读取网格，初始流场信息; 
! 分配内存变量；
! 计算几何量；
!------------------------------------------------------------------------------
   subroutine init
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx1,ny1,nz1,Num_Block1
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B
!------------------------------------------------------------------
!   open(100,file="control.in")
!   read(100,*)
!   read(100,*) Ma, Re, gamma,  AoA 
!   close(100)
    print*, "Please input Ma, Re, gamma, A0A"
	read(*,*) Ma, Re,gamma, AoA

   p00=1.d0/(gamma*Ma*Ma)
! ---------node Coordinates----------------------------------------  
!  网格文件：PLOT3D格式；   
   print*, "read Mesh3d.dat... (PLOT3D Format)"
   open(99,file="Mesh3d.dat")
   read(99,*) Num_Block         ! 总块数
   allocate(Block(Num_Block))             
   allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! 每块的大小
   read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
! 读取每块信息----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   ! nx,ny,nz 每块的大小
     nx1=B%nx ; ny1= B%ny ; nz1=B%nz
! ----------  几何量 -----------------------------------------------
    allocate(B%xc(0:nx1,0:ny1,0:nz1), B%yc(0:nx1,0:ny1,0:nz1), B%zc(0:nx1,0:ny1,0:nz1)) 
    allocate(B%d(0:nx1,0:ny1,0:nz1),B%u(0:nx1,0:ny1,0:nz1),B%v(0:nx1,0:ny1,0:nz1), &
       B%w(0:nx1,0:ny1,0:nz1),B%T(0:nx1,0:ny1,0:nz1),B%p(0:nx1,0:ny1,0:nz1))
    enddo

!---------------从flow3d.dat中读取场作为初值----------------------
        open(99,file="flow3d.dat")
         print*, "Init from 'flow3d.dat' "
         read(99,*)
        do m=1,Num_Block
          B => Block(m)
          read(99,*)
          do k=0,B%nz
	  do j=0,B%ny
	  do i=0,B%nx
	   read(99,*) B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%d(i,j,k),B%u(i,j,k),B%v(i,j,k), &
                      B%w(i,j,k),B%T(i,j,k),B%p(i,j,k)
          enddo
          enddo
          enddo
        enddo
      close(99)
 
 end subroutine init  


