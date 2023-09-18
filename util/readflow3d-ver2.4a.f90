!-----------------------------------------------------------
! Copyright by LiXinliang
! Ver 1.2   将格心坐标转换到格点上
! Ver 1.3 可读取flow3d.dat (tecplot格式文件)及flow.dat (无格式文件)
! Ver 1.4  读取 ver 0.8 以上版本的数据 （包含k,w,mu_t等）
! Ver 1.5 读取ver0.82以上的数据
! Ver 1.6 读取ver 0.84以上版本的数据（.inp文件）
! Ver 1.7 Plot 3D flow
! Ver 1.8 Plot Cp and Cf on the wall
! Ver 1.8a, vt 无需乘以Re
! Ver 2.3, 可读取BC_user
! Ver 2.4a, 修正了read wall_dist() 中的Bug
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Wall_Turbo=201             
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401      
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   integer,parameter::  BC_Zero=0     
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902      
   integer,parameter:: BC_USER_Blow_Suction_Wall=903    
end module Const_Variables

 
  module Global_Variables
  use Const_Variables
  implicit none
 ! global parameter (for all blocks) 
  real(PRE_EC),save:: Ma,Re,gamma,T_inf,p00
  integer,save:: Num_Block, Mesh_File_Format, BC_number,BC_type(100)
!-----------------------------------------------------------------------------------------
! 网格连接信息 
  TYPE BC_MSG_TYPE
   integer::   ib,ie,jb,je,kb,ke,bc
  END TYPE BC_MSG_TYPE

!  核心变量——每块网格存储的信息 （全局变量）
   TYPE Block_TYPE           !  variables for each block 
     integer :: Block_no,nx,ny,nz,subface
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc  ! coordinates of cell center, 网格中心坐标 
     real(PRE_EC),pointer,dimension(:,:,:):: dc,uc,vc,wc,Tc,pc,mut,vt,kt,wt
	 real(PRE_EC),pointer,dimension(:,:,:):: d,u,v,w,T,p,mut1,vt1
	 real(PRE_EC),pointer,dimension(:,:,:):: x,y,z,dw
     TYPE(BC_MSG_TYPE),pointer,dimension(:):: bc_msg
   End TYPE Block_TYPE  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  
  end module Global_Variables
!------------------------------------------------------------------


!-------------------------------------------
     program readflow   
     use Global_Variables
     implicit none
     integer:: If3d,m
	  character(len=50):: filename
    
	   call init
   
    print*, "Plot boundary plane ;  bounary-xxx.dat  Grid value, boundaryc-xxx.dat  Grid center value"
    print*, "xxx:  002 BC_Wall,  003BC_Symmetry, 004 BC_Farfield, 005 BC_Inflow, 006 BC_Outflow"
     
     do m=1,Bc_number
           write(filename,"('boundary-'I3.3'.dat')") Bc_type(m)
           call plot_p(filename,BC_type(m)) 

            write(filename,"('boundaryc-'I3.3'.dat')") Bc_type(m)
            call plot_center(filename,BC_type(m)) 
    enddo
    
    	   
	   print*, "if you want to plot 2D plane flow ?  1 for yes, 0 for no"
	   read(*,*) If3d
       if(If3d .eq. 1) then
	    call Plot_2D
       endif
      
	 
	   print*, "if you want to plot 3D flow ?  1 for yes, 0 for no"
	   read(*,*) If3d
	   if(If3d .eq. 1) then
	    call Plot_3D
	   endif
     end

!=====================================================================================
     subroutine plot_3D
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition
	 open(100,file="flow3d-tec.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx, " j= ", B%ny , " k= ", B%nz  
      do k=1,B%nz
      do j=1,B%ny
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
     enddo
     close(100)
    end



     subroutine plot_2D
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,i0,j0,k0,Iflag

	 print*, "If you want to plot 2D section ? (0/1) "
	 read(*,*)  Iflag
	 if (Iflag ==0) return
	 print*, "Please input i0, j0,k0"
	 read(*,*) i0,j0,k0

     print*, "... flow2d-i.dat..."
     i=i0
	 open(100,file="flow2d-i.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone j=", B%ny , " k= ", B%nz  
      do k=1,B%nz
      do j=1,B%ny
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)

       print*, "... flow2d-j.dat..."
      
     j=j0
  	 open(100,file="flow2d-j.dat")
	 write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx,  " k= ", B%nz  
      do k=1,B%nz
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)

     print*, " ... flow2d-k.dat..."
	 k=k0
	 open(100,file="flow2d-k.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx, " j= ", B%ny   
      do j=1,B%ny
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)
  
    end




   
    subroutine plot_p(filename,bCondition)
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition
     character(len=50):: filename

	
	 open(100,file=trim(filename))
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
     do n=1,B%subface
     Bc=>B%bc_msg(n)
     if(Bc%bc .eq. bCondition) then
      write(100,*) "zone i=", (Bc%ie-Bc%ib+1), " j= ",(Bc%je-Bc%jb+1) , " k= ", Bc%ke-Bc%kb+1  
      do k=Bc%kb,Bc%ke
      do j=Bc%jb,Bc%je
      do i=Bc%ib,Bc%ie
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
     endif
 
     enddo
     enddo
     close(100)
    end

    

 
   
!------------------------------------------------------------------------------------

   subroutine plot_center(filename,bCondition)
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition,ib,ie,jb,je,kb,ke
      real(PRE_EC):: Cpw,Cf,Cfx,mu0,Tsb

	 character(len=50):: filename

     open(100,file=trim(filename))
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     if(BCondition==BC_Wall) then
	 open(101,file="CpCf.dat")
	 write(101,*) "variables=x,y,z,Cp,Cf,Cfx,Pw"
	 endif
  
     Tsb=110.4d0/T_inf
	
	 do m=1,Num_Block
     B=>Block(m)
     do n=1,B%subface
     Bc=>B%bc_msg(n)
     if(Bc%bc .eq. bCondition) then
       ib=Bc%ib; ie=Bc%ie-1; jb=Bc%jb; je=Bc%je-1; kb=Bc%kb; ke=Bc%ke-1
       if(Bc%ib .eq. Bc%ie) then
	     if(Bc%ib .eq. 1) then
		  ie=1
		 else
		  ib=ib-1
		 endif
		endif

       if(Bc%jb .eq. Bc%je) then
	     if(Bc%jb .eq. 1) then
		  je=1
		 else
		  jb=jb-1
		 endif
		endif

       if(Bc%kb .eq. Bc%ke) then
	     if(Bc%kb .eq. 1) then
		  ke=1
		 else
		  kb=kb-1
		 endif
		endif

			    
	  write(100,*) "zone i=", (ie-ib+1), " j= ",(je-jb+1) , " k= ", ke-kb+1  
  	  do k=kb,ke
      do j=jb,je
      do i=ib,ie
       write(100,"(13E18.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%dc(i,j,k),B%uc(i,j,k),B%vc(i,j,k), &
                     B%wc(i,j,k),B%Tc(i,j,k),B%pc(i,j,k),B%mut(i,j,k),B%vt(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
      
	  if(bCondition==BC_Wall) then
 	  write(101,*) "zone i=", (ie-ib+1), " j= ",(je-jb+1) , " k= ", ke-kb+1  
  	  do k=kb,ke
      do j=jb,je
      do i=ib,ie
       mu0=1.d0/Re*(1.d0+Tsb)*sqrt(B%Tc(i,j,k)**3)/(Tsb+B%Tc(i,j,k))
       Cpw=(B%pc(i,j,k)-p00)*2.d0
	   Cf=mu0*(sqrt(B%uc(i,j,k)**2+B%vc(i,j,k)**2+B%wc(i,j,k)**2)/B%dw(i,j,k))*2.d0
	   Cfx=mu0*(B%uc(i,j,k)/B%dw(i,j,k))*2.d0
	   write(101,"(13E18.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),Cpw,Cf,Cfx,B%dc(i,j,k)*B%Tc(i,j,k)
      enddo
      enddo
      enddo
      endif

     endif

     enddo
     enddo
     close(100)
	 close(101)
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
   integer :: i,j,k,m,km,nx,ny,nz,Num_Block1,ksub,tmp,NB1,Iflag
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B
   TYPE(BC_MSG_TYPE),pointer::Bc
   logical:: EX
   integer:: ib,ie,jb,je,kb,ke
!------------------------------------------------------------------
!    print*, "Please input Ma, Re, T_inf, Mesh_File_Format"
!    read(*,*) Ma, Re, T_inf,Mesh_File_Format


    call read_parameter_ec 

    print*, "Ma,Re,gamma,Tinf=", Ma,Re,gamma,T_inf
    print*, "Mesh_File_Format=", Mesh_File_Format

    p00=1.d0/(gamma*Ma*Ma)

! ---------node Coordinates----------------------------------------  
!  网格文件：PLOT3D格式；   
   print*, "read Mesh3d.dat... (PLOT3D Format)"
   if(Mesh_File_Format .eq. 0) then
    open(99,file="Mesh3d.dat",form="unformatted")
    read(99) Num_Block         ! 总块数
   else
    open(99,file="Mesh3d.dat")
    read(99,*) Num_Block         ! 总块数
   endif
   print*, "Num_Block=", Num_Block
   
   
      
    allocate(Block(Num_Block))             
    allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! 每块的大小
  if(Mesh_File_Format .eq. 0) then
   read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  else
   read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  endif

! 读取每块信息----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   ! nx,ny,nz 每块的大小
     nx=B%nx ; ny= B%ny ; nz=B%nz
! ----------  几何量 -----------------------------------------------
    allocate(B%xc(0:nx,0:ny,0:nz), B%yc(0:nx,0:ny,0:nz), B%zc(0:nx,0:ny,0:nz)) 
    allocate(B%dc(0:nx,0:ny,0:nz),B%uc(0:nx,0:ny,0:nz),B%vc(0:nx,0:ny,0:nz), &
       B%wc(0:nx,0:ny,0:nz),B%Tc(0:nx,0:ny,0:nz),B%pc(0:nx,0:ny,0:nz), &
	   B%mut(0:nx,0:ny,0:nz),B%vt(0:nx,0:ny,0:nz),B%Kt(0:nx,0:ny,0:nz),B%Wt(0:nx,0:ny,0:nz))
    allocate(B%x(nx,ny,nz),B%y(nx,ny,nz),B%z(nx,ny,nz))
    allocate(B%d(nx,ny,nz),B%u(nx,ny,nz),B%v(nx,ny,nz), &
             B%w(nx,ny,nz),B%T(nx,ny,nz),B%p(nx,ny,nz),B%mut1(nx,ny,nz),B%vt1(nx,ny,nz))
    allocate(B%dw(nx-1,ny-1,nz-1))

   B%mut(:,:,:)=0.d0
   B%vt(:,:,:)=0.d0
   B%mut1(:,:,:)=0.d0
   B%Kt(:,:,:)=0.d0
   B%Wt(:,:,:)=0.d0
   B%vt1(:,:,:)=0.d0

   print*, "m=", m
   print*, "nx,ny,nz=", nx,ny,nz
   if(Mesh_File_Format .eq. 0) then
	read(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else
  	read(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif
   print*, "read mesh ok "
!----------插值出网格中心点的数据-----------------------------------
   do k=1,B%nz-1
   do j=1,B%ny-1
   do i=1,B%nx-1
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
  enddo
 
  close(99)

    inquire( file="wall_dist.dat", exist=EX)
    if(EX) then
      open(100,file="wall_dist.dat",form="unformatted")
      do m=1,Num_Block
      B => Block(m)
      read(100) (((B%dw(i,j,k),i=1,B%nx-1),j=1,B%ny-1),k=1,B%nz-1)
     enddo
     close(100)
   endif
   
!---------------从flow3d.dat中读取场作为初值----------------------
  print*, "input 1 or 2,  1 read flow3d.dat,  2 read flow3d_average.dat "
  read(*,*) Iflag
  
  if(Iflag == 1) then   
	open(99,file="flow3d.dat",form="unformatted")
  else
	open(99,file="flow3d_average.dat",form="unformatted")
  endif
    
	 do m=1,Num_Block
     B => Block(m)
        read(99)   (((B%dc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%uc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%vc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%wc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
		           (((B%Tc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 
          do k=0, B%nz
		  do j=0, B%ny
		  do i=0, B%nx
		    B%pc(i,j,k)=p00*B%dc(i,j,k)*B%Tc(i,j,k)
		  enddo
		  enddo
		  enddo

	 enddo
	close(99)
    
	 inquire(file="SA3d.dat",exist=EX)
	 if(EX) then
     open(100,file="SA3d.dat",form="unformatted")
      do m=1,Num_Block
        B => Block(m)
        read(100)   (((B%vt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 
 	  enddo
      close(100)
     endif

    inquire(file="SST3d.dat",exist=EX)
	 if(EX) then
     open(100,file="SST3d.dat",form="unformatted")
     do m=1,Num_Block
     B => Block(m)
        read(100)   (((B%Kt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz),  (((B%Wt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz)
 	 enddo
 	 close(100)
    endif
!=================================================================================
   call comput_value_in_mesh


!---------------------------------------------------- 
   BC_number=0
   BC_type=0
   
!--------------------------------------------------
    print*, "read bc3d.inp"
    open(88,file="bc3d.inp")
    read(88,*)
    read(88,*) Num_Block1
    if(Num_Block1 .ne. Num_Block) then
      print*, "Error!  Block number in bc2d.in is not equal to that in Mesh3d.dat !"
      stop
    endif
    
    do m=1,Num_Block
     B => Block(m)
     read(88,*)
     read(88,*)
     read(88,*) B%subface   !number of the subface in the Block m
!-------------------------------------------------------------------------------
     allocate(B%bc_msg(B%subface))
     do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     read(88,*) ib,ie,jb,je,kb,ke,Bc%bc
	 if(Bc%bc .lt. 0) then
  	   read(88,*)
  	  endif
  	  
 ! Find boundary type 	  
   if(Bc%bc .ge. 0) then
      if(BC_number .eq. 0) then
       BC_number=1
       BC_type(1)=Bc%bc
      else
       km=1
       do k=1,BC_number
        if(Bc%bc .eq. BC_type(k)) km=0
       enddo
      if(km .eq. 1) then
       Bc_number=Bc_number+1
       Bc_type(Bc_number)=Bc%bc
      endif
     endif
    endif



	  Bc%ib=min(abs(ib),abs(ie)) ;  Bc%ie=max(abs(ib),abs(ie))
	  Bc%jb=min(abs(jb),abs(je)) ;  Bc%je=max(abs(jb),abs(je))
	  Bc%kb=min(abs(kb),abs(ke)) ;  Bc%ke=max(abs(kb),abs(ke))
     enddo
     enddo
     close(88)
	 print*, "read bc3d.in ok ..."


 end subroutine init  

   subroutine comput_value_in_mesh
   use  Global_Variables
   implicit none
   integer:: nx,ny,nz,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do m=1, NUM_BLOCK
   B=> Block(m)
   nx=B%nx; ny=B%ny; nz=B%nz
!  将格心处的物理量插值到格点上   
   do k=1,B%nz
   do j=1,B%ny
   do i=1,B%nx
     B%d(i,j,k)=(B%dc(i-1,j-1,k-1)+B%dc(i,j-1,k-1)+B%dc(i-1,j,k-1)+B%dc(i,j,k-1)  &
	             +B%dc(i-1,j-1,k)+B%dc(i,j-1,k)+B%dc(i-1,j,k)+B%dc(i,j,k))/8.d0
     B%u(i,j,k)=(B%uc(i-1,j-1,k-1)+B%uc(i,j-1,k-1)+B%uc(i-1,j,k-1)+B%uc(i,j,k-1)  &
	             +B%uc(i-1,j-1,k)+B%uc(i,j-1,k)+B%uc(i-1,j,k)+B%uc(i,j,k))/8.d0
     B%v(i,j,k)=(B%vc(i-1,j-1,k-1)+B%vc(i,j-1,k-1)+B%vc(i-1,j,k-1)+B%vc(i,j,k-1)  &
	             +B%vc(i-1,j-1,k)+B%vc(i,j-1,k)+B%vc(i-1,j,k)+B%vc(i,j,k))/8.d0
     B%w(i,j,k)=(B%wc(i-1,j-1,k-1)+B%wc(i,j-1,k-1)+B%wc(i-1,j,k-1)+B%wc(i,j,k-1)  &
	             +B%wc(i-1,j-1,k)+B%wc(i,j-1,k)+B%wc(i-1,j,k)+B%wc(i,j,k))/8.d0

     B%T(i,j,k)=(B%Tc(i-1,j-1,k-1)+B%Tc(i,j-1,k-1)+B%Tc(i-1,j,k-1)+B%Tc(i,j,k-1)  &
	             +B%Tc(i-1,j-1,k)+B%Tc(i,j-1,k)+B%Tc(i-1,j,k)+B%Tc(i,j,k))/8.d0

     B%p(i,j,k)=(B%pc(i-1,j-1,k-1)+B%pc(i,j-1,k-1)+B%pc(i-1,j,k-1)+B%pc(i,j,k-1)  &
	             +B%pc(i-1,j-1,k)+B%pc(i,j-1,k)+B%pc(i-1,j,k)+B%pc(i,j,k))/8.d0

     B%mut1(i,j,k)=(B%mut(i-1,j-1,k-1)+B%mut(i,j-1,k-1)+B%mut(i-1,j,k-1)+B%mut(i,j,k-1)  &
	             +B%mut(i-1,j-1,k)+B%mut(i,j-1,k)+B%mut(i-1,j,k)+B%mut(i,j,k))/8.d0

     B%vt1(i,j,k)=(B%vt(i-1,j-1,k-1)+B%vt(i,j-1,k-1)+B%vt(i-1,j,k-1)+B%vt(i,j,k-1)  &
	             +B%vt(i-1,j-1,k)+B%vt(i,j-1,k)+B%vt(i-1,j,k)+B%vt(i,j,k))/8.d0

   enddo
   enddo
   enddo
   enddo
   end



!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_ec 
   use Global_Variables
   implicit none
   real(PRE_EC),parameter:: PI=3.14159265358979d0
   real(PRE_EC):: R0, a0, d0, mu0,mu1
   
   real(PRE_EC)::  AoA, AoS, p_outlet, t_end, &
	    PrL, PrT, &
        dt_global,CFL,dtmax,dtmin,Time_Method, &
		w_LU, Twall,Kt_inf,Wt_inf,&
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, &
        Ref_S,Ref_L,Centroid(3), &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
 		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ
  
   integer:: Kstep_save, Kstep_average, Iflag_turbulence_model,Iflag_init,If_viscous,  &
           Iflag_local_dt,  If_Residual_smoothing, If_dtime_mesh,  &
		   Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, Pre_Step_Mesh(3), &
           Kstep_show,Kstep_smooth,Kstep_init_smooth, Bound_Scheme,  &
           Num_Mesh,IF_Debug, NUM_THREADS, Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug(4),  &
           IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, IF_Innerflow,Iflag_savefile

	namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
	    gamma, PrL, PrT, &
	    Kstep_save, &
	    Iflag_turbulence_model,Iflag_init,If_viscous,  &
        Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method, &
		If_Residual_smoothing,w_LU,If_dtime_mesh,  &
        Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, &
		Mesh_File_Format,Kstep_show,Kstep_average,Kstep_smooth,Kstep_init_smooth,  &
        Num_Mesh,T_inf,Twall,Kt_inf,Wt_inf,IF_Debug,NUM_THREADS,  &
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme, &
        Pre_Step_Mesh,Ref_S,Ref_L,Centroid,Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug, &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
        IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ, &
		IF_Innerflow, Iflag_savefile


	

!---- default--------
        Ma=1.d0
        Re=100.d0
    	gamma=1.4d0
        T_inf=288.15
        Mesh_File_Format=0
        
!---------------------------------
	open(99,file="control.ec")
	read(99,nml=control_ec)
    close(99)
 
 !---- convert parameters ----------------------
 ! Ref_medium_usrdef==0 使用默认介质 (Ma=1, 根据总温、总压计算 Re) ； ==1 使用自定义介质 （人为输入Ma, Re等）
    if( (IF_TurboMachinary ==1 .or. IF_Innerflow ==1) .and.  Ref_medium_usrdef == 0) then   ! 默认空气介质，计算Mach数， Reynolds数
      
	  T_inf=Turbo_T0  ! 参考温度 （来流总温）
      gamma=1.4d0    ! 
	  PrL=0.7d0   ! Prandtl数
	  PrT=0.9d0
      R0= 287.06d0   ! 空气的气体常数R
	  a0= sqrt(gamma*R0*Turbo_T0)    ! 参考温度下的声速 
	  mu0=1.179d-5     ! 空气粘性系数 (288.15K)  
      mu1=mu0* sqrt((Turbo_T0/288.15d0)**3)*(288.15d0+110.4d0)/(Turbo_T0+110.4d0)  ! 参考温度下的空气粘性系数
      d0=Turbo_P0/(R0*Turbo_T0)
	  Re=d0*a0*Turbo_L0/mu1    ! 参考温度下，以声速运动的Reynolds数
	  Ma=1.d0     ! Mach数    （以声速作为参考速度，因而参考Mach数为1）
      Turbo_w= 2.d0*PI*Turbo_w/(a0/Turbo_L0)   ! 无量纲角速度 Turbo_W（转/秒）
      P_outlet=P_outlet/Turbo_P0    ! 背压 （无量纲）
	endif

 !  print*, "Re, Ma, gamma, T_inf=", Re, Ma, gamma, T_inf
 !  print*, "Mesh_File_Format=", Mesh_File_Format


 !----------------------------------------------
 end
!--------------------------------------------------
 


