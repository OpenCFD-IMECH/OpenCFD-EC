!-----------------------------------------------------------
! Check mesh ver 1.1;  Code by Li Xinliang
! Plot all boundary mesh (not include the inner boundary) 
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
  integer,parameter:: BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Outflow=6, BC_Periodic=501     ! 与Griggen .inp文件的定义可能有所区别，请注意
 end module Const_Variables

 
!-----------------------------------------------------------------------------------------
 module Global_Variables
  use Const_Variables
  implicit none
  TYPE BC_MSG_TYPE
   integer::   ib,ie,jb,je,kb,ke,bc
  END TYPE BC_MSG_TYPE

  TYPE Block_TYPE           !  variables for each block 
     integer :: Block_no,nx,ny,nz,subface
	 real(PRE_EC),pointer,dimension(:,:,:):: x,y,z
     real(PRE_EC),pointer,dimension(:,:,:):: dw
     TYPE(BC_MSG_TYPE),pointer,dimension(:):: bc_msg
  End TYPE Block_TYPE  
  TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  integer:: Num_Block
 end module Global_Variables
!------------------------------------------------------------------

     program check_mesh   
     use Global_Variables
     implicit none
       call read_mesh
	   call read_inp
	   call check_mesh_size
       call plot_bc
	   call plot_dw
     end

!=====================================================================================
    
   subroutine check_mesh_size
    use  Global_Variables
    implicit none
    integer:: Total_Cell, Total_grid,m
    Type (Block_TYPE),pointer:: B
    
	 Total_Cell=0
	 Total_grid=0
    
	 do m=1,Num_Block
     B => Block(m)
     Total_Cell=Total_Cell+(B%nx-1)*(B%ny-1)*(B%nz-1)
	 Total_grid=Total_grid+B%nx*B%ny*B%nz
	 enddo
	 print*, "Total Cell number=", Total_Cell
	 print*, "Total Grid number=", Total_grid
     end

!------------------------------------------------------------
    subroutine plot_dw
    use  Global_Variables
    implicit none
    integer:: i,j,k,m,nx,ny,nz
    real(PRE_EC),pointer,dimension(:,:,:):: dw
    real(PRE_EC):: xc,yc,zc

	Type (Block_TYPE),pointer:: B
    logical ext
     Inquire(file="wall_dist.dat",Exist=EXT)
 	 if(.not. ext) return
     print*, "Find wall_dist.dat, read it ..."     
	 open(99,file="wall_dist.dat",form="unformatted")
 	 open(100,file="wall_dist_tec.dat")
	 write(100,*) "variables=x,y,z,dw"
	 do m=1,Num_Block
      B => Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
	  allocate(dw(nx-1,ny-1,nz-1))
      read(99) (((dw(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
      write(100,*) "zone i=", nx-1, " j= ", ny-1, " k= ",nz-1
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	     xc= (B%x(i,j,k)+B%x(i,j+1,k)+B%x(i,j,k+1)+B%x(i,j+1,k+1)+ &
              B%x(i+1,j,k)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1)+B%x(i+1,j+1,k+1))*0.125    
	     yc= (B%y(i,j,k)+B%y(i,j+1,k)+B%y(i,j,k+1)+B%y(i,j+1,k+1)+ &
              B%y(i+1,j,k)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1)+B%y(i+1,j+1,k+1))*0.125    
	     zc= (B%z(i,j,k)+B%z(i,j+1,k)+B%z(i,j,k+1)+B%z(i,j+1,k+1)+ &
              B%z(i+1,j,k)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1)+B%z(i+1,j+1,k+1))*0.125    
      write(100,"(4F20.10)") xc,yc,zc,dw(i,j,k)
	  enddo
	  enddo
	  enddo
	  enddo
	  close(99)
	  close(100)
    end

  
    subroutine plot_bc
    use  Global_Variables
    implicit none
    integer:: Nbt, boundtype(100),m,ksub,i,j,k,kt,bct,n
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer::Bc
    character(len=50):: filename
	 Nbt=0
	 boundtype(:)=0

!  统计边界条件的数目
	do m=1,Num_Block
     B => Block(m)
     do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
      
	  if(Bc%bc .ne. -1) then
	   do k=1,Nbt
	    if(Bc%bc .eq. boundtype(k)  ) goto 100
	   enddo
       Nbt=Nbt+1           ! 新的边界条件
       boundtype(Nbt)=Bc%bc
100   continue
      endif
     enddo
    enddo


     print*, "Total types of the boundary is:", Nbt
	 print*, "Boundary types are:", boundtype(1:Nbt)
	 print*, "--------------------------------"

	 
   do kt=1,Nbt
    bct=boundtype(kt)
	 write(filename,"('bc-'I3.3)") bct
	 open(100,file=filename)
     write(100,*) "variables=x,y,z"
  
     do m=1,Num_Block
     B=>Block(m)
     do n=1,B%subface
     Bc=>B%bc_msg(n)
     if(Bc%bc .eq. bct) then
      write(100,*) "zone i=", (Bc%ie-Bc%ib+1), " j= ",(Bc%je-Bc%jb+1) , " k= ", Bc%ke-Bc%kb+1  
      do k=Bc%kb,Bc%ke
      do j=Bc%jb,Bc%je
      do i=Bc%ib,Bc%ie
	   write(100,"(3E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k)
      enddo
      enddo
      enddo
     endif
    enddo
    enddo
   close(100)
  enddo
  end subroutine plot_bc

! -------------------------------------------------  
   subroutine read_mesh
    use  Global_Variables
    implicit none
    integer:: Mesh_File_Format,m,i,j,k,tmp,nx,ny,nz
    integer,allocatable,dimension(:):: NI,NJ,NK
    Type (Block_TYPE),pointer:: B
    print*, "please input the format of Mesh3d.dat, 0 for unformatted, 1 for formatted"
	read(*,*) Mesh_File_Format

    print*, "read Mesh3d.dat... (PLOT3D Format)"
   
   if(Mesh_File_Format .eq. 0) then
    open(99,file="Mesh3d.dat",form="unformatted")
    read(99) Num_Block         ! Total Block Number
   else
    open(99,file="Mesh3d.dat")
    read(99,*) Num_Block         ! 
   endif
   
    allocate(Block(Num_Block))             
    allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! Size of each block
  
   if(Mesh_File_Format .eq. 0) then
    read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
   else
    read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
   endif
    
	print*, "Num_BLock=",Num_Block
	 

    do m=1,Num_Block
     B => Block(m)
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   
     nx=B%nx ; ny= B%ny ; nz=B%nz
     print*, "block ", m, " Cell number", (nx-1)*(ny-1)*(nz-1)

	allocate(B%x(nx,ny,nz),B%y(nx,ny,nz),B%z(nx,ny,nz))
 
    if(Mesh_File_Format .eq. 0) then
	read(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    else
  	read(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif

  enddo

 end subroutine read_mesh 

!---------------------------
   subroutine read_inp
    use  Global_Variables
    implicit none
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer::Bc
    integer:: Num_Block1,m,ksub,ib,ie,jb,je,kb,ke

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

	  Bc%ib=min(abs(ib),abs(ie)) ;  Bc%ie=max(abs(ib),abs(ie))
	  Bc%jb=min(abs(jb),abs(je)) ;  Bc%je=max(abs(jb),abs(je))
	  Bc%kb=min(abs(kb),abs(ke)) ;  Bc%ke=max(abs(kb),abs(ke))
     enddo
     enddo
     close(88)
	 print*, "read bc3d.in ok ..."
   end subroutine read_inp

