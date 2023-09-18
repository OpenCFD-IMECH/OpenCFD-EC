!-----------------------------------------------------------
! Copyright by LiXinliang
! Convert (X,Y,Z) --> (X, Z, -Y)  coordination
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
end module Const_Variables

 
  module Global_Variables
  use Const_Variables
  implicit none
  integer,save:: Num_Block, Mesh_File_Format, IFC
!-----------------------------------------------------------------------------------------

   TYPE Block_TYPE           !  variables for each block 
	 real(PRE_EC),pointer,dimension(:,:,:):: x,y,z
   End TYPE Block_TYPE  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  
  end module Global_Variables
!------------------------------------------------------------------




!-------------------------------------------
     program convert
     use Global_Variables
     implicit none
      Type (Block_TYPE),pointer:: B
      integer,allocatable,dimension(:):: NI,NJ,NK
      integer:: i,j,k,m,nx,ny,nz

   print*, "Convert (X,Y,Z) to (X, Z, -Y)  or to (X, Y, -Z) coordination,  New File is Mesh3d-new.dat"
    print*, "read Mesh3d.dat... (PLOT3D Format)"
    print*, "please input 0 or 1, 1 Format, 0 Unformat"
    read(*,*) Mesh_File_Format
    print*, " please input 1 (counter-clockwise rotate) or -1  (clock-wise rotate) "
	print*, " 1 (X,Y,Z) --> (X,Z,-Y);   -1 (X,Y,Z) --> (X, -Z, Y) "
	read(*,*) IFC
    
	if(IFC .ne. 1  .and. IFC .ne. -1) then
	 print*, "Please input 1 or -1"
	 stop
	endif

	if(Mesh_File_Format .eq. 0) then
     open(99,file="Mesh3d.dat",form="unformatted")
     read(99) Num_Block         
    else
     open(99,file="Mesh3d.dat")
     read(99,*) Num_Block         
    endif
      
    allocate(Block(Num_Block))             
    allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! 每块的大小
  
  if(Mesh_File_Format .eq. 0) then
   read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  else
   read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  endif

! read mesh ----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     nx=NI(m); ny=NJ(m) ; nz=NK(m)   
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
  close(99)
!-------------------------------------------------------
!  write file

   print*, "Convert  coordination ..."
	
	if(Mesh_File_Format .eq. 0) then
     open(99,file="Mesh3d-new.dat",form="unformatted")
     write(99) Num_Block         
    else
     open(99,file="Mesh3d-new.dat")
     write(99,*) Num_Block         
    endif
  if(Mesh_File_Format .eq. 0) then
   write(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  else
   write(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  endif

    do m=1,Num_Block
     B => Block(m)
     nx=NI(m); ny=NJ(m) ; nz=NK(m)   
    if(Mesh_File_Format .eq. 0) then
 	 write(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                 (((IFC*B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                 (((-1*IFC*B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else
  	 write(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((IFC*B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((-1*IFC*B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif
   enddo
  close(99)
   
   end

