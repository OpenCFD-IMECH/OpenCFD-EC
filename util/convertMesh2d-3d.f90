!  PLOT3D mesh file convert :  3D --> 2D  or 2D---> 3D
  module mesh
    implicit none
    integer,parameter:: PRE_EC=8
    integer,save:: Num_Block, Mesh_File_Format
    integer,allocatable,dimension(:):: NI, NJ,NK
    real(PRE_EC),allocatable,dimension(:,:):: x2d, y2d
    real(PRE_EC),allocatable,dimension(:,:,:):: x3d, y3d, z3d
    
    end module 

!-----------------------------------------------------------
    program conver_mesh
    implicit none
    integer::  Iflag_convert
    print*, "convert 2D mesh to 3D mesh;  or convert 3D mesh to 2D mesh.  PLOT3D type "
    print*, " input 1 convert 2D to 3D;   input 2 convert 3D to 2D"
    read(*,*) Iflag_convert
    
    if (Iflag_convert .eq. 1) then
      call convert2d_3d
    else
      call convert3d_2d
    endif
    end
  
 !----------------------------------------------------------------------------------   
    subroutine convert2d_3d
     use mesh
   implicit none
    integer:: nx,ny, nz, m,i,j,k
   real(PRE_EC),allocatable,dimension(:):: zz
    real(PRE_EC):: Lz
    
    print*, "convert 2d to 3d,   input Mesh2d.dat, output Mesh3d.dat"
    print*, "read Mesh2d.dat,  input  1 formatted,  0 unformatted"
    read(*,*) Mesh_File_Format
    print*, "please input nz,  Lz"
    read(*,*) nz, Lz
    allocate(zz(nz))
    
    do k=1,nz
      zz(k)=Lz*(k-1.d0)/(nz-1.d0)
    enddo
    
  ! read Mesh2d,  write Mesh3d  
    if(Mesh_File_Format .eq. 0) then
      open(99,file="Mesh2d.dat",form="unformatted")
      read(99) Num_Block
      allocate(NI(Num_block), NJ(Num_block))
      read(99) ( NI(m),  NJ(m) ,  m=1,Num_Block)
   
      open(100,file="Mesh3d.dat",form="unformatted")
      write(100) Num_Block
      write(100) ( NI(m), NJ(m) , nz ,  m=1,Num_Block)
   
    else
       open(99,file="Mesh2d.dat")
       read(99,*) Num_Block
       allocate(NI(Num_block), NJ(Num_block))
       read(99,*) ( NI(m), NJ(m), nz , m=1,Num_Block)
    
      open(100,file="Mesh3d.dat")
      write(100,*) Num_Block
      write(100,*) ( NI(m), NJ(m) , nz ,  m=1,Num_Block)
  
    endif

     do m=1,Num_Block
        nx=NI(m) ;  ny=NJ(m)
        allocate(x2d(nx,ny), y2d(nx,ny))
            
     if(Mesh_File_Format .eq. 0) then
       read(99)  (( x2d(i,j), i=1,nx),j=1,ny), (( y2d(i,j), i=1,nx),j=1,ny) 
       write(100)   ((( x2d(i,j), i=1,nx),j=1,ny),k=1,nz)  ,  ((( y2d(i,j), i=1,nx),j=1,ny),k=1,nz) , &
                          ((( zz(k), i=1,nx),j=1,ny),k=1,nz)        
    else
       read(99,*)  (( x2d(i,j), i=1,nx),j=1,ny), (( y2d(i,j), i=1,nx),j=1,ny) 
       write(100,*)   ((( x2d(i,j), i=1,nx),j=1,ny),k=1,nz)  ,  ((( y2d(i,j), i=1,nx),j=1,ny),k=1,nz) , &
                            ((( zz(k), i=1,nx),j=1,ny),k=1,nz)        
    endif
 
      deallocate(x2d, y2d)
    enddo
      deallocate(zz)
    end
         
!---------------------------------------------------------------------------------------------------------------       
 !----------------------------------------------------------------------------------   
    subroutine convert3d_2d
    use mesh
    implicit none
    integer:: nx,ny, nz, m,i,j,k
    
    print*, "convert 3d to 2d,   input Mesh3d.dat, output Mesh2d.dat"
    print*, "read Mesh2d.dat,  input  1 formatted,  0 unformatted"
   read(*,*) Mesh_File_Format

   
  ! read Mesh3d,  write Mesh2d  
    if(Mesh_File_Format .eq. 0) then
      open(99,file="Mesh3d.dat",form="unformatted")
      read(99) Num_Block
      allocate(NI(Num_block), NJ(Num_block), NK(Num_Block))
      read(99) ( NI(m), NJ(m), NK(m), m=1,Num_Block)
   
      open(100,file="Mesh2d.dat",form="unformatted")
      write(100) Num_Block
      write(100) (NI(m), NJ(m) ,  m=1,Num_Block)
   
    else
       open(99,file="Mesh3d.dat")
       read(99,*) Num_Block
       allocate(NI(Num_block), NJ(Num_block), NK(Num_Block))
       read(99,*) (NI(m), NJ(m), NK(m),  m=1,Num_Block)
     
      open(100,file="Mesh2d.dat")
      write(100,*) Num_Block
      write(100,*) ( NI(m), NJ(m) ,  m=1,Num_Block)
 
    endif

     do m=1,Num_Block
        nx=NI(m) ;  ny=NJ(m) ; nz=NK(m)
        allocate(x3d(nx,ny,nz), y3d(nx,ny,nz),z3d(nx,ny,nz))
            
     if(Mesh_File_Format .eq. 0) then
        read(99)   ((( x3d(i,j,k), i=1,nx),j=1,ny),k=1,nz)  ,  ((( y3d(i,j,k), i=1,nx),j=1,ny),k=1,nz) , &
                          ((( z3d(i,j,k), i=1,nx),j=1,ny),k=1,nz)        
        write(100)  (( x3d(i,j,1), i=1,nx),j=1,ny), (( y3d(i,j,1), i=1,nx),j=1,ny) 
   
    else
         read(99,*)   ((( x3d(i,j,k), i=1,nx),j=1,ny),k=1,nz)  ,  ((( y3d(i,j,k), i=1,nx),j=1,ny),k=1,nz) , &
                          ((( z3d(i,j,k), i=1,nx),j=1,ny),k=1,nz)        
        write(100,*)  (( x3d(i,j,1), i=1,nx),j=1,ny), (( y3d(i,j,1), i=1,nx),j=1,ny) 
 
    endif
      deallocate(x3d, y3d, z3d)
    enddo
    end
         
!---------------------------------------------------------------------------------------------------------------       
   