! In2inp, Transform BXCFD .in file to Gridgen .inp file  
! （网格连接信息）从BXCFD的 .in格式 转化为Gridgen .inp格式
! Copyright by Li Xinliang, lixl@imech.ac.cn
! Ver 1.0, 2012-7-11
! Ver 1.1, 2013-5-4
!------Types 定义子面和块两种数据结构 ----------------------------------------------------------------------
! 面： 属性有 面号、维数、连接信息
! 块： 属性有 块号、维数、面
!-----------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------
  module Def_block
  implicit none
  Integer:: NB    ! 网格块数

! 边界信息, Gridgen .inp格式
   TYPE BC_MSG_TYPE             ! 边界链接信息 
     integer:: ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
   END TYPE BC_MSG_TYPE
  
   TYPE Block_TYPE           !  块
     integer:: nx,ny,nz
	 integer::  subface  !  子面数目
     TYPE (BC_MSG_TYPE),dimension(:),pointer:: bc_msg   ! 子面 （从属于块）
   End TYPE Block_TYPE  
  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  
  end module Def_block

!============================================================================
  program Inp2In
  use Def_block
  implicit none
  print*, " Transform  from .in file to .inp file ......"
!----------------------------
  call read_bcin
  call trans_in_inp
  call write_inp
  end

  
!--------------------------------------------------------------------

  
 !----Mesh control message (bc2d.in)------------------------------------------
  subroutine read_bcin 
   use Def_block
   implicit none
   integer::m,ksub,Ia,NB1
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc

   print*, "read bc3d.in ......"
   open(88,file="bc3d.in")
   read(88,*)
   read(88,*)
   read(88,*) NB
   allocate (Block(NB))

   do m=1,NB
     B => Block(m)
     read(88,*)
     read(88,*) 
     read(88,*) B%subface   !number of the subface in the Block m
     read(88,*)
     allocate(B%bc_msg(B%subface))
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       read(88,*)  Bc%f_no, Bc%face, Bc%ist, Bc%iend, Bc%jst, Bc%jend,  Bc%kst, Bc%kend, Bc%neighb, Bc%subface, Bc%orient
     enddo
   enddo
   close(88)
   print*, "read bc3d.in OK"
   
   print*, "read Mesh3d.dat for nx,ny,nz"
   print*, "please input the format of Mesh3d.dat, 1 formatted, 2 unformatted, 0 not read"
   read(*,*) Ia
   
   if(Ia .eq. 1 ) then
     open(100,file="Mesh3d.dat")
     read(100,*) NB1
     if(NB1 .ne. NB) then
	  print*, "error ! NB in Mesh3d.dat is not the same as that in bc3d.in "
	  stop
	 endif 
	  read(100,*) ((Block(m)%nx,Block(m)%ny,Block(m)%nz),m=1,NB)
   
   else if (Ia .eq. 2 ) then
     open(100,file="Mesh3d.dat",form="unformatted")
     read(100) NB1
     if(NB1 .ne. NB) then
	  print*, "error ! NB in Mesh3d.dat is not the same as that in bc3d.in "
	  stop
	 endif 
	  read(100) ((Block(m)%nx,Block(m)%ny,Block(m)%nz),m=1,NB)
   else
	 do m=1,NB
	  Block(m)%nx=0; Block(m)%ny=0; Block(m)%nz=0
	 enddo
   endif


  end  subroutine read_bcin
!-------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------- 
! 将 .in文件转化为 .inp文件  
  subroutine trans_in_inp
   use Def_block
   implicit none
   integer,parameter:: BC_Wall_in=-10, BC_Farfield_in=-20, BC_Periodic_in=-30,BC_Symmetry_in=-50,BC_Outlet_in=-22   ! .in 关于边界条件的定义
   integer,parameter:: BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Outlet=6, BC_Periodic=501     ! 与Griggen .inp文件的定义可能有所区别，请注意
   integer:: m,ksub
   integer:: Lp(3),Ls(3),tmp  
   Type (Block_TYPE),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
    do m=1,NB
     B => Block(m)
      do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
        if(Bc%neighb .lt. 0) then    ! 物理边界
          Bc%ib=Bc%ist ; Bc%ie=Bc%iend ;  Bc%jb=Bc%jst ; Bc%je=Bc%jend ;  Bc%kb=Bc%kst ; Bc%ke=Bc%kend
		
		   if(Bc%neighb .eq. BC_Wall_in ) then 
		    Bc%bc=BC_Wall
		   else if (Bc%neighb .eq. BC_Farfield_in ) then 
		    Bc%bc=BC_Farfield
		   else if (Bc%neighb .eq. BC_Periodic_in ) then 
		    Bc%bc=BC_Periodic
		   else if (Bc%neighb .eq. BC_Symmetry_in ) then 
		    Bc%bc=BC_Symmetry
		   else if (Bc%neighb .eq. BC_Outlet_in ) then 
		    Bc%bc=BC_Outlet
 		   else
		    print*, "The boundary condition is not supported!", m, ksub
			print*, "Boundary type is", Bc%neighb
	       endif 
		else
		  Bc%bc=-1
          B1=>Block(Bc%neighb)
		  Bc1=>B1%bc_msg(Bc%subface)
		  Bc%nb1=Bc%neighb
		  Bc%face1=Bc1%face
		  Bc%f_no1=Bc1%f_no

          Bc%ib=Bc%ist ; Bc%ie=Bc%iend ;  Bc%jb=Bc%jst ; Bc%je=Bc%jend ;  Bc%kb=Bc%kst ; Bc%ke=Bc%kend
		  if(Bc%kst .eq. Bc% kend ) then
		    Bc%jb=-Bc%jst ; Bc%je=-Bc%jend
		  else
		    Bc%kb=-Bc%kst ; Bc%ke=-Bc%kend
          endif
          call get_orient(Bc%face,Bc1%face,Bc%orient,Lp,Ls)  ! Lp(k)==-1 交换次序, Ls(k)==-1 改变符号 
		   
		   Bc%ib1= Ls(1)*Bc1%ist ;  Bc%ie1=Ls(1)*Bc1%iend
		   if(Lp(1) .eq. -1) then
		    tmp=Bc%ib1; Bc%ib1=Bc%ie1; Bc%ie1=tmp
           endif
		   Bc%jb1= Ls(2)*Bc1%jst ;  Bc%je1=Ls(2)*Bc1%jend
		   if(Lp(2) .eq. -1) then
		    tmp=Bc%jb1; Bc%jb1=Bc%je1; Bc%je1=tmp
           endif
		   Bc%kb1= Ls(3)*Bc1%kst ;  Bc%ke1=Ls(3)*Bc1%kend
		   if(Lp(3) .eq. -1) then
		    tmp=Bc%kb1; Bc%kb1=Bc%ke1; Bc%ke1=tmp
           endif
		  
		endif
      enddo 
    enddo
 end  subroutine trans_in_inp

!------------------------------------------------------------------------------------------  
   subroutine write_inp
   use Def_block
   implicit none
   integer:: m,ksub
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   open(99,file="bc3d.inp")
    write(99,*) 1
	write(99,*) NB
	do m=1,NB
     B => Block(m)
     write(99,*) B%nx, B%ny, B%nz
	 write(99,*) "Block ", m 
	 write(99,*) B%subface 
	  do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
	   write(99,"(7I6)") Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc
       if(Bc%bc .eq. -1) then
	   write(99,"(7I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1
       endif
	  enddo
    enddo
   close(99)
 end
  


!  根据orient的值，确定.inp文件的连接次序 
!  见 OpenCFD-EC理论手册  
!      subroutine get_ijk_orient(i2,j2,i1,j1,ibegin,iend,jbegin,jend,orient,face1,face2)  ! bug bug but !!! (face1, face2)
      subroutine get_orient(face1,face2,orient,Lp,Ls)
      implicit none 
      integer:: l1,m1,l2,m2,tmp,face1,face2,orient,Lp(3),Ls(3),k0,k1,k2
 
 !           
          if(mod(face1,2) .eq. 1) then    ! i-, j- or k- 面
            l1=2 ; m1= 1                  ! l是第2个下标， m是第1个下标
          else                            ! i+, j+ or k+ 面
            l1=1 ; m1= 2                  ! l是第1个下标， m是第2个下标
          endif
         
          if(orient .eq. 1) then           ! 根据orient来 旋转连接方向 （见《理论手册》）
            l2=l1 ; m2=-m1
          else if (orient .eq. 2) then
            l2=m1 ; m2=l1
          else if (orient .eq. 3) then
            l2=-l1; m2=m1
          else 
            l2=-m1; m2=-l1
          endif
          
          if(mod(face2,2) .eq. 1) then            ! 对于 i-, j- or k- 面 
            tmp=l2; l2=m2; m2=tmp                  ! swap l2 and m2
          endif

! .inp文件的连接描述 为：正<-->正， 负<-->负； 
! Lp(k)==-1 交换次序 (ib,ie)--> (ie,ib)
! Ls(k)==-1 改变符号 (ib,ie) ---> (-ib, -ie)
          
		  if(face2 .eq. 1 .or. face2 .eq. 4) then  ! i- or i+
		    k0=1; k1=2 ; k2=3                  ! k0 单面； k1 -- l;  k2--m
		  else  if(face2 .eq. 2 .or. face2 .eq. 5) then  ! j- or j+
            k0=2; k1=1; k2=3
          else
		    k0=3; k1=1; k2=2
		  endif
		   
		   Lp(k0)=1 ; Ls(k0)=1            ! 无需交换次序，改变符号 
		   Lp(k1)=sign(1,l2)
           if(abs(l2) .eq. 1) then
		     Ls(k1)=1
		   else
             Ls(k1)=-1
		   endif

		   Lp(k2)=sign(1,m2)
           if(abs(m2) .eq. 1) then
		     Ls(k2)=1
		   else
             Ls(k2)=-1
		   endif
       end subroutine get_orient
    


