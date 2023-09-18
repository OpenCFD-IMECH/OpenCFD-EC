!----Boundary message (bc3d.inp, Gridgen general format)----------------------------------------------------------
! 读取边界格式文件(bc3d.inp)，并转化为OpenCFD-EC的内建.inc格式 (见《OpenCFD-EC理论手册》 2.4.3节)
! OpenCFD-EC内建的存储格式比bc3d.inp多了一些冗余信息 （有些类似BXCFD的.in格式），例如多了子面号f_no,
! 面类型face 以及连接的子面号f_no1,连接的面类型face1
! 以及连接次序L1, L2, L3  (例如L1=1表示该维与连接块的第1维正连接， L1=-1表示与连接块的第1为反向连接).
! 这些冗余信息为块-块之间的通信（尤其是MPI并行通信）提供了便利，有利于简化代码
 
 module  Type_def1
    TYPE BC_MSG_TYPE              ! 边界链接信息
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                     ! 子面号，连接顺序描述符
   END TYPE BC_MSG_TYPE


    TYPE Block_TYPE1                          ! 数据结构：仅包含Bc_msg 
	 integer::  nx,ny,nz                      ! 网格数nx,ny,nz
	 integer::  subface                       ! 子面数
 	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
    END TYPE Block_TYPE1   
 End module  Type_def1

!------------------------------------------------------
  subroutine convert_inp_inc 
   use Type_def1
   implicit none
  
   integer:: NB,m,ksub,nx,ny,nz,k,j,k1,ksub1
   integer:: kb(3),ke(3),kb1(3),ke1(3),s(3),p(3),Lp(3)
   TYPE(Block_TYPE1),Pointer,dimension(:):: Block
   Type (Block_TYPE1),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
     
   Interface   
     subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
	   integer,dimension(3):: kb,ke,kb1,ke1
	 end subroutine
    end Interface
   

   print*, "Convert bc3d.inp to bc3d.inc ..."
   open(88,file="bc3d.inp")
   read(88,*)
   read(88,*) NB
   allocate(Block(NB))
   do m=1,NB
    B => Block(m)
    read(88,*) B%nx,B%ny,B%nz
	read(88,*)
    read(88,*) B%subface   !number of the subface in the Block m
     
	allocate(B%bc_msg(B%subface))   ! 边界描述

    do ksub=1, B%subface
      Bc => B%bc_msg(ksub)
      Bc%f_no=ksub                        ! 子面号
	  read(88,*)  kb(1),ke(1),kb(2),ke(2),kb(3),ke(3),Bc%bc
	 
	  if(Bc%bc .lt. 0) then
 !  --------有连接的情况 (内边界)--------------------------------------------------------	    
	   read(88,*) kb1(1),ke1(1),kb1(2),ke1(2),kb1(3),ke1(3),Bc%nb1
      else
 !---------无连接情况 (物理边界)----------------------------------------
            kb1(:)=0; ke1(:)=0; Bc%nb1=0
	  endif
	 call Convert_bc(Bc,kb,ke,kb1,ke1)
  enddo
  enddo
   
   close(88)

!  搜索连接块的块号 f_no1  (便于MPI并行通信是使用)
   do m=1,NB
     B => Block(m)
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .lt. 0) then
         Bc%f_no1=0
		 B1=>Block(Bc%nb1)         ! 指向连接块
         
		 do ksub1=1,B1%subface
		 Bc1=>B1%bc_msg(ksub1)
		 if(Bc%ib1==Bc1%ib .and. Bc%ie1==Bc1%ie .and. Bc%jb1==Bc1%jb .and. Bc%je1==Bc1%je   &
		    .and. Bc%kb1==Bc1%kb .and. Bc%ke1==Bc1%ke  .and. Bc1%nb1==m) then
		    Bc%f_no1=ksub1
		  exit
		 endif 	 
         enddo

         if(Bc%f_no1 ==0) then
		  print*, " Error in find linked block number !!!"
		  print*, "Block, subface=",m, ksub
		  stop
		 endif
	   endif
      enddo
	enddo

   open(99,file="bc3d.inc")
    write(99,*) " Inp-liked file of OpenCFD-EC"
    write(99,*) NB
	do m=1,NB
     B => Block(m)
     write(99,*) B%nx, B%ny, B%nz
	 write(99,*) "Block ", m
	 write(99,*) B%subface
	  do ksub=1,B%subface
        Bc => B%bc_msg(ksub)
	   write(99,"(9I6)") Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
	   write(99,"(12I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
      enddo
	enddo
	close(99)

   print*, "Convert bc3d.inp to bc3d.inc OK"


  end  

!-----------------------------------------------------------------------------------




!   将Gridgen格式 转换为 OpenCFD-EC 的边界连接格式
!     计算面号、 连接次序 (L1,L2,L3)等
      subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
	   integer,dimension(3):: kb,ke,kb1,ke1,s,p,LP
       integer:: k,j,k1

	   if(Bc%bc .ge. 0) then
	     Bc%ib1=0; Bc%ie1=0; Bc%jb1=0; Bc%je1=0; Bc%kb1=0; Bc%ke1=0; Bc%nb1=0
         Bc%L1=0; Bc%L2=0; Bc%L3=0; Bc%face1=0; Bc%f_no1=0
       endif
     

!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
       do k=1,3
  	     if(kb(k) .eq. ke(k) ) then 
	       s(k)=0                           ! 连接维
	     else if (kb(k) .gt. 0) then 
	       s(k)=1                           ! 正
	     else
	       s(k)=-1                          ! 负
	     endif
       enddo

!    边界子面的大小     
	 Bc%ib=min(abs(kb(1)),abs(ke(1))) ;  Bc%ie=max(abs(kb(1)),abs(ke(1)))
     Bc%jb=min(abs(kb(2)),abs(ke(2))) ;  Bc%je=max(abs(kb(2)),abs(ke(2)))
     Bc%kb=min(abs(kb(3)),abs(ke(3))) ;  Bc%ke=max(abs(kb(3)),abs(ke(3))) 


!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
      if(s(1) .eq. 0) then
	     if (Bc%ib .eq. 1) then
	      Bc%face=1               ! i-
	     else
	      Bc%face=4               ! i+
	    endif
      else if(s(2) .eq. 0) then
	    if(Bc%jb .eq. 1) then
	      Bc%face=2                 ! j-
	    else
	      Bc%face=5                 ! j+
	    endif
      else
 	   if(Bc%kb  .eq. 1) then
	    Bc%face=3                 ! k-
	   else
	    Bc%face=6                 ! k+
	   endif
     endif 

!---------------------------------------------------------------------------
!------内边界的情况，建立连接描述
  if( Bc%bc .lt. 0) then            ! 内边界
!      计算连接顺序描述符L1,L2,L3
!      计算各维之间的连接关系      
     do k=1,3  
	   if(kb1(k) .eq. ke1(k) ) then 
	       p(k)=0                      ! 
       else if (kb1(k) .gt. 0) then
	       p(k)=1                      ! .inp 文件的 正数
       else
	       p(k)=-1
       endif
     enddo
 	   

!    对应连接子面的大小     
	 Bc%ib1=min(abs(kb1(1)),abs(ke1(1))) ;  Bc%ie1=max(abs(kb1(1)),abs(ke1(1)))
     Bc%jb1=min(abs(kb1(2)),abs(ke1(2))) ;  Bc%je1=max(abs(kb1(2)),abs(ke1(2)))
     Bc%kb1=min(abs(kb1(3)),abs(ke1(3))) ;  Bc%ke1=max(abs(kb1(3)),abs(ke1(3))) 
    
 	  
!   判断该面连接面的类型 (i-, i+, j-,j+, k-,k+)     
      if(p(1) .eq. 0) then
	     if (Bc%ib1 .eq. 1) then
	      Bc%face1=1               ! i-
	     else
	      Bc%face1=4               ! i+
	    endif
      else if(p(2) .eq. 0) then
	    if(Bc%jb1 .eq. 1) then
	      Bc%face1=2                 ! j-
	    else
	      Bc%face1=5                 ! j+
	    endif
      else
 	   if(Bc%kb1 .eq. 1) then
	    Bc%face1=3
	   else
	    Bc%face1=6
	   endif
     endif  

!  计算“连接对” 描述符  bc%L1, bc%L2, bc%L3 
 	   do k=1,3
	     do j=1,3
	       if(s(k) .eq. p(j)) Lp(k)=j          ! .inp文件的连接格式： 正对正、 负对负、 0对0； 
	     enddo
	   enddo    
	   
!    计算连接次序 （正为顺序；负为拟序）	  
	  do k=1,3
	   if(s(k) .ne. 0) then
	     k1=Lp(k)
	     if( (ke(k)-kb(k))*(ke1(k1)-kb1(k1)) .lt. 0) Lp(k)=-Lp(k)      ! 逆序连接
	   else
         k1=Lp(k)
!		 if( (mod(Bc%face,2)-mod(Bc%face1,2))==0) Lp(k)=-Lp(k)         ! 逆序连接 （单面） ! Bug 2012-7-13
! 正-正连接 Lp为负 (例， i+ 面连接到 j+面， 则为逆序连接)
		 if( (Bc%face-1)/3 .eq. (Bc%face1-1)/3 ) Lp(k)=-Lp(k)        ! (Bc%face=1,2,3为 +面，4,5,6为-面)  ! 逆序连接 （单面）

	   endif
	 enddo 
     
	 Bc%L1=Lp(1); Bc%L2=Lp(2); Bc%L3=Lp(3)   ! 连接次序描述符 （详见《理论手册》）
   endif
  end