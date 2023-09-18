! 构建多块结构网格连接信息 bc_in文件
! Copyright by Li Xinliang, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Ver 1.0, 2010-11-21
! Ver 1.1, 2010-12-3  : 对于方形计算域，可自动识别固壁及外边界， （找出最外侧的6个平面，设定远场或对称边界条件，其他设为固壁）
! Ver 1.2, 2010-12-29: 修正了一个小Bug, 无连接点一律按内点处理
! Ver 1.3, 2010-1-3:   采用“快捷通道法”，提高了搜索算法的效率
! Ver 1.4.1, 2010-1-3:  采用树状多级搜索算法，通过切割面元(分片搜索)，提高了搜索算法的效率
! Ver 1.4.2, 2013-8-13: 将In2inp,f90 合并到本程序中，可生成.in与.inp格式文件
!------------------------------------------------------------------------------------------------------------------------------ 
  module Global_Variables
  implicit none
  integer,parameter:: Nq=7  ! 每一个面分割成若干个面元，每个面元Nq*Nq个点； 用于树状多级快速搜索时使用, 可以调整该数（推荐在5-16之间）以达到最高效率
  integer,parameter:: Max_Link=50   ! 一个点最多有50个连接点 
  integer,parameter:: Max_Subface=100  ! 一个面最多100个子面
  integer,parameter:: BC_Wall=-10, BC_Farfield=-20 ,BC_Symmetry=-50   ! 固壁、远场及对称面
  real*8 :: Dist_limit , Dist_min  ! 最小网格距离 
  Integer:: Num_Block,Form_mesh
!------Types 定义了点、面和块三种数据结构 -------------------------------------------------------------------------------------
! 点： 属性有 坐标，连接点的个数，连接点的信息， 是否为内点标志，是否被搜索过标志
! 面： 属性有 面号、维数，点
! 块： 属性有 块号、维数、面
!-------------------------------------------------------------------------------------------------------------------------------
   TYPE Point_Type    ! 点
   real*8:: x, y, z   ! 坐标
   integer:: num_link ! 连接点的个数
   integer:: search_flag, If_inner, color ! 搜索/染色时使用的标记量    
   integer,dimension(MAX_Link):: block_link, face_link, n1_link, n2_link  ! 连接点的信息（块号、面号、下标）
   End TYPE Point_Type
  
   TYPE Subface_TYPE  ! 子面
   integer:: subface_no, face_no, ib,jb,ie,je   ! 子面号(颜色)，子面所在的面号，起止位置
   integer:: subface_link , block_link , face_link , orient  ! 连接的子面 , 连接的块，连接的面 , 连接指向性（1-4）  
   END TYPE Subface_TYPE
   

   TYPE Face_Type                     ! 面 
   integer:: face_no,n1,n2   ! 面号、数组维
   integer:: nk1,nk2   !面元的数目
   Type (Point_Type),pointer,dimension(:,:) :: point 
   integer,pointer,dimension(:):: n1b,n1e,n2b,n2e
   real*8,pointer,dimension(:,:):: xc,yc,zc,dc  ! 中心坐标，到中心的最大距离
   End TYPE Face_Type
   
   TYPE Block_TYPE           !  块
     integer :: Block_no,nx,ny,nz   ! 块号、数组维
     integer:: color_now, Num_Subface  ! 正在使用的颜色， 子面数目
     TYPE(Face_type),dimension(6):: Face
     TYPE (Subface_TYPE),dimension(Max_subface):: subface   ! 子面 （从属于块）
  End TYPE Block_TYPE  
  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  
  end module Global_Variables

!---------------------------------------------------------------
   program get_bcin
   use Global_Variables
   implicit none
   integer:: zone, face,i,j,m,nf,flag1
   TYPE (Block_Type), Pointer:: B
   TYPE (Face_Type), Pointer:: F

   print*, "----------------------------------------------------------------------------"
   print*, "- This code is used to create the Grid Link file (BXCFD bc_in Format)      -"
   print*, "-               Ver 1.4.1 , 2010-1-3                                       -"  
   print*, "- Copyright by Li Xinliang, Institute of Mechanics, CAS, lixl@imech.ac.cn  -"
   print*, "----------------------------------------------------------------------------"
   write(*,*)
   write(*,*)

   call read_coordinate    ! 读取网格坐标，初始化

   call search_link        ! 搜索连接点，建立各点连接信息 （计算量大）

   call Mark_inner_point   !标记内点

   call creat_subface      ! 创建子面 （核心程序）

 !    print*, "Do you want to search the wall and Far field bounary ? (1 for yes, 0 for no)"
 !    read(*,*)  flag1
 !    if( flag1 .eq. 1) then
 !      call search_wall_and_farfield   ! 搜索固壁、远场及对称面
 !    endif

   call write_bcin
   
   call Inp2In (Form_Mesh)

   print*, "OK, The link message is written into the file bc.in and bc.inp"   
!-----------------------------------
!    print*, "Write the test file test.dat"
!    open(99,file="test.dat")
!    write(99,*) "variables=x,y,z,link,color,inner,blink,flink,i1,i2,block,face"
!    do m=1,Num_Block
!      B => Block(m)
!    do nf=1,6 
!      F=> B%face(nf)
!     write(99,*) "zone i = ", F%n1, " j= ", F%n2
!     do j=1, F%n2
!     do i=1, F%n1
!     write(99,"(3f16.8,9I4)")  F%point(i,j)%x, F%point(i,j)%y,F%point(i,j)%z,F%point(i,j)%num_link,F%point(i,j)%color, &
!                F%point(i,j)%If_inner,F%point(i,j)%block_link(1),F%point(i,j)%face_link(1),&
!                i,j,m,nf
!     enddo
!     enddo
!    enddo
!    enddo
!   close(99)

  end



!---------------------------------------------------------------------------------------------------------------------------
! 标记子面内点  
  subroutine Mark_inner_point
   use Global_Variables
   implicit none
   integer:: m,nf,i1,i2,flag1,flag2,flag3,flag4
   TYPE (Block_Type), Pointer:: B   ! 块指针
   TYPE (Face_Type), Pointer:: F    ! 面指针
   TYPE (Point_Type), Pointer:: P   ! 点指针

  do m=1,Num_Block
    B => Block(m)
  do nf=1,6 
    F=> B%face(nf)

!-----标示 "子面内点"  （ 周围点的连接点均在同一块、同一面）   
      do i2=2,F%n2-1
      do i1=2,F%n1-1
      P=> F%point(i1,i2)
! 修改 （2010-12-29）     
!      if(P%num_link .le. 1) then   ! 1个或0个连接  (多个连接的点必为子面边界点)

! 0个连接必为内点，1个连接需要判断，多个连接必为边界点

       if(P%num_link .le. 0) then   ! 0个连接必为内点
            P%If_inner=1
       else if  (P%num_link .eq. 1) then    
         call compare(m,nf,i1,i2, i1-1,i2, flag1) ! 判断两个点 的连接点是否在同一面
         call compare(m,nf,i1,i2, i1+1,i2, flag2) 
         call compare(m,nf,i1,i2, i1,i2-1, flag3) 
         call compare(m,nf,i1,i2, i1,i2+1, flag4) 
         if( flag1*flag2*flag3*flag4 .eq. 1 ) P%If_inner=1   ! 子面内点 
      endif     
      enddo
      enddo
   enddo
   enddo
  print*, "Mark Inner Point OK"
 end subroutine Mark_inner_point


!------------------------------------------------------------------------------------------------------------------
! 建立子面 （核心算法）  
  subroutine creat_subface
   use Global_Variables
   implicit none
   integer:: m,nf,i1,i2,color_now,ia,ja,ie,je,mt,nft,iat,jat,orient
   integer:: i1t,j1t,i2t,j2t
   TYPE (Block_Type), Pointer:: B,Bt
   TYPE (Face_Type), Pointer:: F,Ft            ! 变量中带t的是连接（面、子面、点）
   TYPE (SubFace_Type), Pointer:: SF,SFt
   TYPE (Point_Type), Pointer:: P,P1

  do m=1,Num_Block
    B => Block(m)
  do nf=1,6 
    F=> B%face(nf)

 !------------------------------------------------------------------
 
    do i2=2,F%n2-1
    do i1=2,F%n1-1
   
    P=> F%point(i1,i2)
    if(P%If_inner .eq. 1 .and. P%color .eq. 0) then  ! 如非内点，或已被染过色，则跳过

! 通过染色算法，找到子块的下边界 (ie,je)     (i1,i2)- (ie,je) 就是子面内点
       call find_subface_boundary(m,nf,i1,i2,ie,je)

!------创立子面，并染色--------------------------------------
       B%color_now=B%color_now+1   ! 设定当前颜色 （创立一个子面，新建一种颜色）
! 本子面各点染色
        do ja=i2,je
        do ia=i1,ie
           F%point(ia,ja)%color=B%Color_now         
        enddo
        enddo

!  登记本子面的属性 (编号、起止点坐标)
      SF=>B%subface(B%Color_now)  ! 指向该子面    
      SF%subface_no=B%Color_now   ! 编号 （颜色）
      SF%face_no=nf               ! 该子面所在的面号
      SF%ib=i1-1 ;  SF%jb=i2-1    ! 子面的左上角点坐标
      SF%ie=ie+1 ;  SF%je=je+1    ! 子面的右下角点坐标


! 如果有连接,创立连接的子面信息，并染色
   if(F%point(i1,i2)%Num_link .ne. 0) then   
          mt=F%point(i1,i2)%block_link(1)     ! 内点连接的块号 
          nft=F%point(i1,i2)%face_link(1)     ! 内点连接的面号
	
   ! 创立连接子面，登记连接子面的连接信息
          Bt=>Block(mt)   ! 连接的块
          Ft=>Bt%face(nft)                    ! 连接的面                    
          Bt%color_now=Bt%color_now+1         ! 设定连接子面的当前颜色  （在链接块中建立一个新的子面）
          SFt=>Bt%subface(Bt%color_now)       ! 连接的子面
          SFt%subface_no=Bt%color_now         ! 连接子面的编号
          SFt%face_no=nft                     ! 连接子面所在的面号
          SFt%block_link=m   ! 连向当前块
          SFt%face_link=nf   ! 连向当前面
          SFt%subface_link=B%Color_now  ! 连向当前子面

  !  登记本子面的连接信息		  
          SF%block_link=mt  ! 子面的连接块号 （=内点的连接块号）
          SF%face_link=nft  ! 子面的连接面号
          SF%subface_link=Bt%color_now        ! 子面的连接子面编号 （颜色,所连接面的当前色）

 ! 确定左上、右下角点的连接点         
          call find_corner_link_point(m,nf,SF%ib,SF%jb,i1t,j1t,mt,nft)
          call find_corner_link_point(m,nf,SF%ie,SF%je,i2t,j2t,mt,nft)

  !  登记连接子面的角点（范围）信息
         SFt%ib=min(i1t,i2t) ; SFt%ie=max(i1t,i2t)
         SFt%jb=min(j1t,j2t) ; SFt%je=max(j1t,j2t)

  ! 连接子面染色
        do ja=i2,je
        do ia=i1,ie
          iat=F%point(ia,ja)%n1_link(1)
          jat=F%point(ia,ja)%n2_link(1)
          Ft%point(iat,jat)%color=Bt%color_now   ! 连接的点也染色 （同时创立子面及连接子面，并染色）
        enddo
        enddo
 
 !  创建连接指向性
      if( mod(F%face_no,2) .eq. 1) then  !  l方向是第2个下标变化的方向 (l*m=n)
           call find_corner_link_point(m,nf,SF%ib,SF%je,i2t,j2t,mt,nft)   ! (ib,je)连接的点
      else                               !  l方向是第1个下标变化的方向
           call find_corner_link_point(m,nf,SF%ie,SF%jb,i2t,j2t,mt,nft)   ! (ie,jb)连接的点
      endif

      call comput_orient(Ft%face_no,i1t,j1t,i2t,j2t,orient )
         SF%orient=orient
         SFt%orient=orient
   else
! 无连接
         SF%block_link=0  
         SF%face_link=0  
         SF%subface_link=0        
         SF%orient=0
   endif
!-----------------------------------------------------------
  endif
   
  enddo
  enddo
 
   B%Num_subface=B%Color_now  ! 子面的数目

 !------------------------------------------
  enddo
    print*, "zone = ", m ,  "Subface Number=", B%Num_Subface
 
  enddo
 
 end subroutine creat_subface

!--------------------------------------------------------------------------------------------
!---通过染色法，找到子块的下边界点(ie,je), 即(i1,i2)-(ie,je)围成的矩形区域内的点，全部是内点
   subroutine find_subface_boundary(m,nf,i1,i2,ie,je)
   use Global_Variables
   implicit none
   integer:: m,nf,i1,i2,ie,je,ia,ja
   TYPE (Face_Type), Pointer:: F
   TYPE (Point_Type), Pointer:: P

    F=> Block(m)%face(nf)
 !  寻找该连续矩形块的终止位置 (ie,je)
 Loop1:  do ia=i1, F%n1
          P=> F%point(ia,i2)
          if(P%If_inner .eq. 0) then 
          ie=ia-1 
          exit Loop1
          endif
         enddo Loop1
         
 Loop2:  do ja=i2, F%n2
           do ia=i1,ie
            P=> F%point(ia,ja)
            if(P%If_inner .eq. 0) then 
            je=ja-1
            exit Loop2
            endif
          enddo
         enddo Loop2
    end subroutine find_subface_boundary

!----------------------------------------------------
! 找到角点(ib,ie)对应的连接点(i1t,j1t) , 由于该点可能不止一个连接，因此需要搜索与内点连接点在同一块、同一面的点
      subroutine find_corner_link_point(m,nf,ib,jb,i1t,j1t,mt,nft)
          use Global_Variables
          implicit none
          integer:: m,nf,ib,jb,i1t,j1t,mt,nft,ns
          TYPE (Point_Type), Pointer:: P
          P=>Block(m)%face(nf)%point(ib,jb)   ! 角点
  Loop3:  do ns=1,P%Num_link        ! 在全部连接点中搜索（连接块=mt, 连接面=nft的点）
            if(P%block_link(ns) .eq. mt .and. P%face_link(ns) .eq. nft) then
             i1t=P%n1_link(ns); j1t=P%n2_link(ns)    ! 连接点
            exit Loop3           
	    endif
	   enddo Loop3
      end subroutine find_corner_link_point


!---------------------------------------------------------------
!  比较两个点 m块,nf面上的两个点(i1,j1) 与 (i2,j2) 是否具有相同的连接信息
    subroutine compare(m,nf, i1,j1, i2,j2,flag1) 
!   subroutine compare(P,P1,flag1)
   use Global_Variables
   implicit none
   integer:: m,nf, i1,j1, i2,j2, flag1 , pb,pf, n
   TYPE (Face_Type), Pointer::F
   TYPE (Point_Type), Pointer:: P1,P2
    F=> Block(m)%face(nf)
    P1=> F%Point(i1,j1)
    P2 => F%point(i2,j2)
    if(P1%num_link .eq. 0)  then   ! P点无连接
       if(P2%num_link .eq. 0) then  
         flag1=1           !两个点均无连接，匹配
       else if (i2 .eq. 1 .or. j2 .eq. 1 .or. i2 .eq. F%n1 .or. j2 .eq. F%n2) then
         flag1=1           ! 面的边界点即使有连接也可以和无连接点匹配
       else   
        flag1=0            ! 有连接点与无连接点，不匹配
       endif
    else     
         pb=P1%block_link(1) ; pf=P1%face_link(1)    ! (第一个点) 连接的块号与面号
         flag1=0
Loop1:   do n=1,P2%num_link                           ! 搜索第2个点的全部连接信息
         if(P2%block_link(n) .eq. Pb .and. P2%face_link(n) .eq. Pf) then   
         flag1=1                      ! 块号、面号均匹配
         exit Loop1
         endif
         enddo Loop1
     endif

   end subroutine compare
!---------------------------------------------------------------

! 搜索坐标相同的点（连接点）
   subroutine search_link
   use Global_Variables
   implicit none
   integer:: m,nf,mt,nft,i1,i2,j1,j2,i1t,i2t,n,ns,mt1,nft1,mt_now,nft_now 
   integer:: zone,face,i,j,flag1,i0,j0,Link_max_point,ni1,nj1,ni2,nj2
   real*8:: m1,m2,mp1,mp2,tmp,d1

   real*8:: dist2
   TYPE (Block_Type), Pointer:: B,Bt
   TYPE (Face_Type), Pointer:: F,Ft
   TYPE (Point_Type), Pointer:: P,Pt


!-------------估算采用 树形多级搜索法 (分片搜索法) 的效率, 该值可作为修改参数Nq的参考-------------------
   print*, "To Estimate the efficient of Multi-stage search program ..."
   m1=0; m2=0;mp1=0;mp2=0
   do m=1,Num_Block
     B => Block(m)
     do nf=1,6 
     F=> B%face(nf)
      do j1=1,F%nk2
 !     print*, F%n2b(j1),F%n2e(j1)
      do i1=1,F%nk1
 !---------------- 被搜索点   
       do mt=1,Num_Block
        Bt=> Block(mt)
        do nft=1,6
         Ft=> Bt%face(nft)
         do j2=1,Ft%nk2
         do i2=1,Ft%nk1
         tmp=1.d0*(F%n2e(j1)-F%n2b(j1)+1)*(F%n1e(i1)-F%n1b(i1)+1)*(Ft%n2e(j2)-Ft%n2b(j2)+1)*(Ft%n1e(i2)-Ft%n1b(i2)+1)
         if(tmp .lt. 0) then
          print*, i1,j1,i2,j2
          print*, F%n2e(j1),F%n2b(j1), F%n1e(i1), F%n1b(i1), Ft%n2e(j2), Ft%n2b(j2), Ft%n1e(i2),F%n1b(i2)
         endif

        if(sqrt((F%xc(i1,j1)-Ft%xc(i2,j2))**2+(F%yc(i1,j1)-Ft%yc(i2,j2))**2+(F%zc(i1,j1)-Ft%zc(i2,j2))**2)  &
                 .gt. F%dc(i1,j1)+Ft%dc(i2,j2)+Dist_limit ) then  ! 两个面距离太远，跳过搜索
         m1=m1+1
         mp1=mp1+1.d0*(F%n2e(j1)-F%n2b(j1)+1)*(F%n1e(i1)-F%n1b(i1)+1)*(Ft%n2e(j2)-Ft%n2b(j2)+1)*(Ft%n1e(i2)-Ft%n1b(i2)+1)
        else
         m2=m2+1
         mp2=mp2+1.d0*(F%n2e(j1)-F%n2b(j1)+1)*(F%n1e(i1)-F%n1b(i1)+1)*(Ft%n2e(j2)-Ft%n2b(j2)+1)*(Ft%n1e(i2)-Ft%n1b(i2)+1)
       endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
 
      print*, "The speed-up by using multi-stage search technique is:",    mp1/mp2
      print*, "If the speed-up is low, you can modify the parameter 'Nq' to increase it" 


!---Search the face point----------------------------------------------------------

  Link_max_point=0   ! 最多连接点的连接数
   print*, "Search the Edge Points ......"
 !------搜索点--------             
   do m=1,Num_Block
    print*, "----------------------------- :  Zone ", m 
     B => Block(m)
     do nf=1,6 
     print*, "face ", nf
     F=> B%face(nf)
     do nj1=1,F%nk2
     do ni1=1,F%nk1
       do i2= F%n2b(nj1),F%n2e(nj1)
       do i1= F%n1b(ni1),F%n1e(ni1)
        P=> F%point(i1,i2)
        P%Search_flag=1     ! 表示已被搜索过
!--------------------------------------------------------------    
!       被搜索点       
       do mt=1,Num_Block
        Bt=> Block(mt)
        do nft=1,6
         Ft=> Bt%face(nft)
         do nj2=1,Ft%nk2
         do ni2=1,Ft%nk1
          d1=sqrt((F%xc(ni1,nj1)-Ft%xc(ni2,nj2))**2+(F%yc(ni1,nj1)-Ft%yc(ni2,nj2))**2+(F%zc(ni1,nj1)-Ft%zc(ni2,nj2))**2)
          if(d1 .gt. F%dc(ni1,nj1)+Ft%dc(ni2,nj2)+Dist_limit )  cycle  ! 两个面距离太远，跳过搜索
 !---------------------------------------
          do i2t= Ft%n2b(nj2),Ft%n2e(nj2)
          do i1t= Ft%n1b(ni2),Ft%n1e(ni2)
            Pt=>Ft%point(i1t,i2t)
           if(Pt%Search_flag .eq. 0 ) then       ! 尚未搜索过的
            dist2=(P%x-Pt%x)**2+(P%y-Pt%y)**2+(P%z-Pt%z)**2

           if(dist2 .lt. Dist_limit**2  ) then       ! 找到坐标相同的点        
!        一对对应点         
           P%num_link=P%num_link+1
           P%block_link(P%num_link)=mt
           P%face_link(P%num_link)=nft
           P%n1_link(P%num_link)=i1t
           P%n2_link(P%num_link)=i2t
           Pt%num_link=Pt%num_link+1
           Pt%block_link(Pt%num_link)=m
           Pt%face_link(Pt%num_link)=nf
           Pt%n1_link(Pt%num_link)=i1
           Pt%n2_link(Pt%num_link)=i2

           if(P%num_link .gt. Link_max_point) Link_max_point=P%num_link
           if(Pt%num_link .gt. Link_max_point) Link_max_point=Pt%num_link
           endif
          endif
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo 
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo

  print*, "Search Edge OK ..."
  print*, "Max Link is ", Link_max_point



  print*, "Search OK ..."

!-----------------------------------
!    print*, "Write the file link.dat"
!    open(99,file="link.dat")

!    do m=1,Num_Block
!      B => Block(m)
!    do nf=1,6 
!      F=> B%face(nf)
!     write(99,*) "Zone= ", m, " Face= ", nf , "-----------------------------------------"
!     do j=1, F%n2
!     do i=1, F%n1
!      P=> F%point(i,j)
!     write(99,*) i,j,P%Num_link
!     do n=1,P%Num_link
!     write(99,*) P%block_link(n),P%face_link(n),P%n1_link(n),P%n2_link(n)
!     enddo
!     enddo
!     enddo
!    enddo
!    enddo
!   close(99)


 
  end subroutine search_link
        
  
!-----------------------------------------------------------------
! 读取网格信息，存储表面网格
! 搜索最小网格间距Dist_min (距离小于0.1倍Dis_min的两个点认为是同一个点） 
  
  subroutine read_coordinate
   use Global_Variables
   implicit none
   integer,allocatable,dimension(:):: NI,NJ,NK
   real*8,allocatable,dimension(:,:,:):: x,y,z
   integer:: nx,ny,nz,i,j,k,m,nf,ib(6),jb(6),kb(6),ie(6),je(6),ke(6)
   integer:: i1,i2,nf1,ni1,nj1,n1,nk1,nk2
    TYPE (Block_Type), Pointer:: B
   TYPE (Face_Type), Pointer:: F
   TYPE (SubFace_Type), Pointer:: SF
   TYPE (Point_Type), Pointer:: P

   real*8:: di,dj,dk,d1,tmp
!-----------------------------------------------------------------------  
   Dist_min=1000.d0  ! 初值
  
 !  print*, "==================================================================="
 !  print*, "Please input the distance threshold Dist_limit"
 !  print*, "two points with distance < Dist_limit will be considered as the same point"
 !  print*, "if Dist_limit<=0, this code will automatically set it"
 !  print*, "????? Input Dist_limit ?????"
 !  read(*,*) Dist_limit
 
 
   print*, "Is Mesh3d.dat formatted file ?  1 for formatted, 0 for unformatted"
   read(*,*) Form_Mesh   
   if(Form_Mesh .eq. 1) then
    open(99,file="Mesh3d.dat")
    read(99,*) Num_Block
   else
    open(99,file="Mesh3d.dat",form="unformatted")
    read(99) Num_Block         ! 总块数
   endif

  
    allocate(Block(Num_Block))               
    allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! 每块的大小
   if(Form_Mesh .eq. 1) then
    read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
   else
    read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
   endif
! 读取每块几何信息, 记录面上的信息----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     B%Block_no=m
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   ! nx,ny,nz 每块的大小
     nx=B%nx ; ny= B%ny ; nz=B%nz
! ----------  几何量 -----------------------------------------------
    allocate(x(nx,ny,nz), y(nx,ny,nz), z(nx,ny,nz))  ! 格点坐标
   if(Form_Mesh .eq. 1) then
    read(99,*) (((x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                (((y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                (((z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else  
    read(99) (((x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                (((y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                (((z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif

!------------
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
      di=sqrt( (x(i+1,j,k)-x(i,j,k))**2+ (y(i+1,j,k)-y(i,j,k))**2 + (z(i+1,j,k)-z(i,j,k))**2 )
      dj=sqrt( (x(i,j+1,k)-x(i,j,k))**2+ (y(i,j+1,k)-y(i,j,k))**2 + (z(i,j+1,k)-z(i,j,k))**2 )
      dk=sqrt( (x(i,j,k+1)-x(i,j,k))**2+ (y(i,j,k+1)-y(i,j,k))**2 + (z(i,j,k+1)-z(i,j,k))**2 )
     Dist_min=min(Dist_min,di,dj,dk)
    enddo
    enddo
    enddo
!-------------

!  记录6个面的信息 
! face 1  (i-)     
     F=> B%face(1)
     F%face_no=1 ; F%n1= ny ; F%n2= nz           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(1,i1,i2);   F%point(i1,i2)%y=y(1,i1,i2);  F%point(i1,i2)%z=z(1,i1,i2)
     enddo
     enddo

! face 2  (j-)     
     F=> B%face(2)
     F%face_no=2 ; F%n1= nx ; F%n2= nz           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(i1,1,i2);   F%point(i1,i2)%y=y(i1,1,i2);  F%point(i1,i2)%z=z(i1,1,i2)
     enddo
     enddo

! face 3  (k-)     
     F=> B%face(3)
     F%face_no=3 ; F%n1= nx ; F%n2= ny           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(i1,i2,1);   F%point(i1,i2)%y=y(i1,i2,1);  F%point(i1,i2)%z=z(i1,i2,1)
     enddo
     enddo

! face 4  (i+)     
     F=> B%face(4)
     F%face_no=4 ; F%n1= ny ; F%n2= nz           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(nx,i1,i2);   F%point(i1,i2)%y=y(nx,i1,i2);  F%point(i1,i2)%z=z(nx,i1,i2)
     enddo
     enddo

! face 5  (j+)     
     F=> B%face(5)
     F%face_no=5 ; F%n1= nx ; F%n2= nz           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(i1,ny,i2);   F%point(i1,i2)%y=y(i1,ny,i2);  F%point(i1,i2)%z=z(i1,ny,i2)
     enddo
     enddo

! face 6  (k+)     
     F=> B%face(6)
     F%face_no=6 ; F%n1= nx ; F%n2= ny           
     allocate(F%point(F%n1,F%n2))
     do i2=1,F%n2
     do i1=1,F%n1
       F%point(i1,i2)%x=x(i1,i2,nz);   F%point(i1,i2)%y=y(i1,i2,nz);  F%point(i1,i2)%z=z(i1,i2,nz)
     enddo
     enddo
!----------------------------------

    deallocate(x,y,z) 

!  初始化面信息    
   do nf=1,6
    F=> B%face(nf)
!   初始化点信息   
    do i2=1,F%n2
    do i1=1,F%n1
    F%point(i1,i2)%Search_flag=0    ! 搜索时使用的临时变量 
    F%Point(i1,i2)%color=0
    F%point(i1,i2)%If_inner=0       ! "子面内点" 标志
    F%point(i1,i2)%num_link=0       ! 连接点总数
    F%point(i1,i2)% block_link(:)=0    ! 
    F%point(i1,i2)% face_link(:)=0
    F%point(i1,i2)% n1_link(:)=0
    F%point(i1,i2)% n2_link(:)=0
    enddo
    enddo
   enddo

    ! 初始化子面信息 
     B%Num_Subface=0  
     B%color_now=0

     do k=1, Max_Subface
     SF=> B%subface(k)
     SF%subface_no=0 ; SF%face_no=0; SF%ib=0 ; SF%jb=0 ; SF%ie=0 ; SF%je=0
     SF%subface_link=0; SF%block_link=0 ; SF%face_link=0 ; SF%orient=0
     enddo
  !--------------------------------------------------------------
  ! 创建面元， 记录每个面元的中心点及半径
     do nf=1,6 
      F=> B%face(nf)
      nk1=int((F%n1-1)/Nq)+1    !面元的数目 （nk1*nk2块）
      nk2=int((F%n2-1)/Nq)+1
      F%nk1=nk1  
      F%nk2=nk2
      allocate(F%n1b(nk1),F%n1e(nk1),F%n2b(nk2),F%n2e(nk2))
      allocate(F%xc(nk1,nk2),F%yc(nk1,nk2),F%zc(nk1,nk2),F%dc(nk1,nk2))

! 把从下标从1到F%n1的点分割成若干段，每段Nq个点; 如果不能整除，最后一段数目少些
!  F%n1b(k) 和F%n1e(k) 是第k段的起、止下标    
      do k=1,nk1
        F%n1b(k)=(k-1)*Nq+1        
        if(k .ne. nk1) then
          F%n1e(k)=k*Nq
        else
          F%n1e(k)=F%n1
        endif   
      enddo

! 把从下标从1到F%n2的点分割成若干段，每段Nq个点; 如果不能整除，最后一段数目少些
      do k=1,nk2
        F%n2b(k)=(k-1)*Nq+1        
        if(k .ne. nk2) then
          F%n2e(k)=k*Nq
        else
          F%n2e(k)=F%n2
        endif   
      enddo


      do nj1=1,nk2
      do ni1=1,nk1

      F%xc(ni1,nj1)=0.d0; F%yc(ni1,nj1)=0.d0 ; F%dc(ni1,nj1)=0.d0
      do j=F%n2b(nj1), F%n2e(nj1)
      do i=F%n1b(ni1), F%n1e(ni1)
       P=> F%point(i,j)
       F%xc(ni1,nj1)=F%xc(ni1,nj1)+P%x
       F%yc(ni1,nj1)=F%yc(ni1,nj1)+P%y
       F%zc(ni1,nj1)=F%zc(ni1,nj1)+P%z
      enddo
      enddo
      tmp=(F%n2e(nj1)-F%n2b(nj1)+1)*(F%n1e(ni1)-F%n1b(ni1)+1)
      F%xc(ni1,nj1)=F%xc(ni1,nj1)/tmp
      F%yc(ni1,nj1)=F%yc(ni1,nj1)/tmp
      F%zc(ni1,nj1)=F%zc(ni1,nj1)/tmp

      do j=F%n2b(nj1), F%n2e(nj1)
      do i=F%n1b(ni1), F%n1e(ni1)
       P=> F%point(i,j)
       d1=sqrt((P%x-F%xc(ni1,nj1))**2+(P%y-F%yc(ni1,nj1))**2 +(P%z-F%zc(ni1,nj1))**2 )
       if(d1 .gt. F%dc(ni1,nj1)) F%dc(ni1,nj1)=d1
      enddo
      enddo
 !     print*, F%xc,F%yc,F%zc,F%dc
    enddo
    enddo

   enddo      ! face

!---------------------------------------------------------------------------------------   
   enddo      ! block
   close(99)
   print*, "read Mesh3d OK ..."
   print*, "Minima Mesh Space is ", Dist_min
    
 !   if(Dist_limit .le. 0) Dist_limit=0.5d0*Dist_min  ! 设置间距门槛 （0.5倍最小网格间距），距离小于该门槛的两个点将认为是同一个点
    Dist_limit=0.5d0*Dist_min  ! 设置间距门槛 （0.5倍最小网格间距），距离小于该门槛的两个点将认为是同一个点
    print*, "The distance threshold is ", Dist_limit
   end subroutine read_coordinate

!---------------------------------------------------------
! 计算连接指向性
!  子面的l方向连接到目标子面从(i1,j1)指向(i2,j2)的一条有向线段
    subroutine comput_orient(fno,i1,j1,i2,j2,orient )
    implicit none
    integer:: fno,i1,j1,i2,j2,orient,p1,p2
! 判断 从(i1,j1) 指向 (i2,j2) 矢量的方向 （1,-1, 2, -2 表示 i, -i, j, -j 方向）     
    if(i1 .eq. i2 ) then
      if(j2 .gt. j1) then
        p1=2            ! j方向
      else
        p1=-2           ! -j方向
      endif
    else
      if(i2 .gt. i1) then
        p1=1
      else
        p1=-1
      endif
    endif
! 根据（连接面）的面号（指向），将P方向(P1以 i, -i, j, -j度量） 转换成以 l, -l, m, -m 度量 的数 P2；
   if( mod(fno,2) .eq. 0) then
     p2=p1
   else
     if(p1 .eq. 1)  p2=2
     if(p1 .eq. -1) p2=-2
     if(p1 .eq. 2)  p2=1
     if(p1 .eq. -2) p2=-1
   endif
!  根据P2，设定指向性orient    
    if(p2 .eq. 1) orient=1
    if(p2 .eq. 2) orient=2
    if(p2 .eq. -1) orient=3
    if(p2 .eq. -2) orient=4
  
   end subroutine comput_orient

! -------写bc_in文件----------------------------------------  
  subroutine write_bcin
   use Global_Variables
   implicit none
   integer:: m,ns,ist,jst,kst,iend,jend,kend,face,block_link
   TYPE (Block_Type), Pointer:: B   ! 块指针
   TYPE (Face_Type), Pointer:: F    ! 面指针
   TYPE (Subface_Type), Pointer:: SF   ! 子块指针
   open(109,file="bc.in")
   write(109,*) " BC file , By Li Xinliang"
   write(109,*) "# Blocks"
   write(109,*) Num_Block
   do m=1,Num_Block
   B => Block(m)
    write(109,*) "Block  ", m
    write(109,*) "Subfaces"
    write(109,*) B%Num_subface, 1, 1, -1
     write(109,*) " f_no, face, istart,iend, jstart, jend, kstart, kend, neighb, subface ori theta"
     do ns=1,B%Num_subface
     SF=>B%subface(ns)
     face=SF%face_no
     call convert_ijk(ist,jst,kst,SF%ib,SF%jb,face,B%nx,B%ny,B%nz) 
     call convert_ijk(iend,jend,kend,SF%ie,SF%je,face,B%nx,B%ny,B%nz) 
     if(SF%block_link .eq. 0) then
       Block_link=-1
     else
       Block_link=SF%block_link
     endif

     write(109,"(11I7,E16.5)") ns,face,ist,iend,jst,jend,kst,kend, block_link, SF%subface_link, SF%orient , 0.d0
     enddo
    enddo
   close(109)
   end subroutine write_bcin


!-----------------------------------------------------------------
! 将面上的局部坐标(i1,i2)转化为块上的全局坐标(i,j,k)  
  subroutine convert_ijk(i,j,k,i1,i2,face_no,nx,ny,nz)
  implicit none
  integer:: i,j,k,i1,i2,face_no,nx,ny,nz
   if(face_no .eq. 1) then
     i=1; j=i1; k=i2
   else if (face_no .eq. 2) then
     i=i1; j=1; k=i2
   else if (face_no .eq. 3) then
     i=i1; j=i2; k=1
   else if (face_no .eq. 4) then
     i=nx; j=i1; k=i2
   else if (face_no .eq. 5) then
     i=i1; j=ny ; k=i2
   else if (face_no .eq. 6) then
     i=i1; j=i2; k=nz
   endif
   end subroutine convert_ijk

!-----------------------------------------------------
! 搜索固壁、远场及对称面
   subroutine search_wall_and_farfield
   use Global_Variables
   implicit none
   integer:: i,j,k,m,ns,nf
   real*8:: xc,yc,zc,tmp,xrms,yrms,zrms
   TYPE (Block_Type), Pointer:: B   ! 块指针
   TYPE (Face_Type), Pointer:: F    ! 面指针
   TYPE (Subface_Type), Pointer:: SF   ! 子块指针
 !---------------------
   open(103,file="surface.dat")
   write(103,*) "variables=x,y,z,bc" 
   do m=1,Num_Block
     B => Block(m)
   do ns=1,B%Num_subface
     SF=>B%subface(ns)
     nf=SF%face_no
     F=>B%face(nf)
    if(SF%block_link .eq. 0 ) then
!  计算各子面的平均(中心)坐标
     xc=0.d0; yc=0.d0; zc=0.d0
     xrms=0.d0;yrms=0.d0;zrms=0.d0
     do j=SF%jb,SF%je
     do i=SF%ib,SF%ie
       xc=xc+F%point(i,j)%x
       yc=yc+F%point(i,j)%y
       zc=zc+F%point(i,j)%z
      enddo
      enddo
       tmp=1.d0*(SF%ie-SF%ib+1)*(SF%je-SF%jb+1)
       xc=xc/tmp; yc=yc/tmp ; zc=zc/tmp

!  计算子面各点与中心点的位置差  
       do j=SF%jb,SF%je
       do i=SF%ib,SF%ie
        xrms=xrms+(F%point(i,j)%x-xc)**2
        yrms=yrms+(F%point(i,j)%y-yc)**2
        zrms=zrms+(F%point(i,j)%z-zc)**2
       enddo
       enddo
        xrms=sqrt(xrms/tmp) ; yrms=sqrt(yrms/tmp) ; zrms=sqrt(zrms/tmp)

!       print*, "------------------"
!       print*, "Block =", B%Block_no, "subface=",SF%subface_no, "x,y,z=", xc, yc,zc
!       print*, "xrms,yrms,zrms=",xrms,yrms,zrms
!        We assume "Y=0" is the symmetry plane !!!!  if Y=0 is not the symmetry plane, please modify the code !
!    如果是平面，则应当是远场或对称面 （只适用于方形计算域）      
      if(xrms .lt. Dist_limit .or. yrms .lt. Dist_limit .or. zrms .lt. Dist_limit ) then
         if(abs(yc) .lt. Dist_limit) then  
          SF%block_link=BC_Symmetry        ! 对称面
         else 
           SF%block_link=BC_Farfield       ! 远场
         endif
      else 
        SF%block_link=BC_Wall
      endif
!------write to surface file-------------------------
      write(103,*) "zone i=", SF%ie-SF%ib+1, " j= ", SF%je-SF%jb+1
      do j=SF%jb,SF%je
      do i=SF%ib,SF%ie
       write(103, "(4f20.10)") F%point(i,j)%x, F%point(i,j)%y,F%point(i,j)%z, SF%block_link*1.d0
      enddo
      enddo
     endif
     enddo
     enddo
 
   close(103)
   end  


!==========================================================================
! In2inp, Transform BXCFD .in file to Gridgen .inp file  
! （网格连接信息）从BXCFD的 .in格式 转化为Gridgen .inp格式
! Copyright by Li Xinliang, lixl@imech.ac.cn
! Ver 1.0, 2012-7-11

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
  subroutine Inp2In (Form_Mesh)
  use Def_block
  implicit none
  integer:: Form_Mesh
  print*, " Transform  from .in file to .inp file ......"
!----------------------------
  call read_bcin(Form_Mesh)
  call trans_in_inp
  call write_inp
  end

  
!--------------------------------------------------------------------

  
 !----Mesh control message (bc2d.in)------------------------------------------
  subroutine read_bcin (Ia) 
   use Def_block
   implicit none
   integer::m,ksub,Ia,NB1
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc

   print*, "read bc.in ......"
   open(88,file="bc.in")
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
   
 !  print*, "read Mesh3d.dat for nx,ny,nz"
 !  print*, "please input the format of Mesh3d.dat, 1 formatted, 2 unformatted, 0 not read"
 !  read(*,*) Ia
   
   if(Ia .eq. 1 ) then
     open(100,file="Mesh3d.dat")
     read(100,*) NB1
     if(NB1 .ne. NB) then
	  print*, "error ! NB in Mesh3d.dat is not the same as that in bc3d.in "
	  stop
	 endif 
	  read(100,*) ((Block(m)%nx,Block(m)%ny,Block(m)%nz),m=1,NB)
   
   else if (Ia .eq. 0 ) then
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
   integer,parameter:: BC_Wall_in=-10, BC_Farfield_in=-20, BC_Periodic_in=-30,BC_Symmetry_in=-40,BC_Outlet_in=-22   ! .in 关于边界条件的定义
   integer,parameter:: BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Outlet=401, BC_Periodic=501     ! 与Griggen .inp文件的定义可能有所区别，请注意
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
		    print*, "The boundary condition is not supported!", NB, ksub
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
   open(99,file="bc.inp")
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
    



