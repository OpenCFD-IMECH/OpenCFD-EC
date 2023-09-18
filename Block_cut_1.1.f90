!----Boundary message (bc3d.inp, Gridgen general format)----------------------------------------------------------
! Cut large Block into pieces
! Ver 1.1, 2013-9-6:  Cut Mesh and bc3d.inp
 
  module const_var1
  implicit none
 
!    integer,parameter:: PRE_SC=4            ! Single precision
     integer,parameter:: PRE_SC=8           ! Double Precision
     integer,parameter::  BC_Wall=2, BC_Symmetry=3,   BC_Inflow=5, BC_Outflow=6 , BC_NonReflection=9, BC_Dilchlet=1 
     integer,parameter::  BC_WallH=20        ! 壁面边界 (使用Ghost Cell)
     integer,parameter::  BC_inZ=-1, BC_In=-2,  BC_Periodic=-3       ! -1硬分区连接； -2软分区连接 ; -3周期边界 （软连接）
     integer,parameter:: Maxsub=100 ,NPMAX=100  ! 每块最大的子面数
  end module

  module  Type_def1
    use const_var1
	implicit none

    TYPE BC_MSG_TYPE              ! 边界链接信息
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                                            ! 子面号，连接顺序描述符
     integer:: b1,e1,b2,e2,b1t,e1t,b2t,e2t           ! 面内两维的描述
   END TYPE BC_MSG_TYPE


    TYPE Block_TYPE1                          ! 数据结构：仅包含Bc_msg 
	 integer::  nx,ny,nz                      ! 网格数nx,ny,nz
	 integer::  subface                       ! 子面数
 	 integer:: Block_No
     real(PRE_SC),pointer,dimension(:,:,:,:):: xyz   ! x, y, z
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
     integer:: Pi,Pj,Pk                                 ! (原块)三个方向的分割数 ; (新块)三个方向的块序号
     integer:: ib,jb,kb                                 ! 子块的起始下标
	 integer:: nb_ori       ! 原块号
	 integer:: nbk0      ! 子块的起始块号
	END TYPE Block_TYPE1   
    
	TYPE bcnp_type              ! 子块
     integer:: idx,blk,idxt,blkt,nn  ! 首地址，子块号 （源、连接）； 数目     
    END TYPE 
   
   End module  Type_def1


    module global_var1
    use Type_def1
    TYPE(Block_TYPE1),Pointer,dimension(:):: Block,Block_new
    integer:: Mesh_form,NB,NBnew
    integer,allocatable,dimension(:):: PI,PJ,PK
    type(bcnp_type),target:: bcn(NPMAX,2)
	end module global_var1

!---------------------------------


   program cut_block
    use global_var1
    implicit none
!----Cut Mesh --------------------
    call read_cutfile
	call cut_mesh
    call write_new_mesh
!---cut bc3d.inp -----------------
     call convert_inp_inc  
     call cut_inp   ! 构造子块的 inp信息
     call set_new_inner  ! 设定新的内联边界
	 call search_fno1(2)  ! 建立 fno1信息
	 call write_inp_new  ! 写入新的边界信息	 


   end


!----------------------------------------------------------
   subroutine cut_inp 
     use global_var1
     implicit none
	 integer:: m,ks
     Type (Block_TYPE1),pointer:: B,B1
     TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
 	 
     call creat_new_bcmsg
!---------------------------
     do m=1,NB
	  B=>Block(m)
	  do ks=1,B%subface
      Bc=>B%bc_msg(ks)
       if(Bc%bc>0) then  
        call search_1d_bd(m,ks) ! 物理边界
       else
	    call search_1d_lk(m,ks)  ! 有连接的边界
       endif
	  enddo
      
     enddo
!----------------------------
	 call set_ibie  ! 根据ib1,ie1等信息，设定ib,ie
		 
   end

!-------------------------------------------------
   subroutine set_ibie
    use global_var1
    Type (Block_TYPE1),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc
    integer:: m,ks
    do m=1,NBnew
     B => Block_new(m)
    do ks=1,B%subface
     Bc=>B%bc_msg(ks)
     select case (Bc%face)
	  case (1)
	    Bc%ib=1; Bc%ie=1; Bc%jb=Bc%b1 ; Bc%je =Bc%e1 ; Bc%kb=Bc%b2; Bc%ke=Bc%e2
      case (4)
	    Bc%ib=B%nx; Bc%ie=B%nx; Bc%jb=Bc%b1 ; Bc%je =Bc%e1 ; Bc%kb=Bc%b2; Bc%ke=Bc%e2
      case(2)
	    Bc%ib=Bc%b1; Bc%ie=Bc%e1; Bc%jb=1 ; Bc%je =1 ; Bc%kb=Bc%b2; Bc%ke=Bc%e2
      case(5)
	    Bc%ib=Bc%b1; Bc%ie=Bc%e1; Bc%jb=B%ny ; Bc%je =B%ny ; Bc%kb=Bc%b2; Bc%ke=Bc%e2
      case(3)
	    Bc%ib=Bc%b1; Bc%ie=Bc%e1;  Bc%jb=Bc%b2; Bc%je=Bc%e2; Bc%kb=1 ; Bc%ke =1 
      case(6)
	    Bc%ib=Bc%b1; Bc%ie=Bc%e1;  Bc%jb=Bc%b2; Bc%je=Bc%e2; Bc%kb=B%nz ; Bc%ke = B%nz
	 end select 


    enddo
	enddo

    end










  subroutine write_inp_new  
    use global_var1
    Type (Block_TYPE1),pointer:: B,B1
    TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1


   open(99,file="bc3d_new.inc")
    write(99,*) " Inp-liked file of OpenCFD-EC"
    write(99,*) NBnew

	do m=1,NBnew
     B => Block_new(m)
     write(99,*) B%nx, B%ny, B%nz
	 write(99,*) "Block ", m
	 write(99,*) B%subface
	  do ksub=1,B%subface
        Bc => B%bc_msg(ksub)
	   write(99,"(9I6)")  Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
	   write(99,"(12I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
      enddo
	enddo
	close(99)
   end



!----creat bc_msg ----------	
	 subroutine creat_new_bcmsg
	 use global_var1
	 implicit none
     Type (Block_TYPE1),pointer:: B
     TYPE (BC_MSG_TYPE),pointer:: Bc
	 integer:: m,ks

	 do m=1,NBnew
      B=>Block_new(m)
	  B%subface=0
	  allocate(B%bc_msg(Maxsub))
	  do ks=1,Maxsub
	  Bc=>B%bc_msg(ks)
	  Bc%ib=0; Bc%ie=0;Bc%jb=0;Bc%je=0;Bc%kb=0;Bc%ke=0; Bc%bc=0; Bc%face=0; Bc%f_no=0
	  Bc%ib1=0; Bc%ie1=0; Bc%jb1=0; Bc%je1=0; Bc%kb1=0; Bc%ke1=0
	  Bc%nb1=0; Bc%face1=0; Bc%f_no1=0; Bc%L1=0; Bc%L2=0; Bc%L3=0
	  Bc%b1=0; Bc%e1=0; Bc%b2=0; Bc%e2=0
	  Bc%b1t=0; Bc%e1t=0; Bc%b2t=0; Bc%e2t=0           
	 enddo
     enddo
	 end

!------------------------------------------------
! 无连接的边界
     subroutine search_1d_bd(m,ks)   
	 use global_var1
	 implicit none
     TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
     TYPE (bcnp_type),pointer ::BP,BP1,BP2
     Type (Block_TYPE1),pointer:: B,B1
	 	  
	  integer,dimension(2)::ib,ie,ns
      integer:: npk(2),blkm(3)  ! npk 每一维的块数； blkm 三维的块地址
      integer:: ks,q,m,k,blk,idx,k1,k2
	  B=>Block(m)
      Bc=>B%bc_msg(ks)
	  npk(:)=1   ! 每一维的块数
	  call get_ib(ib,ie,ns,blkm,m,ks)    ! 每一维的起始、终止地址
      
	  do q=1,2   ! 2 维
      npk(q)=1
	  BP=>bcn(npk(q),q)
	  
	  call get_bk(blk,idx,ib(q),ns(q),m)
	  BP%blk=blk   ! 块地址
	  BP%idx=idx   ! 偏移地址
	  BP%nn=1
	  do k=ib(q)+1,ie(q)
	  call get_bk(blk,idx,k,ns(q),m)
       if(blk .ne. BP%blk ) then
         npk(q)=npk(q)+1  ! 创立新区
         BP=>bcn(npk(q),q)
		 BP%blk=blk  ! 首元素块地址
		 BP%idx=idx  ! 首元素偏移地址
         BP%nn=1
	    else
		 BP%nn=BP%nn+1
		endif
	   enddo
	
	 enddo
 
!-  write bc msg in new block ----
     do k1=1,npk(1)
	 do k2=1,npk(2)
      BP1=>bcn(k1,1)
	  BP2=>bcn(k2,2)

!	  blkm(ns(1))=k1   ! 面上的第1维
!	  blkm(ns(2))=k2
	  blkm(ns(1))=BP1%blk   
	  blkm(ns(2))=BP2%blk


      blk=B%nbk0+B%PI*B%PJ*(blkm(3)-1)+B%PI*(blkm(2)-1)+blkm(1)-1     ! 新块号
      B1=>block_new(blk)
      B1%subface=B1%subface+1   ! 新建边界区域
	  k=B1%subface
	  Bc1=>B1%bc_msg(k)
      Bc1%b1=BP1%idx ; Bc1%e1=BP1%idx+BP1%nn-1
	  Bc1%b2=BP2%idx ; Bc1%e2=BP2%idx+BP2%nn-1
	  Bc1%bc=Bc%bc
	  Bc1%face=Bc%face
	  Bc1%f_no=k
	 enddo
	 enddo

	end
!---------------------------------------------



!------------------------------------------------
! 有连接的边界
     subroutine search_1d_lk(m,ks)   
	 use global_var1
	 implicit none
     TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
     TYPE (bcnp_type),pointer ::BP,BP1,BP2
     Type (Block_TYPE1),pointer:: B,Bt,B1,B1t
	 	  
	  integer,dimension(2)::ib,ie,ns   
      integer,dimension(2):: ibt,nst,Lt  

	  integer:: npk(2),blkm(3),blkmt(3) ! npk 每一维的块数； blkm 三维的块地址
      integer:: ks,q,m,k,blk,idx,k1,k2,mt,kt,blkt,idxt,nst3
	  integer:: kb1t(3),ke1t(3),nkt(6)

	  B=>Block(m)
      Bc=>B%bc_msg(ks)
	  mt=Bc%nb1   ! 连接块号
      Bt=>Block(mt)  !连接块

	  npk(:)=1   ! 每一维的块数

	  call get_ib(ib,ie,ns,blkm,m,ks)           ! 每一维的起始、终止地址
      call get_ibt(ibt,ns,nst,Lt,blkmt, m,ks)   ! 连接维的起始地址，维号


	  do q=1,2   ! 2 维
      npk(q)=1
	  BP=>bcn(npk(q),q)
	  
	  call get_bk(blk,idx,ib(q),ns(q),m)   ! 块地址，索引地址
      call get_bk(blkt,idxt,ibt(q),nst(q),mt)

	  BP%blk=blk   ! 块地址
	  BP%idx=idx   ! 偏移地址
	  BP%nn=1
	  BP%blkt=blkt
	  BP%idxt=idxt


	  do k=ib(q)+1,ie(q)          ! 被连接 index
	     kt=ibt(q)+Lt(q)*(k-ib(q))   ! 连接 index

	  call get_bk(blk,idx,k,ns(q),m)
	  call get_bk(blkt,idxt,kt,nst(q),mt)

       if(blk .ne. BP%blk  .or. blkt .ne. BP%blkt) then
         npk(q)=npk(q)+1  ! 创立新区
         BP=>bcn(npk(q),q)
		 BP%blk=blk  ! 首元素块地址
		 BP%idx=idx  ! 首元素偏移地址
         BP%nn=1
	     BP%blkt=blkt
		 BP%idxt=idxt
	    else
		 BP%nn=BP%nn+1
		endif
	   enddo
	
	 enddo
 
!-  write bc msg in new block ----
     do k1=1,npk(1)
	 do k2=1,npk(2)
      BP1=>bcn(k1,1)
	  BP2=>bcn(k2,2)

      ! 被连接面上 两维的块索引
!	  blkm(ns(1))=k1   
!	  blkm(ns(2))=k2

	  blkm(ns(1))=BP1%blk   
	  blkm(ns(2))=BP2%blk
      blkmt(nst(1))=BP1%blkt
	  blkmt(nst(2))=BP2%blkt

      blk=B%nbk0+B%PI*B%PJ*(blkm(3)-1)+B%PI*(blkm(2)-1)+blkm(1)-1     ! 新块号
      blkt=Bt%nbk0+Bt%PI*Bt%PJ*(blkmt(3)-1)+Bt%PI*(blkmt(2)-1)+blkmt(1)-1     ! 连接块号
         
	 
	 
	  B1=>block_new(blk)
      B1%subface=B1%subface+1   ! 新建边界区域
	  k=B1%subface
	  Bc1=>B1%bc_msg(k)
      Bc1%b1=BP1%idx ; Bc1%e1=BP1%idx+BP1%nn-1
	  Bc1%b2=BP2%idx ; Bc1%e2=BP2%idx+BP2%nn-1
	  Bc1%bc=Bc%bc
	  Bc1%face=Bc%face
	  Bc1%f_no=k
      

	  Bc1%b1t=BP1%idxt
      Bc1%e1t=Bc1%b1t+Lt(1)*(BP1%nn-1)
   	  Bc1%b2t=BP2%idxt
      Bc1%e2t=Bc1%b2t+Lt(2)*(BP2%nn-1)	  


	  Bc1%L1=Bc%L1
	  Bc1%L2=Bc%L2
      Bc1%L3=Bc%L3
	  Bc1%face1=Bc%face1
	  Bc1%nb1=blkt               ! 连接块号
	  B1t=>block_new(Bc1%nb1)  


!-----------------------------------------------------
      kb1t(nst(1))=min(Bc1%b1t,Bc1%e1t)
	  ke1t(nst(1))=max(Bc1%b1t,Bc1%e1t)

      kb1t(nst(2))=min(Bc1%b2t,Bc1%e2t)
	  ke1t(nst(2))=max(Bc1%b2t,Bc1%e2t)
      nst3=6-nst(1)-nst(2)     ! 另外一维 （退化维）
    
	  nkt(1:3)=1; nkt(4)=B1t%nx; nkt(5)=B1t%ny; nkt(6)=B1t%nz
	  kb1t(nst3)=nkt(Bc1%face1)
	  ke1t(nst3)=nkt(Bc1%face1)
 
      Bc1%ib1=kb1t(1) ; Bc1%ie1=ke1t(1) 
	  Bc1%jb1=kb1t(2) ; Bc1%je1=ke1t(2)
	  Bc1%kb1=kb1t(3) ; Bc1%ke1=ke1t(3)

	 enddo
	 enddo

	end
!---------------------------------------------














!------------------------------------------------
! 设定新的内边界

     subroutine set_new_inner
	 use global_var1
	 implicit none
     TYPE (BC_MSG_TYPE),pointer:: Bc
     Type (Block_TYPE1),pointer:: B,B1,B2
	 	  
     integer:: m,m0,m1,ms
	  
	  do m=1,NBnew
	   B1=>Block_new(m)
       m0=B1%nb_ori
	   B=>Block(m0) 

       if(B1%PI .ne. 1) then   ! 左连接
	     m1=m-1   ! 左侧块号
		 B2=>Block_new(m1)  ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=1 ;  Bc%ie=1 ; Bc%jb=1; Bc%je=B1%ny ; Bc%kb=1 ; Bc%ke=B1%nz
		 Bc%bc=BC_In ;     Bc%face=1  ;  Bc%f_no=ms
		 Bc%ib1=B2%nx; Bc%ie1=B2%nx; Bc%jb1=1; Bc%je1=B2%ny; Bc%kb1=1; Bc%ke1=B2%nz
		 Bc%nb1=m1 ;       Bc%face1=4   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif

       if(B1%PI .ne. B%PI) then   ! 右连接
	     m1=m+1   ! 左侧块号
		 B2=>Block_new(m1)  ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=B1%nx ;  Bc%ie=B1%nx ;  Bc%jb=1; Bc%je=B1%ny ; Bc%kb=1 ; Bc%ke=B1%nz
		 Bc%bc=BC_In ;     Bc%face=4  ;  Bc%f_no=ms
		 Bc%ib1=1; Bc%ie1=1; Bc%jb1=1; Bc%je1=B2%ny; Bc%kb1=1; Bc%ke1=B2%nz
		 Bc%nb1=m1 ;       Bc%face1=1   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif


       if(B1%PJ .ne. 1) then   ! j- 连接
	     m1=m-B%PI      ! 下侧块号
		 B2=>Block_new(m1)  ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=1 ;  Bc%ie=B1%nx ;  Bc%jb=1; Bc%je=1 ;           Bc%kb=1 ; Bc%ke=B1%nz
		 Bc%bc=BC_In ;     Bc%face=2  ;  Bc%f_no=ms
		 Bc%ib1=1; Bc%ie1=B2%nx;   Bc%jb1=B2%ny; Bc%je1=B2%ny;  Bc%kb1=1; Bc%ke1=B2%nz
		 Bc%nb1=m1 ;       Bc%face1=5   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif

       if(B1%PJ .ne. B%PJ) then   ! j+ 连接
	     m1=m+B%PI      ! 上侧块号
		 B2=>Block_new(m1)  ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=1 ;  Bc%ie=B1%nx ;  Bc%jb=B1%ny; Bc%je=B1%ny ;      Bc%kb=1 ; Bc%ke=B1%nz
		 Bc%bc=BC_In ;     Bc%face=5  ;  Bc%f_no=ms
		 Bc%ib1=1; Bc%ie1=B2%nx;   Bc%jb1=1; Bc%je1=1;           Bc%kb1=1; Bc%ke1=B2%nz
		 Bc%nb1=m1 ;       Bc%face1=2   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif


       if(B1%PK .ne. 1) then   ! k- 连接
	     m1=m-B%PI*B%PJ      ! 前侧块号
		 B2=>Block_new(m1)   ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=1 ;  Bc%ie=B1%nx ;  Bc%jb=1; Bc%je=B1%ny ;      Bc%kb=1 ; Bc%ke=1
		 Bc%bc=BC_In ;     Bc%face=3  ;  Bc%f_no=ms
		 Bc%ib1=1;  Bc%ie1=B2%nx;   Bc%jb1=1; Bc%je1=B2%ny;    Bc%kb1=B2%nz; Bc%ke1=B2%nz
		 Bc%nb1=m1 ;       Bc%face1=6   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif


       if(B1%PK .ne. B%PK) then   ! k- 连接
	     m1=m+B%PI*B%PJ      ! 前侧块号
		 B2=>Block_new(m1)   ! 连接块

		 B1%subface=B1%subface+1  ! 创建
		 ms=B1%subface
		 Bc=>B1%bc_msg(ms)
	
		 Bc%ib=1 ;  Bc%ie=B1%nx ;  Bc%jb=1; Bc%je=B1%ny ;      Bc%kb=B1%nz ; Bc%ke=B1%nz
		 Bc%bc=BC_In ;     Bc%face=6  ;  Bc%f_no=ms
		 Bc%ib1=1;  Bc%ie1=B2%nx;   Bc%jb1=1; Bc%je1=B2%ny;    Bc%kb1=1; Bc%ke1=1
		 Bc%nb1=m1 ;       Bc%face1=3   
		 Bc%L1=1 ; Bc%L2=2 ; Bc%L3=3
       endif
     enddo
	end
!---------------------------------------------











!-----找出第1,2维的起始、终止地址
	  subroutine get_ib(ib,ie,ns,blkm,m,ks)    
	  use global_var1
	  implicit none
	  integer::ib(2),ie(2),ns(2),blkm(3),m,ks
      Type (Block_TYPE1),pointer:: B
      TYPE (BC_MSG_TYPE),pointer:: Bc
       
	   B=>Block(m)
	   Bc=>B%bc_msg(ks)

       if(Bc%face==1 .or. Bc%face == 4) then
	    ib(1)=Bc%jb; ie(1)=Bc%je ; ns(1)=2       ! 第1维是 j方向 
		ib(2)=Bc%kb; ie(2)=Bc%ke ; ns(2)=3       ! 第2维是 k方向
	   else if(Bc%face== 2 .or. Bc%face == 5) then
        ib(1)=Bc%ib; ie(1)=Bc%ie ; ns(1)=1
		ib(2)=Bc%kb; ie(2)=Bc%ke ; ns(2)=3
	   else
	    ib(1)=Bc%ib; ie(1)=Bc%ie; ns(1)=1
		ib(2)=Bc%jb; ie(2)=Bc%je; ns(2)=2
	   endif


! 退化维的块索引
	  select case (Bc%face)  
       case(1)
		   blkm(1)=1
	   case(4)
		   blkm(1)=B%PI
	   case(2)
		   blkm(2)=1
	   case(5)
		   blkm(2)=B%PJ
	   case(3)
		   blkm(3)=1
	   case(6)
		   blkm(3)=B%PK
	   end select
    end

 
  ! 被连接块的信息 
  ! ibt(1), ibt(2): 第1，第2连接维 
  ! nst(1), nst(2): 第1,2连接维的 维号
  ! Lt(1),Lt(2): 连接次序 (1 or -1)
      subroutine get_ibt(ibt,ns,nst,Lt,blkmt, m,ks)   
	  use global_var1
	  implicit none
	  integer::ibt(2),ns(2),nst(2),Lt(2),blkmt(3),m,ks,Lk(3),bt(3),et(3),k
      Type (Block_TYPE1),pointer:: B,B1
      TYPE (BC_MSG_TYPE),pointer:: Bc
       
	   B=>Block(m)
	   Bc=>B%bc_msg(ks)
       LK(1)=Bc%L1; LK(2)=Bc%L2; LK(3)=Bc%L3 
       bt(1)=Bc%ib1; bt(2)=Bc%jb1 ; bt(3)=Bc%kb1
       et(1)=Bc%ie1; et(2)=Bc%je1; et(3)=Bc%ke1

	   do k=1,2
       nst(k)=abs(LK(ns(k)))    ! 第1,2维的连接
	   Lt(k)=sign(1,LK(ns(k)))  ! 正或负连接
                               ! 连接的起始地址
		 if(Lt(k) > 0) then
          ibt(k)=bt(nst(k))       
         else
	      ibt(k)=et(nst(k))
		 endif
	   enddo
       
	   B1=>block(Bc%nb1)

! （被连接块）退化维的块索引
	  select case (Bc%face1)  
       case(1)
		   blkmt(1)=1
	   case(4)
		   blkmt(1)=B1%PI
	   case(2)
		   blkmt(2)=1
	   case(5)
		   blkmt(2)=B1%PJ
	   case(3)
		   blkmt(3)=1
	   case(6)
		   blkmt(3)=B1%PK
	   end select

    end




! -----------
 !  给定下标k, 计算出是第bk个子块,本地下标bi
  subroutine get_bk(bk,ki,k,ns,mb)    ! k, 下标， ns 维数, nb (原)块号
     use global_var1
     implicit none
     integer:: bk,ki,k,ns,mb,Pn,nn            ! Pn 块数， nn 网格点数
     Type (Block_TYPE1),pointer:: B
	 B=> block(mb)
     if(ns == 1) then  ! 第1维
	  nn=B%nx
	  Pn=B%Pi
	 else if(ns==2) then
	  nn=B%ny
	  Pn=B%Pj
	 else if(ns==3) then
	  nn=B%nz
	  Pn=B%Pk
	 else
	  print*, "Error at get_bk !!"
	  stop
	 endif
     call get_bki(bk,ki,k,nn,Pn)
  end

  subroutine get_bki(bk,ki,k,nn,Pn)
     implicit none
	 integer:: bk,k,nn,Pn,k1,ki
	     bk=0
  	     do k1=1,Pn
          if( k >=int(nn*(k1-1)/Pn)+1  .and. k<=int(nn*k1/Pn) ) then
           bk=k1
	       exit 
 	      endif
		 enddo
   
         if(bk==0) then
		 print*, "error at get_bki !!"
		 stop
		 endif
        
		ki=k-int(nn*(bk-1)/Pn)
  end
   









!--------------------------------------------------------

   subroutine read_cutfile
   use global_var1
   implicit none
   integer:: m
   	
	NBnew=0
	open(99,file="Mesh3d.cut")
	read(99,*)
    read(99,*) Mesh_form, NB
    allocate(PI(NB),PJ(NB),PK(NB))
    do m=1,NB
	 read(99,*)
	 read(99,*) PI(m),PJ(m),PK(m)
	 NBnew=NBnew+PI(m)*PJ(m)*PK(m)
	enddo
    close(99)
    print*, "NB=",NB, " NBnew= ",NBnew
   end
!---------------------------------------------------
   
   subroutine cut_mesh
   use global_var1
   implicit none
   integer:: NB0,NB1,i,j,k,m,nx,ny,nz,nx1,ny1,nz1,i1,j1,k1,m0,m1,n
   integer:: i2,j2,k2,m2
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE1),pointer:: B,B1
    

    allocate(Block(NB))
	allocate(Block_new(NBnew))
    allocate(NI(NB),NJ(NB),NK(NB))

     if(mesh_form==0) then
	  open(100,file="Mesh0.dat",form="unformatted")
	  read(100) NB1
	else
	  open(100,file="Mesh0.dat")
	  read(100,*) NB1
    endif	 


	if(NB1 .ne. NB) then
	 print*, "Block number in Mesh0.dat is not equal to that in Mesh3d.cut! "
	 stop
	endif

    if(mesh_form==0) then
     read(100) ((NI(k),NJ(k),Nk(k)),k=1,NB)
    else
     read(100,*) ((NI(k),NJ(k),Nk(k)),k=1,NB)
	endif

!================================================	 
     m1=1
	 do m=1,NB
 	  B=>Block(m)
      B%Block_no=m
	  B%Pi=PI(m)
	  B%Pj=PJ(m)
	  B%Pk=PK(m)
      B%nbk0=m1    ! 起始子块号
      B%nx=NI(m)
	  B%ny=NJ(m)
	  B%nz=NK(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
      allocate(B%xyz(nx,ny,nz,3))
      
	  if(mesh_form==0) then
	    read(100) ((((B%xyz(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,3)
      else
	   read(100,*) ((((B%xyz(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,3)
	  endif


	   do k1=1,PK(m)
	   do j1=1,PJ(m)
	   do i1=1,PI(m)
	     B1=>Block_new(m1)
		 B1%block_no=m1        ! 块号
		 B1%nb_ori=m           ! 原块号
		 B1%Pi=i1      ! i 方向的块序号
		 B1%Pj=j1
		 B1%Pk=k1
		 
		 B1%nx=int((nx*i1)/PI(m))- int((nx*(i1-1))/PI(m))
		 B1%ny=int((ny*j1)/PJ(m))- int((ny*(j1-1))/PJ(m))
		 B1%nz=int((nz*k1)/PK(m))- int((nz*(k1-1))/PK(m))
         nx1=B1%nx; ny1=B1%ny ; nz1=B1%nz
         B1%ib= int((nx*(i1-1))/PI(m))+1
		 B1%jb= int((ny*(j1-1))/PJ(m))+1
		 B1%kb= int((nz*(k1-1))/PK(m))+1
 		 allocate(B1%xyz(nx1,ny1,nz1,3))
 	 	 B1%block_no=m1
		 
		 do m2=1,3
         do k2=1,nz1
		 do j2=1,ny1
		 do i2=1,nx1
		  i=B1%ib+i2-1
		  j=B1%jb+j2-1
		  k=B1%kb+k2-1
          B1%xyz(i2,j2,k2,m2)=B%xyz(i,j,k,m2)
         enddo
 		 enddo
		 enddo
		 enddo
		 m1=m1+1

	 enddo
	 enddo
	 enddo

   enddo
   close(100)

   end


!--------------------------------------------------------
  subroutine write_new_mesh
   use global_var1
   implicit none
   integer:: i,j,k,m,n
   Type (Block_TYPE1),pointer:: B
   open(99,file="Mesh3d_new.dat",form="unformatted")
   write(99) NBnew
   write(99) (Block_new(m)%nx,Block_new(m)%ny,Block_new(m)%nz,m=1,NBnew)
   do m=1,NBnew
   B=> Block_new(m)
   write(99) ((((B%xyz(i,j,k,n),i=1,B%nx),j=1,B%ny),k=1,B%nz),n=1,3)
   enddo
   close(99)
   end
        


!------------------------------------------------------
  subroutine convert_inp_inc 
   use global_var1
   implicit none
  
   integer:: NB1,m,ksub,nx,ny,nz,k,j,k1,ksub1
   integer:: kb(3),ke(3),kb1(3),ke1(3),s(3),p(3),Lp(3)
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
   read(88,*) NB1
   if(NB1 .ne. NB ) then
    print*, "NB in bc3d.inp is not correct !"
	stop
   endif 

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

  call search_fno1(1)

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



!-----------------------------------------------------
! 搜索块，建立fno1
! Kflag==1 search old block;  2 new block
subroutine search_fno1(Kflag)
   use global_var1
   implicit none
   integer:: NB1,m,kflag,ksub,ksub1
   Type (Block_TYPE1),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
   
   if(Kflag==1) then
    NB1=NB
   else
    NB1=NBnew
   endif

!  搜索连接块的块号 f_no1  (便于MPI并行通信是使用)
   do m=1,NB1
     if(Kflag==1) then
 	  B => Block(m)
     else
	  B=> Block_new(m)
	 endif

	 do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .lt. 0) then
         Bc%f_no1=0
		 
		 if(Kflag==1) then
		   B1=>Block(Bc%nb1)         ! 指向连接块
         else
		   B1=>Block_new(Bc%nb1)
		 endif


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
  end

