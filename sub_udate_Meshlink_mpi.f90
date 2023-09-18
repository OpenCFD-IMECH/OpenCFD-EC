!----多块网格之间的连接--------------------------------------------------------------------------
! 交换LAP层 网格 (格心点) 信息


!---------Continue boundary (inner boundary) ------------------------------------------------------
!  传送相邻面的中心坐标给本块的虚网格  （便于使用周期性边界条件）； 
!  由于坐标存储在节点，物理量存储在中心点， 因而传输坐标与传输物理量的方法（下标对应方式）略有区别
!  MPI并行版； 
!  策略： 如两块(block)位于同一进程(proc),则采用直接传递(不通过MPI),以提高效率；
!  使用MPI传递时，首先使用MPI_Bsend()发送全部信息； 然后使用MPI_recv()接收。 

     subroutine update_Mesh_Center(nMesh)
     use Global_Var
     use interface_defines
     implicit none
     integer:: nMesh,ierr
! --------------------------------------------------------------------------------------- 
! 模块边界通信：  MPI 版本
    call Mesh_send_mpi(nMesh)   ! 使用MPI发送全部信息 ；如目标块也在本进程内，则不通过MPI,直接交换信息.
    call Mesh_recv_mpi(nMesh)   ! 接收全部信息
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call Mesh_Center_Periodic(nMesh)     ! 利用周期条件，修改坐标 （旋转一定角度, 或增加Dx,Dy,Dz）

  end subroutine update_Mesh_Center

!---------------------------------------------------------------

! 同一进程内两块之间的通信 （可直接通信）

     subroutine Mesh_send_mpi(nMesh)  
     use Global_Var
     use interface_defines
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,i1,j1,k1,mb,mBlock,ksub,nMesh,Send_to_ID,Num_Data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),Ks(3)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Usend

 ! --------------------------------------------------------------------------------------- 
 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)

  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

!-----------------------------------------------------------------------------------------------
      if(Bc%bc .ge. 0 ) cycle          ! 非内边界 inner boundary
      
!--------------------------------------------------------
!  B块为源数据； B1块为目标数据；   数据从源数据写入目标数据

!      源数据:临近边界的LAP层 内点;   i=kb(1) to ke(1), j=kb(2) to ke(2),  k= kb(3) to ke(3)
!      目标数据 （临近边界的LAP层 Ghost点,  现版本LAP=2）

!      源数据 
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                 ! k=1,2,3 为i,j,k方向
	   if(Bc%face .gt. 3) kb(k)=kb(k)-LAP   ! i+, j+ or k+ 面
       ke(k)=kb(k)+LAP-1

!      目标数据 （边界外扩充的LAP层Ghost点）
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
       k=mod(Bc%face1-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face1 .le. 3)  kb1(k)=kb1(k)-LAP   ! i-, j- or k- 面
       ke1(k)=kb1(k)+LAP-1


!  L1,L2,L3 : 描述维的连接 ; P(1),P(2),P(3): 描述连接次序(顺序or 逆序)
       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)    ! L1=1 意味着（本维）与目标数据的第1维连接， =2 为与第2维连接, ...
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)

!   ks(k) : 目标数据第k维的起始下标
	   do k=1,3  !   目标数据 起始下标     
	    if(P(k) .gt. 0) then
		 ks(k)=kb(k)     ! 顺序，从kb开始
		else
		 ks(k)=ke(k)     ! 拟序， 从ke开始
	    endif
	  enddo
!
!-------------------------------------------------------------
      allocate(Usend(3,kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
       do k=kb1(3),ke1(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
		     ka(1)=i-kb1(1)            
			 ka(2)=j-kb1(2)
			 ka(3)=k-kb1(3)
			 i1=ks(1)+ka(L1)*P(1)        ! 源数据的下标 (i1,j1,k1), 从kb1到ke1 (需考虑 a.顺序或逆序, b. 维的连接), L1,L2,L3控制维的连接，P(k)控制顺序/逆序
			 j1=ks(2)+ka(L2)*P(2)
			 k1=ks(3)+ka(L3)*P(3)
			 Usend(1,i,j,k)=B%xc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
			 Usend(2,i,j,k)=B%yc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
			 Usend(3,i,j,k)=B%zc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
           
		   enddo
		 enddo
		enddo
  
  !    
	 Send_to_ID=B_proc(Bc%nb1)              ! 发送目标块所在的进程号

	 if( Send_to_ID .ne. my_id) then        ! 目标块不在本进程内
       Num_data=3*(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)      ! 数据量
       tag=  Bc%nb1*1000+Bc%f_no1                                         ! 标记; 块号Bc%nb1, 子面号Bc%f_no1 （向一个进程发送多个数据包时，用于识别）
	   call MPI_Bsend(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
!	   call MPI_send(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
	 else 
!        目标块也在本进程内，直接写入 (不通过MPI)
         mb=B_n(Bc%nb1)          ! 块内部编号
		 B1 =>  Mesh(nMesh)%Block(mb)  ! 相邻块 （目标块）
	     do k=kb1(3),ke1(3)             ! 
	      do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
			  B1%xc(i,j,k)=Usend(1,i,j,k)    
			  B1%yc(i,j,k)=Usend(2,i,j,k)    
			  B1%zc(i,j,k)=Usend(3,i,j,k)    
            enddo
		  enddo
		 enddo
	  endif	
	  deallocate(Usend)

  enddo
  enddo
 
  end
!--------------------------------------------------------------------------

! 将全部信息接收

subroutine Mesh_recv_mpi(nMesh) ! 使用MPI发送全部信息
     use Global_Var
     use interface_defines
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
	 integer:: tag,Num_DATA,ierr,Status(MPI_Status_SIZE)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Urecv

!---------------------------------------------------------------------------------------
! 接收阶段，将全部信息接收
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ge. 0 ) cycle               ! 内边界 inner boundary 或周期性边界
   Recv_from_ID=B_proc(Bc%nb1)           ! 相邻块（接收源块）所在的进程号
   if(Recv_from_ID .eq. my_id) cycle     ! 源块在本进程内，不使用MPI通信 (Umessage_send_mpi()已完成写入操作)  
   
!------------------------------------------------------------------------------------------------
!      目标数据 （边界外扩充的LAP层Ghost点）
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1

       allocate(Urecv(3,kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))    ! 接收数组，按照目标数据的格式
       Num_data=3*(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)      ! 数据量
	   tag=B%Block_no*1000+Bc%f_no                                  ! tag 标记，标记块号+子面号
	   
	   call MPI_Recv(Urecv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,MPI_COMM_WORLD,status,ierr)

        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%xc(i,j,k)=Urecv(1,i,j,k)
             B%yc(i,j,k)=Urecv(2,i,j,k)
             B%zc(i,j,k)=Urecv(3,i,j,k)
           enddo
		  enddo
		enddo
       deallocate(Urecv)
   enddo
   enddo
  end


!------------------------------------------------------------------------
! 根据周期性条件， 对Ghost 点的坐标进行调整
! 左侧(BC_PeriodicL)的点 -Turbo_Seta;  右侧的点 + Turbo_Seta
 
  subroutine Mesh_Center_Periodic(nMesh)
     use Global_Var
     use interface_defines
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,kb(3),ke(3)
     real(PRE_EC):: Seta,SetaP,rr,Xp,Yp,Zp
!---------------------------------------------------------------------------------------
! 
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ne. BC_PeriodicL .and. Bc%bc .ne. BC_PeriodicR ) cycle               ! 非周期性边界
    if(Bc%bc == Bc_PeriodicL ) then    ! 左周期边界
	    setaP= -Turbo_Periodic_Seta    !  叶轮机周向计算域的周角
	    Xp = - Periodic_dX;  Yp=- Periodic_dY; Zp= - Periodic_dZ
	endif

	if(Bc%bc == BC_PeriodicR)  then
	  SetaP=  Turbo_Periodic_Seta
      Xp =  Periodic_dX;  Yp= Periodic_dY; Zp=  Periodic_dZ
	endif

!------------------------------------------------------------------------------------------------
!    边界外扩充的LAP层Ghost点
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1
	   
     if( IF_TurboMachinary == 1) then  ! 叶轮机模式
        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             rr=sqrt(B%yc(i,j,k)**2+B%zc(i,j,k)**2)
			 seta=acos(B%yc(i,j,k)/rr)
			 if(B%zc(i,j,k) < 0) seta=-seta
 			 seta=seta+SetaP     ! 根据周期性条件， 旋转 SetaP 角度
             B%yc(i,j,k)=rr*cos(seta)
             B%zc(i,j,k)=rr*sin(seta)
           enddo
		  enddo
		enddo
	 else             ! 非叶轮机模式
        do k=kb(3),ke(3)             
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%xc(i,j,k)=B%xc(i,j,k)+Xp
             B%yc(i,j,k)=B%yc(i,j,k)+Yp
             B%zc(i,j,k)=B%zc(i,j,k)+Zp
           enddo
		  enddo
		enddo
     endif


   enddo
   enddo
  end
