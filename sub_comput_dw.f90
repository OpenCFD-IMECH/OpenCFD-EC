!-----------------------------------------------------
!  计算到壁面的距离    
!  Code by Li Xinliang

  module Wall_dist
   use precision_EC
   implicit none
   integer,save:: Npw_total                ! 总壁面网格数
   real(PRE_EC),allocatable,dimension(:,:):: Xw   ! 壁面点的坐标
  end module Wall_dist

 subroutine Comput_dist_wall
    use  Global_Var
    use  Wall_dist
    implicit none
    logical EXT
	integer:: Iext,ierr

!  判断 'wall_dist.dat'文件是否存在， 如存在则读取； 如不存在，则计算到壁面的距离

   if(my_id .eq. 0) then
     Inquire(file="wall_dist.dat",Exist=EXT)
 	 if(EXT) then
	   print*, "Find 'wall_dist.dat', read it ..."
	   Iext=1
	 else
	   print*, "Can not find the file 'wall_dist.dat', Comput wall distance ... "
	   Iext=0
     endif
   endif
   call MPI_Bcast(Iext,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)

   if(Iext .eq. 1) then
     call read_dw
   else
     call wall_point         ! 收集壁面网格点
     call comput_wall_dist          ! 计算各点到壁面的距离
     call write_dw           ! 写入文件
   endif

   end subroutine Comput_dist_wall
     
   subroutine comput_wall_dist
    use  Global_Var
    use  Wall_dist
    implicit none
    integer::  mb,i,j,k,k1,nx,ny,nz,ierr
	real(PRE_EC),allocatable,dimension(:,:,:):: dis
	real(PRE_EC):: d1, Dis_max,  Dis_min, Dis_max0,  Dis_min0
	real(PRE_EC),parameter:: d_init=1.d8
	Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer:: BC

     Dis_max=0.d0; Dis_min=d_init
     MP=>Mesh(1)
     do mb=1,MP%Num_Block
       B=>MP%Block(mb)
	    nx=B%nx; ny=B%ny; nz=B%nz
        allocate(dis(nx,ny,nz))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (i,j,k,k1,d1)
		 do k=1,nz
          do j=1,ny
           do i=1,nx
		    dis(i,j,k)=d_init
		    do k1=1,Npw_total 
             d1=sqrt((B%x(i,j,k)-Xw(1,k1))**2+(B%y(i,j,k)-Xw(2,k1))**2+(B%z(i,j,k)-Xw(3,k1))**2)
             if(d1 .lt. dis(i,j,k)) dis(i,j,k)=d1
             enddo
            enddo
           enddo
		  enddo
 !$OMP END PARALLEL DO      
	   
     do k=1,nz-1
      do j=1,ny-1
       do i=1,nx-1
        B%dw(i,j,k)=(dis(i,j,k)+dis(i,j+1,k)+dis(i,j,k+1)+dis(i,j+1,k+1)+ &
            dis(i+1,j,k)+dis(i+1,j+1,k)+dis(i+1,j,k+1)+dis(i+1,j+1,k+1))*0.125    ! 格心点上的值
        Dis_max=max(Dis_max,B%dw(i,j,k))
        Dis_min=min(Dis_min,B%dw(i,j,k))
       enddo
      enddo
     enddo

    deallocate(dis)
   enddo
    call MPI_ALLREDUCE(Dis_max,Dis_max0,1,OCFD_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(Dis_min,Dis_min0,1,OCFD_DATA_TYPE,MPI_MIN,MPI_COMM_WORLD,ierr)
    if(my_id .eq. 0) then
	  print*, "----------comput wall dist OK ----------------"
	  print*, "Maximum distance to the wall is :" , Dis_max0
	  print*, "Minimum distance to the wall is :" , Dis_min0
	  print*, "----------------------------------------------"
	endif

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
   

  end subroutine comput_wall_dist
!-----------------------------------------
  
  
   
   
! 收集壁面上的网格点  
  subroutine wall_point
    use  Global_Var
    use  Wall_dist
    implicit none
    integer:: Npw1, mb,ksub,i,j,k,m,k0,ierr
	Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer:: BC
    integer,allocatable,dimension(:):: Npw,Nc,displs  ! 各进程壁面点的数目
    real(PRE_EC),allocatable,dimension(:,:):: Xw1  ! 壁面点（本进程）
    allocate(Npw(0:Total_proc-1),Nc(0:Total_proc-1),displs(0:Total_proc-1))

    MP=> Mesh(1)
!     统计壁面网格点的数目
      Npw1=0
      do mb=1,MP%Num_Block
          B=>Mp%Block(mb)
           do ksub=1,B%subface
             Bc=>B%bc_msg(ksub)
              if(Bc%bc .eq. BC_Wall) then
                Npw1=Npw1+(Bc%ie-Bc%ib+1)*(Bc%je-Bc%jb+1)*(Bc%ke-Bc%kb+1)
              endif
            enddo
	   enddo
       allocate(Xw1(3,Npw1))

!    读入本进程的壁面网格坐标
      k0=1
      do mb=1,MP%Num_Block
          B=>Mp%Block(mb)
           do ksub=1,B%subface
             Bc=>B%bc_msg(ksub)
              if(Bc%bc .eq. BC_Wall) then
                do k=Bc%kb,Bc%ke
                  do j=Bc%jb,Bc%je
                    do i=Bc%ib,Bc%ie
                      Xw1(1,k0)=B%x(i,j,k)
                      Xw1(2,k0)=B%y(i,j,k)
                      Xw1(3,k0)=B%z(i,j,k)
                      k0=k0+1
                    enddo
                   enddo
                 enddo
              endif
            enddo
	   enddo

! 全部进程的总壁面网格数
       call MPI_ALLREDUCE(Npw1,Npw_Total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
       allocate(Xw(3,Npw_total))
! 进行全搜集操作，得到全部的壁面网格坐标
       call MPI_Allgather(Npw1,1,MPI_Integer,Npw,1,MPI_Integer,MPI_COMM_WORLD,ierr)
	   Nc(:)=3*Npw(:)   !数据量 
	   displs(0)=0
	   do m=1,Total_proc-1
	     displs(m)=displs(m-1)+Nc(m-1)
	   enddo
	   call MPI_Allgatherv(Xw1,3*Npw1,OCFD_DATA_TYPE,Xw,Nc,displs,OCFD_DATA_TYPE,MPI_COMM_WORLD,ierr)
       deallocate(Xw1,Npw,Nc,displs)

     end



