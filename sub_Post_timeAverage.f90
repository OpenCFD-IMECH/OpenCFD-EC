!  后处理模块： 进行时间平均
!--------------------------------------------------------
  subroutine Time_average      
   use Global_Var
   implicit none
   integer:: i,j,k,m,mB,nf,nx,ny,nz
   integer,save:: Iflag=0
   real(PRE_EC):: d1,u1,v1,w1,T1,tmp

   Type (Block_TYPE),pointer:: B

!---仅被运行1次 -----------------
   if(Iflag == 0) then
     Iflag=1             
    do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
 	 allocate(B%U_average(0:nx,0:ny,0:nz,5))       ! 时均量 d,u,v,w,T
     enddo

    call init_average        ! 初始化平均场

   endif
!-------------------------------------
! 时间平均  
   if(my_id .eq. 0) print*, "Time Average ......", Istep_average+1

   tmp=1.d0/(Istep_average+1.d0)
   do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz

!$OMP PARALLEL DO PRIVATE(i,j,k,d1,u1,v1,w1,T1) SHARED(nx,ny,nz,B,Cv,tmp,Istep_average)

     do k=0,nz
	 do j=0,ny
	 do i=0,nx
         d1= B%U(1,i,j,k)
         u1= B%U(2,i,j,k)/d1
         v1= B%U(3,i,j,k)/d1
         w1= B%U(4,i,j,k)/d1
         T1=(B%U(5,i,j,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))/(Cv*d1)
      
	   B%U_average(i,j,k,1)=(Istep_average*B%U_average(i,j,k,1)+d1)*tmp    
	   B%U_average(i,j,k,2)=(Istep_average*B%U_average(i,j,k,2)+u1)*tmp    
	   B%U_average(i,j,k,3)=(Istep_average*B%U_average(i,j,k,3)+v1)*tmp    
	   B%U_average(i,j,k,4)=(Istep_average*B%U_average(i,j,k,4)+w1)*tmp    
	   B%U_average(i,j,k,5)=(Istep_average*B%U_average(i,j,k,5)+T1)*tmp    

	 enddo
	 enddo
	 enddo
   enddo
    Istep_average=Istep_average+1
  end

!-----------------------------------------------------
! 初始化, 目前版本只支持重新开始平均，暂不支持读取flow3d_average.dat
   subroutine init_average      
   use Global_Var
   implicit none
   integer:: i,j,k,mB,nf,nx,ny,nz

   Type (Block_TYPE),pointer:: B
    Istep_average=0
    do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
      do k=0,nz
	  do j=0,ny
	  do i=0,nx
	   B%U_average(i,j,k,1)=0.d0 
	   B%U_average(i,j,k,2)=0.d0    
	   B%U_average(i,j,k,3)=0.d0    
	   B%U_average(i,j,k,4)=0.d0    
	   B%U_average(i,j,k,5)=0.d0    
	  enddo
	  enddo
	  enddo
     enddo
  end
!----------------------------------------------------      

 !  输出平均量 （Plot3d格式）, 最细网格flow3d_average.dat  
  subroutine output_flow_average
   use Global_Var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:,:):: U
   integer:: NB,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)

   character(len=50):: filename   
   
   MP=>Mesh(1)

!---------------------------------------------------------------
 if(my_id .eq. 0) then
  
   print*, "write flow3d_average.dat ......"
   
   open(99,file="flow3d_average.dat",form="unformatted")                    ! d,u,v,w,T

   do m=1, Total_block   ! 全部块
     
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz,5))

	if(B_proc(m) .eq. 0) then             ! 这些块属于根进程
      mt=B_n(m)                           ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	 
	   do m1=1,5
	    do k=0,nz
	    do j=0,ny
	    do i=0,nx
		  U(i,j,k,m1)=B%U_average(i,j,k,m1)                     
        enddo
	    enddo
	    enddo
 	   enddo

    else                        ! 接收该块信息
	   Num_data=5*(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! 在该块中的编号
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
    endif
! write Data ....
    
	write(99) (((( U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=1,5)    
    
	deallocate(U)

   enddo
   
   write(99) Istep_average
   close(99)
 
 else     ! 非0节点

    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,5))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*5
 	  tag=m
	 
      do m1=1,5
	    do k=0,nz
	    do j=0,ny
	    do i=0,nx
		 U(i,j,k,m1)=B%U_average(i,j,k,m1)
        enddo
	    enddo
	    enddo
 	   enddo
	   
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "write flow3d_average.dat OK"

  end subroutine output_flow_average



