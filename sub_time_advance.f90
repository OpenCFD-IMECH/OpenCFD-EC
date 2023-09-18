!----------------------------------------------------------------------
! 在给定的网格上求解N-S方程 （推进1个时间步）
! 对于单重网格，nMesh=1;  对于多重网格，nMesh=1,2,3, ... 分别对应用细网格、粗网格、更粗网格 ...
! 2015-11-26: A bug in Line 299 is removed   (KRK should be a shared data)
 
  subroutine NS_Time_advance(nMesh)
   use Global_var
   implicit none
   integer:: nMesh
   if(Time_Method .eq. Time_Euler1) then
     call NS_Time_advance_1Euler(nMesh)                    ! 1阶Euler
   else if (Time_Method .eq. Time_LU_SGS ) then            !  LU_SGS
     call NS_Time_advance_LU_SGS(nMesh)
   else if (Time_Method .eq. Time_Dual_LU_SGS) then        ! Dual_LU_SGS
     call  NS_Time_Dual_LU_SGS(nMesh)
   else if (Time_Method .eq. Time_RK3) then
     call NS_Time_advance_RK3(nMesh)                       ! 3阶RK
   else
     print*, "This time advance method is not supported!!!"
   endif
    call force_vt_kw(nMesh)     ! 强制 vt, k,w 非负

  end subroutine NS_Time_advance

!---------------------------------------------------------------------------------------------
! 强制vt, k,w非负   
   subroutine force_vt_kw(nMesh)
    use Global_var
    implicit none
    integer:: nMesh,mBlock,nx,ny,nz,i,j,k
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
  
    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny ; nz=B%nz
    if(MP%NVAR == 6) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
      do k=1,nz-1
 	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(6,i,j,k) < 0)  B%U(6,i,j,k)=1.d-10
      enddo
      enddo
	  enddo
!$OMP END PARALLEL DO
    else if (MP%NVAR == 7) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
 	  do k=1,nz-1
	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(6,i,j,k) < 0)  B%U(6,i,j,k)=1.d-10
       if(B%U(7,i,j,k) < 0)  B%U(7,i,j,k)=1.d-10
	  enddo
      enddo
	  enddo
!$OMP END PARALLEL DO
    endif
   enddo    
  end
!--------------------------------------------------------------------------------------











! 采用 LU_SGS方法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
   call Set_Un(nMesh)
   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差 (以及Du)
   if(nMesh .ne. 1) call Add_force_function(nMesh)   !  添加强迫函数（多重网格的粗网格使用）
  
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!   时间推进 

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             B%U(m,i,j,k)=B%Un(m,i,j,k)+B%dU(m,i,j,k)           ! LU_SGS方法
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
  
  enddo

!----------------------------------------------------------------   
    if( IFLAG_LIMIT_Flow == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

!---------------------------------------------------------------------------------------  
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_advance_LU_SGS
!--------------------------------------------------------------------------------------



!  采用双时间步长法 LU_SGS方法进行时间推进一个时间步 
!  目前Dual LU_SGS 方法尚不支持多重网格,因而nMesh只能为1

  subroutine NS_Time_Dual_LU_SGS(nMesh)
    use Global_var
    implicit none
    integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz,Kt_in
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
    real(PRE_EC):: max_res
    
	MP=>Mesh(nMesh)
    NVAR1=MP%NVAR
 do kt_in=1, step_inner_Limit                      ! 内循环迭代

   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差及Du
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             B%U(m,i,j,k)=B%U(m,i,j,k)+B%dU(m,i,j,k)           ! LU_SGS方法
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
   enddo

!----------------------------------------------------------------   
    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 


  call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
  call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
  call comput_max_Res_onemesh(nMesh)                 ! 计算最大残差及均方根残差

     max_res=MP%Res_rms(1)       ! 最大均方根残差 (作为内迭代标准)
     do m=1,NVAR1
 	  max_res=max(max_res,MP%Res_rms(m))
	 enddo
     if( max_res .le. Res_Inner_Limit) exit   ! 达到残差标准，跳出内迭代
 enddo
   
   if(my_id .eq. 0) then 	
	 print*, "Inner step ... ", kt_in
	 print*, "rms residual eq =", MP%Res_rms(1:NVAR1)
   endif

    
 
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
			 B%Un1(m,i,j,k)=B%Un(m,i,j,k)        
             B%Un(m,i,j,k)=B%U(m,i,j,k)  
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
   
  enddo


   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_Dual_LU_SGS
!--------------------------------------------------------------------------------------












! 采用1阶Euler法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_1Euler(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
   call Set_Un(nMesh)
   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差
   if(nMesh .ne. 1) call Add_force_function(nMesh)   !  添加强迫函数（多重网格的粗网格使用）

    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!   时间推进 
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             du=B%Res(m,i,j,k)/B%vol(i,j,k)   
             B%U(m,i,j,k)=B%Un(m,i,j,k)+B%dt(i,j,k)*du
           enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
  
  enddo    


!----------------------------------------------------------------   
    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

!---------------------------------------------------------------------------------------  
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_advance_1Euler
!----------------------------------------------------------------------------------------


! 采用3阶RK方法推进1个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_RK3(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
 
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)

!$OMP PARALLEL DO PRIVATE(i,j,k,m) SHARED(NVAR1,B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Un(m,i,j,k)=B%U(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
 
   enddo

   do KRK=1,3                                          ! 3-step Runge-Kutta Method
	 call Comput_Residual_one_mesh(nMesh)              ! 计算残差
     if(nMesh .ne. 1) call Add_force_function(nMesh)   ! 添加强迫函数（多重网格的粗网格使用）
	 do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)                  ! 第nMesh 重网格的第mBlock块
       nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!    时间推进

!$OMP PARALLEL DO PRIVATE(i,j,k,m,du) SHARED(NVAR1,nx,ny,nz,Ralfa,Rbeta,Rgamma,B,KRK)
       do k=1,nz-1 
         do j=1,ny-1
           do i=1,nx-1
             do m=1,NVAR1
		       du=B%Res(m,i,j,k)/B%Vol(i,j,k)  
               B%U(m,i,j,k)=Ralfa(KRK)*B%Un(m,i,j,k)+Rgamma(KRK)*B%U(m,i,j,k)+B%dt(i,j,k)*Rbeta(KRK)*du        ! 3阶RK
             enddo
           enddo
         enddo
	   enddo
 !$OMP END PARALLEL DO 
   enddo    

!---------------------------------------------------------------------------------------

    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

 
     call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
     call update_buffer_onemesh(nMesh)              ! 同步各块的交界区
   enddo   
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global          ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1            ! 计算步数

  end subroutine NS_Time_advance_RK3


! 计算最大残差和均方根残差（整个网格）
  subroutine comput_max_Res_onemesh(nMesh)
   use Global_var
   implicit none
	integer:: nMesh,mBlock,i,j,k,m,ierr
	logical F_NaN
 	real(PRE_EC):: Res,Res_max(7),Res_rms(7)
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
 
     MP=> Mesh(nMesh)
   
      Res_max(:)=0.d0
	  Res_rms(:)=0.d0
 
   do mBlock=1,MP%NUM_BLOCK 
!     call comput_max_Res_oneblock(nMesh,mBlock)
   	  B => MP%Block(mBlock)                 !第nMesh 重网格的第mBlock块

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(MP,B) REDUCTION(MAX: Res_max) REDUCTION(+: Res_rms)   
	 do k=1,B%nz-1
       do j=1,B%ny-1
         do i=1,B%nx-1
! -------------------------------------------------------------------------------------------
!    时间推进
           do m=1,MP%NVAR
             Res=B%Res(m,i,j,k)
!--------------------------------------------------------------------------------------------------
! detech "NaN", Since Ver 0.72, which is useful for debug
             F_NaN=Isnan(Res)
             if(F_NaN) then
 		       print*, "NaN in Residual is found !, In block",B%block_no
			   print*, "location i,j,k,m=",i,j,k,m
!               B%Res(m,i,j,k)=0.d0    ! 强制为0
			    print*, "Stop"
			   stop
		     endif  
              Res_max(m)=max(Res_max(m),abs(Res))    ! 最大残差
			  Res_rms(m)=Res_rms(m)+Res*Res           ! 均方根残差
!--------------------------------------------------------------------------------------------------       
	       enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO
  enddo
 
    call MPI_ALLREDUCE(Res_max(1),MP%Res_max(1),MP%NVAR,OCFD_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(Res_rms(1),MP%Res_rms(1),MP%NVAR,OCFD_DATA_TYPE,MPI_SUM,MPI_COMM_WORLD,ierr)
    MP%Res_rms(:)=sqrt(MP%Res_rms(:)/(MP%Num_Cell))   !均方根残差

  end  subroutine comput_max_Res_onemesh

!-------------------------------------------------------------

!--------------------------------------------------------------
! 打印残差（最大残差和均方根残差）
  subroutine output_Res(nMesh)
   use Global_var
   implicit none
	integer:: nMesh
    call   comput_max_Res_onemesh(nMesh)
!-----------------------------------
   if(my_id .eq. 0) then
    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
    print*, "----------The Max Residuals are-------- ", " ---Mesh---",nMesh
    write(*, "(7E20.10)") Mesh(nMesh)%Res_max(:)
    print*, "  The R.M.S Residuals are "
    write(*, "(7E20.10)") Mesh(nMesh)%Res_rms(:)
    open(99,file="Residual.dat",position="append")
    write(99,"(I8,15E20.10)") Mesh(nMesh)%Kstep, Mesh(nMesh)%Res_max(:),Mesh(nMesh)%Res_rms(:)
    close(99) 
   endif

  end  subroutine output_Res

!-------------------------------------------------------------


!----------------------------------------------------------
! 对SA,SST方程的物理量进行限制
  subroutine limit_vt(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   Type (Block_TYPE),pointer:: B
   integer nMesh,mBlock,NVAR1,nx,ny,nz,i,j,k
   
   B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
   nx=B%nx; ny=B%ny; nz=B%nz
   NVAR1=Mesh(nMesh)%NVAR
   if(NVAR1 .eq. 6) then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B)
    do k=0,nz
    do j=0,ny
	do i=0,nx
	 if(B%U(6,i,j,k) .lt. 0.d0) B%U(6,i,j,k)=0.d0
	enddo
	enddo
	enddo
!$OMP END PARALLEL DO
   else if (NVAR1 .eq. 7) then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B)
    do k=0,nz
	do j=0,ny
	do i=0,nx
	 if(B%U(6,i,j,k) .lt. 0.d0) B%U(6,i,j,k)=0.d0
	 if(B%U(7,i,j,k) .lt. 0.d0) B%U(7,i,j,k)=0.d0
    enddo
	enddo
	enddo
!$OMP END PARALLEL DO

  endif
  end	 







!----------------------------------------------------------------------
! 多重网格求解N-S方程 （推进1个时间步）
! nMesh=1,2,3 分别对应用细网格、粗网格、更粗网格
! 包括2重网格和3重网格两个子程序；
! Code by Li Xinliang & Leng Yan
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! 两重网格上推进1个时间步 (3阶RK or 1th Euler)
  subroutine NS_2stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   Type (Block_TYPE),pointer:: B
   integer,parameter:: Time_step_coarse_mesh=3       ! 粗网格迭代步数
!---------------------------------------------------
! -------------------------  网格1 -----------------
   if(Time_Method .eq. Time_Euler1) then
	 call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
   else
	 call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
   endif 
   call  Comput_Residual_one_mesh(1)                 ! 计算网格1的残差 R(n+1)  
   call  interpolation2h(1,2,2)                      ! 把残差插值到网格2 (储存在QF里面)
   call  interpolation2h(1,2,1)                      ! 把守恒变量从网格1插值到网格2   （flag=1 插值守恒变量，=2 插值残差）
!------------------------------
   call  Boundary_condition_onemesh(2)               ! 物理边界条件
   call  update_buffer_onemesh(2)                    ! 内边界条件
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差
   call  comput_force_function(2)                    ! 计算强迫函数QF
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(2)                                  ! 记录初始值  （RK方法中已经包含了该步） 
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(2)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(2)                    ! RK方法推进1步 （网格2）
   endif
   call  comput_delt_U(2)                            ! 计算修正量deltU （储存在Un里面）
   call  prolong_U(2,1,2)                            ! 把修正量插值到细网格 (储存在Un里面); flag=2 插值deltU (储存在Un里)
!------------------------------------	 
   call  comput_new_U(1)                             ! 计算新的U  (U=U+deltU)
   call  Boundary_condition_onemesh(1)               ! 物理边界条件
   call  update_buffer_onemesh(1)                    ! 内边界条件

  end subroutine NS_2stge_multigrid
!------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! 三重网格上迭代1个时间步 （V-型迭代） 3阶RK or 1阶Euler
  subroutine NS_3stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   integer,parameter:: Time_step_coarse_mesh=3       ! 粗网格迭代步数 (对1阶Euler有效)
!---------------------------------------------------
! ---- ---------------------------------- 网格1 -----------------
   if(Time_Method .eq. Time_Euler1) then
	 call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
   else
	 call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
   endif
   call  Comput_Residual_one_mesh(1)                 ! 计算网格1的残差 R(n+1)   ! ????? 该步似乎可以省略 ?????  
! -----------------------------  
   call  interpolation2h(1,2,2)                      ! 把残差插值到网格2 (储存在网格2的QF里面)
   call  interpolation2h(1,2,1)                      ! 把守恒变量从网格1插值到网格2 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）
!-------网格2 --------------------
   call  Boundary_condition_onemesh(2)               ! 物理边界条件
   call  update_buffer_onemesh(2)                    ! 内边界条件
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差         Res_2h(0)
   call  comput_force_function(2)                    ! 计算强迫函数QF （网格2）QF_2h=QF_2h-Res_2h(0) 
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(2)                                  ! 记录初始值  （RK方法中已经包含了该步） 
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(2)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(2)                    ! RK方法推进1步 （网格2）
   endif
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差 R_2h(n+1) 
   call  Add_force_function(2)                       ! 添加上强迫残差(储存在Res里面)  RF_2h(n+1)=R_2h(n+1)+QF_2h  ;  目的：插值到网格3上
   call  interpolation2h(2,3,2)                      ! 把残差插值到网格3 (储存在网格3的QF里面)
   call  interpolation2h(2,3,1)                      ! 把守恒变量从网格2插值到网格3 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）
!------网格3----------------------	  
   call  Boundary_condition_onemesh(3)               ! 边界条件: 物理边界 
   call  update_buffer_onemesh(3)                    ! 内边界
   call  Comput_Residual_one_mesh(3)                 ! 计算网格3的残差
   call  comput_force_function(3)                    ! 计算强迫函数QF （网格3）
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(3)
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(3)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(3)                    ! RK方法推进1步 （网格3）
   endif
   call  comput_delt_U(3)                            ! 计算修正量deltU (=U-Un)
   call  prolong_U(3,2,2)                            ! 把修正量插值到网格2 (储存在deltU里面); flag=2 插值deltU 
!------网格2------------------------      
   call  comput_new_U(2)                             ! 网格2计算新的U  (U=U+deltU)
   call  comput_delt_U(2)                            ! 计算修正量deltU =U-Un
   call  prolong_U(2,1,2)                            ! 把修正量插值到细网格 (储存在deltU里面); flag=2 插值deltU 
!------网格1------------------------------
   call  comput_new_U(1)                             ! 计算新的U  (U=U+deltU)
   call Boundary_condition_onemesh(1)                ! 物理边界条件
   call update_buffer_onemesh(1)                     ! 内边界条件

  end subroutine NS_3stge_multigrid

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
  
!  计算强迫函数 QF=Ih_to_2h Res(n-1) - Res(n)        ! QF中储存着细网格插值过来的残差
  subroutine comput_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%QF(m,i,j,k)=B%QF(m,i,j,k)-B%Res(m,i,j,k)            ! QF原先储存着从细网格插值过来的残差
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo

  end  subroutine comput_force_function

!------------------------------------------------------------
!  把强迫函数添加到残差中 RF=R+QF        
  subroutine Add_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Res(m,i,j,k)=B%Res(m,i,j,k)+B%QF(m,i,j,k)            ! 添加强迫函数后的残差仍储存在B%Res里面 （节省内存）
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo
 
  end  subroutine Add_force_function

!----------------------------------------------------------------------  
!  计算修正量 deltU=U-Un
  subroutine comput_delt_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,5
	         B%deltU(m,i,j,k)=B%U(m,i,j,k)-B%Un(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
 !$OMP END PARALLEL DO 
  enddo
  end  subroutine comput_delt_U
!-----------------------------------------------------------------------
! 设定Un=U
  subroutine Set_Un(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Un(m,i,j,k)=B%U(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo
  
  end  subroutine Set_Un
!-------------------------------修正U --------------------------------------
  subroutine comput_new_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,5                                             ! 第6个量是湍流粘性系数vt (SA模型使用),不需要修正
	         B%U(m,i,j,k)=B%U(m,i,j,k)+B%deltU(m,i,j,k)         ! Un里面储存的是U的修正量 （从粗网格插值而来）
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo

  end  subroutine comput_new_U
!---------------------------------------------------------------------------
! 粗网格向细网格的插值(Prolong) 及 细网格向粗网格上插值 (interpolation)
!----------------------------------------------------------------------
! 将网格m1的守恒变量(U) 或U的差插值到网格m2 (上一级细网格)
! flag=1时，将U插值到上一级网格；  (准备初值时使用)
! flag=2时，将deltU插值到上一级网格 （deltU储存着本时间步与上个时间步U的差） 

  Subroutine prolong_U(m1,m2,flag)
   use Global_Var
   use interface_defines
   implicit none
   integer:: m1,m2,mb,flag
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   if(m1 .le. 1 .or. m1-m2 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2)
     do mb=1,MP1%Num_Block
       B1=>MP1%Block(mb)
	   B2=>Mp2%Block(mb)
	   if(flag .eq. 1) then
!	   call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%U(1,-1,-1,-1),B2%U(1,-1,-1,-1))   ! 旧的程序接口，与新版Fortran不兼容
	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%U,B2%U)
  	   else
!	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%deltU(1,-1,-1,-1),B2%deltU(1,-1,-1,-1))
	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%deltU,B2%deltU)
     endif
   enddo

  end Subroutine prolong_U
!-------------------------------------------------------------
!   粗网格向细网格上的插值     
!   U1是粗网格上的值； U2是细网格上的值
  subroutine prolongation(nx1,ny1,nz1,nx2,ny2,nz2,U1,U2)
   use precision_EC
   implicit none
    integer:: i,j,k,m,nx1,ny1,nz1,nx2,ny2,nz2,NV
    real(PRE_EC),dimension(:,:,:,:),pointer:: U1,U2
 !   real(PRE_EC):: U1(NVAR,-1:nx1+1,-1:ny1+1,-1:nz1+1),U2(NVAR,-1:nx2+1,-1:ny2+1,-1:nz2+1)
    integer:: ia(2,0:nx2),ja(2,0:ny2),ka(2,0:nz2),U_bound(4)
 !   integer,parameter::NVAR=6
    real(PRE_EC),parameter:: a1=27.d0/64.d0,a2=9.d0/64.d0,a3=3.d0/64.d0,a4=1.d0/64.d0   ! 插值系数
    
!  寻找插值基架点的下标 
!  ia(1,i) 是距离i点最近的粗网格点的下标；ia(2,i)是次近点的下标	 
    U_bound=UBOUND(U1)   ! 第1维的上界 （NVAR)
    NV=U_bound(1)   ! NVAR= 5, 6 or 7 
    
   do i=0,nx2
	 if(mod(i,2).eq.0) then
	   ia(1,i)=i/2                    !最近点
	   ia(2,i)=i/2+1                  !次近点
	 else  
	   ia(1,i)=i/2+1                  !最近点
	   ia(2,i)=i/2                    !次近点
	 endif
   enddo
   do j=0,ny2
	 if( mod(j,2).eq. 0) then
	   ja(1,j)=j/2
	   ja(2,j)=j/2+1
	 else
	   ja(1,j)=j/2+1
	   ja(2,j)=j/2
	 endif
   enddo
   do k=0,nz2
	 if(mod(k,2).eq.0) then
	   ka(1,k)=k/2                    !最近点
	   ka(2,k)=k/2+1                  !次近点
	 else  
	   ka(1,k)=k/2+1                  !最近点
	   ka(2,k)=k/2                    !次近点
	 endif
   enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,m)

   do k=0,nz2
	 do j=0,ny2
	   do i=0,nx2
	     do m=1,NV
!               插值，最近点的权重a1, 次近点的权重a2, 最远点的权重a3	 
	       U2(m,i,j,k)=a1*U1(m,ia(1,i),ja(1,j),ka(1,k))+a2*(U1(m,ia(2,i),ja(1,j),ka(1,k))+U1(m,ia(1,i),ja(2,j),ka(1,k)) &
	                   +U1(m,ia(1,i),ja(1,j),ka(2,k)))+a3*(U1(m,ia(2,i),ja(1,j),ka(2,k))+U1(m,ia(1,i),ja(2,j),ka(2,k)) &
	                   +U1(m,ia(2,i),ja(2,j),ka(1,k)))+a4*U1(m,ia(2,i),ja(2,j),ka(2,k))
         enddo
	   enddo
     enddo
   enddo
!$OMP END PARALLEL DO 
  end subroutine prolongation
!---------------------------------------------------------------------------------
! 将网格m1的守恒变量U插值到网格m2 (细网格->粗网格) 
  Subroutine interpolation2h(m1,m2,flag)
   use Global_Var
   implicit none
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   real(PRE_EC),dimension(:,:,:,:),pointer:: P1,P2
   integer:: NVAR1,flag,m1,m2,mb,i,j,k,m,i1,i2,j1,j2,k1,k2
!  flag==1 插值守恒变量； flag==2 插值残差
   if( m2-m1 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2) 
	 NVAR1=5              ! 只插值5个守恒变量  
     do mb=1,MP1%Num_Block
       B1=>MP1%Block(mb)
  	   B2=>Mp2%Block(mb)
       if(flag .eq. 1) then  ! 插值守恒变量
	     P1=>B1%U
	     P2=>B2%U
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B1,B2,NVAR1,P1,P2)
	     do k=1,B2%nz-1
           do j=1,B2%ny-1
	         do i=1,B2%nx-1
	           i1=2*i-1 ; i2=2*i
	           j1=2*j-1 ; j2=2*j
	           k1=2*k-1 ; k2=2*k
               do m=1,NVAR1
!     以控制体体积为权重的加权平均 	 
	             P2(m,i,j,k)=(P1(m,i1,j1,k1)*B1%Vol(i1,j1,k1)+P1(m,i1,j2,k1)*B1%Vol(i1,j2,k1)   &
		                      +P1(m,i1,j2,k2)*B1%Vol(i1,j2,k2)+P1(m,i2,j1,k1)*B1%Vol(i2,j1,k1)   &
	                          +P1(m,i2,j1,k2)*B1%Vol(i2,j1,k2)+P1(m,i2,j2,k1)*B1%Vol(i2,j2,k1)   &
					          +P1(m,i2,j2,k2)*B1%Vol(i2,j2,k2)+P1(m,i1,j1,k2)*B1%Vol(i1,j1,k2))/B2%Vol(i,j,k)
  	           enddo
	         enddo
	       enddo
	     enddo
!$OMP END PARALLEL DO 

       else     ! 插值残差  （把m1网格上的残差B%Res 插值到m2网格上B%QF (然后减去本m2网格上的残差，形成强迫函数)）
   	     P1=>B1%Res
	     P2=>B2%QF

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B2,P1,P2,NVAR1)
	     do k=0,B2%nz
           do j=0,B2%ny
	         do i=0,B2%nx
	           i1=2*i-1 ; i2=2*i
	           j1=2*j-1 ; j2=2*j
	           k1=2*k-1 ; k2=2*k
               do m=1,NVAR1
	             P2(m,i,j,k)=P1(m,i1,j1,k1)+P1(m,i1,j2,k1)+P1(m,i1,j2,k2)+P1(m,i2,j1,k1)   &    ! 残差的插值： 简单相加
		                     +P1(m,i2,j1,k2)+P1(m,i2,j2,k1)+P1(m,i2,j2,k2)+P1(m,i1,j1,k2)
  	           enddo
	         enddo
	       enddo
	     enddo
 !$OMP END PARALLEL DO 

 	   endif
     enddo

  end Subroutine interpolation2h
