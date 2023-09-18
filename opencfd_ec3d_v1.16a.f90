!--------------------------OpenCFD-EC 3D--------------------------------------------  
!----OpenCFD-EC 3D: A multi-block Finite Volume  Navier-Stokes solver---------------
!  Copyright by Li Xinliang, LHD, Institute of Mechanics, lixl@imech.ac.cn  
!  Code developed by Li Xinliang and Leng Yan               
!  Ref. J. Blazek's Book: Computational Fluid Dynamics: principle and application 
!----------------------------------------------------------------------------------
! Ver 0.2:   2010-10 the first version
! Ver 0.21:  Jocabian coefficients are restored
! Ver 0.22:  All Goemetry variables are restored
! Ver 0.3 :  BL model ; Goemetry variables are NOT restored to save memory
! Ver 0.31:  3rd Runge-Kutta and 1st order Euler time-step method
! Ver 0.32:  Minor revised in find negative temperature 
! Ver 0.40:  Aritifical Viscous method
! Ver 0.41:  Modify the defination of ori (in bc.in) to suitable to BXCFD 
! Ver 0.42:  Residual Smoothing method
! Ver 0.43:  Append Van Leer, AUSM (by Leng Yan) 2010-11-18
! Ver 0.43a: Two bugs are corrected (by LengYan)
! Ver 0.43b: New method for Corner Point  2010-11-24
! Ver 0.44:  New analysis code for computing the aerodynamic force and momentium
! Ver 0.45:  3rd order MUSCL scheme  
! Ver 0.46:  Boundary conditions are revised
! Ver 0.46a: 2010-12-7   A Bug in get_ijk_orient() is removed
! Ver 0.46b: 2010-12-9   A Bug in comput_force() is removed
! Ver 0.50:  2010-12-9   Two layers of Ghost Cells are used; AUSM-PW ;  A bug in updata_buffer (line 98) is removed
! Ver 0.51:  2010-12-17  SA model is supported
! Ver 0.51a: 2010-12-22  Artificial viscous method supportting is removed; some bug in SA model is removed
! Ver 0.52: 2010-12-23   support OpenMP 
! Ver 0.52a: 2011-3-21:  A bug in Flux_Roe_1D() is removed
! Ver 0.52b: 2011-4-13:  A bug in  subroutine boundary_Farfield() is removed
! Ver 0.6: 2011-4-22:    Support Multi-grid  (modified by Leng Yan)
! Ver 0.6a: 2011-5-5:    A bug in Roe flux is removed
! Ver 0.6b: 2011-6-1:    A bug in the initial of multi-grid is removed (read data )
! Ver 0.6c: 2011-6-20:   A bug in read Mesh3d.dat is removed (unformatted data)
! Ver 0.6d: 2011-10-8:   A bug in AUSM-PW() is removed
! Ver 0.7:  2011-12-14:  Can read Gridgen (CFL3D) .inp file (.inp file can be tranformed into .in file );
!                        Isothermal wall boundary condition is supported
!                        auto-detect distance-to-the-wall (for SA model) 
! Ver 0.71: 2011-12-20:  Support OpenMP 
! Ver 0.72: 2011-12-25:  More Debug function is adding (such as detect the location of "NaN")
! Ver 0.72a: 2011-12-26: Can be initial with zero flow (u=v=w=0)
! Ver 0.72b: 2011-12-26: Some incompatible code (such as  using  pointers) is modified
! Ver 0.73:  2011-12-27: Some bugs (such as in BL model) is removed
! Ver 0.74:  2011-12-30: LU-SGS is supported
! Ver 0.75:  2012-1-3:   Robustness is enhanced by initial smoothing (filtering) and smoothing (filtering) during every Kstep_smooth steps 
! Ver 0.75b: 2012-1-4:   Bugs in LU-SGS are fixed; LU-SGS support viscous flow
! Ver 0.75c: 2012-1-11:  Turbulent viscous is limited; Forces are reported by each block
! Ver 0.75d: 2012-1-16:  New post-analysis code
! Ver 0.75e: 2012-1-19:  Modify in Open-MP feature, Limit for vt/v increased to 300
! Ver 0.76: 2012-1-31:   New feature of UIRS (Upwind Implicit Residual Smoothing); 
! Ver 0.76a: 2012-2-1:   Bug in output( ) (flow2d-surf.dat) is fixed;  iso_thermal_wall( ) is modified
! Ver 0.76b: 2012-2-3:   Viscous terms in du_LU_SGS is modified, DU(1) (density) need not viscous term.
! Ver 0.76c: 2012-2-4:   Improve Robustness by lower local time-step when Negative Temperature points appear 
! Ver 0.76d: 2012-2-6:   Support outflow boundary condition 
! Ver 0.76e: 2012-2-29:  Viscous terms in du_LU_SGS is modified, DU(1) (density) use viscous term; Local average d and T when them < Limit 
! Ver 0.77:  2012-3-2:   New feature of Data file 
! Ver 0.80:  2012-5-5:   New feature of Turbulence models,Support K-W SST model;
! Ver 0.81:  2012-5-8:   New feature for Open-MP
! Ver 0.81a: 2012-5-13:  Code in SA model is modified
! Ver 0.81b: 2012-5-14:  Code in SA model is modeified, Sorce term is limited (see : http://turbmodels.larc.nasa.gov/spalart.html )
! Ver 0.81c: 2012-5-16:  Code in SST is modified
! Ver 0.82:  2012-5-17:  User can set OpenMP threads
! Ver 0.82a: 2012-5-21:  Codes in boundary condition are modified, outflow boundary 
! Ver 0.83:  2012-6-1:   Check Mesh quality before computing, Lower local dt at poor quality mesh
! Ver 0.84:  2012-6-26:  Built-in support Gridgen (.inp) format boundary file  (origin version built-in support BXCFD .in format file)
! Ver 0.90:  2012-8-30:  MPI version
! Ver 0.91:  2012-9-14:  MPI_Bsend -> MPI_send
! Ver 0.92:  2012-9-28:  Support Distributed input file (such as Mesh3d.dat.00000 ...)
! Ver 0.93:  2012-10-9:  Support Dual-time-step LU-SGS method; a bug in read_flow() is removed
! Ver 0.94:  2012-10-31: 5th-order WENO scheme supported; boundary ghost Cell (LAP)=3
! Ver 0.95:  2012-11-7:  automatic comput wall-distance  (former version read data file) ;  automatic partation for MPI    
! Ver 0.95a: 2012-11-30:  A bug in allocate_mem_blocks is removed (LAP)
! Ver 0.95b: 2012-12-23:  A bug with LAP for d, uu,v,w,T,p,cc  is removed
! Ver 0.95c: 2012-12-25:  code in Far-field boundary condition is modified (more robust)
! Ver 0.96: 2013-4-26: Extrapolate boundary is supported;  
!                      Support CGNS format mesh by using convert tool convert-cgns-inp-1.0.f90
!                      Comput_Goemetric_var() is modified to support degenerate line (plane with 0 aera)
! Ver 0.97: 2013-5-3:  New feature in viscous term (new Jocabian for viscous terms)   
! Ver 0.97a: 2013-5-4: new Jocabian in the corner points, New directive in viscous term for the corner points;
! Ver 0.98: 2013-6-3: Bugs in SST k-w model are removed (nondimensional)
! Ver 0.98b: 2013-8-19: mu_t is limited
! Ver 0.98c: 2013-9-14: code in SA model is modified (nondimensional)
! Ver 0.99: 2013-10-6: Code in Filtering is modified
! Ver 1.00: 2013-11-10: Support user-defined boundary
! Ver 1.01: 2013-11-13: Support boundary scheme ; 1.01a :  about read wall_dist.dat
! Ver 1.02: 2013-11-25: WENO7 scheme is supported
! Ver 1.03: 2013-12-24: Support NAMELIST type input file (control.ec)
! Ver 1.04: 2014-1-3:  New feature in output (error.log) /input 
! Ver 1.05: 2014-3-14: code for Numerical scheme are modified; WENO3 scheme is deleted; MUSCL with k=-1 (upwind MUSCL) is added;  OMUSCL2 is added
! Ver 1.05a: 2014-3-25: Limit for pressure increase
! Ver 1.05b: 2014-3-26: New treatment of Corner point
! Ver 1.05c: 2014-3-31: A bug in comput_duvwpc() (limit d, T) is removed
! Ver 1.06: 2014-4-3: Limit flow (d,p,u)
! Ver 1.06a: 2014-4-5: Limit SA 
! Ver 1.07a: 2014-4-20: New SA model
! Ver 1.1: 2015-2-4: Turbo-Machinary Solver
! Ver 1.1a: 2015-3-11: A bug for the nondimensionlize of pressure is fixed
! Ver 1.11: 2015-4-12: A bug in computing CL, CD is fixed (for 	Cood_Y_UP ==0 )
! Ver 1.12: 2015-5-4:  Support hybrid Finite-Valume/Finite-Difference method
! Ver 1.13: 2015-5-12: Support Inner-flow mode
! Ver 1.13b: 2015-8-31: A bug in Roe split is removed
! Ver 1.13c: 2015-11-29: A OpenMP bug in Runge-Kutta is removed; 
!            Time-serial inlet boundary is supported;
!            Time average is supported; 
! Ver 1.14: 2015-12-29:  output flile flow3d-xxxxxx.dat
! Ver 1.14a: 2016-1-5:  Wall boundary with blow-and-suction perturbation 
! Ver 1.15: 2017-3-13: A Bug in comput_force( ) is removed (in Computing Viscous force)
! Ver 1.15a: 2017-5-11: Bug in output_vt removed;  boundary condition for BC_Inlet (not inner flow) is regarded as fixed boundary condition;
! Ver 1.16: 2017-7-11:  IF_OverLimit is used,  IF_OverLimit=1 when overflow in the block, and then UD1 scheme is used for inviscous flow;  
！Ver 1.16a: 2018-9-29:  A bug in SST model is removed;
!-------------------------------------------------------------------------------------------------------------------------------------------------

! 流场物理量 （计算每块时申请内存，该块计算结束后释放；属于临时变量） 
!  include "sub_modules.f90"

  module Flow_Var
   use precision_EC
   real(PRE_EC), save,pointer,dimension(:,:,:)::  d,uu,v,w,T,p,cc ! 密度、x-速度、y-速度、z-速度、压力、声速
   real(PRE_EC), save,pointer,dimension(:,:,:,:):: Flux                 ! i- ,j-及k-方向的通量
   real(PRE_EC), save,pointer,dimension(:,:,:):: Lvi,Lvj,Lvk,Lci,Lcj,Lck   ! 无粘项及粘性项Jocabian的谱半径 

  end module Flow_Var

!------------------------------------------------------------------------------------------
! 主程序 主程序 主程序
!-----------------------------------------------------------------------------------------
  program main
   use Global_Var
   implicit none
   integer:: ierr

   call Init_mpi
  
   if(my_id .eq. 0) then
    print*,  "----------------- OpenCFD-EC3D ver 1.14a (MPI-OpenMP version)------------------"
    print*,  "        Copyright by Li Xinliang, lixl@imech.ac.cn                            "
    print*,  "        Programming by Li Xinliang  2016-1                                   "
    print*,  "----------------------------------------------------------------------------- " 
   endif

   call read_parameter                     ! 读取流动参数及控制信息
!$ call omp_set_num_threads(NUM_THREADS)   ! 设置OpenMP的运行线程数 （并行数目）， 本语句对openmp编译器不是注释!

!$ if(my_id ==0) then         ! 测试一下运行的进程 （openmp编译时，不是注释）
!$OMP Parallel
!$  print*, "omp run ..."
!$OMP END parallel
!$ endif 

   allocate( Mesh(Num_Mesh) )                                 ! 主数据结构： “网格” （其成员是“网格块”）
   
   if(my_id .eq. 0)  call check_mesh_multigrid               ! 检查网格配置所允许的最大重数,并设定多重网格的重数
   
   call Init                               ! 初始化变量（分配内存，读取网格）
   call set_control_para                   ! 设定各重网格上的控制信息（数值方法、通量技术、湍流模型、时间推进方式）
   call check_mesh_quality                 ! 检查网格质量,在网格质量差的区域降低局部时间步长
   call Init_flow                          ! 初始化流场 （初值）
   

   if(my_id .eq. 0) print*, " Start ......"

!------------------------------------------------------------------------
! 时间推进，采用单重网格、二重网格或三重网格； 采用1阶Euler或3阶RK
   do while(Mesh(1)%tt .lt. t_end )
     call show_Wall_time()

     if(Num_Mesh .eq. 1)  then                          ! 单重网格推进1个时间步
       call NS_Time_advance(1)
	 else  if(Num_Mesh .eq. 2)  then                    ! 2重网格推进1个时间步
  	   call NS_2stge_multigrid
     else                                               ! 3重网格推进1个时间步
  	   call NS_3stge_multigrid 
     endif	  
 
 !  滤波 ,可以增强稳定性. 如Kstep_Filter=0则不使用滤波   
	if(Kstep_smooth .gt. 0) then
 	  if(mod(Mesh(1)%Kstep, Kstep_smooth).eq.0)   call Filtering_oneMesh(1)                      ! 滤波            
    endif


!  每隔一定步数输出气动力及残差（输出到屏幕及文件: force.log, Residual.dat）
     if(mod(Mesh(1)%Kstep, Kstep_show).eq.0) then
      call comput_force
      call output_Res(1)
     endif
! 每隔一定步数输出数据文件(flow3d.dat, PLOT3D 格式)
      if(mod(Mesh(1)%Kstep, Kstep_Save).eq.0) then
	     call output_flow 
!         if(If_debug == 1 ) call output_vt                ! 输出湍流粘性系数，供debug使用   ! Bug 2017-5-11
          if(If_debug == 1 .and.  If_viscous==1 .and.  Iflag_turbulence_model .ne. 0) call output_vt                ! 输出湍流粘性系数，供debug使用
	  endif
 
 ! 进行时间平均
    if(Kstep_average > 0) then     
      if(mod(Mesh(1)%Kstep, Kstep_average) .eq.0) then
         call Time_average           ! 时间平均
	  endif
      if(mod(Mesh(1)%Kstep, Kstep_Save).eq.0) then
	      call output_flow_average   ! 输出时均场,PLOT3D格式
	  endif    
    endif 
   
   enddo

  end
!----------------------------------------------------------------------------------------------
  subroutine Init_mpi
   use Global_var
   implicit none
      integer,parameter:: IBuffer_Size=10000000    
!      real*8,allocatable,dimension(:)::  Buffer_mpi    ! Buffer for MPI  message transfer (used by MPI_Bsend)
       real*8:: Buffer_mpi(IBuffer_Size)
      integer   ierr, status(MPI_status_size)

!------------------------------------------------
       call mpi_init(ierr)                                     ! 初始化MPI
       call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)           ! 获取本进程编号
       call mpi_comm_size(MPI_COMM_WORLD,Total_proc,ierr)      
!       allocate(Buffer_mpi(IBuffer_Size))
	   call MPI_BUFFER_ATTACH(Buffer_mpi,8*IBuffer_Size,ierr)   ! 创建消息发送缓冲区，供MPI_Bsend()使用
   end subroutine Init_mpi


!  显示（墙钟）时间，用于统计MPI并行效率	  
    subroutine show_Wall_time()
      use Global_var
	  real*8:: wtime
	  real*8,save:: wtime0,wtime1    ! 初始时间，上一步的时间
	  integer,save:: KP=0  ! 计算步
      if(my_id .eq. 0) then
	    wtime=MPI_Wtime()
        if(KP .eq. 0) then   
		  wtime0=wtime   ! 初始CPU时间
		else
          if(mod(Mesh(1)%Kstep, Kstep_Show).eq.0) then
		  print*, "CPU wall time in this step:", wtime-wtime1 
          print*, "Averaged CPU wall time is:", (wtime-wtime0)/KP 
          endif
        endif
		 wtime1=wtime  
         KP=KP+1    ! 统计计算步
      endif

    end subroutine show_wall_time
       

!----------------------------------------------------------------------------
!  include "sub_interfaces.f90"
!  include "sub_init.f90"
!  include "sub_Residual.f90"
!  include "sub_flux_split.f90"
!  include "sub_boundary.f90"
!  include "sub_turbulence_BL.f90"
!  include "sub_turbulence_SA.f90"
!  include "sub_turbulence_SST.f90" 
!  include "sub_time_advance.f90"
!  include "sub_updata_buffer_mpi.f90"
!  include "sub_Post.f90"
!  include "sub_time_acceleraction.f90"
!  include "sub_LU_SGS.f90"
!  include "sub_IO.f90"
!  include "sub_geometry.f90"
!  include "sub_partation_mpi.f90"
!  include "sub_convert_inp.f90"
!  include "sub_read_parameter.f90"
