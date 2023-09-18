! 全局变量及类定义
!-----------------------------------------------------------------------------------
! Consts  
  module precision_EC
  implicit none
  include "mpif.h"
!     integer,parameter:: PRE_EC=4            ! Single precision
     integer,parameter:: PRE_EC=8           ! Double Precision
!    integer,parameter:: OCFD_DATA_TYPE=MPI_REAL  
	integer,parameter:: OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION ! Double precision
  end module  precision_EC 

!------------------------------------------------------------------------------------

! 常数
  module const_var
   use precision_EC
   implicit none
   real(PRE_EC),parameter::  PI=3.1415926535897932d0, Lim_Zero=1.d-20   ! 小于该值认为0
   integer,parameter:: LAP=4                          ! 虚网格的数目 （使用3阶格式该值不小于2；如使用5阶WENO, 该值不小于3; 如果使用WENO7, 则该值不小于4）
   integer,parameter:: Scheme_UD1=0, Scheme_NND2=1, Scheme_UD3=2,Scheme_MUSCL2U=3,Scheme_MUSCL2C=4,   &
                       Scheme_MUSCL3=5,Scheme_OMUSCL2=6,Scheme_WENO5=7,Scheme_UD5=8, Scheme_WENO7=9
   integer,parameter:: Scheme_CD2=20, Scheme_none=-1         ! 不使用（边界）格式
   integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_Van_Leer=5,Flux_Ausm=6
   integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2


!---------------边界条件-------------------------
!  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Periodic=-30,BC_Symmetry=-40,BC_Outlet=-22
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Wall_Turbo=201             ! 叶轮机械中机匣的边界条件 （绝对速度为0， 相对坐标系有旋转） 
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401      ! (扩展) 与Griggen .inp文件的定义可能有所区别，请注意
!                     周期性边界条件并非设置成BC_Peridodic=501, 而是设置BC_PeriodicL=-2, BC_PeriodicR=-3

!   周期边界按照内边界处理， 但几何变量须专门处理。 叶轮机模式下，周期边界也须特殊处理   

!  -1 内边界 （非物理边界）；  -2  左周期边界； -3 右周期边界 （区分左、右便于处理几何坐标）
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   
   integer,parameter::  BC_Zero=0       ! 无边界条件

!                 用户自定义的边界条件，要求代码 >=900  
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902       !给定入口流动； 给定入口时间序列
   integer,parameter:: BC_USER_Blow_Suction_Wall=903    ! 吹吸扰动壁面

   integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
   integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2, Turbulence_SST=3,Turbulence_NewSA=21
   integer,parameter:: Init_continue=1, Init_By_FreeStream=0, Init_By_Zeroflow=-1,  Smooth_2nd=0,Smooth_4th=1 
!   real(PRE_EC), parameter::  Density_LIMIT=1.d-4,Temperature_LIMIT=1.d-4,Pressure_LIMIT=1.d-4

   integer,parameter:: Method_FVM=0, Method_FDM=1        ! 差分、有限体积
   integer,parameter:: FD_WENO5=1,FD_WENO7=2,FD_OMP6=3  ! 差分法采用的数值格式
   integer,parameter:: FD_Steger_Warming=1,FD_Van_Leer=2


  end module const_var

! 定义类 (边界连接， 网格块， "网格")
 module Mod_Type_Def
   use precision_EC
   implicit none

    TYPE BC_MSG_TYPE              ! 边界链接信息
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                     ! 子面号，连接顺序描述符
   END TYPE BC_MSG_TYPE

!------------------------------------网格块--------------------------------------
   TYPE Block_TYPE                                 ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
     integer::  Block_no,mpi_id           ! 块号；所属的进程号
	 integer::  nx,ny,nz                   ! 网格数nx,ny,nz
	 integer::  subface                   ! 子面数
 !   几何量  
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z     ! coordinates of vortex, 网格节点坐标
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc  ! coordinates of cell center, 网格中心坐标 
     real(PRE_EC),pointer,dimension(:,:,:):: Vol,Si,Sj,Sk ! Volume and surface area, 控制体的体积，i,j,k方向控制体边界面的面积 
	 real(PRE_EC),pointer,dimension(:,:,:):: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3  !  i,j,k方向三个控制面的法方向
!      Jocabian 变换系数	 
	 real(PRE_EC),pointer,dimension(:,:,:):: ix1,iy1,iz1,jx1,jy1,jz1,kx1,ky1,kz1
	 real(PRE_EC),pointer,dimension(:,:,:):: ix2,iy2,iz2,jx2,jy2,jz2,kx2,ky2,kz2
	 real(PRE_EC),pointer,dimension(:,:,:):: ix3,iy3,iz3,jx3,jy3,jz3,kx3,ky3,kz3
	 real(PRE_EC),pointer,dimension(:,:,:):: ix0,iy0,iz0,jx0,jy0,jz0,kx0,ky0,kz0
!     物理量
	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U,Un,Un1    ! 守恒变量 (本时间步及前一、二个时间步的值), conversation variables 
     real(PRE_EC),pointer,dimension(:,:,:,:) :: Res         ! 残差 （净通量）
     real(PRE_EC),pointer,dimension(:,:,:):: dt          ! (局部)时间步长
     real(PRE_EC),pointer,dimension(:,:,:,:) :: QF      ! 强迫函数 (多重网格法中粗网格使用)
     real(PRE_EC),pointer,dimension(:,:,:,:) :: deltU   ! 守恒变量的差值, dU=U(n+1)-U(n)  多重网格使用
     real(PRE_EC),pointer,dimension(:,:,:,:) :: DU      ! U(n+1)-U(n)    LU-SGS中使用
	 real(PRE_EC),pointer,dimension(:,:,:):: dw         ! turbulent viscous ; distance to the wall  (used in SA model)
	 real(PRE_EC),pointer,dimension(:,:,:):: surf1,surf2,surf3,surf4,surf5,surf6  ! 边界处的通量，供计算气动力、热使用；
	 real(PRE_EC),pointer,dimension(:,:,:):: mu,mu_t    ! 层流粘性系数和湍流粘性系数
     real(PRE_EC),pointer,dimension(:,:,:):: dtime_mesh   ! 时间步长因子 （根据网格质量情况）

	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U_average ! 时间平均场 （d,u,v,w,T）
	 
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
     integer,pointer,dimension(:,:,:):: BcI,BcJ,BcK     ! 边界指示符 （物理边界 or 内边界）
   	 
	 integer:: IFLAG_FVM_FDM               ! 差分 or 有限体积
	 integer:: IF_OverLimit                          ! 物理量超限（如负温度）， 需要降低精度（1阶）
	End TYPE Block_TYPE  

!---------------------------网格 -------------------------------------------------------- 
!  (如单重网格，只有1套；如多重网格，可以有多套) 
  
   TYPE Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep          ! 网格编号 (1号为最细网格，2号为粗网格， 3号为更粗网格...)，网格块数，网格数目(本进程中), 时间步 
     integer:: NVAR       ! 变量的数目；基本变量5个+ 0，1或2个附件变量 （BL模型0个，SA模型1个，SST模型2个）； 粗网格不使用湍流模型
	 real(PRE_EC)::  tt                   !  推进的时间
	 real(PRE_EC),pointer,dimension(:)::  Res_max,Res_rms                 ! 最大残差，均方根残差, 推进的时间
	 TYPE (Block_TYPE),pointer,dimension(:):: Block       ! “网格块”  （从属于“网格”）

!                                                       控制参数，用于控制数值方法、通量技术、湍流模型等    
!             这些控制参数从属于“网格”，不同“网格”可以采用不同的计算方法、湍流模型等。	 （例如，粗网格用低精度方法，粗网格不使用湍流模型,...）
!  If_dtime_mesh     是否根据网格情况降低局部时间步长（局部时间步长法有效）	
    integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction, Bound_Scheme
   End TYPE Mesh_TYPE
  
  end module Mod_Type_Def


!-------------------------------------------------------------------------------------- 
! Global Variables:
! Ma: Mach number ; Re: Reynolds number; gamma: Specific rato (=Cp/Cv); Pr: Prandtl number; 
! AoA: Angle of Attack; p00=1/(gamma*Ma*Ma), p=p00*d*T; 
! t_end: end time; 
! Num_Block: Total Block number 
! Kstep_save: Save data every Kstep_save step
! Iflag_turbulence_model: 0 no model, 1 BL model
! Iflag_flux: type of Splitting (or Riemann solver) : 0 Steger-Warming 1 HLL 2 HLLC 3 Roe
! Iflag_local_dt:  0 global time step, 1 local time step
! Iflag_Reconstruction: Schemes 1 NND 2 3rd Upwind 3 WENO3  4 MUSCL 
!----------------------------------------------------------------------------------------
! 全局变量，包括参数、几何量及物理量
!========================================================================================

  module Global_Var    
   use const_var        ! 常量
   use mod_type_def     ! 边界连接
   implicit none


!---------------------------------------------------------------------------------------------
! global variables                                       各子程序均可见的全局变量
!----------------------------------------------------------------------------

   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh                          ! 主数据 “网格”
   integer,save:: Num_Mesh,NVAR , Total_block, Num_block                      ! 网格的套数 ，变量数， 总网格块数, 本mpi进程的网格块数  
   integer,pointer,dimension(:):: bNi,bNj,bNk                            ! 各块的维数（全局）， bNi(k)为（全局）第k块的nx
   integer,save::  Kstep_save, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method, &
	  Kstep_show,If_viscous,If_Residual_smoothing,Mesh_File_Format,IF_Debug, &
	  Kstep_smooth,Kstep_init_smooth,NUM_THREADS,If_dtime_mesh, Step_Inner_Limit, &  
      Bound_Scheme, &                         ! 边界格式
      IF_Walldist, IFLAG_LIMIT_FLOW, &             ! 是否需要读取 到壁面距离; 是否需要限制压力增长率
      IF_Scheme_Positivity,          &    ! 检查插值过程中压力、密度是否非负，否则使用1阶迎风；
      Kstep_average,                 &          ! 时均统计的步数间隔， 0 不统计
      Iflag_savefile                       ! 0 保存到flow3d.dat, 1 保存到flow3d-xxxxxxx.dat
   integer,save:: IF_TurboMachinary , Ref_medium_usrdef   !  是否启用叶轮机械求解器模式（0 or 1）； 是否使用用户自定义介质 （默认空气介质）
   integer,save:: IF_InnerFlow   ! 内流模式 (入口边界条件与叶轮机模式类似)

   integer,save:: FD_Flux,FD_scheme   ! 内嵌差分法采用的通量方式、数值格式
   integer,save:: KRK=0            ! Runge-Kutta方法中的 子步
   integer,save:: Istep_average=0  ! 统计步数
 ! global parameter (for all Meshes )                     流动参数, 对全体“网格”都适用

   real(PRE_EC),save:: Ma,Re,gamma,Cp,Cv,t_end,P_OUTLET,&
                       A_alfa,A_beta,PrL,PrT,T_inf,Twall,w_LU,Kt_inf,Wt_inf , Res_Inner_Limit, MUT_MAX,AoA,Aos
   real(PRE_EC),save:: Turbo_Periodic_seta, Turbo_w, Turbo_P0, Turbo_T0 , Turbo_L0   ! 周向跨度（角）， 转速，总压，总温, 参考长度
   real(PRE_EC),save:: Periodic_dX,Periodic_dY,Periodic_dZ

 
 ! 全局控制参数，控制数值方法、通量技术及湍流模型等 （有些只对最细网格有效）
   real(PRE_EC),save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin    ! RK方法中的常量，与时间步长有关的量
   integer,save:: Pre_Step_Mesh(3)                                                 ! 构建初值时，粗网格预迭代步数  
   integer,save:: Cood_Y_UP             ! 1 Y轴垂直向上， 0 Z轴垂直向上
   integer,save:: Pdebug(4)             ! debug 使用， 输出某块的某一点的值
   real(PRE_EC),save:: Ref_S,Ref_L, Centroid(3)                           ! （计算六分量使用）参考面积, 参考长度, 矩心坐标
   real(PRE_EC),save:: Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax                ! 对密度、压力、速度、及SA模型中变量的限制条件 (最小与最大值)
   real(PRE_EC),save:: CP1_NSA,CP2_NSA       ! parameters in New SA model
 !-----------mpi data ----------------------------------------------------------- 
   integer:: my_id,Total_proc                   ! my_id (本进程号), Total_proc 总进程数
   integer,pointer,dimension(:):: B_Proc, B_n    ! B_proc(m) m块所在的进程号; B_n(m) m块所在进程中的内部编号
   integer,pointer,dimension(:):: my_Blocks       ! 本进程包含的块号
  end module Global_Var  
!----------------------------------------------------------------------------



! 差分-有限体积混合方法使用
!------------SEC part---------------------------------------------------------------------

   module FDM_data
    use precision_EC
    implicit none
    TYPE FDM_Block_TYPE                              ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
    real(PRE_EC), pointer,dimension(:,:,:):: ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac   ! Jocabian变换系数
    End TYPE FDM_Block_TYPE  
  
   TYPE FDM_Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
	 TYPE (FDM_Block_TYPE),pointer,dimension(:):: Block       ! “网格块”  （从属于“网格”）
   End TYPE FDM_Mesh_TYPE 
  
   TYPE (FDM_Mesh_TYPE),pointer,dimension(:):: FDM_Mesh       ! 主数据 “网格”
  
   end  module FDM_data
