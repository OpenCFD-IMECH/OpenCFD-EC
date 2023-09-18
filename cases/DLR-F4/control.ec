$control_ec
  Iflag_init=0  
  Ma=0.75d0
  Re=2.125d4
  AoA=0.d0              ! Angle of Attack
  If_viscous=1
  Iflag_turbulence_model=2    ! SA   
  Kstep_save=5000
  t_end=100.d0     
  CFL=5.d0           
  Time_Method=0              ! LU-SGS 
  T_inf=288.15d0        
  Twall=-1.d0                !  Adabitic wall
  Iflag_Scheme= 5            !  5 MUSCL3 
  Iflag_Flux=5               !  5 Van_Leer
  IFlag_Reconstruction=0     !  0 Original
  Ref_S=0.0728d6                ! Ref. area  (1/2 of S0)
  Ref_L=141.2d0                 ! Ref. length
  Centroid= 504.9d0, 0.d0, 0.d0    ! Centroid 
  Cood_Y_UP=1          ! 1 Y- is up; 0 Z- is up
  NUM_THREADS=1         
$end

!-----Following is comment----------------------------------------------
! control file for   OpenCFD-EC ver 1.11
!
!======C C C ============================================
! Followings are comments
! Default  values for opencfd-ec ver 1.11  

    Ma=1.d0        ! Mach number     
    Re= 1000.d0    ! Reynolds number
    gamma=1.4d0    ! 

	AOA=0.d0       ! Angle of attack    
	AOS=0.d0       ! Angle of Slide
	P_outlet=-1.d0   ! Outlet pressure (<0 extrapolation)  
	t_end=100.d0     ! End time (non-dimensional)
	Kstep_save=1000  ! Save data per xxx steps
    Iflag_turbulence_model=0   ! turbulence model (0 none, 1 BL, 2 SA, 3 SST)
	Iflag_init=0  ! 0 从初始值（均匀来流）开始计算； 1  续算； -1 从0 流场开始计算
    If_viscous=1  ! 0 无粘； 1 有粘
    Iflag_local_dt=1   ! 0 全局步长；  1 局部时间步长
    dt_global=0.01     ! Global time step
    CFL=1.d0           ! CFL number 
    dtmax=10.d0       !Limit of maximum time step 
	dtmin=1.d-9       !Limit of minimum time step
	Time_Method=0     ! Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
    If_Residual_smoothing=0   ! 0 Do not need smoothing
    w_LU=1.d0         ! factor in LU-SGS
    If_dtime_mesh=1   ! decrease time step when grid quality is not good
    Iflag_Scheme= 5    ! 0 UD1;  1 NND2 ; 2 UD3 ; 3 MUSCL2U ;  4 MUSCL2C; 5 MUSCL3; 6 OMUSCL2; 7 WENO5 ; 8 UD5; 9 WENO7 
    Iflag_Flux=5      ! 1 Steger_Warming; 2 HLL, 3 HLLC,  4 Roe, 5  Van_Leer, 6 Ausm
    IFlag_Reconstruction=0   !  0 Original, 1 Conservative, 2 Characteristic
    Mesh_File_Format=0    ! 0 unformatted, 1 formatted
    Kstep_show=1          ! show per xxx steps
    Kstep_smooth=-1       ! Smoothing flow per xxx steps (<0 donot smooth)
    Kstep_init_smooth=0   ! Smoothing flow at initial time
    Num_Mesh=1            ! Number of mesh (1 single-grid), 2,3 multi-grid
    T_inf=288.15d0        ! Reference temperature  (in viscous coefficient)
    Twall=-1.d0           ! Wall temperature (in K degree) (<0 adabitic)
    Kt_inf=1.d-5          ! initial of Kt for SST model
    Wt_inf=0.01d0         ! Initial of Wt for SST model
    IF_Debug=0            ! 1 Debug mode
    NUM_THREADS=1         ! Threads for OpenMP
    Step_Inner_Limit=20   ! Inner time advance limit (Dual-time)
    Res_Inner_Limit=1.d-10 ! Inner Resdial limit (Dual-time)
    MUT_MAX=-1.d0          ! Limit for vt/vs (<0 no limit)
	Bound_Scheme= Scheme_MUSCL2C         ! Boundary scheme (Default: MUSCL2C)
    Pre_Step_Mesh(1:3)=0   ! Pre step for multi-grid
    Ref_S=1.d0             ! Ref. area
	Ref_L=1.d0             ! Ref. length
	Centroid(1:3)=0.d0  ! Centroid coordinate  
    Cood_Y_UP=1   !       默认Y轴垂直向上
	IFLAG_LIMIT_FLOW=0         ! 限制流场（密度、速度、压力）
	Pdebug(1:4)=1
    PrL=0.7d0          ! Linear Prandtl number
	PrT=0.9d0          ! Turbulent Prandtl number
    Ldmin=1.d-6
	Ldmax=1000.d0
	Lpmin=1.d-6
	Lpmax=1000.d0
	Lumax=1000.d0
	LSAmax=1000.d0
	CP1_NSA=0.2d0   ! for New SA
	CP2_NSA=100.d0
	Periodic_dX=0.d0   ! 周期边界的几何增量
	Periodic_dY=0.d0 
	Periodic_dZ=0.d0

!----for Turbomachinary solver------------
    IF_TurboMachinary=0    ! 启用叶轮机计算模式
	Ref_medium_usrdef=0    ! 启用自定义介质 （为默认空气）
    IF_Scheme_Positivity=1     ! 检查插值过程中压力、密度是否非负，否则使用1阶迎风；
    Turbo_P0= 101330.d0    ! 总压 （默认为1个大气压）
	Turbo_T0= 288.15d0     ! 总温 （默认288.15K)
    Turbo_L0= 1.d0         ! 参考长度 （默认为1m)
    Turbo_w=0.d0           ! 转速  ( 转/秒 ， 默认0)
    Turbo_Periodic_seta=0.d0   ! 周期方向计算域，角

  
!==========================================  
  Ma          Re        A_alfa   A_beta     P_outlet  t_end   Kstep_save Turbulence_model Iflag_init If_viscous   
 0.75d0       2.125d4   0.d0    0.0d0     -1.d0    100.d0    5000                 2        0        1 
  Iflag_local_dt,  dt_global,  CFL,    dtmax,   dtmin  Time_Method  If_Residual_smoothing   w_LU     If_dtime_mesh
      1            1.d-4       5.d0    10.d0    1.d-9     0              0                  1.0      1
  Iflag_scheme Iflag_flux  IFlag_Reconstruction Mesh_File_Format    Kstep_show Kstep_smooth,Kstep_init_smooth
      5         5           0                        0                         1           -1           0
  Num_Mesh   T_inf     T_Wall   Kt_inf   Wt_inf     IF_Debug    NUM_THREADS  Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme
     1        288.15   -1.0      1.d-5    0.01          0            1           0            0                  -1      -1
Pre_Step_Mesh(1)   Pre_Step_Mesh(2)   Pre_Step_Mesh(3)
     0                  0                0
#------------------------------------------------------
 # OpenCFD-EC Ver 1.02 
     open(99,file="control.in")
     read(99,*)
     read(99,*) Ma, Re, A_alfa, A_beta, p_outlet, t_end,Kstep_save, Iflag_turbulence_model,Iflag_init,If_viscous
     read(99,*)
     read(99,*) Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method,If_Residual_smoothing,w_LU,If_dtime_mesh
     read(99,*)
     read(99,*) Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction,Mesh_File_Format,Kstep_show,Kstep_smooth,Kstep_init_smooth
     read(99,*)
     read(99,*) Num_Mesh,T_inf,Twall,Kt_inf,Wt_inf,IF_Debug,NUM_THREADS, Step_Inner_Limit, Res_Inner_Limit
     read(99,*)
     read(99,*) (Pre_Step_Mesh(k), k=1,Num_Mesh)
     close(99)

   integer,parameter:: Scheme_UD1=0,Scheme_NND2=1, Scheme_UD3=2,Scheme_WENO3=3,Scheme_MUSCL2=4,Scheme_MUSCL3=5,Scheme_OMUSCL2=6,Scheme_WENO5=7,Scheme_UD5=8 
   integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_Van_Leer=5,Flux_Ausm=6
   integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2
   integer,parameter:: BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6, BC_Periodic=501     

   integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
   integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2, Turbulence_SST=3
   integer,parameter:: Init_By_FreeStream=0,Smooth_2nd=0,Smooth_4th=1 
   real(PRE_EC), parameter::  Temperature_LIMIT=1.d-5,Density_LIMIT=1.d-5
 
