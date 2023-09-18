!  -----------------------读取流动参数及控制变量----------------------
  subroutine read_parameter
   use Global_var
   implicit none
   logical ext1
!-----------------------------------------------------------------------------
   call set_default_parameter ! 设置参数默认值

   if(my_id .eq. 0) then
     inquire(file="control.ec",exist=ext1)
      if(ext1) then           ! 优先读取control.ec (Namelist 格式控制文件) 
       call read_parameter_ec
      else
       print*, "Can not find 'control.ec', stop !"
	   stop
	  endif

 
   endif

   call bcast_para      ! 广播至全部进程
   
   call set_const_para
!--------------------------------------------------------------------

  end subroutine read_parameter
   
!--------------------------------------------------------------------    
! 设置参数的默认值  
  subroutine set_default_parameter 
   use Global_var
   implicit none
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
    Kstep_average=0       ! do not average
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
	  IFLAG_LIMIT_FLOW=1         ! 限制流场（密度、速度、压力），设定为1 ， 2016-10-21
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

    Iflag_savefile=0       ! 默认写入flow3d.dat
!----for Turbomachinary solver------------
    IF_TurboMachinary=0    ! 启用叶轮机计算模式
	Ref_medium_usrdef=0    ! 启用自定义介质 （为默认空气）
    IF_Scheme_Positivity=1     ! 检查插值过程中压力、密度是否非负，否则使用1阶迎风；
    Turbo_P0= 101330.d0    ! 总压 （默认为1个大气压）
	Turbo_T0= 288.15d0     ! 总温 （默认288.15K)
    Turbo_L0= 1.d0         ! 参考长度 （默认为1m)
    Turbo_w=0.d0           ! 转速  ( 转/秒 ， 默认0)
    Turbo_Periodic_seta=0.d0   ! 周期方向计算域，角
!-------------------------------------
    IF_InnerFlow=0   ! 内流模式
end

!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_ec 
   use Global_var
   implicit none
   real(PRE_EC):: R0, a0, d0, mu0,mu1

 	namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
	    gamma, PrL, PrT, &
	    Kstep_save, &
	    Iflag_turbulence_model,Iflag_init,If_viscous,  &
        Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method, &
		If_Residual_smoothing,w_LU,If_dtime_mesh,  &
        Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, &
		Mesh_File_Format,Kstep_show,Kstep_average,Kstep_smooth,Kstep_init_smooth,  &
        Num_Mesh,T_inf,Twall,Kt_inf,Wt_inf,IF_Debug,NUM_THREADS,  &
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme, &
        Pre_Step_Mesh,Ref_S,Ref_L,Centroid,Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug, &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
        IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ, &
		IF_Innerflow, Iflag_savefile


	open(99,file="control.ec")
	read(99,nml=control_ec)
    close(99)
 
 !---- convert parameters ----------------------
 ! Ref_medium_usrdef==0 使用默认介质 (Ma=1, 根据总温、总压计算 Re) ； ==1 使用自定义介质 （人为输入Ma, Re等）
 ! 仅适用于叶轮机或内流模式；

    if( (IF_TurboMachinary ==1 .or. IF_Innerflow==1 )    &  
	  .and.  Ref_medium_usrdef == 0) then   ! 默认空气介质，计算Mach数， Reynolds数
      
	  T_inf=Turbo_T0  ! 参考温度 （来流总温）
      gamma=1.4d0    ! 
	  PrL=0.7d0   ! Prandtl数
	  PrT=0.9d0
      R0= 287.06d0   ! 空气的气体常数R
	  a0= sqrt(gamma*R0*Turbo_T0)    ! 参考温度下的声速 
	  mu0=1.179d-5     ! 空气粘性系数 (288.15K)  
      mu1=mu0* sqrt((Turbo_T0/288.15d0)**3)*(288.15d0+110.4d0)/(Turbo_T0+110.4d0)  ! 参考温度下的空气粘性系数
      d0=Turbo_P0/(R0*Turbo_T0)
	  Re=d0*a0*Turbo_L0/mu1    ! 参考温度下，以声速运动的Reynolds数
	  Ma=1.d0     ! Mach数    （以声速作为参考速度，因而参考Mach数为1）
      Turbo_w= 2.d0*PI*Turbo_w/(a0/Turbo_L0)   ! 无量纲角速度 Turbo_W（转/秒）
!	  P_outlet=P_outlet/Turbo_P0    ! 背压,  Bug !!
	  P_outlet=P_outlet/(d0*a0*a0)    ! 背压 （用动压 无量纲）
	endif
      Turbo_Periodic_seta=Turbo_Periodic_seta*PI/180.d0                 ! Turbo_Periodic_seta 角度


 !---output paramters----------------------------
    open(99,file="output_para.out")
	write(99,*) "-------OpenCFD-EC (Ver 1.1t), (c) Li Xinliang, lixl@imech.ac.cn--"
	write(99,*) "-------------------------------------------"
	write(99,*) "Ma=", Ma,  "  Re= ", Re , "gamma=", gamma
	write(99,*) "PrL=", PrL, "PrT=",PrT

	write(99,*) "A_alfa=", A_alfa, " A_beta=", A_beta 
	write(99,*) "P_outlet=", P_outlet
	write(99,*) "t_end=", t_end, " Iflag_local_dt=", Iflag_local_dt
	write(99,*) "dt_global=",dt_global 
	write(99,*) "Time_Method=", Time_Method,  " CFL= ", CFL 
	write(99,*) "dtmax=", dtmax, " dtmin=", dtmin 
    write(99,*) "Iflag_init=", Iflag_init," If_viscous=", If_viscous
	write(99,*) "Iflag_turbulence_model=",Iflag_turbulence_model
    write(99,*) "Iflag_Scheme= ",Iflag_Scheme, " Iflag_Flux=", Iflag_Flux
	write(99,*) "IFlag_Reconstruction=", IFlag_Reconstruction, "Bound_Scheme= ", Bound_Scheme
    write(99,*) "Kstep_save= ",Kstep_save,  " Kstep_show=", Kstep_show, "Kstep_average=",Kstep_average
	write(99,*) "Kstep_smooth= ", Kstep_smooth, " Kstep_init_smooth=",Kstep_init_smooth
    write(99,*) "If_Residual_smoothing=", If_Residual_smoothing, " If_dtime_mesh=",If_dtime_mesh
	write(99,*) "w_LU=",w_LU
    write(99,*) "Mesh_File_Format=",Mesh_File_Format, " Num_Mesh=",Num_Mesh
    write(99,*) "T_inf=", T_inf, " Twall=", Twall
	write(99,*) "Kt_inf=",Kt_inf,  " Wt_inf=",Wt_inf
    write(99,*) "Step_Inner_Limit=",Step_Inner_Limit, " MUT_MAX=", MUT_MAX
	write(99,*) "Pre_Step_Mesh(:)=", Pre_Step_Mesh(1:Num_Mesh)
    write(99,*) "IF_Debug=",IF_Debug
	write(99,*) "Ref_S=", Ref_S, "Ref_L=", Ref_L
	write(99,*) "Centroid=", Centroid(1:3)
	write(99,*) "Cood_Y_UP=",Cood_Y_UP
    write(99,*) "Periodic_dX, dY, dZ=", Periodic_dX,Periodic_dY,Periodic_dZ
	write(99,*) "IFLAG_LIMIT_FLOW=",IFLAG_LIMIT_FLOW
	write(99,*) "Ldmin,Ldmax,Lpmin,Lpmax,Lumax=", Ldmin,Ldmax,Lpmin,Lpmax,Lumax
	write(99,*) "CP1_NSA,CP2_NSA=",CP1_NSA,CP2_NSA
    write(99,*) "IF_TurboMachinary=",  IF_TurboMachinary
	write(99,*) "Turbo_Periodic_seta=", Turbo_Periodic_seta
    write(99,*) "Ref_medium_usrdef=", Ref_medium_usrdef
	write(99,*) "Turbo_w= ", Turbo_w, "Turbo_L0=", Turbo_L0
	write(99,*) "a0 (Reference velocity)=", a0
	write(99,*) "d0, T0, p0 (Reference values)=", d0,Turbo_T0, Turbo_P0
	write(99,*) "IF_Scheme_Positivity =", IF_Scheme_Positivity
	write(99,*) "Iflag_savefile=", Iflag_savefile
	write(99,*) "NUM_THREADS=",NUM_THREADS
    write(99,*) "--------------------------------------------"

    close(99)

 !----------------------------------------------
 end
!--------------------------------------------------
 

!------------------------------------------------------------------
  
  
   subroutine bcast_para
   use Global_var
   implicit none
    integer:: Ipara(100),ierr
    real(PRE_EC):: rpara(100)
    Ipara=0
	rpara=0.d0
!----
	rpara(1)=Ma 
	rpara(2)=Re
	rpara(3)=AoA
	rpara(4)=AoS
	rpara(5)=p_outlet
	rpara(6)=t_end
	rpara(7)=dt_global
	rpara(8)=CFL
	rpara(9)=dtmax
	rpara(10)=dtmin
	rpara(11)=w_LU
	rpara(12)=T_inf
	rpara(13)=Twall
	rpara(14)=Kt_inf
	rpara(15)=Wt_inf
    rpara(16)=Res_Inner_Limit
    rpara(17)=MUT_MAX
    rpara(18)=Ref_S
	rpara(19)=Ref_L
	rpara(20:22)=Centroid(1:3)
	rpara(23)=Ldmin
	rpara(24)=Ldmax
	rpara(25)=Lpmin
	rpara(26)=Lpmax
	rpara(27)=Lumax
	rpara(28)=LSAmax
	rpara(29)=CP1_NSA
	rpara(30)=CP2_NSA
	rpara(31)=gamma
	rpara(32)=PrL
	rpara(33)=PrT
    rpara(34)=Turbo_w
    rpara(35)=Turbo_Periodic_seta
    rpara(36)=Periodic_dX
	rpara(37)=Periodic_dY
	rpara(38)=Periodic_dZ




    Ipara(1)=Kstep_save
	Ipara(2)=Iflag_turbulence_model
	Ipara(3)=Iflag_init
	Ipara(4)=If_viscous
    Ipara(5)=Iflag_local_dt
	Ipara(6)=Time_Method
	Ipara(7)=If_Residual_smoothing
	Ipara(8)=If_dtime_mesh
    Ipara(9)=Iflag_Scheme
	Ipara(10)=Iflag_Flux
	Ipara(11)=IFlag_Reconstruction
	Ipara(12)=Mesh_File_Format
	Ipara(13)=Kstep_show
	Ipara(14)=Kstep_smooth
	Ipara(15)=Kstep_init_smooth
    Ipara(16)=Num_Mesh
	Ipara(17)=IF_Debug
	Ipara(18)=NUM_THREADS
    Ipara(19)=Step_Inner_Limit
    Ipara(20)=Bound_Scheme
    Ipara(21)=Cood_Y_UP
	Ipara(22)=IFLAG_LIMIT_Flow
    Ipara(23:26)=Pdebug(1:4)
    Ipara(27)=IF_TurboMachinary
	Ipara(28)=IF_Scheme_Positivity
   	Ipara(29)=IF_Innerflow
   	Ipara(30)=Kstep_average
    Ipara(31)=Iflag_savefile

	 call MPI_bcast(rpara,100,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)
	 call MPI_bcast(Ipara,100,MPI_Integer,0,  MPI_COMM_WORLD,ierr)

	Ma=rpara(1) 
	Re=rpara(2)
	AoA=rpara(3)
	AoS=rpara(4)
	p_outlet=rpara(5)
	t_end=rpara(6)
	dt_global=rpara(7)
	CFL=rpara(8)
	dtmax=rpara(9)
	dtmin=rpara(10)
	w_LU=rpara(11)
	T_inf=rpara(12)
	Twall=rpara(13)
	Kt_inf=rpara(14)
	Wt_inf=rpara(15)
    Res_Inner_Limit=rpara(16)
    MUT_MAX=rpara(17)
    Ref_S=rpara(18)
	Ref_L=rpara(19)
	Centroid(1:3)=rpara(20:22)
	Ldmin=rpara(23)
	Ldmax=rpara(24)
	Lpmin=rpara(25)
	Lpmax=rpara(26)
	Lumax=rpara(27)
	LSAmax=rpara(28)
	CP1_NSA=rpara(29)
	CP2_NSA=rpara(30)
	gamma=rpara(31)
	PrL=rpara(32)
	PrT=rpara(33)
    Turbo_w=rpara(34)
    Turbo_Periodic_seta=rpara(35)
    Periodic_dX=rpara(36)
	Periodic_dY=rpara(37)
	Periodic_dZ=rpara(38)



    Kstep_save=Ipara(1)
	Iflag_turbulence_model=Ipara(2)
	Iflag_init=Ipara(3)
	If_viscous=Ipara(4)
    Iflag_local_dt=Ipara(5)
	Time_Method=Ipara(6)
	If_Residual_smoothing=Ipara(7)
	If_dtime_mesh=Ipara(8)
    Iflag_Scheme=Ipara(9)
	Iflag_Flux=Ipara(10)
	IFlag_Reconstruction=Ipara(11)
	Mesh_File_Format=Ipara(12)
	Kstep_show=Ipara(13)
	Kstep_smooth=Ipara(14)
	Kstep_init_smooth=Ipara(15)
    Num_Mesh=Ipara(16)
	IF_Debug=Ipara(17)
	NUM_THREADS=Ipara(18)
    Step_Inner_Limit=Ipara(19)
    Bound_Scheme=Ipara(20)
    Cood_Y_UP=Ipara(21)
	IFLAG_LIMIT_FLOW=Ipara(22)
    Pdebug(1:4)=Ipara(23:26)
    IF_TurboMachinary=Ipara(27)
	IF_Scheme_Positivity=Ipara(28)
   	IF_Innerflow=Ipara(29)
   	Kstep_average=Ipara(30)
    Iflag_savefile=Ipara(31)


    call MPI_bcast(Pre_Step_Mesh,Num_Mesh,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
    end  subroutine bcast_para
  
  
      
! 设定常数
   subroutine set_const_para 
   use Global_var
   implicit none

   Twall=Twall/T_inf                         ! wall temperature
   Lpmin=Lpmin/(gamma*Ma*Ma)                 ! rato of free-stream pressure
   Lpmax=Lpmax/(gamma*Ma*Ma)



   if(Bound_Scheme== Scheme_none)  Bound_Scheme=Iflag_Scheme    ! 如不使用边界格式，则与内点格式一致   
   if(If_viscous .eq. 0) Iflag_turbulence_model=Turbulence_NONE !    求解无粘方程，不采用湍流模型
    
   if(Iflag_turbulence_model .eq. Turbulence_SA .or. & 
      Iflag_turbulence_model .eq. Turbulence_NewSA) then
     NVAR=6                       ! SA模型，总共6个变量
   else if (Iflag_turbulence_model .eq. Turbulence_SST) then
     NVAR=7                       ! SST 模型，总共7个变量
   else
     NVAR=5                       ! 5个变量
   endif

  if(Iflag_turbulence_model .eq. 0) then
   IF_Walldist=0                   ! 无需该数据
  else
   IF_Walldist=1
  endif
 
!--------------------------------------------------------------------------------
  AoA=AoA*PI/180.d0   ! Angle of attack
  AoS=AoS*PI/180.d0   ! Angle of Slide
 if(Cood_Y_UP ==1) then   ! Y 轴垂直向上 or Z轴垂直向上 
   A_alfa=AoA             ! Y轴向上， A_alfa为攻角
   A_beta=AoS
 else
   A_alfa=AoS  
   A_beta=AoA
 endif
   
  
   Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma) 
   Cp=Cv*gamma 
!--------------------------------------------------------------------

!  计算到壁面的距离
!   if(Iflag_turbulence_model .eq. Turbulence_SA .or. Iflag_turbulence_model .eq. Turbulence_SST) then
!     inquire(file="wall_dist.dat",exist=file_exist)
!	 if( .not. file_exist) then
!	   call  comput_wall_dist  
!     endif
!   endif
    
	if(Num_Mesh .ne. 1) then
	   if(Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_Dual_LU_SGS) then
	    print*, "In this version, LU_SGS (or Dual_time_LU_SGS) DO NOT Support Multi-Grid !!!"
		stop
	   endif
	endif
  end



