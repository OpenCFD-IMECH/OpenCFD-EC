-------OpenCFD-EC  (c) by Li Xinliang, lixl@imech.ac.cn-------------------

编译方法： 
         主程序编译，运行make即可 （如有问题，请修改makefile).
            注： 需要MPI库及Fortran 90 编译器支持；
         
         CGNS网格转换器编译:  ifort -O3 -o convert-cgns-inp-1.1a.out convert-cgns-inp-1.1a.f90  -lcgns
            注： 需要CGNS库， 下载地址 http://cgns.sourceforge.net/
                安装方法（以cgnslib_3.1.4为例）： 
                 1) 下载原始文件 （例如cgnslib_3.1.4.tar.gz), 解压
                 2） 进入 ./cgnslib_3.1.4/src
                 3) ./configure --with-fortran
                 4) make  (或make all, 编译全部工具)
                 5) sudo make install (需要超级用户口令）
                   （该命令将头文件(*.h)拷贝到./usr/local/include 文件夹下面； 将库文件 （libcgns.a）拷贝到/usr/local/lib文件夹下）
                   如无超级用户口令，可手工将 *.h及 libcgns.a 拷贝到用户指定的目录。

运行方法：
          1) 准备好网格文件Mesh3d.dat (PLOT3D格式， format与unformat 均可； 即Gridgen的 Generic格式）
          2）准备好边界条件文件 bc3d.inp  (Gridgen 的Generic 格式）
          3) 如使用CGNS格式的网格文件 (grid.cgns)，可用convert-cgns-1.1a.out 将其转换为Mesh3d.dat及bc3d.inp 
          4）编辑修改控制文件control.ec 设置 流动参数及计算方法等；
             其中Ma,Re 分别为Mach数及Reynolds数；  
             AoA 为攻角 ， AoS为侧滑角 ; (单位为度)

          5）运行 mpirun -np xx ./opencfd-1.04.out （或用bsub 等方式提交）
            xx 可为1,2,3,4, ... 但不能超过网格的总块数。
          
          注1：
           如无wall_dist.dat (记录每个网格点到最近壁面的距离)， 则程序开始时先生成该文件，需要消耗一些时间；
           运行时输出残差到Residual.dat （各方程的最大残差即均方根残差）
           运行时输出整个物体承受的气动力6分量到force.log 
              各列含义： 计算步， CL, CD, CS, CMx, CMy, CMz 
			  (升力系数、阻力系数、侧向力系数，三个方向的力矩系数， Ver1.04起，均按照通常定义无量纲化）。

           对于定常问题，请根据残差（或气动力）的收敛情况，判断是否可结束计算。

         注2：
            如使用OpenMP, 则需要：
                1）在makefile中打开openmp开关（设定 opt=-O3 -openmp），然后重新编译。         
                2）在控制文件中设定NUM_THREADS的数目 （例如NUM_THREADS=4， 设定为1为单线程串行）。


后处理方法： 运行 readflow3d-ver1.7.f90； 该程序将以tecplot格式输出壁面及对称面上的二维场，以及全流场的三维场。 

CGNS网格处理方法：
      本版通过转换器本支持CGNS网格，具体方法：
       1） 准备好CGNS格式的网格文件 grid.cgns
       2） 运行 ./convert-cgns-1.1a.out (读入grid.cgns, 生成Mesh3d.dat 及bc3d.inp)


---------------------相关输入数据-----------------------------------------------------
 以下文件夹中包含了一些示例计算，为了减少数据量，文件夹中并不包含网格文件，网格文件请到指定位置下载。
 
 1） ./case/DLR-F4:   DLR-F4 翼身组合体 （粗网格）
     具体使用方法见该目录下的readme.txt
   
 2） ./case/M6-wing:   ONEAR M6-wing 三维翼 （粗网格）
      control.ec 主控文件
      bc3d.inp 边界条件文件  
      网格文件(Mesh3d.dat) 请到作者的网盘下载 （https://skydrive.live.com/?cid=1cc0dcbff560c149, 公开->OpenCFD-EC->算例->M6-Wing）
      该网格为4块，因而MPI并行运行，最多不超过4进程 （推荐2进程）。



!-----------------------------
 版本更新记录：
  0.98c: 改写了SA模型部分，采用类似无量纲计算（含Re）， 参考《CFL3D手册》 （与CFL3D无量纲化有所区别，因而不包含Mach数）；
  0.99: 改写了滤波部分
  1.00: 支持自定义边界条件
  1.01: 增加了边界格式功能（允许边界使用低阶精度格式）
  1.02: 增加了WENO7格式
  1.03: 对控制文件的格式进行了修改支持 NAMELIST格式
  1.04: 对输出进行了修改（报错信息，气动力、矩信息）
  1.05: 对差分格式进行了修改； 去掉了WENO3格式，补充了 MUSCL2 with k=-1 (迎风型); 原MUSCL2 是k=1的中心型MUSSCL
  1.05a: 对压力、密度的增量进行了限制(类似CFL3D的做法)，以增强稳定性；
  1.05b: 修改了角点处理方式 (角点薄层近似)
  1.05d: 在压力、密度限制方面进行了修改，（局部平均）
  1.06:  采用限制措施（限制最低、最高密度、压力及速度）
  1.06a: 对SA等扩展方程也进行了限制；
  1.07:  对OMUSCL2进行了修改， 添加了New-SA模型；
  1.07a: new SA模型的系数可在控制文件中修改；
  1.1:  增加了叶轮机械（压气机、涡轮等）计算功能，使用旋转坐标系 （含源项）
  1.1a: 修正了背压定义的Bug  （压力应当用动压无量纲化，而不是用来流总压）
  1.11：修正了与X,Y,Z坐标定义及计算CL,CD的Bug, 使得对于Z-坐标轴向上的网格也能输出正确的气动力系数 （需设置Cood_Y_UP=0)
  1.12: 恢复对差分-有限体积混合方法的支持；
  1.13: 增加了内流模式；与叶轮机械模式类似(但不含旋转);
  1.13b: 修改了Roe 格式中的Bug;
  1.13c：添加了时间相关入口边界条件 （读取时间序列）， 修正了Runge-Kutta 3格式中OpenMP bug， 支持时间平均
  1.14:  增加了数据文件output方式（可写入不同的文件名）；
  1.14a: 添加了含吹吸气扰动的壁面边界条件 （作为用户自定义边界条件）
  1.15： A Bug in comput_force( ) is removed;  修正了后处理粘性力计算的Bug （i-, j-, k- 面 cf应为负！）
  1.15a: 修正了call output_vt() 中的Bug;  修正了readflow3d-ver2.4.f90中的Bug (read wall_dist.dat)；
         对边界条件进行了修改，  Bc%bc .eq. BC_Inflow .and. IF_InnerFlow .eq. 0 时， 采用完全给定来流的边界条件 （不区分超、亚声速）；
  1.16： 提升程序鲁棒性。  如发现某块出现物理量超限（如负温度），则强制该块采用1阶迎风格式计算，该块时间步降为原先的1/10， 该块不计算粘性项；
   1.16a:  sub_turbulence_SST 中的Bug 被修正；

   