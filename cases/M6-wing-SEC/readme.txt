  ONEAR-M6 Wing 计算 （差分-有限体积混合方法）使用方法:

  1.  准备好网格文件Mesh3d.dat (可到../M6-wing/ 目录下拷贝），边界条件文件bc3d.inp, 主控制文件control.ec, 以及差分方法控制文件FDM.in;
  2.  准备好可执行代码opencfd-ec1.12.out
  3.  运行 mpirun -np 4 ./opencfd-ec1.12.out
  4.  待残差收敛后 （例如20000步）停止计算；
  5. 运行 readflow3d-ver2.1.out 后处理 （在../../util/ 目录下）。

 注：
    自1.12版本起， OpenCFD-EC 恢复对内嵌差分方法的支持。 有限体积-差分混合方法的理论说明见《OpenCFD-EC理论手册》;
    内嵌差分方法由控制文件FDM.in 控制。 如运行目录中找不到该文件，则程序不启用内嵌差分法。
    FDM.in 第3 行的参数FD_Flux 为内嵌差分法所采用的通量方法(目前版本只支持Steger-Warming 及Van Leer 分裂）；
                  参数 FD_Scheme 为内嵌差分法所采用的差分格式 （目前版本支持WENO5, WENO7 及OMP6);
    第5行的参数nbk 为采用内嵌差分法的块数；第7行为内嵌差分块的列表。 如 2， 3  则表示第2块、第3块采用内嵌差分计算。
