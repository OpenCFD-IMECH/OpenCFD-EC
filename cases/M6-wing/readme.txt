M6 Wing 算例， 参考文献： http://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html
使用方法：
  1. 将可执行文件 opencfd-ec-1.04.out 拷贝入该文件夹；
  2. 运行 mpirun -np 2 ./opencfd-ec-1.04.out
  3. 监控残差 (屏幕及Residual.dat中输出）， 当残差下降到一定程度时，终止计算；
  4. 利用后处理程序readflow3d-ver1.8a.f90 或 readflow3d-wing.f90 (程序/util目录下)，进行后处理 （读取flow3d.dat)。

 