!  后处理模块,  计算力和力矩
!  A bug removed, 2017-3-13 
!--------------------------------------------------------
  subroutine comput_force      
   use Global_Var
   implicit none
   integer:: i,j,k,m,mB,nf,nx,ny,nz,NM,ierr
   real(PRE_EC):: Fx,Fy,Fz,Mx,My,Mz   ! 6分量
   real(PRE_EC):: Px,Py,Pz,Cfx,Cfy,Cfz   ! 压力和摩擦力
   real(PRE_EC):: fx0,fy0,fz0,xc,yc,zc,p1,p2
   real(PRE_EC):: P_inf,Pw  
   real(PRE_EC):: CL,CD,CS      ! 升力系数、阻力系数、侧向力系数
   real(PRE_EC),dimension(:),allocatable:: Mx1,My1,Mz1,Px1,Py1,Pz1,Cfx1,Cfy1,Cfz1    ! 各块的气动力、力矩
   real(PRE_EC),dimension(9):: Ft,Ft0  ! mpi reduce 汇总
   
   integer:: nMesh
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   character(len=50):: filename

!   print*, "comput force ..."    ! 如没有该语句在SW上运行出错 ???

   p_inf=1.d0/(gamma*Ma*Ma)
! 搜索细网格所有块的所有子面，如果发现固壁边界条件，则统计气动力及力矩
   NM=Mesh(1)%Num_Block
   allocate(Mx1(NM),My1(NM),Mz1(NM),Px1(NM),Py1(NM),Pz1(NM),Cfx1(NM),Cfy1(NM),Cfz1(NM))
   
   Mx1=0.d0; My1=0.d0; Mz1=0.d0
   Px1=0.d0; Py1=0.d0; Pz1=0.d0; Cfx1=0.d0; Cfy1=0.d0; Cfz1=0.d0
   
   do mB=1,NM
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
     do nf=1,B%subface  
       Bc=> B%bc_msg(nf)

       if(Bc%bc .eq. BC_WALL) then
         if(Bc%face .eq. 1 ) then              ! i- 面
           do k=Bc%kb,Bc%ke-1
             do j=Bc%jb,Bc%je-1
!-------------------------------表面压力----------------------------------------------------  
               p1=(gamma-1.d0)*(B%U(5,1,j,k)-(B%U(2,1,j,k)**2+B%U(3,1,j,k)**2+B%U(4,1,j,k)**2)/B%U(1,1,j,k))
               p2=(gamma-1.d0)*(B%U(5,0,j,k)-(B%U(2,0,j,k)**2+B%U(3,0,j,k)**2+B%U(4,0,j,k)**2)/B%U(1,0,j,k))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%si(1,j,k)    ! 外法向
! -------------------------积分表面压力及表面摩擦阻力 -----------------------------------          
               Px1(mB)=Px1(mB)+Pw*B%ni1(1,j,k) ;   Py1(mB)=Py1(mB)+Pw*B%ni2(1,j,k) ;  Pz1(mB)=Pz1(mB)+Pw*B%ni3(1,j,k)   
!              粘性力 （i-, j-, k- 为正； i+ , j+, k+ 为负）
               Cfx1(mB)=Cfx1(mB)+B%Surf1(j,k,1) ;  Cfy1(mB)=Cfy1(mB)+B%Surf1(j,k,2) ;  Cfz1(mB)=Cfz1(mB)+B%Surf1(j,k,3)     

! --------------------------计算力矩 (矩心坐标 centroid(1:3) )------------------------            
               fx0=Pw*B%ni1(1,j,k)+B%Surf1(j,k,1)             
               fy0=Pw*B%ni2(1,j,k)+B%Surf1(j,k,2)
               fz0=Pw*B%ni3(1,j,k)+B%Surf1(j,k,3)
               xc=(B%xc(1,j,k)+B%xc(0,j,k))*0.5d0 -centroid(1)            
               yc=(B%yc(1,j,k)+B%yc(0,j,k))*0.5d0 -centroid(2)             
               zc=(B%zc(1,j,k)+B%zc(0,j,k))*0.5d0 -centroid(3)            
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo


         else if(Bc%face .eq. 2 ) then       ! j- 面
           do k=Bc%kb,Bc%ke-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,1,k)-(B%U(2,i,1,k)**2+B%U(3,i,1,k)**2+B%U(4,i,1,k)**2)/B%U(1,i,1,k))
               p2=(gamma-1.d0)*(B%U(5,i,0,k)-(B%U(2,i,0,k)**2+B%U(3,i,0,k)**2+B%U(4,i,0,k)**2)/B%U(1,i,0,k))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%sj(i,1,k)    
               Px1(mB)=Px1(mB)+Pw*B%nj1(i,1,k);   Py1(mB)=Py1(mB)+Pw*B%nj2(i,1,k) ;   Pz1(mB)=Pz1(mB)+Pw*B%nj3(i,1,k)   
               Cfx1(mB)=Cfx1(mB)+B%Surf2(i,k,1); Cfy1(mB)=Cfy1(mB)+B%Surf2(i,k,2);  Cfz1(mB)=Cfz1(mB)+B%Surf2(i,k,3)

 ! --------------------------计算力矩 (以坐标原点为中心) ---------------------------           
               fx0=Pw*B%nj1(i,1,k)+B%Surf2(i,k,1)                
               fy0=Pw*B%nj2(i,1,k)+B%Surf2(i,k,2)
               fz0=Pw*B%nj3(i,1,k)+B%Surf2(i,k,3)
               xc=(B%xc(i,1,k)+B%xc(i,0,k))*0.5d0 -centroid(1)            
               yc=(B%yc(i,1,k)+B%yc(i,0,k))*0.5d0 -centroid(2)            
               zc=(B%zc(i,1,k)+B%zc(i,0,k))*0.5d0 -centroid(3)            
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo        


         else if(Bc%face .eq. 3 ) then       ! k- 面
           do j=Bc%jb,Bc%je-1
             do i=Bc%ib,Bc%ie-1   
               p1=(gamma-1.d0)*(B%U(5,i,j,1)-(B%U(2,i,j,1)**2+B%U(3,i,j,1)**2+B%U(4,i,j,1)**2)/B%U(1,i,j,1))
               p2=(gamma-1.d0)*(B%U(5,i,j,0)-(B%U(2,i,j,0)**2+B%U(3,i,j,0)**2+B%U(4,i,j,0)**2)/B%U(1,i,j,0))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%sk(i,j,1)   
               Px1(mB)=Px1(mB)+Pw*B%nk1(i,j,1);  Py1(mB)=Py1(mB)+Pw*B%nk2(i,j,1) ;  Pz1(mB)=Pz1(mB)+Pw*B%nk3(i,j,1)   
               Cfx1(mB)=Cfx1(mB)+B%Surf3(i,j,1);  Cfy1(mB)=Cfy1(mB)+B%Surf3(i,j,2);  Cfz1(mB)=Cfz1(mB)+B%Surf3(i,j,3)         
			     
               fx0=Pw*B%nk1(i,j,1)+B%Surf3(i,j,1)          
               fy0=Pw*B%nk2(i,j,1)+B%Surf3(i,j,2)
               fz0=Pw*B%nk3(i,j,1)+B%Surf3(i,j,3)

               xc=(B%xc(i,j,1)+B%xc(i,j,0))*0.5d0  - centroid(1)
               yc=(B%yc(i,j,1)+B%yc(i,j,0))*0.5d0  - centroid(2)
               zc=(B%zc(i,j,1)+B%zc(i,j,0))*0.5d0  - centroid(3)
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo 

         else if(Bc%face .eq. 4 ) then              ! i+ 面
           do k=Bc%kb,Bc%ke-1
             do j=Bc%jb,Bc%je-1
 !-----------------------------------------表面压力 ---------------------------------------------------               
               p1=(gamma-1.d0)*(B%U(5,nx-1,j,k)-(B%U(2,nx-1,j,k)**2+B%U(3,nx-1,j,k)**2+B%U(4,nx-1,j,k)**2)/B%U(1,nx-1,j,k))
               p2=(gamma-1.d0)*(B%U(5,nx,j,k)-(B%U(2,nx,j,k)**2+B%U(3,nx,j,k)**2+B%U(4,nx,j,k)**2)/B%U(1,nx,j,k))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%si(nx,j,k)  ! 外法向
! ------------------------------积分表面压力及表面摩擦阻力------------------------------------------           
               Px1(mB)=Px1(mB)+Pw*B%ni1(nx,j,k) ;   Py1(mB)=Py1(mB)+Pw*B%ni2(nx,j,k) ;  Pz1(mB)=Pz1(mB)+Pw*B%ni3(nx,j,k)   
!            粘性力 （i+, j+, k+ 面为负， 壁面所受力）
!               Cfx1(mB)=Cfx1(mB)+B%Surf4(j,k,1) ;  Cfy1(mB)=Cfy1(mB)+B%Surf4(j,k,2) ;  Cfz1(mB)=Cfz1(mB)+B%Surf4(j,k,3)         ! Bug removed
               Cfx1(mB)=Cfx1(mB)-B%Surf4(j,k,1) ;  Cfy1(mB)=Cfy1(mB)-B%Surf4(j,k,2) ;  Cfz1(mB)=Cfz1(mB)-B%Surf4(j,k,3)

! -------------------------------计算力矩 (以坐标原点为中心) --------------------------------------           
               fx0=Pw*B%ni1(nx,j,k)-B%Surf4(j,k,1)              ! Bug removed 
               fy0=Pw*B%ni2(nx,j,k)-B%Surf4(j,k,2)
               fz0=Pw*B%ni3(nx,j,k)-B%Surf4(j,k,3)
               xc=(B%xc(nx-1,j,k)+B%xc(nx,j,k))*0.5d0   - centroid(1)           
               yc=(B%yc(nx-1,j,k)+B%yc(nx,j,k))*0.5d0   - centroid(2)          
               zc=(B%zc(nx-1,j,k)+B%zc(nx,j,k))*0.5d0   - centroid(3)         
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo


         else if(Bc%face .eq. 5 ) then       ! j+ 面
           do k=Bc%kb,Bc%ke-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,ny-1,k)-(B%U(2,i,ny-1,k)**2+B%U(3,i,ny-1,k)**2+B%U(4,i,ny-1,k)**2)/B%U(1,i,ny-1,k))
               p2=(gamma-1.d0)*(B%U(5,i,ny,k)-(B%U(2,i,ny,k)**2+B%U(3,i,ny,k)**2+B%U(4,i,ny,k)**2)/B%U(1,i,ny,k))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%sj(i,ny,k)
               Px1(mB)=Px1(mB)+Pw*B%nj1(i,ny,k);   Py1(mB)=Py1(mB)+Pw*B%nj2(i,ny,k) ;   Pz1(mB)=Pz1(mB)+Pw*B%nj3(i,ny,k)   
!              Cfx1(mB)=Cfx1(mB)+B%Surf5(i,k,1); Cfy1(mB)=Cfy1(mB)+B%Surf5(i,k,2);  Cfz1(mB)=Cfz1(mB)+B%Surf5(i,k,3)
               Cfx1(mB)=Cfx1(mB)-B%Surf5(i,k,1); Cfy1(mB)=Cfy1(mB)-B%Surf5(i,k,2);  Cfz1(mB)=Cfz1(mB)-B%Surf5(i,k,3)

! --------------------------------计算力矩 (以坐标原点为中心)  ---------------------------------------          
               fx0=Pw*B%nj1(i,ny,k)-B%Surf5(i,k,1)                 ! Bug removed
               fy0=Pw*B%nj2(i,ny,k)-B%Surf5(i,k,2)
               fz0=Pw*B%nj3(i,ny,k)-B%Surf5(i,k,3)
               xc=(B%xc(i,ny-1,k)+B%xc(i,ny,k))*0.5d0  - centroid(1)           
               yc=(B%yc(i,ny-1,k)+B%yc(i,ny,k))*0.5d0  - centroid(2)          
               zc=(B%zc(i,ny-1,k)+B%zc(i,ny,k))*0.5d0  - centroid(3)          
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo        


         else if(Bc%face .eq. 6 ) then       ! k+ 面
           do j=Bc%jb,Bc%je-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,j,nz-1)-(B%U(2,i,j,nz-1)**2+B%U(3,i,j,nz-1)**2+B%U(4,i,j,nz-1)**2)/B%U(1,i,j,nz-1))
               p2=(gamma-1.d0)*(B%U(5,i,j,nz)-(B%U(2,i,j,nz)**2+B%U(3,i,j,nz)**2+B%U(4,i,j,nz)**2)/B%U(1,i,j,nz))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%sk(i,j,nz)  
               Px1(mB)=Px1(mB)+Pw*B%nk1(i,j,nz);  Py1(mB)=Py1(mB)+Pw*B%nk2(i,j,nz) ;  Pz1(mB)=Pz1(mB)+Pw*B%nk3(i,j,nz)   
!               Cfx1(mB)=Cfx1(mB)+B%Surf6(i,j,1);  Cfy1(mB)=Cfy1(mB)+B%Surf6(i,j,2);  Cfz1(mB)=Cfz1(mB)+B%Surf6(i,j,3)            ! Bug removed
               Cfx1(mB)=Cfx1(mB)-B%Surf6(i,j,1);  Cfy1(mB)=Cfy1(mB)-B%Surf6(i,j,2);  Cfz1(mB)=Cfz1(mB)-B%Surf6(i,j,3)

               fx0=Pw*B%nk1(i,j,nz)-B%Surf6(i,j,1)        ! Bug removed
               fy0=Pw*B%nk2(i,j,nz)-B%Surf6(i,j,2)
               fz0=Pw*B%nk3(i,j,nz)-B%Surf6(i,j,3)

               xc=(B%xc(i,j,nz-1)+B%xc(i,j,nz))*0.5d0 - centroid(1)
               yc=(B%yc(i,j,nz-1)+B%yc(i,j,nz))*0.5d0 - centroid(2)
               zc=(B%zc(i,j,nz-1)+B%zc(i,j,nz))*0.5d0 - centroid(3)
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo
         endif
       endif
     enddo
   enddo
 
   Fx=0.d0; Fy=0.d0; Fz=0.d0; Mx=0.d0; My=0.d0; Mz=0.d0
   Px=0.d0; Py=0.d0; Pz=0.d0; Cfx=0.d0; Cfy=0.d0; Cfz=0.d0
 
 !  把各块的气动力、力矩加起来
   do mB=1,NM
   Px=Px+Px1(mB); Py=Py+Py1(mB); Pz=Pz+Pz1(mB)
   Cfx=Cfx+Cfx1(mB); Cfy=Cfy+Cfy1(mB); Cfz=Cfz+Cfz1(mB)
   Mx=Mx+Mx1(mB); My=My+My1(mB); Mz=Mz+Mz1(mB)
   enddo

   Ft(1)=Px ; Ft(2)=Py ; Ft(3)=Pz
   Ft(4)=Cfx; Ft(5)=Cfy; Ft(6)=Cfz
   Ft(7)=Mx;  Ft(8)=My; Ft(9)=Mz

!  各进程归约求和
   call MPI_ALLREDUCE(Ft,Ft0,9,OCFD_DATA_TYPE,MPI_SUM,MPI_COMM_WORLD,ierr)

!   Fx=Px+Cfx; Fy=Py+Cfy; Fz=Pz+Cfz
    Ft0(1:6)=2.d0*Ft0(1:6)/Ref_S
	Ft0(7:9)=2.d0*Ft0(7:9)/(Ref_S*Ref_L)

	Fx=Ft0(1)+Ft0(4)              ! 气动力系数
	Fy=Ft0(2)+Ft0(5)
	Fz=Ft0(3)+Ft0(6)


!------   
 if(Cood_Y_UP ==1) then   ! Y 轴垂直向上 
	CL=Fy*cos(AoA)-Fx*sin(AoA)
	CD=Fx*cos(AoA)+Fy*sin(AoA)
    Cs=Fz
 else                      ! Z轴垂直向上 
	CL=Fz*cos(AoA)-Fx*sin(AoA)
	CD=Fx*cos(AoA)+Fz*sin(AoA)
    Cs=Fy
 endif


   if(my_id .eq. 0) then
   print*, "--------------Force and Moment -----------------------------"
   print*, "Total Force Coefficient (CL, CD, CS) ="
   print*, CL, CD, CS
   print*, "Total Moment Coefficient(CMx, CMy, CMz)="
   print*, Ft0(7),Ft0(8),Ft0(9)
   if(If_Debug==1) then
    print*, "Inviscous force: Fx, Fy, Fz="
    print*, Ft0(1),Ft0(2),Ft0(3)
    print*, "viscous force: Fx, Fy, Fz="
    print*, Ft0(4),Ft0(5),Ft0(6)
   endif



   open(99,file="force.log",position="append")
      write(99,"(I7,6E20.8)") Mesh(1)%Kstep, CL,CD,CS, Ft0(7),Ft0(8),Ft0(9)
   close(99)

!   if(If_Debug==1) then
      open(99,file="force-invis-vis.log",position="append")
      write(99,"(I7,6E20.8)") Mesh(1)%Kstep, Ft0(1), Ft0(2), Ft0(3), Ft0(4), Ft0(5), Ft0(6)
   close(99)
!   endif

   endif

   if(If_Debug == 1) then
      write(filename,"('force-debug-'3I3)") my_id
	  open(88,file=filename)
      do mB=1,NM
	     B => Mesh(1)%Block(mB)                                        
	    write(88,"(I4,6F20.10)")  B%Block_no,   Px1(mB), Py1(mB), Pz1(mB), Cfx1(mB), Cfy1(mB), Cfz1(mB)
	  enddo
      close(88)
   endif

	deallocate(Mx1,My1,Mz1,Px1,Py1,Pz1,Cfx1,Cfy1,Cfz1)


  end  subroutine comput_force   
        
!-----------------------------------------------------------------------------------
!  光顺（滤波） 操作， 耗散很大的滤波操作。 用于对初值的光顺，或者计算异常（如负温度）时的光顺
!  滤波运算可以消除高频振荡，提高计算的稳定性；但也会增加耗散，降低精度
!  2阶精度滤波耗散非常大，只能在处理初值或异常是使用，不可在常规的计算中使用。
!  4阶精度滤波也有一定耗散，计算过程中需谨慎使用 
  subroutine smoothing_oneMesh(nMesh,Smooth_method)     
   use Global_Var
   implicit none
   integer:: nMesh,mBlock,Smooth_method
!   print*, "Filtering ......", nMesh
   do mBlock=1,Mesh(nMesh)%Num_Block
   if(Smooth_method .eq. Smooth_2nd)then
    call   smoothing_oneBlock_2nd(nMesh,mBlock)     
   else
    call   smoothing_oneBlock_4th(nMesh,mBlock)     
   endif
   enddo
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   end subroutine smoothing_oneMesh

!----------------------------------------------------------
! 低精度滤波（2阶精度）
  subroutine smoothing_oneBlock_2nd(nMesh,mBlock)     
   use Global_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,istat
   Type (Block_TYPE),pointer:: B
   real(PRE_EC)::tmpa(0:2000),tmpb(0:2000),tmpc(0:2000)
   NVAR1=Mesh(nMesh)%NVAR
   B=>Mesh(nMesh)%block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
!----------------------------------------------------------------------------
! Warning, allocatable不能作为私有变量 !!!   
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)

!$OMP DO   
   do k=1,nz-1
   do j=1,ny-1
   do m=1, NVAR1
    do i=0,nx
    tmpa(i)=B%U(m,i,j,k)
    enddo
   do i=1,nx-1
     B%U(m,i,j,k)=0.25d0*(tmpa(i-1)+tmpa(i+1))+0.5d0*tmpa(i)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do k=1,nz-1
   do i=1,nx-1
   do m=1, NVAR1
   do j=0,ny
    tmpb(j)=B%U(m,i,j,k)
   enddo
   do j=1,ny-1
     B%U(m,i,j,k)=0.25d0*(tmpb(j-1)+tmpb(j+1))+0.5d0*tmpb(j)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do j=1,ny-1
   do i=1,nx-1
   do m=1, NVAR1
   do k=0,nz
    tmpc(k)=B%U(m,i,j,k)
   enddo
   do k=1,nz-1
     B%U(m,i,j,k)=0.25d0*(tmpc(k-1)+tmpc(k+1))+0.5d0*tmpc(k)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!-----------------------------------------------------------------------------------------------
   end  subroutine smoothing_oneBlock_2nd
   
   
!------------------------------------------------------------------------------------------------
! 高精度滤波（4阶精度）
  subroutine smoothing_oneBlock_4th(nMesh,mBlock)     
   use Global_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1
   Type (Block_TYPE),pointer:: B
   real(PRE_EC)::tmpa(0:2000),tmpb(0:2000),tmpc(0:2000)

   NVAR1=Mesh(nMesh)%NVAR
   B=>Mesh(nMesh)%block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   
!---------------------------------------------------------------------------- 
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)

!$OMP DO   
   do k=1,nz-1
   do j=1,ny-1
   do m=1, NVAR1
    do i=0,nx
    tmpa(i)=B%U(m,i,j,k)
    enddo
    do i=2,nx-2
     B%U(m,i,j,k)=(-tmpa(i-2)+4.d0*tmpa(i-1)+10.d0*tmpa(i)+4.d0*tmpa(i+1)-tmpa(i+2))/16.d0
    enddo
     B%U(m,1,j,k)=0.25d0*(tmpa(0)+tmpa(2))+0.5d0*tmpa(1)
     B%U(m,nx-1,j,k)=0.25d0*(tmpa(nx-2)+tmpa(nx))+0.5d0*tmpa(nx-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do k=1,nz-1
   do i=1,nx-1
   do m=1, NVAR1
   do j=0,ny
    tmpb(j)=B%U(m,i,j,k)
   enddo
    do j=2,ny-2
     B%U(m,i,j,k)=(-tmpb(j-2)+4.d0*tmpb(j-1)+10.d0*tmpb(j)+4.d0*tmpb(j+1)-tmpb(j+2))/16.d0
    enddo
     B%U(m,i,1,k)=0.25d0*(tmpb(0)+tmpb(2))+0.5d0*tmpb(1)
     B%U(m,i,ny-1,k)=0.25d0*(tmpb(ny-2)+tmpb(ny))+0.5d0*tmpb(ny-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do j=1,ny-1
   do i=1,nx-1
   do m=1, NVAR1
   do k=0,nz
    tmpc(k)=B%U(m,i,j,k)
   enddo
    do k=2,nz-2
     B%U(m,i,j,k)=(-tmpc(k-2)+4.d0*tmpc(k-1)+10.d0*tmpc(k)+4.d0*tmpc(k+1)-tmpc(k+2))/16.d0
    enddo
     B%U(m,i,j,1)=0.25d0*(tmpc(0)+tmpc(2))+0.5d0*tmpc(1)
     B%U(m,i,j,nz-1)=0.25d0*(tmpc(nz-2)+tmpc(nz))+0.5d0*tmpc(nz-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL
!-----------------------------------------------------------------------
   end  subroutine smoothing_oneBlock_4th



