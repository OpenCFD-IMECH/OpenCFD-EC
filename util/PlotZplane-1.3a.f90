!-----------------------------------------------------------
! Copyright by LiXinliang
! Ver 1.2; Plot 3D data in a given range
! Ver 1.3a: read flow3d.dat or flow3d-average.dat
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
  end module Const_Variables

 
  module Global_Variables
  use Const_Variables
  implicit none

   TYPE Block_TYPE           !  variables for each block 
     integer :: Block_no,nx,ny,nz
	 real(PRE_EC),pointer,dimension(:,:,:)::  x,y,z
	 real(PRE_EC),pointer,dimension(:,:,:)::  xc,yc,zc
	 real(PRE_EC),pointer,dimension(:,:,:) :: dc,uc,vc,wc,Tc,pc
	 real(PRE_EC),pointer,dimension(:,:,:) :: d,u,v,w,T,p

   End TYPE Block_TYPE  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
   integer:: Num_Block, Mesh_File_Format
   integer:: Zn, Zplane(1000,2),If_X,If_K  
  end module Global_Variables
!------------------------------------------------------------------




!-------------------------------------------
     program choose_slit   
     use Global_Variables
     implicit none
     integer:: If3d   
       call read_flow
       call comput_value_in_mesh
	   call choose_k
       call Plot_Zplane	   	 
       call plot_cf
       call Plot_Cf2d	   	 

	   call Plot_3D	   	 

	 end

!=====================================================================================
   subroutine read_flow
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,Num_Block1,ksub,tmp,NB1,Iflag
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B

! ---------node Coordinates----------------------------------------  
   print*, "read Mesh3d.dat... (PLOT3D Format)"
   print*, " please input the format, 0: unformatted, 1 formatted"
   read(*,*) Mesh_File_Format

   if(Mesh_File_Format .eq. 0) then
    open(99,file="Mesh3d.dat",form="unformatted")
    read(99) Num_Block         
   else
    open(99,file="Mesh3d.dat")
    read(99,*) Num_Block        
   endif

      
  allocate(Block(Num_Block))             
  allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   
  if(Mesh_File_Format .eq. 0) then
   read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  else
   read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  endif

! ------read mesh----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   
     nx=B%nx ; ny= B%ny ; nz=B%nz

    allocate(B%x(nx,ny,nz),B%y(nx,ny,nz),B%z(nx,ny,nz))
    allocate(B%xc(0:nx,0:ny,0:nz), B%yc(0:nx,0:ny,0:nz), B%zc(0:nx,0:ny,0:nz)) 

    allocate(B%dc(0:nx,0:ny,0:nz),B%uc(0:nx,0:ny,0:nz),B%vc(0:nx,0:ny,0:nz), &
             B%wc(0:nx,0:ny,0:nz),B%Tc(0:nx,0:ny,0:nz),B%pc(0:nx,0:ny,0:nz))
    allocate(B%d(nx,ny,nz),B%u(nx,ny,nz),B%v(nx,ny,nz), &
             B%w(nx,ny,nz),B%T(nx,ny,nz),B%p(nx,ny,nz))

!------read mesh ---------------------------------------
   if(Mesh_File_Format .eq. 0) then
	read(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else
  	read(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif
 
!----------插值出网格中心点的数据-----------------------------------
   do k=1,B%nz-1
   do j=1,B%ny-1
   do i=1,B%nx-1
     B%xc(i,j,k)=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j+1,k)+        &
	              B%x(i,j,k+1)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1) + &
                  B%x(i,j+1,k+1)+B%x(i+1,j+1,k+1))/8.d0
     B%yc(i,j,k)=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j+1,k)+        &
	              B%y(i,j,k+1)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1) + &
                  B%y(i,j+1,k+1)+B%y(i+1,j+1,k+1))/8.d0
     B%zc(i,j,k)=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j+1,k)+        &
	              B%z(i,j,k+1)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1) + &
                  B%z(i,j+1,k+1)+B%z(i+1,j+1,k+1))/8.d0
   
   enddo
   enddo
   enddo

  enddo

  close(99)
  deallocate(NI,NJ,NK)

!-------Read flow ----------------------------
    print*, "read flow3d.dat  or flow3d_average.dat ..."
	print*, " input 1 for read flow3d.dat, 2 for flow3d-agerage.dat"
    read(*,*) Iflag
	if(Iflag == 1) then
     open(99,file="flow3d.dat",form="unformatted")
    else
     open(99,file="flow3d_average.dat",form="unformatted")
    endif

	 do m=1,Num_Block
     B => Block(m)
        read(99)   (((B%dc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%uc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%vc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%wc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
		           (((B%Tc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 
          do k=0, B%nz
		  do j=0, B%ny
		  do i=0, B%nx
		    B%pc(i,j,k)=B%dc(i,j,k)*B%Tc(i,j,k)
		  enddo
		  enddo
		  enddo

	 enddo
	close(99)

 end 

!----------------------------------------------------
   subroutine choose_k
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nn,k0
   real(PRE_EC):: X0, XT , dx, dxmin
   real(PRE_EC),allocatable,dimension(:):: Xmean
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),pointer,dimension(:,:,:):: Cx   ! x, y or z 

   print*, "choose a plane with specified x, y, or z cordiation ..."
   print*, " Please input:  1 for x,  2 for y, 3 for z"
   read(*,*) If_x
   print*, "Please input X (or Y, Z) location X0, and the toleration XT"
   read(*,*) X0, XT
   print*, "please input the search plane: 1 for i-; 2 for j-; 3 for k- plane"
   read(*,*) If_k
   print*, " -------search mesh3d ......"
   
   open(99,file="Zplane.dat")
   write(99,*) If_K
   write(99,*) "---------------------------"
    Zn=0
   do m=1,Num_Block
     B => Block(m)
     nx=B%nx ; ny= B%ny ; nz=B%nz

     if(If_x == 1) then
	  Cx => B%x
	 else if( If_x ==2 ) then
	  Cx => B%y
	 else if (If_x == 3) then
	  Cx => B%z
     endif


    if(If_k ==1 ) then   ! i-
	   nn=nx
	   allocate(Xmean(nn))   
      Xmean=0.d0
   
	  do k=1, nz
	  do j=1, ny
	  do i=1,nx
	   Xmean(i)=Xmean(i)+Cx(i,j,k)
	  enddo
	  enddo
	  enddo
	  Xmean=Xmean/(1.d0*ny*nz)
    
	else if(If_k ==2 ) then   ! i-
	   nn=ny
	  allocate(Xmean(nn))   
      Xmean=0.d0
   
	  do k=1, nz
	  do j=1, ny
	  do i=1, nx
	   Xmean(j)=Xmean(j)+Cx(i,j,k)
	  enddo
	  enddo
	  enddo
	  Xmean=Xmean/(1.d0*nx*nz)

	else   ! k-
	  nn=nz
	  allocate(Xmean(nn))   
      Xmean=0.d0
   
	  do k=1, nz
	  do j=1, ny
	  do i=1, nx
	   Xmean(k)=Xmean(k)+Cx(i,j,k)
	  enddo
	  enddo
	  enddo
	  Xmean=Xmean/(1.d0*nx*ny)
    endif
        
		dxmin=abs(Xmean(1)-X0)
        k0=1
	   do k=1,nn
        dx=abs(Xmean(k)-X0)	   
        if(dx < dxmin) then
		 dxmin=dx
		 k0=k
		endif
       enddo

	   print*, " ------- Block : ", m
	   print*, " Dxmin=", dxmin, " k= ", k0
	   if(dxmin <= XT) then
	   Zn=Zn+1
	   Zplane(Zn,1)=m
	   Zplane(Zn,2)=k0

	   write(99,*) m,k0
	   endif

	   deallocate(Xmean)
	 
	 enddo
   close(99)
  end


   subroutine comput_value_in_mesh
   use  Global_Variables
   implicit none
   integer:: nx,ny,nz,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do m=1, NUM_BLOCK
   B=> Block(m)
   nx=B%nx; ny=B%ny; nz=B%nz
!  将格心处的物理量插值到格点上   
   do k=1,B%nz
   do j=1,B%ny
   do i=1,B%nx
     B%d(i,j,k)=(B%dc(i-1,j-1,k-1)+B%dc(i,j-1,k-1)+B%dc(i-1,j,k-1)+B%dc(i,j,k-1)  &
	             +B%dc(i-1,j-1,k)+B%dc(i,j-1,k)+B%dc(i-1,j,k)+B%dc(i,j,k))/8.d0
     B%u(i,j,k)=(B%uc(i-1,j-1,k-1)+B%uc(i,j-1,k-1)+B%uc(i-1,j,k-1)+B%uc(i,j,k-1)  &
	             +B%uc(i-1,j-1,k)+B%uc(i,j-1,k)+B%uc(i-1,j,k)+B%uc(i,j,k))/8.d0
     B%v(i,j,k)=(B%vc(i-1,j-1,k-1)+B%vc(i,j-1,k-1)+B%vc(i-1,j,k-1)+B%vc(i,j,k-1)  &
	             +B%vc(i-1,j-1,k)+B%vc(i,j-1,k)+B%vc(i-1,j,k)+B%vc(i,j,k))/8.d0
     B%w(i,j,k)=(B%wc(i-1,j-1,k-1)+B%wc(i,j-1,k-1)+B%wc(i-1,j,k-1)+B%wc(i,j,k-1)  &
	             +B%wc(i-1,j-1,k)+B%wc(i,j-1,k)+B%wc(i-1,j,k)+B%wc(i,j,k))/8.d0

     B%T(i,j,k)=(B%Tc(i-1,j-1,k-1)+B%Tc(i,j-1,k-1)+B%Tc(i-1,j,k-1)+B%Tc(i,j,k-1)  &
	             +B%Tc(i-1,j-1,k)+B%Tc(i,j-1,k)+B%Tc(i-1,j,k)+B%Tc(i,j,k))/8.d0

     B%p(i,j,k)=(B%pc(i-1,j-1,k-1)+B%pc(i,j-1,k-1)+B%pc(i-1,j,k-1)+B%pc(i,j,k-1)  &
	             +B%pc(i-1,j-1,k)+B%pc(i,j-1,k)+B%pc(i-1,j,k)+B%pc(i,j,k))/8.d0


   enddo
   enddo
   enddo
   enddo
   end


   subroutine   Plot_Zplane	   	 

   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nk
   Type (Block_TYPE),pointer:: B
   open(100,file="Zplane-flow.dat")
   write(100,*) "variables=x,y,z,d,u,v,w,T,p"
   do nk=1,Zn
    m=Zplane(nk,1)
	B=>Block(m)
	nx=B%nx; ny=B%ny; nz=B%nz
     if(If_k == 1) then    !i- plane
	 write(100,*) "zone  i= ", ny, " j= ", nz
	 i=Zplane(nk,2)
     
	  do k=1,nz
	  do j=1,ny
	   write(100,"(9E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                         B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), B%p(i,j,k)
      enddo
	  enddo
     else if(If_k == 2) then    ! j- plane
	 write(100,*) "zone  i= ", nx, " j= ", nz

	  j=Zplane(nk,2)
	  do k=1,nz
	  do i=1,nx
	   write(100,"(9E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                         B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), B%p(i,j,k)
      enddo
	  enddo
    else
	 write(100,*) "zone  i= ", nx, " j= ", ny

	  k=Zplane(nk,2)
	  do j=1,ny
	  do i=1,nx
	   write(100,"(9E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                         B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), B%p(i,j,k)
      enddo
	  enddo
    endif
  enddo
   close(100)
   end


! 计算k-截面的摩阻，假设壁面位于y=0处
   subroutine   Plot_Cf	   	 

   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nk
   real*8,parameter:: ylim=0.1d0
   Type (Block_TYPE),pointer:: B
   real*8:: Re=4.16d3,Tinf=222.65d0,Tw=600.d0
   real*8:: Tsb,mu
       Tw=Tw/Tinf
       Tsb=110.4d0/Tinf
       mu=1.d0/Re*(1.d0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
       print*, "mu=",mu
   If(If_k .ne. 3) return

   open(100,file="Cf.dat")
   write(100,*) "variables=x,Cf,u"
   do nk=1,Zn
    m=Zplane(nk,1)
	B=>Block(m)
	nx=B%nx; ny=B%ny; nz=B%nz
     if(If_k == 3 .and. B%y(1,1,1) .le. ylim) then    !k- plane
 	  write(100,*) "zone  i= ", nx
	  k=Zplane(nk,2)
	  do i=1,nx
	   write(100,"(3E18.9)") B%x(i,1,k),2.d0*mu*B%u(i,2,k)/B%y(i,2,k),B%u(i,2,k)
	  enddo
    endif
  enddo
   close(100)
   end



   subroutine   Plot_Cf2d	   	 
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nk
   real*8,parameter:: ylim=0.1d0
   Type (Block_TYPE),pointer:: B
   real*8:: Re=4.16d3,Tinf=222.65d0,Tw=600.d0
   real*8:: Tsb,mu,Cf,ya, Fw1,Fw2,Dx
   real*8,parameter:: x1=500, x2=503.25
      

       Tw=Tw/Tinf
       Tsb=110.4d0/Tinf
       mu=1.d0/Re*(1.d0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
       print*, "mu=",mu
       Fw1=0.d0
       Fw2=0.d0
       print*, "Comput total fraction force in the forward or backward of the block"
	   print*, "Please input the length DX "
       read(*,*) Dx

   open(100,file="Cf2d.dat")
   write(100,*) "variables=x,y,z,Cf"
   do m=1,Num_Block
     B => Block(m)
     nx=B%nx ; ny= B%ny ; nz=B%nz
     ya=0.d0
	 do k=1,nz
	 do i=1,nx
	   ya=ya+B%y(i,1,k)
	 enddo
	 enddo
       ya=ya/(nx*nz)
	   if(ya .le. ylim) then
 	    write(100,*) "zone  i= ", nx , " j= ", nz
		do k=1,nz
		do i=1,nx
		  Cf=2.d0*mu*B%u(i,2,k)/B%y(i,2,k)
		  write(100,"(4E18.9)") B%x(i,1,k),B%y(i,1,k),B%z(i,1,k),Cf
		enddo
		enddo
       do k=1,nz-1
	   do i=1,nx-1
	   Cf=2.d0*mu*B%u(i,2,k)/B%y(i,2,k)
	   if(B%x(i,1,k) >= x1-Dx .and. B%x(i,1,k) <= x1) then
	     Fw1=Fw1+Cf*(B%x(i+1,1,k)-B%x(i,1,k))*(B%z(i,1,k+1)-B%z(i,1,k))
	   endif
       if(B%x(i,1,k) >=x2 .and. B%x(i,1,k)<=x2+Dx) then
	     Fw2=Fw2+Cf*(B%x(i+1,1,k)-B%x(i,1,k))*(B%z(i,1,k+1)-B%z(i,1,k))
	   endif 
       enddo
	   enddo
      endif
    enddo
   close(100)
    print*, "Fw1=", Fw1, "  Fw2=", Fw2
   end




   subroutine   Plot_3D	   	 
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nk
   integer:: If_3d, If_zone
   real*8:: xb,xe,yb,ye,zb,ze
   Type (Block_TYPE),pointer:: B
   print*, "Do you want Plot 3D flow (1 Yes, 0 No)?"
   read(*,*) If_3D
   if(If_3d ==0 ) return
   print*, "Please input the 3D range:  xb,xe,yb,ye,zb,ze "
   read(*,*) xb,xe,yb,ye,zb,ze

   open(100,file="flow3d-tec.dat")
   write(100,*) "variables=x,y,z,d,u,v,w,T,p"

   do m=1, NUM_BLOCK
   B=> Block(m)
   nx=B%nx; ny=B%ny; nz=B%nz
   if_zone=0
! test whether Block m is in the range
   do k=1,B%nz
   do j=1,B%ny
   do i=1,B%nx
       if( B%x(i,j,k) .ge. xb  .and. B%x(i,j,k) .le. xe .and. &
           B%y(i,j,k) .ge. yb  .and. B%y(i,j,k) .le. ye .and. &
           B%z(i,j,k) .ge. zb  .and. B%z(i,j,k) .le. ze  ) then
		  If_zone=1
		  goto 120
	   endif
   enddo
   enddo
   enddo
 120 continue
  if( If_zone ==1) then
   write(100,*) "zone ", " i= ",B%nx, " j= ", B%ny, " k= ", B%nz
   do k=1,B%nz
   do j=1,B%ny
   do i=1,B%nx
    write(100,"(9E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
                         B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), B%p(i,j,k) 
   enddo
   enddo
   enddo
  endif
  enddo
  end
