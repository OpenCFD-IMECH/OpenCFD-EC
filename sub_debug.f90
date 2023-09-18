   subroutine Output_mesh_debug
   use Global_var
   implicit none
   integer:: m,nx,ny,nz,i,j,k
   Type (Block_TYPE),pointer:: B
   character(len=50):: filename
   write(filename,"('Meshdeb-'I5.5'.dat')") my_id
   open(50,file=filename)
    
!  输出网格文件 tecplot格式；   
   write(50,*) "variables=x,y,z"
   do m=1,Mesh(1)%Num_Block
     B=>Mesh(1)%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
     write(50,*) "zone i=",nx+2, " j= ", ny+2, " k= ", nz+2
	 do k=0,nz+1
	 do j=0,ny+1
	 do i=0,nx+1
	 write(50,"(3f20.10)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k)
	 enddo
	 enddo
	 enddo
	enddo
	close(50)
   end
 

 
   subroutine debug1
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,m,i,j,k
   real*8:: u1,u2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    u1=0.d0
	u2=0.d0
	do m=1,MP%Num_Block
    B=>MP%Block(m)
	nx=B%nx; ny=B%ny; nz=B%nz
	do k=1-LAP,nz+LAP-1
	do j=1-LAP,ny+LAP-1
	do i=1-LAP,nx+LAP-1
	 u1=u1+B%U(1,i,j,k)**2
	 u2=u2+B%U(5,i,j,k)**2
	enddo
	enddo
	enddo
	enddo
	print*, "************ u1,u2=",u1,u2
	end



   subroutine debug2(nx,ny,nz)
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,i,j,k
   real*8:: f1,f2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    f1=0.d0
	f2=0.d0
	do k=1,nz
	do j=1,ny
	do i=1,nx
	 f1=f1+Flux(1,i,j,k)**2
     f2=f2+Flux(5,i,j,k)**2
	enddo
	enddo
	enddo

	print*, "************ f1,f2=",f1,f2
	end
