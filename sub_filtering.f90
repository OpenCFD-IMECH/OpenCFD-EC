!  filtering data 
! Revised by Li Xinliang, 2013-10-4
  module filting_Var
   use precision_EC
   real(PRE_EC), save,pointer,dimension(:,:,:,:)::  f,f0 ! 变量
  end module Filting_Var
   
! 高精度滤波（4阶精度）
  subroutine Filtering_oneMesh(nMesh)     
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1
   Type (Block_TYPE),pointer:: B

   if(my_id .eq. 0) print*, "filtering ......"
 
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
    allocate(f(NVAR1,nx,ny,nz),f0(NVAR1,nx,ny,nz))
    
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
         f(1,i,j,k)= B%U(1,i,j,k)                 ! d
         f(2,i,j,k)= B%U(2,i,j,k)/B%U(1,i,j,k)    ! u
         f(3,i,j,k)= B%U(3,i,j,k)/B%U(1,i,j,k)    ! v
         f(4,i,j,k)= B%U(4,i,j,k)/B%U(1,i,j,k)    ! w
         f(5,i,j,k)=(B%U(5,i,j,k)-0.5d0*B%U(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2))*(gamma-1.d0)  !p
 !        f(5,i,j,k)=(B%U(5,i,j,k)-0.5d0*B%U(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2))/(Cv*f(1,i,j,k))  !T
 
        do m=6,NVAR1
         f(m,i,j,k)=B%U(m,i,j,k)
        enddo
	   
	   enddo
	   enddo
	   enddo
       


       call filter_x3d(nMesh,mBlock)
       call filter_y3d(nMesh,mBlock)
       call filter_z3d(nMesh,mBlock)

      do k=1,nz
      do j=1,ny
	  do i=1,nx
      B%U(1,i,j,k)=f(1,i,j,k)
      B%U(2,i,j,k)=f(1,i,j,k)*f(2,i,j,k)
      B%U(3,i,j,k)=f(1,i,j,k)*f(3,i,j,k)
      B%U(4,i,j,k)=f(1,i,j,k)*f(4,i,j,k)
      B%U(5,i,j,k)=f(5,i,j,k)/(gamma-1.d0)+0.5d0*f(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2)
!      B%U(5,i,j,k)=Cv*f(1,i,j,k)*f(5,i,j,k)+0.5d0*f(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2)
       do m=6,NVAR1
	    B%U(m,i,j,k)=f(m,i,j,k)
	   enddo
	  enddo
      enddo
      enddo

      deallocate(f,f0)
    enddo

  end 

!---------------------------------------------------

  subroutine filter_x3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=nx-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=1,nz
        do j=1,ny
        do i=i1,i2
            p1=dabs(f0(5,i+1,j,k)-f0(5,i  ,j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i  ,j,k)-f0(5,i-1,j,k))
		    
			alpha=0.1d0*min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i+1,j,k)+f0(m,i-1,j,k))+a3*(f0(m,i+3,j,k)+f0(m,i-3,j,k)))
             f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i+1,j,k)+f0(m,i-1,j,k))-1.d0/64.d0*(f0(m,i+2,j,k)+f0(m,i-2,j,k)))

           enddo
		enddo
        enddo
        enddo
  
      end



!--------------------------------------------------------------------------

  subroutine filter_y3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=ny-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=1,nz
        do j=i1,i2
        do i=1,nx
            p1=dabs(f0(5,i,j+1,k)-f0(5,i,  j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i,j,  k)-f0(5,i,  j-1,k))
		 
		    alpha=0.1d0*min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i,j+1,k)+f0(m,i,j-1,k))+a3*(f0(m,i,j+3,k)+f0(m,i,j-3,k)))
             f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i,j+1,k)+f0(m,i,j-1,k))-1.d0/64.d0*(f0(m,i,j+2,k)+f0(m,i,j-2,k)))

		   enddo
		enddo
        enddo
        enddo
  
      end


!--------------------------------------------------------------------------

  subroutine filter_z3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=nz-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=i1,i2    
        do j=1,ny
        do i=1,nx
            p1=dabs(f0(5,i,j,k+1)-f0(5,i,  j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i,j,  k)-f0(5,i,  j,k-1))
		    alpha=0.1d0* min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i,j,k+1)+f0(m,i,j,k-1))+a3*(f0(m,i,j,k+3)+f0(m,i,j,k-3)))
		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i,j,k+1)+f0(m,i,j,k-1))-1.d0/64.d0*(f0(m,i,j,k+2)+f0(m,i,j,k-2)))

           enddo
		enddo
        enddo
        enddo
  
      end


