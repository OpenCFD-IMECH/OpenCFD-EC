!  2017-7-11  对流场进行限制 （防止物理量超界）
!-------------------------------------------
  subroutine limit_flow(nMesh)
   use Global_var
   implicit none
   integer::nMesh
   call   limit_flow_U(nMesh)             ! Limit d,u,v,w,p
   if(Mesh(nMesh)%NVAR == 6) then 
    call limit_flow_SA (nMesh)       ! limit U(6) (usuall for SA)
   endif

   if(Mesh(nMesh)%NVAR == 7) then 
!    call limit_flow_U7       ! limit U(6) (usuall for SA)
     ! you can add your code here  !!!!!
   endif
  
  
  end

!----------------------------------------
!  限定压力及密度的变化幅度 (see CFL3D User's manual: p236, Time Advancement)
!  变化幅度超过阈值(如, -20%) 则进行特殊处理

   subroutine limit_flow_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,nx,ny,nz,i,j,k
   real(PRE_EC):: dn,un,vn,wn,pn,dd,dp
   real(PRE_EC),parameter:: alfac=-0.2d0, phic=2.d0     ! (Phic > 1)
   real(PRE_EC),dimension(:,:,:),pointer:: d1,u1,v1,w1,p1
   integer:: ia,ja,ka,i1,j1,k1,kn,Nneg,Kneg(3,10)
   real(PRE_EC):: d0,u0,v0,w0,p0, sn1
   ! alfac 限定值； phic 目标值 (p1 -> pn/phic)
   
   Type (Block_TYPE),pointer:: B

   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

     allocate(d1(nx-1,ny-1,nz-1),u1(nx-1,ny-1,nz-1),v1(nx-1,ny-1,nz-1), &
	          w1(nx-1,ny-1,nz-1),p1(nx-1,ny-1,nz-1))

!--------------------------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,gamma,d1,u1,v1,w1,p1)
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
 !        dn=B%Un(1,i,j,k)
 !        un=B%Un(2,i,j,k)/dn
 !        vn=B%Un(3,i,j,k)/dn
 !        wn=B%Un(4,i,j,k)/dn
 !        pn=(B%Un(5,i,j,k)-0.5d0*dn*(un*un+vn*vn+wn*wn)) *(gamma-1.d0)

         d1(i,j,k)=B%U(1,i,j,k)
		 u1(i,j,k)=B%U(2,i,j,k)/d1(i,j,k)
         v1(i,j,k)=B%U(3,i,j,k)/d1(i,j,k)
         w1(i,j,k)=B%U(4,i,j,k)/d1(i,j,k)
         p1(i,j,k)=(B%U(5,i,j,k)-0.5d0*d1(i,j,k)    &
		           *(u1(i,j,k)*u1(i,j,k)+v1(i,j,k)*v1(i,j,k)+w1(i,j,k)*w1(i,j,k))) *(gamma-1.d0)

!        dd=d1(i,j,k)-dn
!        dp=p1(i,j,k)-pn
     
!        if( dd/dn < alfac .or. dp/pn < alfac) then
!	      if(dd/dn < alfac)   d1(i,j,k)=dn+dd/(1.d0+phic*(alfac+abs(dd/dn)))
!		  if(dp/pn < alfac)   p1(i,j,k)=pn+dp/(1.d0+phic*(alfac+abs(dp/pn)))
! 		   B%U(1,i,j,k)=d1(i,j,k)
!		   B%U(2,i,j,k)=d1(i,j,k)*u1(i,j,k)
!		   B%U(3,i,j,k)=d1(i,j,k)*v1(i,j,k)
!		   B%U(4,i,j,k)=d1(i,j,k)*w1(i,j,k)
!		   B%U(5,i,j,k)=p1(i,j,k)/(gamma-1.d0)+   &
!		          0.5d0*d1(i,j,k)*(u1(i,j,k)**2+v1(i,j,k)**2+w1(i,j,k)**2)
!  	     endif
	 
	  enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO   
	  
	   Nneg=0  		   
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,gamma,d1,u1,v1,w1,p1,Nneg,Kneg,Ldmin,Lpmin,Ldmax,Lpmax,Lumax)
       do k=1,nz-1
       do j=1,ny-1
       do i=1,nx-1
        if(d1(i,j,k) <Ldmin .or. p1(i,j,k) < Lpmin .or. d1(i,j,k) > Ldmax .or. p1(i,j,k) > Lpmax    &
		  .or. abs(u1(i,j,k)) > Lumax .or. abs(v1(i,j,k)) > Lumax  .or. abs(w1(i,j,k)) > Lumax   ) then                 ! 物理量超限
	
	
		  Nneg=Nneg+1
          if(Nneg <=10) then        ! The location of Limited point
		   Kneg(1,Nneg)=i; Kneg(2,Nneg)=j; Kneg(3,Nneg)=k
		  endif
		  
		  d0=0.d0; u0=0.d0; v0=0.d0 ; w0=0.d0; p0=0.d0; kn=0		   		   
          do ka=-1,1
		  do ja=-1,1
		  do ia=-1,1
		    k1=k+ka ; j1=j+ja ; i1=i+ia
			if(i1 > 0 .and. i1< nx .and. j1 > 0 .and. j1< ny .and. k1>0 .and. k1<nz ) then
            if( .not. (d1(i1,j1,k1) <Ldmin .or. p1(i1,j1,k1) < Lpmin .or. d1(i1,j1,k1) > Ldmax .or. p1(i1,j1,k1) > Lpmax   &
		       .or. abs(u1(i1,j1,k1)) > Lumax .or. abs(v1(i1,j1,k1)) > Lumax  .or. abs(w1(i1,j1,k1)) > Lumax )   ) then

			d0=d0+d1(i1,j1,k1)
			u0=u0+u1(i1,j1,k1)
			v0=v0+v1(i1,j1,k1)
			w0=w0+w1(i1,j1,k1)
			p0=p0+p1(i1,j1,k1)
			kn=kn+1
		    endif
			endif
		   enddo
		   enddo
		   enddo

		   if( kn ==0 ) then         ! 周围全部为“坏点”
			 d0=d1(i,j,k)
			 u0=u1(i,j,k)
			 v0=v1(i,j,k)
			 w0=w1(i,j,k)
			 p0=p1(i,j,k)
			 
			 if( d0 < Ldmin)  d0=Ldmin
			 if( d0 > Ldmax)  d0=Ldmax
			 if( p0 < Lpmin)  p0=Lpmin
			 if( p0 > Lpmax)  p0=Lpmax
			 if(abs(u0) > Lumax ) u0=sign(Lumax,u0)
			 if(abs(v0) > Lumax ) v0=sign(Lumax,v0)
			 if(abs(w0) > Lumax ) w0=sign(Lumax,w0)
		 
		   else
             sn1=1.d0/(kn)
             d0=d0*sn1
			 u0=u0*sn1
			 v0=v0*sn1
			 w0=w0*sn1
			 p0=p0*sn1
            endif

		   B%U(1,i,j,k)=d0
		   B%U(2,i,j,k)=d0*u0
		   B%U(3,i,j,k)=d0*v0
		   B%U(4,i,j,k)=d0*w0
		   B%U(5,i,j,k)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0+w0*w0)
       endif
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO   
    
   if(Nneg > 0) then
     print*, "------------------------------------------------------"
	 print*, "Limters ..."
	 print*, Ldmin,Ldmax,Lpmin,Lpmax,Lumax
	 print*, "Warning !!! Mesh, Block=", nMesh, B%block_no ,  "has", Nneg, " Limitted points !!!!" 
	 do k=1,Min(Nneg,10)
	 print*, (Kneg(i,k),i=1,3)
     i1=Kneg(1,k); j1=Kneg(2,k); k1=Kneg(3,k)
	 print*, d1(i1,j1,k1),u1(i1,j1,k1),v1(i1,j1,k1),w1(i1,j1,k1),p1(i1,j1,k1)
	 enddo

     open(101,file="error.log",position="append" )
     write(101,*) "------------------------------------------------------"
	 write(101,*) "Warning !!! Mesh, Block=", nMesh, B%block_no ,  "has", Nneg, " Limitted points !!!!" 
	 do k=1,Min(Nneg,10)
	 write(101,*) (Kneg(i,k),i=1,3)
	 enddo
     close(101)

     B%IF_OverLimit=1                ! 设定物理量超限标志
   endif
    
   deallocate(d1,u1,v1,w1,p1)
  enddo
  end subroutine limit_flow_U





   subroutine limit_flow_SA(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,nx,ny,nz,i,j,k,kn,ia,ja,ka,i1,j1,k1
   real(PRE_EC)::  s0
   real(PRE_EC),parameter:: SAmin=1.d-8
    
   Type (Block_TYPE),pointer:: B

   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

!--------------------------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,LSAmax)
	
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       if(B%U(6,i,j,k) < SAmin) B%U(6,i,j,k)=SAmin
       if(B%U(6,i,j,k) > LSAmax ) then
	  	 
   	    s0=0.d0; kn=0		   		   
         do ka=-1,1
	     do ja=-1,1
	     do ia=-1,1
	       k1=k+ka ; j1=j+ja ; i1=i+ia
	       if(i1 > 0 .and. i1< nx .and. j1 > 0 .and. j1< ny .and. k1>0 .and. k1<nz ) then
           if( B%U(6,i1,j1,k1) >=SAmin .and. B%U(6,i1,j1,k1) <= LSAmax  ) then
			s0=s0+B%U(6,i1,j1,k1)
			kn=kn+1
		    endif
			endif
		  enddo
		  enddo
		  enddo

		  if( kn ==0 ) then         ! 周围全部为“坏点”
            B%U(6,i,j,k)=LSAmax			 
		  else
		    B%U(6,i,j,k)=s0/kn
          endif
        endif
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO   
  enddo
  end subroutine limit_flow_SA
