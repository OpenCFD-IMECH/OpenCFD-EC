!-----------------------------------------------------
    subroutine boundary_user(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
	 integer:: mBlock,ksub,nMesh
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
   Select case  (Bc%bc)
     case ( BC_USER_FixedInlet )    ! 901, 给定入口流动
       call boundary_user_Inlet(nMesh,mBlock,ksub)
     case ( BC_USER_Inlet_time )    ! 902, 给定入口流动时间序列
       call  boundary_user_Inlet_time(nMesh,mBlock,ksub)
     case ( BC_USER_Blow_Suction_Wall )    ! 903, 吹吸壁面
       call boundary_user_blow_suction_Wall(nMesh,mBlock,ksub)

	 case default
	  
	  print*, "This USER defined boundary is not supported !!!"
	  stop

   end Select
   end




!------------------用户自定义的边界条件-----------------------------------
! 给定入口分布，文件名 inlet.xxx  (xxx为块号)
! 仅使用1层虚网格

    subroutine boundary_user_Inlet(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc

     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1
     integer:: n1,n2,n3
	 integer,save:: read_flag=0
     real(PRE_EC),allocatable,save,dimension(:,:,:,:):: uc,fc       ! 边界值 
     character(len=30):: fname

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke    
     n1=max(ie-ib,1)
	 n2=max(je-jb,1)
	 n3=max(ke-kb,1)

 !-----读取数据--------------------------------------   
	 if(read_flag==0) then
	  read_flag=1
      allocate(uc(n1,n2,n3,5),fc(5,n1,n2,n3))
      write(fname,"('inlet.'I3.3)") B%Block_no
	  open(99,file=fname,form="unformatted")
	  print*, "read file : ", fname
	  print*, "n1,n2,n3=",n1,n2,n3
	  do m=1,5
	   read(99,*) (((uc(i,j,k,m),i=1,n1),j=1,n2),k=1,n3)
      enddo
	  close(99)
	 endif

        do k=1,n3
		do j=1,n2
		do i=1,n1
          fc(1,i,j,k)=uc(i,j,k,1)
          fc(2,i,j,k)=uc(i,j,k,1)*uc(i,j,k,2)
          fc(3,i,j,k)=uc(i,j,k,1)*uc(i,j,k,3)
          fc(4,i,j,k)=uc(i,j,k,1)*uc(i,j,k,4)
          fc(5,i,j,k)=uc(i,j,k,1)*(Cv*uc(i,j,k,5)+0.5d0*(uc(i,j,k,2)**2+uc(i,j,k,3)**2+uc(i,j,k,4)**2))
        enddo
	    enddo
	    enddo

!-----------------------------------------------------------------

 ! 仅使用1层Ghost Cell    
             if(Bc%face .eq. 1) then                 ! i- 面
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5 
				   B%U(m,0,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 2) then
 				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,0,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 3) then             
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,0)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
  
             else if(Bc%face .eq. 4) then 
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5
				   B%U(m,ie,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo

             else if(Bc%face .eq. 5) then
				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,je,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
   
             else if(Bc%face .eq. 6) then
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,ke)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
             endif

     end 
!----------------------------------------------------------------------------------------

!------------------用户自定义的边界条件-----------------------------------
! 给定入口时间序列，文件名 inlet-time.xxx  (xxx为块号)
! 仅使用1层虚网格

    subroutine boundary_user_Inlet_time(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc

     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1
     integer:: n1,n2,n3,Istep_inlet
	 integer,save:: read_flag=0
     real(PRE_EC),allocatable,save,dimension(:,:,:,:):: uc,fc       ! 边界值 
     real(PRE_EC):: tt_inlet
	 character(len=30):: fname

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke    
     n1=max(ie-ib,1)
	 n2=max(je-jb,1)
	 n3=max(ke-kb,1)

 !-----读取数据--------------------------------------   
	 if(read_flag==0) then
	  read_flag=1
      allocate(uc(n1,n2,n3,5),fc(5,n1,n2,n3))
      write(fname,"('inlet-time.'I3.3)") B%Block_no
	  open(999,file=fname,form="unformatted")
	  print*, "read file : ", fname
	  print*, "n1,n2,n3=",n1,n2,n3
    endif

! 初始时刻 (KRK=0) 以及 更新时间步  (KRK=3) 时， 读取 uc数据
! 如采用1步 Euler法， 则KRK保持0

   if(KRK ==0 .or. KRK==3) then
      read(999) Istep_inlet,tt_inlet
      print*, "Istep_inlet,tt_inlet=", Istep_inlet, tt_inlet
	  do m=1,5
	   read(999) (((uc(i,j,k,m),i=1,n1),j=1,n2),k=1,n3)
      enddo

        do k=1,n3
		do j=1,n2
		do i=1,n1
          fc(1,i,j,k)=uc(i,j,k,1)
          fc(2,i,j,k)=uc(i,j,k,1)*uc(i,j,k,2)
          fc(3,i,j,k)=uc(i,j,k,1)*uc(i,j,k,3)
          fc(4,i,j,k)=uc(i,j,k,1)*uc(i,j,k,4)
          fc(5,i,j,k)=uc(i,j,k,1)*(Cv*uc(i,j,k,5)+0.5d0*(uc(i,j,k,2)**2+uc(i,j,k,3)**2+uc(i,j,k,4)**2))
        enddo
	    enddo
	    enddo
    endif

!-----------------------------------------------------------------

 ! 仅使用1层Ghost Cell    
             if(Bc%face .eq. 1) then                 ! i- 面
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5 
				   B%U(m,0,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 2) then
 				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,0,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 3) then             
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,0)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
  
             else if(Bc%face .eq. 4) then 
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5
				   B%U(m,ie,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo

             else if(Bc%face .eq. 5) then
				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,je,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
   
             else if(Bc%face .eq. 6) then
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,ke)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
             endif

     end 
!----------------------------------------------------------------------------------------


!-------------------------------------------------------------------  
! 包含吹吸扰动的壁面 （扰动范围x-z方向）
! 仅使用1层虚网格


    subroutine boundary_user_blow_suction_Wall(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1,i1,j1,k1,i2,j2,k2,nk,ierr
     real(PRE_EC):: n1,n2,n3,Vn
	 integer,save:: Iflag=0
     real(PRE_EC):: d1,u1,v1,w1,T1,p1,d2,u2,v2,w2,T2,p2,xb,xe,zb,ze,wx,wz,Aw
     real(PRE_EC),save:: rpara(100)
! --------- 仅在第1次调用本子程序时执行 ------------- 
!   读取辅助参数 (供自定义边界条件使用)
   if(Iflag==0) then
     Iflag=1
        open(99,file="bc_user.ec")
		read(99,*) 
		read(99,*)
		read(99,*) nk
		read(99,*)
		read(99,*) (rpara(k),k=1,nk)
        close(99)
      print*, "my_id=",  my_id, "rpara=", (rpara(k),k=1,nk)
   endif
! --------------------------------------------
	 
	 NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      

    if(NVAR1 > 5) then
	  print*, "Warning, NVAR1 should =5 !!!"
	  print*, "if you using blow-and-suction wall, you should not use Turbulence model"
	endif

! 壁面吹吸扰动区域、波数及振幅    
	xb=rpara(1)
	xe=rpara(2)
	zb=rpara(3)
	ze=rpara(4)
	wx=rpara(5)
	wz=rpara(6)
	Aw=rpara(7)
!-------------------------

!   仅使用1层Ghost Cell
!   i2 Ghost Cell;  i1 内点
     if(Bc%face .eq. 1 .or. Bc%face .eq. 4) then   ! i- or i+
       if(Bc%face .eq. 1) then
         i=ib; i1=ib; i2=ib-1
       else
         i=ie; i1=ie-1; i2=ie 
       endif

       do k=kb,ke-1
         do j=jb,je-1
    

		    d1=B%U(1,i1,j,k); u1=B%U(2,i1,j,k)/d1; v1=B%U(3,i1,j,k)/d1; w1=B%U(4,i1,j,k)/d1
            p1=(B%U(5,i1,j,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
			p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i1,j,k) >= xb  .and. B%xc(i1,j,k) <=xe .and. B%zc(i1,j,k) >= zb .and. B%zc(i1,j,k) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i1,j,k)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i1,j,k)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

		    n1=B%ni1(i,j,k) ; n2=B%ni2(i,j,k); n3=B%ni3(i,j,k)   ! 归一化法方向  
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i2,j,k)= d2
           B%U(2,i2,j,k)= d2*u2
           B%U(3,i2,j,k)= d2*v2
           B%U(4,i2,j,k)= d2*w2
           B%U(5,i2,j,k)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
		 
		 enddo
       enddo

     else if(Bc%face .eq. 2  .or. Bc%face .eq. 5) then   ! j- or j+
       if(Bc%face .eq. 2) then
         j=jb; j1=jb; j2=jb-1 
       else
         j=je; j1=je-1; j2=je 
       endif

       do k=kb,ke-1
         do i=ib,ie-1
 
 		    d1=B%U(1,i,j1,k); u1=B%U(2,i,j1,k)/d1; v1=B%U(3,i,j1,k)/d1; w1=B%U(4,i,j1,k)/d1
            p1=(B%U(5,i,j1,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
    		p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i,j1,k) >= xb  .and. B%xc(i,j1,k) <=xe .and. B%zc(i,j1,k) >= zb .and. B%zc(i,j1,k) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i,j1,k)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i,j1,k)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

            n1=B%nj1(i,j,k) ; n2=B%nj2(i,j,k); n3=B%nj3(i,j,k)     
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i,j2,k)= d2
           B%U(2,i,j2,k)= d2*u2
           B%U(3,i,j2,k)= d2*v2
           B%U(4,i,j2,k)= d2*w2
           B%U(5,i,j2,k)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)

          enddo
       enddo

     else if(Bc%face .eq. 3 .or. Bc%face .eq. 6) then   ! k- or k+
       if(Bc%face .eq. 3) then
         k=kb; k1=kb; k2=kb-1 
       else
         k=ke; k1=ke-1; k2=ke
       endif

       do j=jb,je-1
         do i=ib,ie-1

 		    d1=B%U(1,i,j,k1); u1=B%U(2,i,j,k1)/d1; v1=B%U(3,i,j,k1)/d1; w1=B%U(4,i,j,k1)/d1
            p1=(B%U(5,i,j,k1)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
    		p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i,j,k1) >= xb  .and. B%xc(i,j,k1) <=xe .and. B%zc(i,j,k1) >= zb .and. B%zc(i,j,k1) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i,j,k1)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i,j,k1)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

		    n1=B%nk1(i,j,k) ; n2=B%nk2(i,j,k); n3=B%nk3(i,j,k)     
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i,j,k2)= d2
           B%U(2,i,j,k2)= d2*u2
           B%U(3,i,j,k2)= d2*v2
           B%U(4,i,j,k2)= d2*w2
           B%U(5,i,j,k2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
		 enddo
       enddo

     endif
        
    end 




