!  采用LU-SGS方法，计算DU=U(n+1)-U(n)
!  Code by Li Xinliang, 2011-12-29
!-----------------------------------------------------------------------------------------
    subroutine  du_LU_SGS(nMesh,mBlock,Sfac1)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
    use Global_Var
    use Flow_Var 
    implicit none
	integer:: nMesh,mBlock,NV,nx,ny,nz,plane,i,j,k,m
    real(PRE_EC),dimension(7)::alfa,dui,duj,duk,DF
    Type (Block_TYPE),pointer:: B
    real(PRE_EC):: Sfac1    ! 双时间步时使用

     NV=Mesh(nMesh)%NVAR
	 B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny; nz=B%nz

! LU-SGS的两次扫描
!----------------------------------
!   从i=1,j=1,k=1 到i=nx-1,j=ny-1,k=nz-1的扫描过程  (向上扫描过程)
!   扫描 i+j+k=plane 的平面
!   w_LU是松弛因子（1到2之间），增大w_LU会提高稳定性，但会降低收敛速度
   do plane=3,nx+ny+nz-3            

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(plane,nx,ny,nz,NV,B,Lci,Lcj,Lck,Lvi,Lvj,Lvk,gamma,w_LU,If_viscous)
     do k=1,nz-1
	 do j=1,ny-1
	 i=plane-k-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle    ! 超出了这个平面
       alfa(1:NV)=B%vol(i,j,k)/B%dt(i,j,k)+w_LU*(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k))     ! 对角线项
     if(If_viscous .eq. 1)  then
	   alfa(1:NV)=alfa(1:NV)+2.d0*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k))          
       if(NV .eq. 6) then
	     alfa(6)=alfa(6)+(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k))      ! 再增大些对角线 （考虑到sigma_SA=2.d0/3.d0）
!	   else if(NV .eq. 7) then
!	      alfa(6)=alfa(6)+0.09*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)                          !处理源项刚性
!         alfa(7)=alfa(7)+2.d0*0.0828*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
	   endif
     endif
	 	 alfa(1:NV)=alfa(1:NV)+Sfac1*B%vol(i,j,k)         ! 单时间步长Sfac1=0
		  		   
	 if(i.ne. 1) then
!                                              通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i-1,j,k),B%DU(1:NV,i-1,j,k),B%ni1(i,j,k),B%ni2(i,j,k),B%ni3(i,j,k),gamma)  
       dui(1:NV)=0.5d0*(DF(1:NV)*B%si(i,j,k)+w_LU*Lci(i-1,j,k)*B%DU(1:NV,i-1,j,k))
       if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+Lvi(i-1,j,k)*B%DU(1:NV,i-1,j,k)         
      else
	   dui(1:NV)=0.d0                             ! 左侧没有点
      endif
	 
	 if(j.ne.1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j-1,k),B%DU(1:NV,i,j-1,k),B%nj1(i,j,k),B%nj2(i,j,k),B%nj3(i,j,k),gamma)  
        duj(1:NV)=0.5d0*(DF(1:NV)*B%sj(i,j,k)+w_LU*Lcj(i,j-1,k)*B%DU(1:NV,i,j-1,k))
        if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+Lvj(i,j-1,k)*B%DU(1:NV,i,j-1,k)   ! 2012-2-29
	 else
	   duj(1:NV)=0.d0
	 endif

    if(k .ne. 1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j,k-1),B%DU(1:NV,i,j,k-1),B%nk1(i,j,k),B%nk2(i,j,k),B%nk3(i,j,k),gamma)  
        duk(1:NV)=0.5d0*(DF(1:NV)*B%sk(i,j,k)+w_LU*Lck(i,j,k-1)*B%DU(1:NV,i,j,k-1))
        if(If_viscous .eq. 1)    duk(1:NV)=duk(1:NV)+Lvk(i,j,k-1)*B%DU(1:NV,i,j,k-1)   ! 2012-2-29
	else
	   duk(1:NV)=0.d0
	endif
	do m=1,NV
    B%DU(m,i,j,k)=(B%Res(m,i,j,k)+dui(m)+duj(m)+duk(m))/alfa(m)
	enddo
    enddo
    enddo
!$OMP END PARALLEL DO 
   enddo
!----------------------------------------------------------
!  从 (nx-1,ny-1,nz-1)到(1,1,1)的扫描过程 （向下扫描过程）
!  plane=i+j+k
   do plane=nx+ny+nz-3,3,-1   

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(plane,nx,ny,nz,NV,B,Lci,Lcj,Lck,Lvi,Lvj,Lvk,gamma,w_LU,If_viscous)
     do k=nz-1,1,-1
	 do j=ny-1,1,-1
	 i=plane-k-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle            ! 超出了这个平面
      alfa(1:NV)=B%vol(i,j,k)/B%dt(i,j,k)+w_LU*(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k))
      if(If_viscous .eq. 1) then
	    alfa(1:NV)=alfa(1:NV)+2.d0*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k) )         
	   if(NV .eq. 6) then
	     alfa(6)=alfa(6)+Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k)             ! Sigma_SA=2./3.
!	    elseif(NV .eq. 7) then
!	      alfa(6)=alfa(6)+0.09*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
!         alfa(7)=alfa(7)+2.d0*0.0828*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
	   endif
      endif
	 	 alfa(1:NV)=alfa(1:NV)+Sfac1*B%vol(i,j,k)         ! 单时间步长Sfac1=0

	 if(i.ne. nx-1) then
!                                              通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i+1,j,k),B%DU(1:NV,i+1,j,k),B%ni1(i,j,k),B%ni2(i,j,k),B%ni3(i,j,k),gamma)  
         dui(1:NV)=-0.5d0*(DF(1:NV)*B%si(i+1,j,k)-w_LU*Lci(i+1,j,k)*B%DU(1:NV,i+1,j,k))
         if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+Lvi(i+1,j,k)*B%DU(1:NV,i+1,j,k)
  	   
	   else
	   dui(1:NV)=0.d0
       endif
	 
	 if(j.ne. ny-1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j+1,k),B%DU(1:NV,i,j+1,k),B%nj1(i,j,k),B%nj2(i,j,k),B%nj3(i,j,k),gamma)  
       duj(1:NV)=-0.5d0*(DF(1:NV)*B%sj(i,j+1,k)-w_LU*Lcj(i,j+1,k)*B%DU(1:NV,i,j+1,k))
       if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+Lvj(i,j+1,k)*B%DU(1:NV,i,j+1,k)

	 else
	   duj(1:NV)=0.d0
	 endif

    if(k .ne. nz-1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j,k+1),B%DU(1:NV,i,j,k+1),B%nk1(i,j,k),B%nk2(i,j,k),B%nk3(i,j,k),gamma)  
       duk(1:NV)=-0.5d0*(DF(1:NV)*B%sk(i,j,k+1)-w_LU*Lck(i,j,k+1)*B%DU(1:NV,i,j,k+1))
       if(If_viscous .eq. 1)    duk(1:NV)=duk(1:NV)+Lvk(i,j,k+1)*B%DU(1:NV,i,j,k+1)
	else
	   duk(1:NV)=0.d0
	endif
	do m=1,NV
      B%DU(m,i,j,k)=B%DU(m,i,j,k)+(dui(m)+duj(m)+duk(m))/alfa(m)
	enddo
   enddo
   enddo
!$OMP END PARALLEL DO 

   enddo

   end subroutine  du_LU_SGS

!-----------------------------------------------------------------------


!  计算通量的差量 DF=F(Unew)-F(Uold),  LU-SGS方法中使用，用来近似A*DU  
    subroutine comput_DFn(NVAR1,DF,U,DU,n1,n2,n3,gamma)
    use precision_EC
	implicit none
    integer:: NVAR1
    real(PRE_EC),dimension(NVAR1):: DF,U,DU,U2
	real(PRE_EC):: n1,n2,n3,un1,un2,gamma,p1,p2
    U2=U+DU
	un1=(U(2)*n1+U(3)*n2+U(4)*n3)/U(1)              !un
	p1=(gamma-1.d0)*(U(5)-0.5d0*(U(2)**2+U(3)**2+U(4)**2)/U(1))
	un2=(U2(2)*n1+U2(3)*n2+U2(4)*n3)/U2(1)          !un
	p2=(gamma-1.d0)*(U2(5)-0.5d0*(U2(2)**2+U2(3)**2+U2(4)**2)/U2(1))

    DF(1)=U2(1)*un2-U(1)*un1                  ! d*un
	DF(2)=(U2(2)*un2+p2*n1)-(U(2)*un1+p1*n1)
	DF(3)=(U2(3)*un2+p2*n2)-(U(3)*un1+p1*n2)
	DF(4)=(U2(4)*un2+p2*n3)-(U(4)*un1+p1*n3)
	DF(5)=(U2(5)+p2)*un2-(U(5)+p1)*un1
	if(NVAR1 .eq. 6) then
	  DF(6)=U2(6)*un2-U(6)*un1
	else if(NVAR1 .eq. 7) then
      DF(6)=U2(6)*un2-U(6)*un1   ! k方程的对流通量
	  DF(7)=U2(7)*un2-U(7)*un1   ! w方程的对流通量
    endif
	end subroutine comput_dFn