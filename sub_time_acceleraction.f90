 !  加速收敛技术
 !   1) 局部时间步长
 !   2) 残差光顺

 !--------------------------------------------------------------------------------------------------
 ! 计算谱半径 see: Blazek's book, p189-190    
     subroutine comput_Lijk(nMesh,mBlock) 
     use Global_Var
	 use Flow_Var
     implicit none
     integer:: mBlock,nx,ny,nz,i,j,k,nMesh
     real(PRE_EC) D0,un,S0,vol1
     Type (Block_TYPE),pointer:: B

     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny; nz=B%nz

! $OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,gamma,PrL,Prt,uu,v,w,cc,B,Lci,Lvi,Lcj,Lvj,Lck,Lvk)
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
       if(If_viscous .eq. 1) then
	    D0=max(gamma,4./3.)/d(i,j,k)*(B%mu(i,j,k)/PrL+B%mu_t(i,j,k)/Prt)
       else
	    D0=0.d0
	   endif
	   vol1=1.d0/B%vol(i,j,k)
       
       un=0.5d0*( uu(i,j,k)*(B%ni1(i,j,k)+B%ni1(i+1,j,k))+ &
                   v(i,j,k)*(B%ni2(i,j,k)+B%ni2(i+1,j,k))+ &
                   w(i,j,k)*(B%ni3(i,j,k)+B%ni3(i+1,j,k)) )
       S0=(B%si(i,j,k)+B%si(i+1,j,k))*0.5
       Lci(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvi(i,j,k)=D0*S0*S0*vol1

       un=0.5d0*( uu(i,j,k)*(B%nj1(i,j,k)+B%nj1(i,j+1,k))+ &
                   v(i,j,k)*(B%nj2(i,j,k)+B%nj2(i,j+1,k))+ &
                   w(i,j,k)*(B%nj3(i,j,k)+B%nj3(i,j+1,k)) )
       S0=(B%sj(i,j,k)+B%sj(i,j+1,k))*0.5
       Lcj(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvj(i,j,k)=D0*S0*S0*vol1


       un=0.5d0*( uu(i,j,k)*(B%nk1(i,j,k)+B%nk1(i,j,k+1))+ &
                   v(i,j,k)*(B%nk2(i,j,k)+B%nk2(i,j,k+1))+ &
                   w(i,j,k)*(B%nk3(i,j,k)+B%nk3(i,j,k+1)) )
       S0=(B%sk(i,j,k)+B%sk(i,j,k+1))*0.5
       Lck(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvk(i,j,k)=D0*S0*S0*vol1

       enddo
       enddo
       enddo
! $OMP END PARALLEL DO 

    end subroutine comput_Lijk


!---------------------------------------------------------------------------------
! 计算（当地）时间步长  ! J. Blazek, P.190
  subroutine comput_dt(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   integer  nMesh,mBlock,nx,ny,nz,i,j,k
   real(PRE_EC) D0,un,S0,vol1,dt_fac
   real(PRE_EC):: C

   Type (Block_TYPE),pointer:: B   
   C=1.d0
   B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
   nx=B%nx; ny=B%ny; nz=B%nz

   if( B%IF_OverLimit .eq. 0) then                ! 物理量超限
     dt_fac=1.0
   else
     dt_fac= 0.1d0              ! 时间步长降低10倍
	 print*, " ------- In Block No. ", B%Block_no,  "flow OverLimit, time step 1/10 ----"
   endif
 

   if(Iflag_local_dt .eq. 1)    then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
       do k=1,nz-1
       do j=1,ny-1
       do i=1,nx-1
         B%dt(i,j,k)=dt_fac*CFL*B%Vol(i,j,k)  &
		           /(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k)+C*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k)))

         if(If_dtime_mesh .eq. 1) B%dt(i,j,k)=B%dt(i,j,k)*B%dtime_mesh(i,j,k)           ! 根据网格质量，修正时间步长

		 if(B%dt(i,j,k) .gt. dtmax) B%dt(i,j,k)=dtmax
         if(B%dt(i,j,k) .lt. dtmin) B%dt(i,j,k)=dtmin
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO 

   else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
          B%dt(i,j,k)=dt_global ! 全局时间步长法
      enddo
	  enddo
      enddo
!$OMP END PARALLEL DO 
  endif


  end subroutine comput_dt

 
 
! 迎风隐式残差光顺
! Ver 0.76  Upwind Implicit Residual smoothing
! 见J. Blazek's Book "Computational Fluid Dynamic: Principle and Application" 
   Subroutine Residual_smoothing(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,m,nMesh
   integer,parameter::Nmax=4000
   real(PRE_EC):: as(Nmax),bs(Nmax),cs(Nmax),R(Nmax),Rs(Nmax),ei,Mn
   real(PRE_EC),parameter:: epsl=1.d0
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)                                         ! 指向其一块
   nx=B%nx;  ny= B%ny;  nz=B%nz

!--i-direction---------------------------
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
      ei=epsl*min(1.0,Lci(i,j,k)/Lcj(i,j,k),Lci(i,j,k)/Lck(i,j,k))
 !     法向Mach数     
	  Mn=0.5d0*( uu(i,j,k)*(B%ni1(i,j,k)+B%ni1(i+1,j,k))+ &
                 v(i,j,k)*(B%ni2(i,j,k)+B%ni2(i+1,j,k))+ &
                 w(i,j,k)*(B%ni3(i,j,k)+B%ni3(i+1,j,k)) )/cc(i,j,k)    
      if(Mn .gt. 1.0) then
       as(i)=-ei;  bs(i)=1.0+ei ;  cs(i)=0.0
      else if (Mn .lt. -1.0) then
	   as(i)=0.0; bs(i)=1.0+ei; cs(i)=-ei
	  else
	   as(i)=-ei; bs(i)=1.+2.*ei ; cs(i)=-ei
	  endif  
    enddo
     do m=1,5
       do i=1,nx-1
        R(i)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(nx-1,as,bs,cs,R,Rs)    ! R残差； Rs光滑后的残差
       do i=1,nx-1
        B%Res(m,i,j,k)=Rs(i)
       enddo
     enddo
   enddo
   enddo

!---j- direction ------------------
   do k=1,nz-1
   do i=1,nx-1
   do j=1,ny-1
    ei=epsl*min(1.0,Lcj(i,j,k)/Lci(i,j,k),Lcj(i,j,k)/Lck(i,j,k))

!     法向Mach数     
     Mn=0.5d0*( uu(i,j,k)*(B%nj1(i,j,k)+B%nj1(i,j+1,k))+ &
                v(i,j,k)*(B%nj2(i,j,k)+B%nj2(i,j+1,k))+ &
                w(i,j,k)*(B%nj3(i,j,k)+B%nj3(i,j+1,k)) )/cc(i,j,k)
     if(Mn .gt. 1.0) then
	  as(j)=-ei ; bs(j)=1.0+ei; cs(j)=0.0
     else if(Mn .lt. -1.0) then
	  as(j)=0.0; bs(j)=1.0+ei; cs(j)=-ei
	 else
	  as(j)=-ei; bs(j)=1.0+2.0*ei; cs(j)=-ei
	 endif
    enddo
      do m=1,5
       do j=1,ny-1
        R(j)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(ny-1,as,bs,cs,R,Rs)
       do j=1,ny-1
        B%Res(m,i,j,k)=Rs(j)
       enddo
      enddo
    enddo
    enddo
!---k- direction ----------------
   do j=1,ny-1
   do i=1,nx-1
      do k=1,nz-1
       ei=epsl*min(1.0,Lck(i,j,k)/Lci(i,j,k),Lck(i,j,k)/Lcj(i,j,k))
        Mn=0.5d0*( uu(i,j,k)*(B%nk1(i,j,k)+B%nk1(i,j,k+1))+ &
                   v(i,j,k)*(B%nk2(i,j,k)+B%nk2(i,j,k+1))+ &
                   w(i,j,k)*(B%nk3(i,j,k)+B%nk3(i,j,k+1)) )/cc(i,j,k)
      if(Mn .gt. 1.0) then
	   as(k)=-ei ; bs(k)=1.0+ei; cs(k)=0.0
      else if(Mn .lt. -1.0) then
	   as(k)=0.0; bs(k)=1.0+ei; cs(k)=-ei
	  else
	   as(k)=-ei; bs(k)=1.0+2.0*ei; cs(k)=-ei
	  endif
     enddo
      do m=1,5
       do k=1,nz-1
         R(k)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(nz-1,as,bs,cs,R,Rs)
       do k=1,nz-1
        B%Res(m,i,j,k)=Rs(k)
       enddo
      enddo
    enddo
    enddo

   end subroutine Residual_smoothing
 

! 三对角方程组求解
! 求解 a(i)*Rs(i-1)+b(i)*Rs(i)+c(i)*Rs(i+1)=R(i)
   subroutine tridiagonal(n,a,b,c,R,Rs)
   use precision_EC
   implicit none
   integer:: n,i
   integer,parameter::Nmax=4000
   real(PRE_EC):: a(n),b(n),c(n),R(n),Rs(n),P(Nmax),Q(Nmax),tmp
   P(1)=0.d0; Q(1)=R(1)
   P(n)=0.d0; Q(n)=R(n)

   do i=2,n
    tmp=1.d0/(a(i)*P(i-1)+b(i))
    P(i)=-c(i)*tmp
    Q(i)=(R(i)-a(i)*Q(i-1))*tmp
   enddo
   Rs(n)=R(n)
   Rs(1)=R(1)
   do i=n-1,2,-1
   Rs(i)=P(i)*Rs(i+1)+Q(i)
   enddo
   end subroutine tridiagonal

