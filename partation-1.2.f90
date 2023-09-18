!-----------------------------------------------------------
! Partation ver 1.1
! by Li Xinliang, 2012-11-13
! 多块结构网格的多进程区域分割，使得各进程的总网格数尽量均匀 （本版本暂不考虑通信量是否均匀）
! 算法： 各块从大到小排列，依次优先投入当前总网格数目最小的进程中
!------------------------------------------------------------
    implicit none
	integer:: formt,Num_Proc
	print*, "Please input Mesh3d.dat format (0 for unformatted /  1 for formatted)"
	read(*,*) formt
	print*, "please input Proc number: Num_proc "
	read(*,*) Num_Proc
	call part(formt,Num_proc)
	end
!------------------------------
  
     subroutine part(formt,Num_Proc)
     implicit none
     integer,allocatable,dimension(:,:):: B_grid
     integer,allocatable,dimension(:):: NI,NJ,NK,Pgrid,Bproc
	 integer:: formt,Num_Proc,Num_Block,m,m0,m1,m2,G0,mb,mg,mp,n,t1,t2,k,Total_cell,Total_grid
	 real:: grid_av
     print*, "Partation ......"
	 if(formt .eq. 0) then
       open(99,file="Mesh3d.dat",form="unformatted")
       read(99) Num_Block         
       allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! Size of each block
       read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
	 else
       open(99,file="Mesh3d.dat")
       read(99,*) Num_Block         ! 
       allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! Size of each block
       read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
     endif
     close(99)

 !    读入每块网格数，按从多到少次序排列
     allocate(B_grid(Num_Block,2))
	 allocate(Pgrid(Num_proc),Bproc(Num_Block))

     Total_Cell=0
	 Total_grid=0
	  do m=1,Num_Block
       B_grid(m,1)=NI(m)*NJ(m)*NK(m) 
	   B_grid(m,2)=m
	   Total_Cell=Total_Cell+(NI(m)-1)*(NJ(m)-1)*(NK(m)-1)
	   Total_grid=Total_grid+NI(m)*NJ(m)*NK(m)
      enddo
	  print*, "Total Cell number=", Total_Cell
	  print*, "Total Grid number=", Total_grid
 
!    按网格点从多到少的次序排序
     do m=1,Num_block
	   G0=B_grid(m,1)     ! 点数
	    mg=m
	    do n=m+1,Num_block   ! 找出数目最大的
         if(B_grid(n,1) .gt. G0 ) then
		    G0=B_grid(n,1)
			mg=n
         endif
       enddo
!        mg块与m块交换
         t1=B_grid(mg,1)
		 t2=B_grid(mg,2)
		 B_grid(mg,1)=B_grid(m,1)
		 B_grid(mg,2)=B_grid(m,2)
         B_grid(m,1)=t1
		 B_grid(m,2)=t2
	  enddo
!----------------------------
      Pgrid(:)=0
	 do m=1,Num_Block
        mb=B_grid(m,2) 
  !     寻找网格数目最小的进程
        mg=Pgrid(1)
		m0=1
		do mp=1,Num_Proc
         if(Pgrid(mp) .lt. mg) then
		   mg=Pgrid(mp)
		   m0=mp
		 endif
		enddo
       Pgrid(m0)=Pgrid(m0)+B_grid(m,1)
	   Bproc(mb)=m0-1
     enddo
!-------------输出--------------------------------------------------
     open(99,file="partation-auto.dat")
	 write(99,*) " Total_block_number    Total_Proc_number "
	 write(99,*)  Num_block, Num_proc
	 write(99,*) " Block_no, Proc_no "
	 do m=1,Num_block
	 write(99,*) m, Bproc(m)
	 enddo
	 close(99)

     open(100,file="part_grid.dat")
	 do m=1,Num_proc
	  write(100,*) m-1, Pgrid(m)
	 enddo
     close(100)

     m1=Pgrid(1)
	 m2=m1
     do mp=1,Num_Proc
      m1=max(m1,Pgrid(mp))
	  m2=min(m2,Pgrid(mp))
	 enddo
	 
	 grid_av=1.d0*Total_grid/Num_Proc
	 print*, "Max grid number is ", m1,  " ... rato to mean grid is", m1/grid_av
	 print*, "Min grid number is ", m2,  " ... rato to mean grid is", m2/grid_av
	 print*, "rato max to min is ", 1.0*m1/m2
   deallocate(NI,NJ,NK,Pgrid,Bproc,B_grid)
   end subroutine part


