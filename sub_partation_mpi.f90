!---------------------------------------------------------
!  分区： 确定 “进程”与“块”之间的联系
!  B_proc(m), 给出m块所属的进程号 
!  B_n(m), 给出m块在该进程中的内部编号
!  my_blocks(:), 给出本进程所包含的块号
!    
   subroutine partation
   use Global_var
   implicit none
   integer:: ierr,m,mt,k,TP
   integer,dimension(:),pointer:: Bproc
   logical EXT

   interface
     subroutine part(formt,Num_Proc,Num_Block,Bproc)
     implicit none
     integer:: formt,Num_Proc,Num_Block
	 integer,dimension(:),pointer:: Bproc
     end subroutine
   end interface

   if(my_id .eq. 0) then
    Inquire(file="partation.dat",EXIST=EXT)
     if(EXT) then
	   print*, "Find partation.dat, read it"
	   open(100,file="partation.dat")
	   read(100,*)
       read(100,*)  Total_block, TP          ! 总块数
	   read(100,*) 
       if(TP .ne. Total_proc) then
        print*, "Error ! Total_proc is not the same as that in  'partation.dat' "
        print*, "Tp, Total_proc=", Tp, Total_Proc
        stop
       endif
	 else 
	    print*, "can not find partation.dat,  create partation ..."
        call part(Mesh_File_Format,Total_Proc,Total_Block,Bproc)
	 endif
    endif

	 call MPI_bcast(Total_block,1,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
     allocate (B_proc(Total_block),B_n(Total_block))   
!   读入B_proc(:)  , 块-->进程 对应关系     
	 
	 if(my_id .eq. 0) then
     if(EXT) then
	   do m=1,  Total_block 
	   read(100,*) mt, B_proc(m)
	   enddo
	   close(100)
	  else
        B_proc(1:Total_Block)=Bproc(1:Total_Block)
!		print*, "------------------"
!		print*, B_proc(:)
	  endif
	  
	  endif

     call MPI_bcast(B_proc(1),Total_block,MPI_Integer,0,  MPI_COMM_WORLD,ierr)

!  计算B_n(m), 第m块在该进程的内部编号
    do m=1, Total_block
	 B_n(m)=0
	 do k=1,m                 ! m块前面有多数个块在B_proc(m)进程
	  if(B_proc(k) .eq. B_proc(m))  B_n(m)=B_n(m)+1
	 enddo
	enddo 


!  统计每个进程包含的块数
     Num_block=0
	 do m=1, Total_block
	 if(B_proc(m) .eq. my_id )   Num_block=Num_block+1    ! my_id进程包含的块数
     enddo
     
	 allocate (my_blocks(Num_block))       ! 本进程包含的块号（数组）
	 k=1
	 do m=1, Total_block
  	  if(B_proc(m) .eq. my_id )  then
	    my_blocks(k)=m
		k=k+1
	  endif
     enddo

   end subroutine 


     subroutine part(formt,Num_Proc,Num_Block,Bproc)
     implicit none
     integer:: formt,Num_Proc,Num_Block
	 integer,dimension(:),pointer:: Bproc
     integer,allocatable,dimension(:,:):: B_grid
     integer,allocatable,dimension(:):: NI,NJ,NK,Pgrid
	 integer:: m,m0,m1,m2,G0,mb,mg,mp,n,t1,t2,k,Total_cell,Total_grid
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
   deallocate(NI,NJ,NK,Pgrid,B_grid)
   end subroutine part



