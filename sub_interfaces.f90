  module interface_defines
     INTERFACE   
     subroutine get_U_conner(nx,ny,nz,NVAR,U)
     use precision_EC
     implicit none
     integer:: nx,ny,nz,NVAR
     real(PRE_EC),dimension(:,:,:,:),Pointer::U
     end subroutine 
      
     subroutine get_xyz_conner(nx,ny,nz,x)
     use precision_EC
     implicit none
     integer:: nx,ny,nz
     real(PRE_EC),dimension(:,:,:),pointer::x
     end subroutine
     
     subroutine prolongation(nx1,ny1,nz1,nx2,ny2,nz2,U1,U2)
     use precision_EC
     implicit none
     integer:: nx1,ny1,nz1,nx2,ny2,nz2
     real(PRE_EC),dimension(:,:,:,:),pointer:: U1,U2
     end subroutine    
     
    subroutine allocate_mem_Blocks(B)
     use global_var
     implicit none
     Type (Block_TYPE),pointer:: B
    end subroutine


     ENDINTERFACE
end module interface_defines 
