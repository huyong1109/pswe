!
MODULE module_array
!
    use module_para
    use communicate, only : nprocs
    use distribution
!
    implicit none
!
    real*8,  allocatable  ::    u (:,:)              ! zonal wind
    real*8,  allocatable  ::    v (:,:)              ! meridional wind
    real*8,  allocatable  ::    wh(:,:)              ! geopotential height
!
    real*8,  allocatable  ::    wu(:,:)              ! u*sqrt(wh)
    real*8,  allocatable  ::    wv(:,:)              ! v*sqrt(wh)(
    
    
!   caculate often used vectors 
    real*8,  dimension(:), allocatable  ::    c1(:), s1(:)		           ! c1(j)=cos(theta(j)), s1(j)=sin(theta(j)), where theta is latitude
    real*8,  dimension(:), allocatable  ::    c2, c11, c12, c13, c14  ! 
    real*8,  dimension(:), allocatable  ::    f1, f2		           ! 
!
    public :: init_array, &
	      destroy_array
contains 

    subroutine init_array

    call create_distribution(nprocs, nproc_x, nproc_y, p,n) 
    
    if(.not. allocated (u )) allocate(u (0:nloc_x+1,0:nloc_y+1))
    if(.not. allocated (v )) allocate(v (0:nloc_x+1,0:nloc_y+1))
    if(.not. allocated (wh)) allocate(wh(0:nloc_x+1,0:nloc_y+1))
    if(.not. allocated (wu)) allocate(wu(0:nloc_x+1,0:nloc_y+1))
    if(.not. allocated (wv)) allocate(wv(0:nloc_x+1,0:nloc_y+1))
    
    if(.not. allocated (c1 )) allocate(c1 (0:n+1))
    if(.not. allocated (s1 )) allocate(s1 (0:n+1))
    if(.not. allocated (c2 )) allocate(c2 (0:n+1))
    if(.not. allocated (c11)) allocate(c11(0:n+1))
    if(.not. allocated (c12)) allocate(c12(0:n+1))
    if(.not. allocated (f1 )) allocate(f1 (0:n+1))
    if(.not. allocated (f2 )) allocate(f2 (0:n+1))



    end subroutine init_array
    
    subroutine destroy_array
    
    deallocate(u )
    deallocate(v )
    deallocate(wh)
    deallocate(wu)
    deallocate(wv)
    
    deallocate(c1 )
    deallocate(s1 )
    deallocate(c2 )
    deallocate(c11)
    deallocate(c12)
    deallocate(f1 )
    deallocate(f2 )



    end subroutine destroy_array
END MODULE module_array
!
