module module_io
    use module_para
    use distribution
    implicit none
contains
    
    subroutine mpi_print(mat)

	use mpi
	!include 'mpif.h'
	integer :: i,j,ierr
	real*8 :: mat(0:nloc_x+1,0:nloc_y+1)
	character(len=5) :: s1,s2
	write(s1,"(I4)") jproc
	write(s2,"(I4)") iproc
	
	open(12,file='out/'//s1//'_'//s2)
	if (jproc /= 0 .and. jproc /= nproc_y -1) then
	    do j = 1,nloc_y
		do i = 1,nloc_x
		    write(12,*),mat(i,j)
		enddo
	    enddo
	endif

	if (jproc == 0 ) then
	    do j = 0,nloc_y
		do i = 1,nloc_x
		    write(12,*),mat(i,j)
		enddo
	    enddo
	endif

	if (jproc == nproc_y-1 ) then
	    do j = 1,nloc_y+1
		do i = 1,nloc_x
		    write(12,*),mat(i,j)
		enddo
	    enddo
	endif



	close(12)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	stop
    end subroutine
    subroutine check_nan(mat, str)
	integer :: i,j
	real*8 :: mat(0:nloc_x+1,0:nloc_y+1)
	character (*) :: str	
	    do j = 0,nloc_y+1
		do i = 0,nloc_x+1
		    if (abs(mat(i,j)) > 1.0E+10) then 
			write(*,*) str, my_task, i, j, mat(i,j)
		    endif 
		enddo
	    enddo

    end subroutine
    subroutine check_nan2(mat, str)
	integer :: i,j
	real*8 :: mat(1:p,0:nloc_y+1)
	character (*) :: str	
	if (iproc == 0 ) then
	    do j = 0,nloc_y+1
		do i = 1,p
		    if (abs(mat(i,j)) > 1.0E+10) then 
			write(*,*) str, my_task, i, j, mat(i,j)
		    endif 
		enddo
	    enddo
	endif 



    end subroutine
end module
