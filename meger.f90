program main
    implicit none 
    integer :: i,k,j,f
    real*8 :: mat(1:80,1:41)
    real*8 :: mat2(1:80,1:41)
    real*8 :: b
    integer :: ed, jed
    character*5 :: s1,s2
    do j=0,4
	do i=0,3
	    write(s1,"(I4)") j
	    write(s2,"(I4)") i
	    open(j*5+i,file = 'out/'//s1//'_'//s2)
	    ed = 8
	    if(j==4) ed = ed + 1
	    do k = 1,ed
		do f = 1,20
		    read(j*5+i,*) mat(i*20+f,j*8+k)
		enddo
	    enddo
	    close(j*5+i)
	enddo
    enddo
    
    jed = 0 
    do j=0,4
	    if(j==0) then
		ed = 9
	    else 
		ed = 8
	    endif 
	do i=0,3
	    write(s1,"(I4)") j
	    write(s2,"(I4)") i
	    open(j*5+i,file = '/home/hy/hy/pswe/out/'//s1//'_'//s2)


	    do k = 1,ed
		do f = 1,20
		    read(j*5+i,*) mat2(i*20+f,jed+k)
		    !write(*,*) j, i, k, f,i*20+f,jed+k,
		enddo
	    enddo
	    close(j*5+i)
	enddo
	jed = jed +ed
    enddo
    !open(111,file='~/nyfshallow/debug.dat')
    
    do k = 1,41
	do f = 1,80
	    if(mat(f,k) /= mat2(f,k)) write(*,*),f,k,mat(f,k), mat2(f,k)
	enddo
    enddo
    !do k = 1,41
    !    do f = 1,80
    !        if(mat(f,k) == mat2(f,k)) write(*,*),f,k,mat2(f,k)
    !    enddo
    !enddo
    
130 format(200f16.8)
end
