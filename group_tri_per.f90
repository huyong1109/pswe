program main
    implicit none
    include 'mpif.h'

    integer rank,np,err

    real(8),allocatable :: a(:,:)
    real(8),allocatable :: b(:,:)
    real(8),allocatable :: c(:,:)
    real(8),allocatable :: x(:,:)
    real(8),allocatable :: r(:,:)
    real(8),allocatable :: ca(:),cb(:),cc(:),cx(:),cr(:)
    integer :: idx,idy
    integer :: LN !local number for each process
    integer :: M
    integer :: N
    real(8),allocatable :: ra(:,:),rb(:,:),rc(:,:),rr(:,:),rx(:,:)

    call MPI_INIT(err)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    M=4
    N=360
    LN=N/np

    allocate(a(LN,M))
    allocate(b(LN,M))
    allocate(c(LN,M))
    allocate(x(LN,M))
    allocate(r(LN,M))

    allocate(ra(N,M))
    allocate(rb(N,M))
    allocate(rc(N,M))
    allocate(rx(N,M))
    allocate(rr(N,M))

    allocate(ca(LN))
    allocate(cb(LN))
    allocate(cc(LN))
    allocate(cx(LN))
    allocate(cr(LN))

    !call  prepare_check(solver,LN,M,N,loop,a,b,c,x,r,ra,rb,rc,rx,rr,ca,cb,cc,cx,cr)
    !call  prepare_check(4,LN,M,N,1,a,b,c,x,r,ra,rb,rc,rx,rr,ca,cb,cc,cx,cr)
    !call  prepare_check(3,LN,M,N,1,a,b,c,x,r,ra,rb,rc,rx,rr,ca,cb,cc,cx,cr)
    call  prepare_check(5,LN,M,N,1,a,b,c,x,r,ra,rb,rc,rx,rr,ca,cb,cc,cx,cr)
    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(x)
    deallocate(r)
    deallocate(ca)
    deallocate(cb)
    deallocate(cc)
    deallocate(cx)
    deallocate(cr)
    deallocate(ra)
    deallocate(rb)
    deallocate(rc)
    deallocate(rr)
    deallocate(rx)

    call MPI_FINALIZE(err)

end program

subroutine prepare_check(solver,LN,M,N,loop,a,b,c,x,r,ra,rb,rc,rx,rr,ca,cb,cc,cx,cr)
    implicit none
    include 'mpif.h'

    integer :: solver,loop
    integer :: rank,np,err
    real(8) :: t,maxt

    real(8):: a(LN,M)
    real(8):: b(LN,M)
    real(8):: c(LN,M)
    real(8):: x(LN,M)
    real(8):: r(LN,M)
    real(8) :: ca(LN),cb(LN),cc(LN),cx(LN),cr(LN)
    integer :: idx,idy
    integer :: LN !local number for each process
    integer :: M
    integer :: N
    real(8) :: ra(N,M),rb(N,M),rc(N,M),rr(N,M),rx(N,M)

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    if(rank==0)then
        open(123,file="./check/a.dat")
        open(124,file="./check/b.dat")
        open(125,file="./check/c.dat")
        open(126,file="./check/r.dat")
        open(127,file="./check/x.dat")
        do idy = 1,M
            do idx = 1,N
        	read(123,*)ra(idx,idy)
        	read(124,*)rb(idx,idy)
        	read(125,*)rc(idx,idy)
        	read(126,*)rr(idx,idy)
        	read(127,*)rx(idx,idy)
            enddo
        enddo
        close(123)
        close(124)
        close(125)
        close(126)
        close(127)
    endif

    do idx=1,M
        call MPI_SCATTER(ra(:,idx),LN,MPI_DOUBLE_PRECISION,a(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rb(:,idx),LN,MPI_DOUBLE_PRECISION,b(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rc(:,idx),LN,MPI_DOUBLE_PRECISION,c(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rr(:,idx),LN,MPI_DOUBLE_PRECISION,r(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rx(:,idx),LN,MPI_DOUBLE_PRECISION,x(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    enddo

    !do idx=1,M
    !    ca(:)=a(:,idx)
    !    cb(:)=b(:,idx)
    !    cc(:)=c(:,idx)
    !    cr(:)=r(:,idx)
    !    cx(:)=x(:,idx)
    !    call check(ca,cb,cc,cr,cx,LN,N,idx)
    !end do

    !a(:,:)  =1666.0
    !b(:,:)  =-832.0
    !c(:,:)  =-832.0
    !r(:,:)  =78400.0

    t = MPI_WTIME()
    if (solver==4) then
	call tridiagnol_solver4(a,b,c,r,LN,M)
    else if(solver==3)then
	call tridiagnol_solver3(a,b,c,r,LN,M)
    else if(solver==2)then
	call tridiagnol_solver2(a,b,c,r,LN,M)
    else if(solver==5)then
	call tridiagnol_solver5(a,b,c,r,LN,M,np)
    else if(solver==0)then
	call LU(b,a,c,r,M,LN)
    end if
    !call tridiagnol_solver3(a,b,c,r,LN,M)
    !call tridiagnol_solver4(a,b,c,r,LN,M)
    !call tridiagnol_solver2(a,b,c,r,LN,M)
    t = MPI_WTIME()-t
    call MPI_REDUCE(t,maxt,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,err)
    if(rank==0)print *, maxt,'is the running time of solver ',solver

    x=r
    !a(:,:)  =1666.0
    !b(:,:)  =-832.0
    !c(:,:)  =-832.0
    !r(:,:)  =78400.0
    do idx=1,M
        call MPI_SCATTER(ra(:,idx),LN,MPI_DOUBLE_PRECISION,a(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rb(:,idx),LN,MPI_DOUBLE_PRECISION,b(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rc(:,idx),LN,MPI_DOUBLE_PRECISION,c(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_SCATTER(rr(:,idx),LN,MPI_DOUBLE_PRECISION,r(:,idx),LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    enddo
    do idx=1,M
        ca(:)=a(:,idx)
        cb(:)=b(:,idx)
        cc(:)=c(:,idx)
        cr(:)=r(:,idx)
        cx(:)=x(:,idx)
        call check(ca,cb,cc,cr,cx,LN,N,idx)
    end do

end subroutine

subroutine check(a,b,c,r,x,LN,N,g)
    implicit none
    include 'mpif.h'

    integer :: idx
    integer :: rank,np,err
    integer :: LN,N,g
    real(8) :: a(LN),b(LN),c(LN),r(LN),x(LN)
    real(8) :: aa(N),bb(N),cc(N),rr(N),xx(N)
    real(8) :: res,xn,xp !x_next x_previous
    logical :: flag

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    !!! gather calculate send
    call MPI_GATHER(a,LN,MPI_DOUBLE_PRECISION,aa,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_GATHER(b,LN,MPI_DOUBLE_PRECISION,bb,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_GATHER(c,LN,MPI_DOUBLE_PRECISION,cc,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_GATHER(r,LN,MPI_DOUBLE_PRECISION,rr,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_GATHER(x,LN,MPI_DOUBLE_PRECISION,xx,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)

    flag=.FALSE.
    if (rank==0) then
	do idx=1,N
	    if (idx==1) then
		xp=xx(N)
	    else
		xp=xx(idx-1)
	    end if
	    if (idx==N) then
		xn=xx(1)
	    else
		xn=xx(idx+1)
	    end if
	    res=abs(rr(idx)-bb(idx)*xp-aa(idx)*xx(idx)-cc(idx)*xn)
	    if (res>1.0d-7) then
		print *,"error:","N:",N,"idx:",idx,"res:",res,"xx:", xx(idx),"b:", bb(idx),"c", cc(idx),"a", aa(idx),"r",rr(idx),"xp", xp,"xn", xn
		flag=.TRUE.
	    end if
	end do
	if (flag) then
	    print *,"group",g,"check failed"
	    do idx=1,N
		    print *,idx,xx(idx)
	    end do
	else
	    print *, "group",g,"check pass!"
	    !do idx=1,N
	    !     print *,idx,xx(idx)
	    !end do
	end if
    end if

    return
end subroutine

! note:
! 1, the matrix should be like
! a b c
!   a b c
!     a b c
! 2, all arrays maybe changed in this subroutine
! 3, output is r
! 4, we can change the i,j order to get better memory access pattern
! 5, LN must > 1
subroutine tridiagnol_solver4(b,a,c,r,LN,M)
    implicit none
    include 'mpif.h'

    integer :: i
    integer :: rank,np,err,p_down,p_up,stat(MPI_STATUS_SIZE)
    ! rank of process; num of processes; process down/up 
    integer :: LN,M
    ! local vector length; group number
    real(8),dimension(1:LN,1:M) :: a,b,c,r
    real(8),dimension(1:M) :: p,q,au,bu,cu,ru,cn,bn,rn,xd,xn
    ! p/q is temp parameter; 
    ! au/bu/cu/ru is coefficient from up process;
    ! cn/bn/rn is the Nth equation passing from top to bottom; 
    ! xd/xn is for the chasing back
    real(8),dimension(1:7,1:M) :: sendbuf1,recvbuf1
    ! send/recv buffer for chasing down: au/bu/cu/ru/cn/bn/rn
    real(8),dimension(1:2,1:M) :: sendbuf2,recvbuf2
    ! send/recv buffer for chasing back: xd/xn

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,  err)
    p_up   = rank-1
    p_down = rank+1
    if(rank==0)    p_up  = MPI_PROC_NULL
    if(rank==np-1) then
	p_down= MPI_PROC_NULL
	LN = LN - 1
    end if

    !!! P[np-1] send the last equation to p[0]
    if(rank==np-1)then
	sendbuf1(5,:)=b(LN+1,:)
	sendbuf1(6,:)=c(LN+1,:)
	sendbuf1(7,:)=r(LN+1,:)
	call MPI_SEND(sendbuf1,7*M,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,err)
    endif

    !!! Chase down and pass the last equation down
    if(rank/=0)then
	call MPI_RECV(recvbuf1,7*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,stat,err)
    else
	call MPI_RECV(recvbuf1,7*M,MPI_DOUBLE_PRECISION,np-1,1,MPI_COMM_WORLD,stat,err)
    endif
    au=recvbuf1(1,:)
    bu=recvbuf1(2,:)
    cu=recvbuf1(3,:)
    ru=recvbuf1(4,:)
    bn=recvbuf1(5,:)
    cn=recvbuf1(6,:)
    rn=recvbuf1(7,:)

    ! note : LN > 1
    if(rank/=0)then
	p      = a(1,:)/bu
	a(1,:) =	-p*au
	b(1,:) = b(1,:) -p*cu
	r(1,:) = r(1,:) -p*ru
	q      = cn/bu
	cn     =	-q*cu
	bn     = bn     -q*au
	rn     = rn     -q*ru
    end if

    do i=2,LN
	p      = a(i,:)/b(i-1,:)
	a(i,:) =	-p*a(i-1,:)
	b(i,:) = b(i,:) -p*c(i-1,:)
	r(i,:) = r(i,:) -p*r(i-1,:)
	q      = cn/b(i-1,:)
	cn     =	-q*c(i-1,:)
	bn     = bn	-q*a(i-1,:)
	rn     = rn	-q*r(i-1,:)
    enddo

    sendbuf1(1,:)=a(LN,:)
    sendbuf1(2,:)=b(LN,:)
    sendbuf1(3,:)=c(LN,:)
    sendbuf1(4,:)=r(LN,:)
    sendbuf1(5,:)=bn
    sendbuf1(6,:)=cn
    sendbuf1(7,:)=rn
    call MPI_SEND(sendbuf1,7*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,err)
    !!! Handle the last equation, get xn
    if(rank==np-1)then
	a(LN+1,:)=a(LN+1,:)+cn
	q        =a(LN+1,:)/b(LN,:)
	b(LN+1,:)=bn-q*(c(LN,:)+a(LN,:))
	r(LN+1,:)=rn-q*r(LN,:)
	r(LN+1,:)=r(LN+1,:)/b(LN+1,:)
    endif
    !!! Chase back
    call MPI_RECV(recvbuf2,2*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,stat,err)
    xd=recvbuf2(1,:)
    xn=recvbuf2(2,:)

    if(rank==np-1)then
	xd=r(LN+1,:)
	xn=r(LN+1,:)
    endif

    r(LN,:)=(r(LN,:)-c(LN,:)*xd-a(LN,:)*xn)/b(LN,:)
    do i=LN-1,1,-1
	r(i,:)=(r(i,:)-c(i,:)*r(i+1,:)-a(i,:)*xn)/b(i,:)
    enddo

    sendbuf2(1,:)=r(1,:)
    sendbuf2(2,:)=xn
    call MPI_SEND(sendbuf2,2*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,err)

    if(rank==np-1) LN=LN+1
    return
end subroutine
! This is a solver for periodic tridiagonal linear system.
! The main idea is Wang's method.
! We solve two linear systems simutaniously,
!              A * x'  = r'
!              A * xc = rc
!
! input: a, coefficient of x(n)
!	 b, coefficient of x(n-1)
!	 c, coefficietn of x(n+1)
!	 r, right hand vector 
!	 LN,local vector r's length
!	 N, the total length of r
!
! output: x,local solution vector x
!
! note:
! 1, the matrix A should be look like
!    b a c  
!      b a c   
!        b a c
!    
! 2, the input a,b,c,r will be changed in this function

!subroutine tridiagnol_solver2(a,b,c,r,x,LN,N)
subroutine tridiagnol_solver2(a,b,c,r,LN,M)
    implicit none
    include 'mpif.h'

    integer :: idx
    integer :: rank,np,err ! MPI parameters
    integer :: p_down,p_up ! neighbor process number
    integer :: LN,M        
    integer :: stat(MPI_STATUS_SIZE)
    real(8) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M),x(LN,M)
    real(8) :: f(LN,M),g(LN,M) ! Wang's method parameters
    real(8) :: mc(M)
    real(8) :: ad(M),gd(M),fd(M),rd(M) ! communication variables, d is for "from first line in DOWN process"
    real(8) :: au(M),gu(M),ru(M)    ! u is for "from last line in UP process"
    real(8) :: sendbuf(5,M),recvbuf(5,M)

    real(8) :: rc(LN,M),xc(LN,M) !c is for "CYCLE bound". We solve two linear system simutaniously
    real(8) :: rcd(M),rcu(M)       !communication variables for cycle vectors

    rc(:,:)=0
    xc(:,:)=0

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    p_up   = rank-1
    p_down = rank+1
    if(rank==0) then
	p_up  =MPI_PROC_NULL
	rc(1,:)=-1*b(1,:)
	b(1,:)=0  ! for the later operation
    end if
    if(rank==np-1) then
	p_down=MPI_PROC_NULL
	LN=LN-1 ! save x(LN+1) 
	rc(LN,:)=-1*c(LN,:)
    end if

    !! Wang's method
    f(1,:)=b(1,:) ! notice rank=0 b(1)=0
    do idx=2,LN
	mc     =b(idx,:)/a(idx-1,:)
	f(idx,:)=0     -f(idx-1,:)*mc
	a(idx,:)=a(idx,:)-c(idx-1,:)*mc
	r(idx,:)=r(idx,:)-r(idx-1,:)*mc
	rc(idx,:)=rc(idx,:)-rc(idx-1,:)*mc
    end do

    g(LN-1,:)=c(LN-1,:)
    do idx=LN-2,1,-1
	mc     =c(idx,:)/a(idx+1,:)
	g(idx,:)=0     -g(idx+1,:)*mc
	f(idx,:)=f(idx,:)-f(idx+1,:)*mc
	r(idx,:)=r(idx,:)-r(idx+1,:)*mc
	rc(idx,:)=rc(idx,:)-rc(idx+1,:)*mc
    end do

    !!! here we can change

    sendbuf(1,:)=a(1,:)
    sendbuf(2,:)=g(1,:)
    sendbuf(3,:)=f(1,:)
    sendbuf(4,:)=r(1,:)
    sendbuf(5,:)=rc(1,:)
    call MPI_SENDRECV(sendbuf,5*M,MPI_DOUBLE_PRECISION,p_up,1,recvbuf,5*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,stat,err)
    ad=recvbuf(1,:)
    gd=recvbuf(2,:)
    fd=recvbuf(3,:)
    rd=recvbuf(4,:)
    rcd=recvbuf(5,:)
    
    if(rank<np-1)then
	mc    =b(LN,:)/ad !notice rank=np-1
	g(LN,:)=0    -gd*mc
	a(LN,:)=a(LN,:)-fd*mc
	r(LN,:)=r(LN,:)-rd*mc
	rc(LN,:)=rc(LN,:)-rcd*mc
    end if

    call MPI_RECV(recvbuf,5*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,stat,err)
    au=recvbuf(1,:)
    gu=recvbuf(2,:)
    ru=recvbuf(3,:)
    rcu=recvbuf(4,:)

    if(rank>0) then
	mc    =f(LN,:)/au
	a(LN,:)=a(LN,:)-gu*mc
	r(LN,:)=r(LN,:)-ru*mc
	rc(LN,:)=rc(LN,:)-rcu*mc
    end if

    if(rank>0) then !simd
	do idx=LN-1,1,-1
	    mc     =f(idx,:)/au
	    g(idx,:)=g(idx,:)-gu*mc
	    r(idx,:)=r(idx,:)-ru*mc
	    rc(idx,:)=rc(idx,:)-rcu*mc
	end do
    end if

    sendbuf(1,:)=a(LN,:)
    sendbuf(2,:)=g(LN,:)
    sendbuf(3,:)=r(LN,:)
    sendbuf(4,:)=rc(LN,:)
    call MPI_SEND(sendbuf,5*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,err)

    call MPI_RECV(recvbuf,5*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,stat,err)
    rcd=recvbuf(3,:)
    rd=recvbuf(2,:)
    ad=recvbuf(1,:)
    if(rank<np-1)then
	r(LN,:)=r(LN,:)-rd*g(LN,:)/ad
	rc(LN,:)=rc(LN,:)-rcd*g(LN,:)/ad
    end if

    x(LN,:)=r(LN,:)/a(LN,:)
    xc(LN,:)=rc(LN,:)/a(LN,:)
    do idx=LN-1,1,-1 !simd
	r(idx,:)=r(idx,:)-r(LN,:)*g(idx,:)/a(LN,:)
	rc(idx,:)=rc(idx,:)-rc(LN,:)*g(idx,:)/a(LN,:)
	x(idx,:)=r(idx,:)/a(idx,:)
	xc(idx,:)=rc(idx,:)/a(idx,:)
    end do

    sendbuf(1,:)=a(LN,:)
    sendbuf(2,:)=r(LN,:)
    sendbuf(3,:)=rc(LN,:)
    call MPI_SEND(sendbuf,5*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,err)

    !! For periodic bound
    if(rank==0)then
	sendbuf(1,:)=x(1,:)
	sendbuf(2,:)=xc(1,:)
	call MPI_SEND(sendbuf,5*M,MPI_DOUBLE_PRECISION,np-1,1,MPI_COMM_WORLD,err)
    end if
    if(rank==np-1)then
	call MPI_RECV(recvbuf,5*M,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,stat,err)
	x(LN+1,:)=(r(LN+1,:)-c(LN+1,:)*recvbuf(1,:)-b(LN+1,:)*x(LN,:))/(a(LN+1,:)+c(LN+1,:)*recvbuf(2,:)+b(LN+1,:)*xc(LN,:))
	sendbuf(1,:)=x(LN+1,:)
    end if

    call MPI_BCAST(sendbuf,5*M,MPI_DOUBLE_PRECISION,np-1,MPI_COMM_WORLD,err)
    
    do idx=1,LN
	x(idx,:)=x(idx,:)+sendbuf(1,:)*xc(idx,:)
    end do

    r=x

    if(rank==np-1) LN=LN+1

    return
end subroutine

! 1, the equation should be like:
!
!	     |b c     a| 
!  	     |a b c    |
!  	     |  a b c  | * x = r
!  	     |    a b c|
!  	     |c     a b|
subroutine tridiagnol_solver5(b,a,c,r,LN,M,np)
    implicit none
    include 'mpif.h'

    integer :: i,j
    integer :: rank,np,err ! MPI parameters
    integer :: p_down,p_up ! neighbor process number
    integer :: LN,M        
    integer :: stat(MPI_STATUS_SIZE)
    real(8) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M)
    real(8) :: f(LN,M),g(LN,M) ! Wang's method parameters
    real(8) :: mc(M)

    real(8) :: sendbuf(M,8),recvbuf(M,8,np)
    real(8) :: sendbuf2(M,2*np),recvbuf2(M,2)
    real(8), dimension(2*np,M):: ga,gb,gc,gr ! extract equation

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    !call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    f(1,:)=a(1,:)
    f(2,:)=a(2,:)
    do i=3,LN
	mc     = a(i,:)/b(i-1,:)
	f(i,:) = -f(i-1,:)*mc
	b(i,:) = b(i,:) - c(i-1,:)*mc
	r(i,:) = r(i,:) - r(i-1,:)*mc
    enddo ! prove

    g(LN,:)=c(LN,:)
    g(LN-1,:)=c(LN-1,:)
    do i=LN-2,2,-1
	mc     = c(i,:)/b(i+1,:)
	g(i,:) = -g(i+1,:)*mc
	f(i,:) = f(i,:) - f(i+1,:)*mc
	r(i,:) = r(i,:) - r(i+1,:)*mc
    enddo
    mc = c(1,:)/b(2,:)
    g(1,:) = -g(2,:)*mc
    b(1,:) = b(1,:) - f(2,:)*mc
    r(1,:) = r(1,:) - r(2,:)*mc
    
    sendbuf(:,1)=f(1,:)
    sendbuf(:,2)=b(1,:)
    sendbuf(:,3)=g(1,:)
    sendbuf(:,4)=r(1,:)
    sendbuf(:,5)=f(LN,:)
    sendbuf(:,6)=b(LN,:)
    sendbuf(:,7)=g(LN,:)
    sendbuf(:,8)=r(LN,:)
    !if (rank==0) then
    ! print *, "rank",rank,"M",M,"np",np
    ! print *, "====="
    ! print *, sendbuf
    ! print *, "-----"
    ! print *, recvbuf
    !endif
    call MPI_GATHER(sendbuf,8*M,MPI_DOUBLE_PRECISION,recvbuf,8*M,MPI_DOUBLE_PRECISION,np-1,MPI_COMM_WORLD,err)
    if(rank==np-1)then
	do i=1,np
	    ! f>a b>b g>c r>r
	    j=2*(i-1)+1
	    ga(j,:)=recvbuf(:,1,i)
	    gb(j,:)=recvbuf(:,2,i)
	    gc(j,:)=recvbuf(:,3,i)
	    gr(j,:)=recvbuf(:,4,i)
	    j=j+1
	    ga(j,:)=recvbuf(:,5,i)
	    gb(j,:)=recvbuf(:,6,i)
	    gc(j,:)=recvbuf(:,7,i)
	    gr(j,:)=recvbuf(:,8,i)
	enddo
	!subroutine LU(a,b,c,r,m,n)
	call LU(ga,gb,gc,gr,M,2*np)
	do i=1,2*np
	    sendbuf2(:,i)=gr(i,:)
	enddo
    endif
    
    call MPI_SCATTER(sendbuf2,2*M,MPI_DOUBLE_PRECISION,recvbuf2,2*M,MPI_DOUBLE_PRECISION,np-1,MPI_COMM_WORLD,err)
    r(1,:)=recvbuf2(:,1)
    r(LN,:)=recvbuf2(:,2)

    do i=2,LN-1
	r(i,:)=(r(i,:)-g(i,:)*r(LN,:)-f(i,:)*r(1,:))/b(i,:)
    enddo

    return
end subroutine
! This is a solver for periodic tridiagonal linear system.
! The main idea is chasing method.
! We solve two linear systems simutaniously,
!              A * x'  = r'
!              A * xc = rc
!
! input: a, coefficient of x(n)
!	 b, coefficient of x(n-1)
!	 c, coefficietn of x(n+1)
!	 r, right hand vector 
!	 LN,local vector r's length
!	 N, the total length of r
!
! output: x,local solution vector x
!
! note:
! 1, the matrix A should be look like
!    b a c  
!      b a c   
!        b a c
!    
! 2, the input a,b,c,r will be changed in this function

subroutine tridiagnol_solver3(a,b,c,r,LN,M)
!subroutine tridiagnol_solver3(a,b,c,r,x,LN,N)
    implicit none
    include 'mpif.h'

    integer :: idx
    integer :: rank,np,err ! MPI paramenters
    integer :: p_down,p_up ! Neighbor processs down and up 
    
    integer :: LN,M
    integer :: stat(MPI_STATUS_SIZE)
    real(8) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M)
    real(8) :: mc(M)

    real(8) :: xd(M)	    ! the first x in Down process, for MPI communication
    real(8) :: au(M),cu(M),ru(M)     ! the last a,c,r in Up process, for MPI communication
    real(8) :: sendbuf(4,M),recvbuf(4,M)

    real(8) :: rcu(M),xcd(M)      ! last rc in Up process, first xc in Down process, for MPI communication
    real(8) :: rc(LN,M),xc(LN,M)! c is for "CYCLE bound"

    rc(:,:)=0
    xc(:,:)=0

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)

    p_up   = rank-1
    p_down = rank+1
    if(rank==0) then
	p_up  =MPI_PROC_NULL
	rc(1,:)=-1*b(1,:)
	b(1,:)=0 ! for the later operation
    end if
    if(rank==np-1) then
	p_down=MPI_PROC_NULL
	!we change LN for the last proc 
	LN=LN-1
	rc(LN,:)=-1*c(LN,:)
    end if

    !! Chasing method
    call MPI_RECV(recvbuf,4*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,stat,err)
    au(:) =recvbuf(1,:)
    cu(:) =recvbuf(2,:)
    ru(:) =recvbuf(3,:)
    rcu(:)=recvbuf(4,:)

    if(rank==0)then
	au(:)=1
	cu(:)=0
	ru(:)=0
	rcu(:)=0
    end if

    mc   =b(1,:)/au
    a(1,:)=a(1,:)-cu*mc
    r(1,:)=r(1,:)-ru*mc
    rc(1,:)=rc(1,:)-rcu*mc

    do idx = 2,LN
	mc     =b(idx,:)/a(idx-1,:)
	a(idx,:) =a(idx,:)-c(idx-1,:)*mc
	r(idx,:) =r(idx,:)-r(idx-1,:)*mc
	rc(idx,:)=rc(idx,:)-rc(idx-1,:)*mc
    end do

    sendbuf(1,:)=a(LN,:)
    sendbuf(2,:)=c(LN,:)
    sendbuf(3,:)=r(LN,:)
    sendbuf(4,:)=rc(LN,:)
    call MPI_SEND(sendbuf,4*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,err)

    ! We send the full send buf, for coding convinient
    call MPI_RECV(recvbuf,4*M,MPI_DOUBLE_PRECISION,p_down,1,MPI_COMM_WORLD,stat,err)
    xd(:)=recvbuf(1,:)
    xcd(:)=recvbuf(2,:)

    if(rank==np-1) then
	xd=0
	xcd=0
    end if

    r(LN,:)=(r(LN,:)-c(LN,:)*xd)/a(LN,:)
    rc(LN,:)=(rc(LN,:)-c(LN,:)*xcd)/a(LN,:)

    do idx = LN-1,1,-1
	r(idx,:)=(r(idx,:)-c(idx,:)*r(idx+1,:))/a(idx,:)
	rc(idx,:)=(rc(idx,:)-c(idx,:)*rc(idx+1,:))/a(idx,:)
    end do

    sendbuf(1,:)=r(1,:)
    sendbuf(2,:)=rc(1,:)
    call MPI_SEND(sendbuf,4*M,MPI_DOUBLE_PRECISION,p_up,1,MPI_COMM_WORLD,err)
    
    !! For periodic bound, exchange x(N)
    if(rank==0)then
	sendbuf(1,:)=r(1,:)
	sendbuf(2,:)=rc(1,:)
	call MPI_SEND(sendbuf,4*M,MPI_DOUBLE_PRECISION,np-1,1,MPI_COMM_WORLD,err)
    end if
    if(rank==np-1)then
	call MPI_RECV(recvbuf,4*M,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,stat,err)
	r(LN+1,:)=(r(LN+1,:)-c(LN+1,:)*recvbuf(1,:)-b(LN+1,:)*r(LN,:))/(a(LN+1,:)+c(LN+1,:)*recvbuf(2,:)+b(LN+1,:)*rc(LN,:))
	sendbuf(1,:)=r(LN+1,:)
    end if

    call MPI_BCAST(sendbuf,4*M,MPI_DOUBLE_PRECISION,np-1,MPI_COMM_WORLD,err)
    
    do idx=1,LN
	r(idx,:)=r(idx,:)+sendbuf(1,:)*rc(idx,:)
    end do

    if(rank==np-1) LN=LN+1
    return
end subroutine

!subroutine tridiagnol_solver1(a,b,c,r,LN,N,M)
!!   triagle
!!   b a c
!!     b a c
!!       b a c
!    implicit none
!    include 'mpif.h'
!
!    integer :: idx
!    integer :: rank,np,err
!    integer :: LN,M,N
!    real(8) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M),x(LN,M)
!    real(8) :: aa(N,M),bb(N,M),cc(N,M),rr(N,M),xx(N,M)
!    real(8) :: mc(M)
!
!    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
!    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,err)
!
!    !!! gather calculate send
!    call MPI_GATHER(a,LN*M,MPI_DOUBLE_PRECISION,aa,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!    call MPI_GATHER(b,LN,MPI_DOUBLE_PRECISION,bb,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!    call MPI_GATHER(c,LN,MPI_DOUBLE_PRECISION,cc,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!    call MPI_GATHER(r,LN,MPI_DOUBLE_PRECISION,rr,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!
!    if (rank==0) then
!	
!	!do idx=1,N
!    	!    print *,aa(idx),bb(idx),cc(idx),rr(idx)
!    	!end do
!
!	do idx=2,N
!    	    mc = bb(idx)/aa(idx-1)
!    	    aa(idx) = aa(idx)-mc*cc(idx-1)
!    	    rr(idx) = rr(idx)-mc*rr(idx-1)
!    	end do
!
!    	xx(N) = rr(N)/aa(N)
!    	do idx=N-1,1,-1
!    	    xx(idx) = (rr(idx)-cc(idx)*xx(idx+1))/aa(idx)
!    	end do
!    
!	!do idx=1,N
!    	!    print *,xx(idx)
!    	!end do
!	!call CHECK(aa,bb,cc,rr,xx,N,0,1)
!    end if
!
!    call MPI_SCATTER(xx,LN,MPI_DOUBLE_PRECISION,x,LN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!
!    !do idx=1,LN
!    !    print *,rank,idx,x(idx)
!    !end do
!
!    return
!end

!subroutine CHECK(a,b,c,r,x,N,rank,np,xn)
!    real(8) :: a(N),b(N),c(N),r(N),x(N)
!    integer :: i
!    real(8) :: res
!    real(8) :: xn
!    do i = 2,N-1
!        res = r(i)-b(i)*x(i-1)-a(i)*x(i)-c(i)*x(i+1)
!        if (res > 0.000001) print *,"Error!",res
!    end do
!
!    print *,"CHECK PASS!"
!end

! 1, the equation should be like:
!
!	     |b c     a| 
!  	     |a b c    |
!  	     |  a b c  | * x = r
!  	     |    a b c|
!  	     |c     a b|
!
! 2, we handle m equations at the same time
! 3, size of equation is n
! subroutine LU(a,b,c,r,m,n)
!    !Gauss elimination method/Chasing method  with periodic boundary
!    ! on Wang's book P159
!    !use kinds_mod
!    implicit none
!    !
!    integer                      :: i,j, n, m
!    real(8) ,dimension(1:m)      :: sn,rn
!    real(8)		         :: ai
!    !
!    real(8),dimension(1:n, 1:m) :: a,b,c,r
!    real(8),dimension(1:n, 1:m) :: s,t
!    !
!    s(1, 1:m)=a(1, 1:m)  ! s = \tilde(a)
!    t(1, 1:m)=c(n, 1:m)  ! t = \tilde(c)
!    !
!    sn(:)=0.0	   ! \tilde(b)_N
!    rn(:)=0.0	   ! \tilde(r)_N
!    
!    do j = 1, m
!	do i=2,n-1
!	    ai=a(i, j)/b(i-1, j)
!	    b(i, j)=b(i, j)-ai*c(i-1,j)
!	    r(i, j)=r(i, j)-ai*r(i-1,j)
!	    s(i, j)=-ai*s(i-1, j)
!	    !
!	    ai=t(i-1, j)/b(i-1, j)
!	    t(i, j)=-ai*c(i-1, j)
!	    sn(j)  =sn(j)-ai*s(i-1, j)
!	    rn(j)  =rn(j)-ai*r(i-1, j)
!	enddo
!    enddo
!    
!    
!    do j = 1, m
!	a(n, j)=a(n, j)+t(n-1, j)
!    	b(n, j)=b(n, j)+sn(j)
!    	c(n-1, j)=c(n-1, j)+s(n-1, j)
!    	r(n, j)=r(n, j)+rn(j)
!
!	ai=a(n, j)/b(n-1, j)
!    	b(n, j)=b(n, j)-ai*c(n-1, j)
!    	r(n, j)=r(n, j)-ai*r(n-1, j)
!    	!
!    	r(n, j)=r(n, j)/b(n, j)
!    	r(n-1, j)=(r(n-1, j)-c(n-1, j)*r(n, j))/b(n-1, j)
!
!	do i=n-2,1,-1
!	    ai=r(i, j)-s(i, j)*r(n, j)-c(i, j)*r(i+1, j)
!	    r(i, j)=ai/b(i, j)
!	enddo
!    enddo 
!
!    return
!    
! end subroutine LU
subroutine LU(a,b,c,r,m,n)
    !Gauss elimination method/Chasing method  with periodic boundary
    ! on Wang's book P159
    !use kinds_mod
    implicit none
    !
    integer                 :: i,j, n, m
    real(8) ,dimension(1:m)     :: sn,rn
    real(8)		         :: ai
    !
    real(8),dimension(1:n, 1:m) :: a,b,c,r
    real(8),dimension(1:n, 1:m) :: s,t
    !
    s(1, 1:m)=a(1, 1:m)  ! s = \tilde(a)
    t(1, 1:m)=c(n, 1:m)  ! t = \tilde(c)
    !
    sn(:)=0.0	   ! \tilde(b)_N
    rn(:)=0.0	   ! \tilde(r)_N
    
    do j = 1, m
	do i=2,n-1
	    ai=a(i, j)/b(i-1, j)
	    b(i, j)=b(i, j)-ai*c(i-1,j)
	    r(i, j)=r(i, j)-ai*r(i-1,j)
	    s(i, j)=-ai*s(i-1, j)
	    !
	    ai=t(i-1, j)/b(i-1, j)
	    t(i, j)=-ai*c(i-1, j)
	    sn(j)  =sn(j)-ai*s(i-1, j)
	    rn(j)  =rn(j)-ai*r(i-1, j)
	enddo
    enddo
    
    
    do j = 1, m
	a(n, j)=a(n, j)+t(n-1, j)
    	b(n, j)=b(n, j)+sn(j)
    	c(n-1, j)=c(n-1, j)+s(n-1, j)
    	r(n, j)=r(n, j)+rn(j)

	ai=a(n, j)/b(n-1, j)
    	b(n, j)=b(n, j)-ai*c(n-1, j)
    	r(n, j)=r(n, j)-ai*r(n-1, j)
    	!
    	r(n, j)=r(n, j)/b(n, j)
    	r(n-1, j)=(r(n-1, j)-c(n-1, j)*r(n, j))/b(n-1, j)

	do i=n-2,1,-1
	    ai=r(i, j)-s(i, j)*r(n, j)-c(i, j)*r(i+1, j)
	    r(i, j)=ai/b(i, j)
	enddo
    enddo

    do j = 1,m
	do i=1,n
	    print *,i,j,r(i,j)
	enddo
    enddo

    return
    
 end subroutine LU
  
