!
subroutine euler(dt,iter)

    use module_para
    use kinds_mod
    use communicate
    use module_array
    use module_io
    use boundary
    use global_reductions
    use mpi    

    implicit none

    integer i,i1,i2,j,k,iter, ierr, break_cir

    real(ld8),dimension(0:nloc_x+1,0:nloc_y+1)  :: tu,tv,th,h 
    real(ld8),dimension(0:nloc_x+1,0:nloc_y+1)  :: du,dv,dh 

    real(ld8),dimension(0:nloc_x+1,0:nloc_y+1)      :: a1,a2,rh,ru 

    real(ld8)                      :: ai,aj,dt,dt2,en,en0,den

    real(ld8),external             :: inner

    !	half timestep
    dt2=dt*0.5

    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    tu(i,j)=wu(i,j)
	    tv(i,j)=wv(i,j)
	    th(i,j)=wh(i,j)
	end do
    end do

    !	Initial energy
    en0=inner(wu,wv,wh,wu,wv,wh)
    call master_print_message(en0, 'en0')

    !	Iterate until convergence ( |E-E0| < 1.0E-15 E0 )
    break_cir = 0
    do k=1,1000
	!call master_print_message(k,'iter')

	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		h(i,j)=sqrt(th(i,j))
	    end do
	end do
	
	!	consider U, phi, compute advection and coriolios terms 

	call DIFUH(tu,tv,du,dh,h)

	call update_boundary(du)
	call update_boundary(dh)
	
	!	get the tri-diagonal matrix
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		tu(i,j)=wu(i,j)-dt2*du(i,j)
		th(i,j)=wh(i,j)-dt2*dh(i,j)
	    end do
	end do
	
	do j=0,nloc_y+1
	    ai=dt2*c11(j)
	    do i=0,nloc_x+1
		aj=ai*h(i,j)
		a1(i,j)=aj
		ru(i,j)=tu(i,j)*aj
		a2(i,j)=aj*aj
	    end do
	end do 
	
		
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		rh(i, j)=th(i,j)-ru(i+1, j)+ru(i-1,j )
	    enddo
	enddo
	
	call update_boundary(rh) 
    
	call Tri_diag(a2, rh, th)

	!call mpi_print(th)
	call update_boundary(th) 


	call mpi_print(th)
	do j = 1, nloc_y
	    do i = 1, nloc_x
		tu(i,j)=tu(i,j)-a1(i,j)*(th(i+1,j)-th(i-1,j))
	    enddo
	enddo 

	call update_boundary(tu) 
	
	call DIFV(tu,tv,th,dv,h)
	
	call update_boundary(dv) 

	do j=0,nloc_y + 1
	    do i=0, nloc_x+1
		tv(i,j)=wv(i,j)-dt2*dv(i,j)
	    enddo
	enddo
	
	en=inner(tu,tv,th,tu,tv,th)
	den=abs(en-en0)*2.0/(en+en0)
	if(my_task .eq. 0) write(*,*),k,en0,en
	en0=en
	if (den.lt.1.0d-10) break_cir = 1
	call MPI_Bcast(break_cir,1,MPI_INT,0,comm_group,ierr)
	if (break_cir .eq. 1) goto 10
    enddo

    10   continue
    iter=k
    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    wu(i,j)=tu(i,j)*2.0d0-wu(i,j)
	    wv(i,j)=tv(i,j)*2.0d0-wv(i,j)
	    wh(i,j)=th(i,j)*2.0d0-wh(i,j)
	end do
    enddo
    return
 end subroutine euler

subroutine tridiagnol_solver3(aa,bb,cc,rr,LN,M)
    use communicate
    use kinds_mod
    implicit none
    include 'mpif.h'

    integer :: idx
    integer :: rank,np,err ! MPI paramenters
    integer :: p_down,p_up ! Neighbor processs down and up 
    
    integer :: LN,M
    integer :: stat(MPI_STATUS_SIZE)
    real(ld8) :: aa(LN,M),bb(LN,M),cc(LN,M),rr(LN,M)
    real(16) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M)
    real(16) :: mc(M)

    real(16) :: xd(M)	    ! the first x in Down process, for MPI communication
    real(16) :: au(M),cu(M),ru(M)     ! the last a,c,r in Up process, for MPI communication
    real(16) :: sendbuf(4,M),recvbuf(4,M)

    real(16) :: rcu(M),xcd(M)      ! last rc in Up process, first xc in Down process, for MPI communication
    real(16) :: rc(LN,M),xc(LN,M)! c is for "CYCLE bound"

    rc(:,:)=0.0d0
    xc(:,:)=0.0d0
    
    a(:,:) = aa(:,:)
    b(:,:) = bb(:,:)
    c(:,:) = cc(:,:)
    r(:,:) = rr(:,:)
    call MPI_COMM_RANK(comm_group,rank,err)
    call MPI_COMM_SIZE(comm_group,np,err)

    p_up   = rank-1
    p_down = rank+1
    if(rank==0) then
	p_up  =MPI_PROC_NULL
	rc(1,:)=-1.0d0*b(1,:)
	b(1,:)=0.0d0! for the later operation
    end if
    if(rank==np-1) then
	p_down=MPI_PROC_NULL
	!we change LN for the last proc 
	LN=LN-1
	rc(LN,:)=-1.0d0*c(LN,:)
    end if

    !! Chasing method
    call MPI_RECV(recvbuf,4*M,MPI_LONG_DOUBLE,p_up,1,comm_group,stat,err)
    au(:) =recvbuf(1,:)
    cu(:) =recvbuf(2,:)
    ru(:) =recvbuf(3,:)
    rcu(:)=recvbuf(4,:)

    if(rank==0)then
	au(:)=1.0d0
	cu(:)=0.0d0
	ru(:)=0.0d0
	rcu(:)=0.0d0
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
    call MPI_SEND(sendbuf,4*M,MPI_LONG_DOUBLE,p_down,1,comm_group,err)

    ! We send the full send buf, for coding convinient
    call MPI_RECV(recvbuf,4*M,MPI_LONG_DOUBLE,p_down,1,comm_group,stat,err)
    xd(:)=recvbuf(1,:)
    xcd(:)=recvbuf(2,:)

    if(rank==np-1) then
	xd=0.0d0
	xcd=0.0d0
    end if

    r(LN,:)=(r(LN,:)-c(LN,:)*xd)/a(LN,:)
    rc(LN,:)=(rc(LN,:)-c(LN,:)*xcd)/a(LN,:)

    do idx = LN-1,1,-1
	r(idx,:)=(r(idx,:)-c(idx,:)*r(idx+1,:))/a(idx,:)
	rc(idx,:)=(rc(idx,:)-c(idx,:)*rc(idx+1,:))/a(idx,:)
    end do

    sendbuf(1,:)=r(1,:)
    sendbuf(2,:)=rc(1,:)
    call MPI_SEND(sendbuf,4*M,MPI_LONG_DOUBLE,p_up,1,comm_group,err)
    
    !! For periodic bound, exchange x(N)
    if(rank==0)then
	sendbuf(1,:)=r(1,:)
	sendbuf(2,:)=rc(1,:)
	call MPI_SEND(sendbuf,4*M,MPI_LONG_DOUBLE,np-1,1,comm_group,err)
    end if
    if(rank==np-1)then
	call MPI_RECV(recvbuf,4*M,MPI_LONG_DOUBLE,0,1,comm_group,stat,err)
	r(LN+1,:)=(r(LN+1,:)-c(LN+1,:)*recvbuf(1,:)-b(LN+1,:)*r(LN,:))/(a(LN+1,:)+c(LN+1,:)*recvbuf(2,:)+b(LN+1,:)*rc(LN,:))
	sendbuf(1,:)=r(LN+1,:)
    end if

    call MPI_BCAST(sendbuf,4*M,MPI_LONG_DOUBLE,np-1,comm_group,err)
    
    do idx=1,LN
	r(idx,:)=r(idx,:)+sendbuf(1,:)*rc(idx,:)
    end do

    rr(:,:) = r(:,:)
    if(rank==np-1) LN=LN+1
    return
end subroutine



! subroutine tridiagnol_solver2(a,b,c,r,LN,M)
!    use distribution,only:iproc  
!    use module_para,only:nproc_x
!    implicit none
!    include 'mpif.h'
!    
!    integer :: idx
!    integer :: rank,np,err ! MPI parameters
!    integer :: p_down,p_up ! neighbor process number
!    integer :: LN,M        
!    integer :: stat(MPI_STATUS_SIZE)
!    real(8) :: a(LN,M),b(LN,M),c(LN,M),r(LN,M),x(LN,M)
!    real(8) :: f(LN,M),g(LN,M) ! Wang's method parameters
!    real(8) :: mc(M)
!    real(8) :: ad(M),gd(M),fd(M),rd(M) ! communication variables, d is for "from first line in DOWN process"
!    real(8) :: au(M),gu(M),ru(M)    ! u is for "from last line in UP process"
!    real(8) :: sendbuf(5,M),recvbuf(5,M)
!
!    real(8) :: rc(LN,M),xc(LN,M) !c is for "CYCLE bound". We solve two linear system simutaniously
!    real(8) :: rcd(M),rcu(M)       !communication variables for cycle vectors
!
!    rc(:,:)=0
!    xc(:,:)=0
!
!    call MPI_COMM_RANK(comm_group,rank,err)
!    call MPI_COMM_SIZE(comm_group,np,err)
!
!    p_up   = rank-1
!    p_down = rank+1
!    if(iproc==0) then
!	p_up  =MPI_PROC_NULL
!	rc(1,:)=-1*b(1,:)
!	b(1,:)=0  ! for the later operation
!    end if
!    if(iproc==nproc_x-1) then
!	p_down=MPI_PROC_NULL
!	LN=LN-1 ! save x(LN+1) 
!	rc(LN,:)=-1*c(LN,:)
!    end if
!
!    !! Wang's method
!    f(1,:)=b(1,:) ! notice rank=0 b(1)=0
!    do idx=2,LN
!	mc     =b(idx,:)/a(idx-1,:)
!	f(idx,:)=0     -f(idx-1,:)*mc
!	a(idx,:)=a(idx,:)-c(idx-1,:)*mc
!	r(idx,:)=r(idx,:)-r(idx-1,:)*mc
!	rc(idx,:)=rc(idx,:)-rc(idx-1,:)*mc
!    end do
!
!    g(LN-1,:)=c(LN-1,:)
!    do idx=LN-2,1,-1
!	mc     =c(idx,:)/a(idx+1,:)
!	g(idx,:)=0     -g(idx+1,:)*mc
!	f(idx,:)=f(idx,:)-f(idx+1,:)*mc
!	r(idx,:)=r(idx,:)-r(idx+1,:)*mc
!	rc(idx,:)=rc(idx,:)-rc(idx+1,:)*mc
!    end do
!
!    sendbuf(1,:)=a(1,:)
!    sendbuf(2,:)=g(1,:)
!    sendbuf(3,:)=f(1,:)
!    sendbuf(4,:)=r(1,:)
!    sendbuf(5,:)=rc(1,:)
!    call MPI_SENDRECV(sendbuf,5*M,MPI_DBL,p_up,1,recvbuf,5*M,MPI_DBL,p_down,1,comm_group,stat,err)
!    ad=recvbuf(1,:)
!    gd=recvbuf(2,:)
!    fd=recvbuf(3,:)
!    rd=recvbuf(4,:)
!    rcd=recvbuf(5,:)
!    
!    if(iproc<nproc_x-1)then
!	mc    =b(LN,:)/ad !notice rank=np-1
!	g(LN,:)=0    -gd*mc
!	a(LN,:)=a(LN,:)-fd*mc
!	r(LN,:)=r(LN,:)-rd*mc
!	rc(LN,:)=rc(LN,:)-rcd*mc
!    end if
!    
!    write(*,*) rank, "1" 
!    call MPI_RECV(recvbuf,5*M,MPI_DBL,p_up,1,comm_group,stat,err)
!    au=recvbuf(1,:)
!    gu=recvbuf(2,:)
!    ru=recvbuf(3,:)
!    rcu=recvbuf(4,:)
!
!    if(iproc>0) then
!	mc    =f(LN,:)/au
!	a(LN,:)=a(LN,:)-gu*mc
!	r(LN,:)=r(LN,:)-ru*mc
!	rc(LN,:)=rc(LN,:)-rcu*mc
!    end if
!
!    if(iproc>0) then !simd
!	do idx=LN-1,1,-1
!	    mc     =f(idx,:)/au
!	    g(idx,:)=g(idx,:)-gu*mc
!	    r(idx,:)=r(idx,:)-ru*mc
!	    rc(idx,:)=rc(idx,:)-rcu*mc
!	end do
!    end if
!
!    sendbuf(1,:)=a(LN,:)
!    sendbuf(2,:)=g(LN,:)
!    sendbuf(3,:)=r(LN,:)
!    sendbuf(4,:)=rc(LN,:)
!    call MPI_SEND(sendbuf,5*M,MPI_DBL,p_down,1,comm_group,err)
!    write(*,*) rank, "2" 
!    call MPI_RECV(recvbuf,5*M,MPI_DBL,p_down,1,comm_group,stat,err)
!    rcd=recvbuf(3,:)
!    rd=recvbuf(2,:)
!    ad=recvbuf(1,:)
!    if(iproc<nproc_x-1)then
!	r(LN,:)=r(LN,:)-rd*g(LN,:)/ad
!	rc(LN,:)=rc(LN,:)-rcd*g(LN,:)/ad
!    end if
!
!    x(LN,:)=r(LN,:)/a(LN,:)
!    xc(LN,:)=rc(LN,:)/a(LN,:)
!    do idx=LN-1,1,-1 !simd
!	r(idx,:)=r(idx,:)-r(LN,:)*g(idx,:)/a(LN,:)
!	rc(idx,:)=rc(idx,:)-rc(LN,:)*g(idx,:)/a(LN,:)
!	x(idx,:)=r(idx,:)/a(idx,:)
!	xc(idx,:)=rc(idx,:)/a(idx,:)
!    end do
!
!    sendbuf(1,:)=a(LN,:)
!    sendbuf(2,:)=r(LN,:)
!    sendbuf(3,:)=rc(LN,:)
!    call MPI_SEND(sendbuf,5*M,MPI_DBL,p_up,1,comm_group,err)
!
!    !! For periodic bound
!    if(iproc==0)then
!	sendbuf(1,:)=x(1,:)
!	sendbuf(2,:)=xc(1,:)
!	call MPI_SEND(sendbuf,5*M,MPI_DBL,np-1,1,comm_group,err)
!    end if
!    if(iproc==nproc_x-1)then
!	write(*,*) rank, "3" 
!	call MPI_RECV(recvbuf,5*M,MPI_DBL,0,1,comm_group,stat,err)
!	x(LN+1,:)=(r(LN+1,:)-c(LN+1,:)*recvbuf(1,:)-b(LN+1,:)*x(LN,:))/(a(LN+1,:)+c(LN+1,:)*recvbuf(2,:)+b(LN+1,:)*xc(LN,:))
!	sendbuf(1,:)=x(LN+1,:)
!    end if
!
!    call MPI_BCAST(sendbuf,5*M,MPI_DBL,np-1,comm_group,err)
!    
!    do idx=1,LN
!	x(idx,:)=x(idx,:)+sendbuf(1,:)*xc(idx,:)
!    end do
!
!    r=x
!
!    if(rank==np-1) LN=LN+1
!
!    return
!end subroutine
!
!
 subroutine Tri_diag(ra, rh, th2)
    use module_para
    use communicate
    use boundary
    use mpi
    use module_io
    use kinds_mod
    implicit none
    
    real(ld8),dimension(0:nloc_x+1,0:nloc_y+1), intent(in)    :: ra,rh
    real(ld8),dimension(0:nloc_x+1,0:nloc_y+1), intent(inout)   :: th2
    real(ld8),dimension(nloc_x/2,nloc_y)  ::  &
    fp,f0,fm,rg,gp,gm,g0,rf,rf_tmp,fp_tmp,fm_tmp,f0_tmp
    real(ld8),dimension(1:nloc_y) :: sendbuf,recvbuf 
    integer :: i,j,kn,i1,i2,ierr
    integer(int_kind) :: snd_req,rcv_req
    integer:: stat(MPI_STATUS_SIZE)
    kn = nloc_x / 2

 do j = 1,nloc_y
    do  i = 1, kn
	i1 = i*2 - 1
	i2 = i1+1
	fp(i,j) = -ra(i2,j)
	rf(i,j) = rh(i1,j)

	gm(i,j) = -ra(i1,j)
	rg(i,j) = rh(i2,j)
    enddo 

 enddo

 do j = 1,nloc_y
      sendbuf(j) = fp(kn,j)
 enddo
    call MPI_ISEND(sendbuf(1),nloc_y,MPI_DBL,e_proc,1,comm,snd_req,ierr)
    call MPI_IRECV(recvbuf(1),nloc_y,MPI_DBL,w_proc,1,comm,rcv_req,ierr)
    call MPI_WAIT(rcv_req,stat,ierr)
 do j =1,nloc_y 
    do i=2,kn
       fm(i,j)=fp(i-1,j)
    enddo
       fm(1,j)=recvbuf(j)
 enddo

 do j = 1,nloc_y
      sendbuf(j) = gm(1,j)
 enddo
    call MPI_ISEND(sendbuf(1),nloc_y,MPI_DBL,w_proc,2,comm,snd_req,ierr)
    call MPI_IRECV(recvbuf(1),nloc_y,MPI_DBL,e_proc,2,comm,rcv_req,ierr)
    call MPI_WAIT(rcv_req,stat,ierr)
 do j =1,nloc_y 
    do i=1,kn-1
	gp(i,j)=gm(i+1,j)
    enddo
	gp(kn,j)=recvbuf(j)
 enddo

 do j =1,nloc_y
    do i=1,kn
	f0(i,j)=1.0d0-fm(i,j)-fp(i,j)
	g0(i,j)=1.0d0-gm(i,j)-gp(i,j)
    enddo
 enddo

 f0_tmp(:,:) = g0(:,:)
 fm_tmp(:,:) = gm(:,:)
 fp_tmp(:,:) = gp(:,:)
 rf_tmp(:,:) = rg(:,:)
 !write(*,*),"====",my_task,"previous",rf(90,13),"======"
 call tridiagnol_solver3(f0,fm,fp,rf,kn,nloc_y)
 call tridiagnol_solver3(g0,gm,gp,rg,kn,nloc_y)
 
! write(*,*),"====",my_task,"later",rf_tmp(90,13),"====="
! do j = 1,nloc_y
! call check(f0_tmp(:,j),fm_tmp(:,j),fp_tmp(:,j),rf_tmp(:,j),rg(:,j),kn,kn*nproc_x,j)
! enddo

! call mpi_print2(rf)

 do j = 1,nloc_y
    do  i = 1, kn
	i1 = i*2 - 1
	i2 = i1+1
	th2(i1,j)=rf(i,j)
	th2(i2,j)=rg(i,j)
!	if(my_task .eq. 0) write(*,*) i1,j,th2(i1,j)
    enddo
 enddo
 
end subroutine

subroutine check(a,b,c,r,x,LN,N,g)
    use communicate
    implicit none
    include 'mpif.h'

    integer :: idx
    integer :: rank,np,err
    integer :: LN,N,g
    real(8) :: a(LN),b(LN),c(LN),r(LN),x(LN)
    real(8) :: aa(N),bb(N),cc(N),rr(N),xx(N)
    real(8) :: res,xn,xp !x_next x_previous
    logical :: flag

    call MPI_COMM_RANK(comm_group,rank,err)
    call MPI_COMM_SIZE(comm_group,np,err)

    !!! gather calculate send
    call MPI_GATHER(a,LN,MPI_DBL,aa,LN,MPI_DBL,0,comm_group,err)
    call MPI_GATHER(b,LN,MPI_DBL,bb,LN,MPI_DBL,0,comm_group,err)
    call MPI_GATHER(c,LN,MPI_DBL,cc,LN,MPI_DBL,0,comm_group,err)
    call MPI_GATHER(r,LN,MPI_DBL,rr,LN,MPI_DBL,0,comm_group,err)
    call MPI_GATHER(x,LN,MPI_DBL,xx,LN,MPI_DBL,0,comm_group,err)

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
	    if (res>1.0d-12) then
		print *,"error:","N:",N,"idx:",idx,"res:",res,"xx:", xx(idx),"b:", bb(idx),"c", cc(idx),"a", aa(idx),"r",rr(idx),"xp", xp,"xn", xn
		flag=.TRUE.
	    end if
	end do
	if (flag) then
	    print *,"group",g,"check failed"
!	    do idx=1,N
!		    print *,idx,xx(idx)
!	    end do
	else
	    print *, "group",g,"check pass!"
	    !do idx=1,N
	    !     print *,idx,xx(idx)
	    !end do
	end if
    end if

    return
end subroutine

subroutine tridiagnol_solver4(b,a,c,r,LN,M)
    use communicate
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

    call MPI_COMM_RANK(comm_group,rank,err)
    call MPI_COMM_SIZE(comm_group,np,  err)
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
	call MPI_SEND(sendbuf1,7*M,MPI_DOUBLE_PRECISION,0,1,comm_group,err)
    endif

    !!! Chase down and pass the last equation down
    if(rank/=0)then
	call MPI_RECV(recvbuf1,7*M,MPI_DOUBLE_PRECISION,p_up,1,comm_group,stat,err)
    else
	call MPI_RECV(recvbuf1,7*M,MPI_DOUBLE_PRECISION,np-1,1,comm_group,stat,err)
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
    call MPI_SEND(sendbuf1,7*M,MPI_DOUBLE_PRECISION,p_down,1,comm_group,err)
    !!! Handle the last equation, get xn
    if(rank==np-1)then
	a(LN+1,:)=a(LN+1,:)+cn
	q        =a(LN+1,:)/b(LN,:)
	b(LN+1,:)=bn-q*(c(LN,:)+a(LN,:))
	r(LN+1,:)=rn-q*r(LN,:)
	r(LN+1,:)=r(LN+1,:)/b(LN+1,:)
    endif
    !!! Chase back
    call MPI_RECV(recvbuf2,2*M,MPI_DOUBLE_PRECISION,p_down,1,comm_group,stat,err)
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
    call MPI_SEND(sendbuf2,2*M,MPI_DOUBLE_PRECISION,p_up,1,comm_group,err)

    if(rank==np-1) LN=LN+1
    return
end subroutine

