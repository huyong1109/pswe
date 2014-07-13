!
subroutine euler(dt,iter)

    use module_para
    use communicate, only:master_print_message
    use module_array, only: wu, wv, wh, c11
    use module_io
    use boundary, only: update_boundary, update_latitude
    use global_reductions, only: global_sum

    implicit none

    integer i,i1,i2,j,k,iter, ierr, break_cir

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: tu,tv,th,h 
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: du,dv,dh 

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)      :: a1,a2,rh,ru 

    real(r8)                      :: ai,aj,dt,dt2,en,en0,den

    real(r8),external             :: inner

    !	half timestep
    dt2=dt*0.5d0
    call update_boundary(wu)
    call update_boundary(wv)
    call update_boundary(wh)
    do j=0,nloc_y+1
	do i=0,nloc_x+1
	    tu(i,j)=wu(i,j)
	    tv(i,j)=wv(i,j)
	    th(i,j)=wh(i,j)
	end do
    end do

    !	Initial energy
    en0=inner(wu,wv,wh,wu,wv,wh)

    !	Iterate until convergence ( |E-E0| < 1.0E-15 E0 )
    break_cir = 0
    do k=1,1000

	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		h(i,j)=dsqrt(th(i,j))
	    end do
	end do
	
	!	consider U, phi, compute advection and coriolios terms 

	call DIFUH(tu,tv,du,dh,h)
	!	get the tri-diagonal matrix
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		tu(i,j)=wu(i,j)-dt2*du(i,j)
		th(i,j)=wh(i,j)-dt2*dh(i,j)
	    end do
	end do
	
	a1(:,:) = 0.0
	a2(:,:) = 0.0
	do j=1,nloc_y
	    ai=dt2*c11(j)

	    do i=0,nloc_x+1
		aj=ai*h(i,j)
		a1(i,j)=aj
		ru(i,j)=tu(i,j)*aj
		a2(i,j)=aj*aj
	    end do
	end do 
	call update_boundary(ru) 
		
	do j=0,nloc_y+1
	    do i=0,nloc_x+1
		rh(i, j)=th(i,j)-ru(i+1, j)+ru(i-1,j )
	    enddo
	enddo
	call  Tri_diag(a2, rh, th)

	call update_boundary(th) 

	do j = 1, nloc_y
	    do i = 1, nloc_x
		tu(i,j)=tu(i,j)-a1(i,j)*(th(i+1,j)-th(i-1,j))
	    enddo
	enddo 

	call update_boundary(tu) 
	
	call DIFV(tu,tv,th,dv,h)
	
	do j=1,nloc_y
	    do i=1, nloc_x
		tv(i,j)=wv(i,j)-dt2*dv(i,j)
	    enddo
	enddo
	
	call update_boundary(tv) 
	
	en=inner(tu,tv,th,tu,tv,th)
	
	den=dabs(en-en0)*2.0/(en+en0)
	!if (my_task .eq. master_task ) write(*,*) k ,  en
	en0=en
	if (den.lt.1.0d-15) break_cir = 1
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

 end subroutine euler

 subroutine Tri_diag(ra, rh, th )
    use kinds_mod, only: r8, int_kind
    use distribution, only: nloc_x, nloc_y, jproc, iproc
    use module_para, only: nproc_x
    use communicate
    use module_io, only: check_nan
    use boundary , only: update_boundary, update_latitude, e_proc, w_proc
    use global_reductions, only:  lat_gather, lat_scatter
    use mpi


    implicit none

    integer i,i1,i2,j,k,hn,iter,ierr

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(in)    :: ra,rh
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1), intent(inout)   :: th
    real(r8),dimension(1:nloc_x/2,1:nloc_y)	   :: fm,fp,f0,gm,gp,g0,rf,rg 
    real(r8),dimension(1:nloc_y) :: sendbuf,recvbuf 
    integer(int_kind) :: snd_req,rcv_req
    integer:: stat(MPI_STATUS_SIZE)
    real(r8)                      :: ai,aj

    hn = nloc_x/2

    !	grid divided into two groups. 
    do j = 1, nloc_y
	do i=1,hn
	    i1=i*2-1
	    i2=i1+1

	    fp(i,j) = -ra(i2,j)
	    rf(i,j) = rh(i1,j)

	    gm(i,j) = -ra(i1,j)
	    rg(i,j) = rh(i2,j)
	enddo
    enddo 


 do j =1,nloc_y 
    do i=2,hn-1
       fm(i,j)=fp(i-1,j)
       gp(i,j)=gm(i+1,j)
    enddo
 enddo

 do j =1,nloc_y 
       fm(hn,j)=fp(hn-1,j)
       fm(1,j)=fp(hn,j)
       gp(1,j)=gm(2,j)
       gp(hn,j)=gm(1,j)
 enddo

 call update_latitude(fm, gp)

    do j =1,nloc_y 
	do i=1,hn
	    f0(i,j)=1.0-fm(i,j)-fp(i,j)
	    g0(i,j)=1.0-gm(i,j)-gp(i,j)
	enddo
    enddo


    call tridiagnol_solver5(f0(:,:),fm(:,:),fp(:,:),rf(:,:), hn, nloc_y, nproc_x)
    call tridiagnol_solver5(g0(:,:),gm(:,:),gp(:,:),rg(:,:), hn, nloc_y, nproc_x)

    
    do j=1,nloc_y
	do i=1,hn
	    i1=i*2-1
	    i2=i1+1

	    th(i1,j)=rf(i,j)
	    th(i2,j)=rg(i,j)
 	!if(my_task .eq. 0) write(*,*) i1,j,th(i1,j)
	enddo
    enddo
 
 end subroutine Tri_diag
 
 subroutine LU(a,b,c,r,m,n)
    !Gauss elimination method/Chasing method  with periodic boundary
    ! on Wang's book P159
    use kinds_mod
    implicit none
    !
    integer                 :: i,j, n, m
    real(r8) ,dimension(1:m)     :: sn,rn
    real(r8)		         :: ai
    !
    real(r8),dimension(1:n, 1:m) :: a,b,c,r
    real(r8),dimension(1:n, 1:m) :: s,t
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

    return
    
 end subroutine LU
    
 subroutine LU0(a,b,c,r,n)
    !Gauss elimination method/Chasing method  with periodic boundary
    ! on Wang's book P159
    use kinds_mod
    implicit none
    !
    integer                 :: i,n
    real(r8)                :: ai,sn,rn
    !
    real(r8),dimension(1:n) :: a,b,c,r
    real(r8),dimension(1:n) :: s,t
    !
    s(1)=a(1)  ! s = \tilde(a)
    t(1)=c(n)  ! t = \tilde(c)
    !
    sn=0.0	   ! \tilde(b)_N
    rn=0.0	   ! \tilde(r)_N
    !
    do i=2,n-1
	ai=a(i)/b(i-1)
	b(i)=b(i)-ai*c(i-1)
	r(i)=r(i)-ai*r(i-1)
	s(i)=-ai*s(i-1)
	!
	ai=t(i-1)/b(i-1)
	t(i)=-ai*c(i-1)
	sn  =sn-ai*s(i-1)
	rn  =rn-ai*r(i-1)
    enddo
    !
    a(n)=a(n)+t(n-1)
    b(n)=b(n)+sn
    c(n-1)=c(n-1)+s(n-1)
    r(n)=r(n)+rn
    !
    ai=a(n)/b(n-1)
    b(n)=b(n)-ai*c(n-1)
    r(n)=r(n)-ai*r(n-1)
    !
    r(n)=r(n)/b(n)
    r(n-1)=(r(n-1)-c(n-1)*r(n))/b(n-1)
    !
    do i=n-2,1,-1
	ai=r(i)-s(i)*r(n)-c(i)*r(i+1)
	r(i)=ai/b(i)
    enddo
    !
    return
    !
 end subroutine LU0


! 1, the equation should be like:
!
!	     |b c     a| 
!  	     |a b c    |
!  	     |  a b c  | * x = r
!  	     |    a b c|
!  	     |c     a b|
subroutine tridiagnol_solver5(b,a,c,r,LN,M,np)
    use communicate
    use mpi
    
    implicit none
    !include 'mpif.h'

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

    call MPI_COMM_RANK(comm_group,rank,err)
    !call MPI_COMM_SIZE(comm_group,np,err)

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
    call MPI_GATHER(sendbuf,8*M,MPI_DOUBLE_PRECISION,recvbuf,8*M,MPI_DOUBLE_PRECISION,np-1,comm_group,err)
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
    
    call MPI_SCATTER(sendbuf2,2*M,MPI_DOUBLE_PRECISION,recvbuf2,2*M,MPI_DOUBLE_PRECISION,np-1,comm_group,err)
    r(1,:)=recvbuf2(:,1)
    r(LN,:)=recvbuf2(:,2)

    do i=2,LN-1
	r(i,:)=(r(i,:)-g(i,:)*r(LN,:)-f(i,:)*r(1,:))/b(i,:)
    enddo

    return
end subroutine


























































































































