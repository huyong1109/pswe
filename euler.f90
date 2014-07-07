!
subroutine euler(dt,iter)

    use module_para
    use module_array
    use module_io
    use boundary
    use global_reductions
    

    implicit none

    integer i,i1,i2,j,k,iter, ierr

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: tu,tv,th,h 
    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)  :: du,dv,dh 

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)      :: a1,a2,rh,ru 
    real(r8),dimension(:), allocatable        :: fm,fp,f0,gm,gp,g0,rf,rg 
    real(r8),dimension(:,:),allocatable   ::  ga2, grh, gth

    real(r8)                      :: ai,aj,dt,dt2,en,en0,den

    real(r8),external             :: inner
    !	half timestep
    dt2=dt*0.5d0

    do j=1,nloc_y
	do i=1,nloc_x
	    tu(i,j)=wu(i,j)
	    tv(i,j)=wv(i,j)
	    th(i,j)=wh(i,j)
	end do
    end do

    !	Initial energy
    en0=inner(wu,wv,wh,wu,wv,wh)
    if (my_task == master_task) then
	write(*,*) 'en0 ', en0
    endif 

    !	Iterate until convergence ( |E-E0| < 1.0E-15 E0 )

    do k=1,1!1000
    if (my_task == master_task) then
	write(*,*) '===================== ',k, '============='
    endif
    call MPI_BARRIER(comm, ierr)


	do j=1,nloc_y
	    do i=1,nloc_x
		h(i,j)=dsqrt(th(i,j))
	    end do
	end do
	
	!	consider U, phi, compute advection and coriolios terms 

	call DIFUH(tu,tv,du,dh,h)
	call check_nan(du, 'du')
	call check_nan(dh, 'dh')
	!	get the tri-diagonal matrix
	do j=1,nloc_y
	    do i=1,nloc_x
		tu(i,j)=wu(i,j)-dt2*du(i,j)
		th(i,j)=wh(i,j)-dt2*dh(i,j)
	    end do
	end do
	a1(:,:) = 0.0
	a2(:,:) = 0.0
	do j=1,nloc_y
	    ai=dt2*c11(j)

	    do i=1,nloc_x
		aj=ai*h(i,j)
		a1(i,j)=aj
		ru(i,j)=tu(i,j)*aj

		a2(i,j)=aj*aj
	    end do
	end do 
	call check_nan(a1, 'a1')
		
	call mpi_print(a1)

	call update_boundary(ru)

	do j=1,nloc_y
	    do i=1,nloc_x
		rh(i, j)=th(i,j)-ru(i+1, j)+ru(i-1,j )
	    enddo
	enddo
	call Tri_diag(a2, rh, th)

	call update_boundary(th) 

	do j = 1, nloc_y
	    do i = 1, nloc_x
		tu(i,j)=tu(i,j)-a1(i,j)*(th(i+1,j)-th(i-1,j))
	    enddo
	enddo 
	!
	call DIFV(tu,tv,th,dv,h)
	!
	do j=1,nloc_y
	    do i=1, nloc_x
		tv(i,j)=wv(i,j)-dt2*dv(i,j)
	    enddo
	enddo
	!
	en=inner(tu,tv,th,tu,tv,th)
	write(*,*) 'en', k, en
	!
	den=dabs(en-en0)*2.0/(en+en0)
	en0=en
	if (den.lt.1.0d-15) goto 10

    enddo

    10   continue
    iter=k

    do j=1,nloc_x
	do i=1,nloc_y
	    wu(i,j)=tu(i,j)*2.0d0-wu(i,j)
	    wv(i,j)=tv(i,j)*2.0d0-wv(i,j)
	    wh(i,j)=th(i,j)*2.0d0-wh(i,j)
	end do
    enddo

    return
 end subroutine euler
 subroutine Tri_diag(ra, rh, th )
    use module_para
    use module_io
    use boundary
    use global_reductions
    

    implicit none

    integer i,i1,i2,j,k,iter, ierr

    real(r8),dimension(0:nloc_x+1,0:nloc_y+1)      :: a1,a2,rh,ru 
    real(r8),dimension(:), allocatable        :: fm,fp,f0,gm,gp,g0,rf,rg 
    real(r8),dimension(:,:),allocatable   ::  gra, grh, gth

    real(r8)                      :: ai,aj


	if (iproc == 0 ) then
	    allocate(gth(1:nloc_x, 1:nloc_y*nproc_x))
	    allocate(gra(1:nloc_x, 1:nloc_y*nproc_x))
	    allocate(grh(1:nloc_x, 1:nloc_y*nproc_x))
	    allocate(fm(1:kn))
	    allocate(fp(1:kn))
	    allocate(f0(1:kn))
	    allocate(gm(1:kn))
	    allocate(gp(1:kn))
	    allocate(g0(1:kn))
	    allocate(rf(1:kn))
	    allocate(rg(1:kn))


	    gth(:,:) = 0.0
	    gra(:,:) = 0.0
	    grh(:,:) = 0.0
	endif 

	!if (iproc ==0 ) then 
	!    do j=1,nloc_y
	!	do i=1,p
	!	   write(*,*) my_task, i, j, ga1(i,j) 
	!	enddo
	!    enddo
	!endif 
	call lat_gather(gra, ra)	
	call lat_gather(grh, rh)	

	!	grid divided into two groups. 
	if (iproc == 0 ) then
	    fm(1:kn)= 0.0 
	    fp(1:kn)= 0.0 
	    f0(1:kn)= 0.0 
	    gm(1:kn)= 0.0 
	    gp(1:kn)= 0.0 
	    g0(1:kn)= 0.0 
	    rf(1:kn)= 0.0 
	    rg(1:kn)= 0.0 

	    do j = 1, nloc_y
		do i=1,kn
		    i1=i*2
		    i2=i1+1

		    fp(i)=-ga2(mod(i2, nloc_x), j+i2/nloc_x*nloc_y)
		    rf(i)=grh(mod(i1, nloc_x),j+i1/nloc_x*nloc_y)

		    gm(i)=-ga2(mod(i1, nloc_x),j+i1/nloc_x*nloc_y)
		    rg(i)=grh(mod(i2, nloc_x), j+i2/nloc_x*nloc_y)
		enddo

		do i=2,kn
		    fm(i)=fp(i-1)
		enddo
		fm(1)=fp(kn)

		do i=1,kn-1
		    gp(i)=gm(i+1)
		enddo
		gp(kn)=gm(1)

		do i=1,kn
		    f0(i)=1.0-fm(i)-fp(i)
		    g0(i)=1.0-gm(i)-gp(i)
		enddo

		!	    solve two tri-diagonal linear system
		call LU0(fm,f0,fp,rf,kn)
		call LU0(gm,g0,gp,rg,kn)


		do i=1,kn
		    i1=i*2
		    i2=i1+1

		    gth(mod(i1, nloc_x),j+i1/nloc_x*nloc_y)=rf(i)
		    gth(mod(i2, nloc_x), j+i2/nloc_x*nloc_y)=rg(i)
		enddo
	    enddo

	    deallocate(fm)
	    deallocate(fp)
	    deallocate(f0)
	    deallocate(gm)
	    deallocate(gp)
	    deallocate(g0)
	    deallocate(rf)
	    deallocate(rg)

	endif

	!write(*,*) 'lat_gather gth'
	call lat_scatter(gth, th)	
	if (iproc == 0 ) then
	    deallocate(gth)
	    deallocate(ga2)
	    deallocate(grh)
	endif 
 
 end subroutine Tri_diag
    
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


