!! TO COMPILE : 
!! gfortran -o calc_extremes.x -I/opt/local/include libs/ncio.f90 libs/index.f90 calc_extremes.f90 -L/opt/local/lib -lnetcdff -lnetcdf

program calc_extremes 

    use ncio 
    use index 

    implicit none 


    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(wp), parameter :: ERR_DIST = real(1E8,wp) 
    integer,    parameter :: ERR_IND  = -1 
    real(wp), parameter :: tol_underflow = real(1E-15,wp)

    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    

    type dataset_class 
        real(wp), allocatable :: lon(:)
        real(wp), allocatable :: lat(:)
        real(wp), allocatable :: tas(:,:,:,:)
    end type 

    type(dataset_class) :: dat 

    character(len=512) :: filename_in 
    character(len=512) :: filename_out 

    ! Test time series calculations 
    call test_timeseries("test.nc",n=100,mu=0.0_wp,sigma=2.0_wp,alpha=0.0_wp)
    stop 


    filename_in  = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1.nc"
    filename_out = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats.nc"

    ! Load data 
    call load_best(dat,filename_in,year0=1850,year1=2020,mv=mv)

    ! Write output
    call write_extremes_init(filename_out,dat%lon,dat%lat,year0=1850,year1=2020)
    call nc_write(filename_out,"tas",dat%tas,dim1="lon",dim2="lat",dim3="month",dim4="year")

contains 

subroutine test_timeseries(filename,n,mu,sigma,alpha)

    implicit none 

    character(len=*), intent(IN) :: filename 
    integer, intent(IN) :: n 

    real(wp), intent(IN) :: mu 
    real(wp), intent(IN) :: sigma 
    real(wp), intent(IN) :: alpha 

    ! Local variables 
    integer :: i, j
    real(wp) :: tmp  
    real(wp) :: x(n)
    real(wp) :: y(n)
    real(wp) :: ysm(n)
    real(wp) :: mean, stdev 

    call random_seed
    
    ! Get a random time series
    do i = 1, n
        x(i) = real(i,wp) 
        !call random_number(y(i))
        call random_number_normal(y(i),mu,sigma) 
        y(i) = y(i) + (x(i)-1.0)*alpha
    end do 

    ! Add missing values 
    do i = 1, int(n*0.2)
        call random_number(tmp)     ! Random number from 0-1
        j = floor(tmp*(n-1))+1      ! Random number from 1 to n
        y(j) = mv
    end do 

    ! Test running mean smoothing 
    call smooth_runmean(ysm,y,x,L=15.0_wp,mv=mv)

    ! Test standard deviation
    call calc_mean(mean,y,mv)
    call calc_stdev(stdev,y,mv)

    write(*, "(a, i0)") "sample size = ", n
    write(*, "(a, f17.15)") "Mean :   ", mean
    write(*, "(a, f17.15)") "Stddev : ", stdev  

    call nc_create(filename)
    call nc_write_dim(filename,"x",x)
    call nc_write(filename,"y",y,dim1="x")
    call nc_write(filename,"ysm",ysm,dim1="x")

    stop 

    return 

end subroutine test_timeseries


subroutine load_best(dat,filename,year0,year1,mv)

    implicit none 

    type(dataset_class), intent(OUT) :: dat 
    character(len=*), intent(IN) :: filename 
    integer, intent(IN) :: year0 
    integer, intent(IN) :: year1 
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: nx, ny, nt, nyr, nm 
    real(wp), allocatable :: tas3D(:,:,:)
    
    write(*,"(a)",advance="no") "load_best..."

    nx = nc_size(filename_in,"longitude")
    ny = nc_size(filename_in,"latitude")
    nt = nc_size(filename_in,"time")

    allocate(dat%lon(nx))
    allocate(dat%lat(ny)) 
    allocate(tas3D(nx,ny,nt))

    call nc_read(filename,"longitude",dat%lon)
    call nc_read(filename,"latitude", dat%lat)
    
    call nc_read(filename,"temperature",tas3D,missing_value=mv)

    write(*,*) "done."
    
    ! Reshape best data to 4D array
    call reshape_to_4D(dat%tas,tas3D,month0=1,mv=mv)


    return 

end subroutine load_best

subroutine reshape_to_4D(var4D,var3D,month0,mv)

    implicit none 

    real(wp), allocatable, intent(OUT) :: var4D(:,:,:,:) 
    real(wp), intent(IN)  :: var3D(:,:,:) 
    integer,  intent(IN)  :: month0 
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: i, j, y, m, k 
    integer :: nx, ny, nt, nyr, nm 

    nx = size(var3D,1)
    ny = size(var3D,2)
    nt = size(var3D,3)
    nm = 12 

    nyr = ceiling(real(nt,wp) / real(nm,wp))

    if (allocated(var4D)) deallocate(var4D)
    allocate(var4D(nx,ny,nm,nyr))

    var4D = mv 

    k = 0 

    write(*,"(a)",advance="no") "reshape_to_4D..."

    do y = 1, nyr 
    do m = 1, nm 

        if ( (y .eq. 1 .and. m .ge. month0) .or. y .gt. 1 ) then 
            
            k = k+1 
            var4D(:,:,m,y) = var3D(:,:,k)

            if (k .eq. nt) exit 
        end if 

    end do 
    end do 

    write(*,*) "done."

    return 

end subroutine reshape_to_4D


subroutine calc_mean(mean,var,mv)
    ! Calculat the standard deviation of a
    ! vector of values, excluding missing values.

    implicit none 

    real(wp), intent(OUT) :: mean 
    real(wp), intent(IN)  :: var(:)
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: i, n
    real(wp) :: nsamples 
    real(wp) :: var_sum

    n = size(var,1) 

    nsamples   = 0.0
    var_sum    = 0.0 

    do i = 1, n 
        if (var(i) .ne. mv) then 
            nsamples = nsamples + 1
            var_sum   = var_sum + var(i) 
        end if 
    end do 

    if (nsamples .gt. 0.0) then 
        ! If samples exist, calculate standard deviation 

        mean   = var_sum / nsamples

    else 
        ! Set standard deviation equal to missing value 

        mean = mv 

    end if 

    return 

end subroutine calc_mean

subroutine calc_stdev(stdev,var,mv)
    ! Calculat the standard deviation of a
    ! vector of values, excluding missing values.

    implicit none 

    real(wp), intent(OUT) :: stdev 
    real(wp), intent(IN)  :: var(:)
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: i, n
    real(wp) :: nsamples 
    real(wp) :: var_sum, var_sum_sq
    real(wp) :: mean 


    n = size(var,1) 

    nsamples   = 0.0
    var_sum    = 0.0 
    var_sum_sq = 0.0 

    do i = 1, n 
        if (var(i) .ne. mv) then 
            nsamples = nsamples + 1
            var_sum   = var_sum + var(i) 
            var_sum_sq = var_sum_sq + var(i)*var(i) 
        end if 
    end do 

    if (nsamples .gt. 0.0) then 
        ! If samples exist, calculate standard deviation 

        mean   = var_sum / nsamples
        stdev = sqrt(var_sum_sq/nsamples - mean*mean)

    else 
        ! Set standard deviation equal to missing value 

        stdev = mv 

    end if 

    return 

end subroutine calc_stdev

subroutine random_number_normal(val,mu,sigma)
    ! Given a random number sampled from a uniform distribution
    ! (ie, using `random_number`), transform number to a sample
    ! from a normal distribution with mean mu and stdev sigma.
    ! Uses the Box-Muller transform:
    ! https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

    implicit none 

    real(wp), intent(OUT) :: val 
    real(wp), intent(IN)  :: mu 
    real(wp), intent(IN)  :: sigma 
    
    ! Local variables 
    real(wp) :: u1, u2, z0 
    real(wp), parameter :: two_pi = 2.0_wp * pi 
    
    ! Note, random seed should be initialized externally
    !call random_seed

    call random_number(u1)
    call random_number(u2)

    ! Normal distribution with (0,1)
    z0 = sqrt(-2.0_wp * log(u1)) * cos(two_pi * u2)

    ! Apply our mean and stdev
    val = z0 * sigma + mu

    return 

end subroutine random_number_normal

subroutine smooth_runmean(ysm,y,x,L,mv)
    ! Calculate smooth vector using the running
    ! mean, with a window half-length of L 
    ! (smoothing window is L*2+1), accounting for 
    ! missing values

    implicit none 

    real(wp), intent(OUT) :: ysm(:) 
    real(wp), intent(IN)  :: y(:)
    real(wp), intent(IN)  :: x(:) 
    real(wp), intent(IN)  :: L 
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: i, n, n_now 
    integer, allocatable :: idx(:) 

    n = size(y,1) 

    ! Initialize output vector with missing values 
    ysm = mv 

    do i = 1, n
        call which(x .ge. x(i)-L .and. x .le. x(i)+L .and. y .ne. mv,idx,n_now)
        if (n_now .gt. 0) then 
            ysm(i) = sum(y(idx)) / real(n_now,wp)
        end if 
    end do  

    return 

end subroutine smooth_runmean

subroutine define_time_vectors()

    implicit none 



    return 

end subroutine define_time_vectors



subroutine write_extremes_init(filename,lon,lat,year0,year1)

    implicit none 

    character(len=*), intent(IN) :: filename 
    real(wp), intent(IN) :: lon(:)
    real(wp), intent(IN) :: lat(:)
    integer,  intent(IN) :: year0 
    integer,  intent(IN) :: year1 
    
    call nc_create(filename)
    call nc_write_dim(filename,"lon",lon,units="degrees_east", &
                        long_name="Longitude")
    call nc_write_dim(filename,"lat",lat,units="degrees_north", &
                        long_name="Latitude")
    
    call nc_write_dim(filename,"month",x=1,dx=1,nx=12, & 
                        units="month",long_name="Month of the year")
    call nc_write_dim(filename,"seas",x=1,dx=1,nx=5, & 
                        units="season",long_name="Season number: 1=[12,1,2], &
                        &2=[3,4,5], 3=[6,7,8], 4=[9,10,11], 5=[1:12]")
    call nc_write_dim(filename,"year",x=year0,dx=1,nx=(year1-year0+1), & 
                        units="year",long_name="Calendar year")
    
    return 

end subroutine write_extremes_init

    
end program calc_extremes

