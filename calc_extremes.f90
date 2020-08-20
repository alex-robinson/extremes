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
        integer :: nx, ny, nm, ns, nyr, nt, nsig 
        real(wp), allocatable :: lon(:)
        real(wp), allocatable :: lat(:)
        real(wp), allocatable :: month(:)
        real(wp), allocatable :: seas(:) 
        real(wp), allocatable :: year(:) 
        real(wp), allocatable :: sigma(:) 
        real(wp), allocatable :: tas3D(:,:,:)
        real(wp), allocatable :: tas(:,:,:,:)
        real(wp), allocatable :: tas_ref(:,:,:)
        real(wp), allocatable :: tas_sm(:,:,:,:)
        real(wp), allocatable :: tas_detrnd(:,:,:,:)
        real(wp), allocatable :: tas_sd(:,:,:)
        real(wp), allocatable :: tas_sigma(:,:,:,:)
    end type 

    type(dataset_class) :: dat 

    character(len=512) :: filename_in 
    character(len=512) :: filename_out 

    ! Test time series calculations 
    !call test_timeseries("test.nc",n=100,mu=0.0_wp,sigma=2.0_wp,alpha=0.0_wp)
    !stop 


    filename_in  = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1.nc"
    filename_out = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats.nc"

    ! Load data including dimension info (nx,ny)
    call load_best(dat,filename_in,year0=1850,year1=2020,mv=mv)

    ! Perform stats calculations
    call stats_calc_1(dat,L=15.0_wp,mv=mv)
    

    ! Write output
    call write_extremes_init(filename_out,dat%lon,dat%lat,year0=1850,year1=2020)
    call nc_write(filename_out,"tas",dat%tas,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out,"tas_sm",dat%tas_sm,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    !call nc_write(filename_out,"tas_detrnd",dat%tas_detrnd,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out,"tas_sd",dat%tas_sd,dim1="lon",dim2="lat",dim3="month",missing_value=mv)
    call nc_write(filename_out,"tas_sigma",dat%tas_sigma,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)

contains

subroutine stats_calc_1(dat,L,mv)
    ! Given a 4D dataset (nx,ny,nm,nyr), perform 
    ! initial steps of statistics calculations. 

    implicit none 

    type(dataset_class), intent(INOUT) :: dat 
    real(wp), intent(IN) :: L 
    real(wp), intent(IN) :: mv 

    ! Local variables 
    integer :: i, j, m, y 
    integer :: nx, ny, nm, nyr, n 
    integer, allocatable :: idx_ref(:) 
    integer, allocatable :: idx_sd(:) 
    integer, allocatable :: idx_numeric(:) 

    integer :: idx_mm(12,3)
    integer :: m1, m2, m3 

    nx  = size(dat%tas,1)
    ny  = size(dat%tas,2)
    nm  = size(dat%tas,3) 
    nyr = size(dat%tas,4)
    
    ! Populate month indices to get 
    ! three months surrounding a given month 
    do m = 1, nm 
        idx_mm(m,:) = [-1,0,1] + m
    end do 
    where(idx_mm .eq. 0)  idx_mm = 12 
    where(idx_mm .eq. 13) idx_mm = 1 
    
    write(*,*) "stats_calc_1..."

    ! Get indices of year ranges
    call which(dat%year .ge. 1951 .and. dat%year .le. 2010,idx_sd)

    ! Initialize all variables with missing values 
    dat%tas_sm     = mv 
    dat%tas_detrnd = mv 
    dat%tas_sd     = mv 

    ! do i = int(nx/4), 2*int(nx/4)
    !     write(*,*) "stats_calc_1...", i, "/", nx 
    ! do j = int(ny/4), 2*int(ny/2)

    do i = 1, nx 
        write(*,*) "stats_calc_1...", i, "/", nx 
    do j = 1, ny 

        ! Perform some initial monthly calculations
        do m = 1, nm 

            ! Calculate reference climate 
            call which(dat%year .ge. 1951 .and. dat%year .le. 1980 .and. &
                        dat%tas(i,j,m,:) .ne. mv, idx_ref,n)
            dat%tas_ref(i,j,m) = mv
            if (n .gt. 0) then 
                dat%tas_ref(i,j,m) = sum(dat%tas(i,j,m,idx_ref)) / real(n,wp)
            end if
            
            ! Calculate smooth time series 
            call smooth_runmean(dat%tas_sm(i,j,m,:),dat%tas(i,j,m,:),dat%year,L,mv)
            
            ! Calculate detrended data
            call which(dat%tas(i,j,m,:) .ne. mv,idx_numeric,n)
            dat%tas_detrnd(i,j,m,:) = mv
            if (n .gt. 0) then 
                dat%tas_detrnd(i,j,m,idx_numeric) = dat%tas(i,j,m,idx_numeric) - dat%tas_sm(i,j,m,idx_numeric) 
            end if
              
        end do 

        ! With reference, smoothed and detrended data available for all months,
        ! perform additional monthly calculations 
        do m = 1, nm 
            m1 = idx_mm(m,1)
            m2 = idx_mm(m,2)
            m3 = idx_mm(m,3)
            
            ! Perform standard deviation calculations over multiple months 
            ! for the year range of interest
            call calc_stdev(dat%tas_sd(i,j,m),[ dat%tas_detrnd(i,j,m1,idx_sd), &
                                                dat%tas_detrnd(i,j,m2,idx_sd), &
                                                dat%tas_detrnd(i,j,m3,idx_sd) ], mv=mv)

            ! Calculate the normalized temp anomaly
            call calc_tas_sigma(dat%tas_sigma(i,j,m,:),dat%tas(i,j,m,:), &
                                        dat%tas_ref(i,j,m),dat%tas_sd(i,j,m),mv)

        end do

    end do 
    end do 

    write(*,*) "stats_calc_1... done."

    return 

end subroutine stats_calc_1

subroutine calc_tas_sigma(tas_sigma,tas,tas_ref,tas_sd,mv)

    implicit none 

    real(wp), intent(OUT) :: tas_sigma(:) 
    real(wp), intent(IN)  :: tas(:) 
    real(wp), intent(IN)  :: tas_ref 
    real(wp), intent(IN)  :: tas_sd 
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: i 

    tas_sigma = mv 

    if (tas_ref .ne. mv .and. tas_sd .ne. mv) then 
        where(tas .ne. mv) tas_sigma = (tas - tas_ref) / tas_sd 
    end if

    return 

end subroutine calc_tas_sigma

subroutine load_best(dat,filename,year0,year1,mv)

    implicit none 

    type(dataset_class), intent(OUT) :: dat 
    character(len=*), intent(IN) :: filename 
    integer,  intent(IN)  :: year0 
    integer,  intent(IN)  :: year1
    real(wp), intent(IN)  :: mv 

    ! Local variables  
    real(wp), allocatable :: tas3D(:,:,:)
    
    write(*,"(a)",advance="no") "load_best..."

    dat%nx = nc_size(filename_in,"longitude")
    dat%ny = nc_size(filename_in,"latitude")
    dat%nt = nc_size(filename_in,"time")

    allocate(dat%lon(dat%nx))
    allocate(dat%lat(dat%ny)) 
    allocate(tas3D(dat%nx,dat%ny,dat%nt))

    call nc_read(filename,"longitude",dat%lon)
    call nc_read(filename,"latitude", dat%lat)
    
    call nc_read(filename,"temperature",tas3D,missing_value=mv)

    write(*,*) "done."
    
    ! Allocate dataset variables to prepare for populating them.
    call dataset_alloc(dat,year0=1850,year1=2020)

    ! Reshape data to 4D array
    call reshape_to_4D(dat%tas,tas3D,month0=1,mv=mv)

    return 

end subroutine load_best

subroutine dataset_alloc(dat,year0,year1)

    implicit none 

    type(dataset_class), intent(INOUT) :: dat 
    integer, intent(IN) :: year0 
    integer, intent(IN) :: year1 
    
    ! Local variables
    integer :: k 
    integer :: nyr, nm, ns  

    dat%nm   = 12
    dat%ns   = 5 
    dat%nyr  = year1 - year0 + 1
    dat%nsig = 5

    ! First ensure variables are deallocated 
    call dataset_dealloc(dat)

    ! Allocate and populate dimensions
    allocate(dat%month(dat%nm))
    allocate(dat%seas(dat%ns))
    allocate(dat%year(dat%nyr))
    allocate(dat%sigma(dat%nsig))

    do k = 1, dat%nm 
        dat%month(k) = k 
    end do 

    do k = 1, dat%ns 
        dat%seas(k) = k 
    end do 
    
    do k = 1, dat%nyr 
        dat%year(k) = year0 + (k-1)
    end do 
    
    dat%sigma = [1,2,3,4,5] 

    ! Allocate variables to correct size
    allocate(dat%tas(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_ref(dat%nx,dat%ny,dat%nm))
    allocate(dat%tas_sm(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_detrnd(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_sd(dat%nx,dat%ny,dat%nm))
    allocate(dat%tas_sigma(dat%nx,dat%ny,dat%nm,dat%nyr))

    return 

end subroutine dataset_alloc

subroutine dataset_dealloc(dat)

    implicit none 

    type(dataset_class), intent(INOUT) :: dat 

    ! deallocate variables 
    if (allocated(dat%tas))         deallocate(dat%tas)
    if (allocated(dat%tas_ref))     deallocate(dat%tas_ref)
    if (allocated(dat%tas_sm))      deallocate(dat%tas_sm)
    if (allocated(dat%tas_detrnd))  deallocate(dat%tas_detrnd)
    if (allocated(dat%tas_sd))      deallocate(dat%tas_sd)
    if (allocated(dat%tas_sigma))   deallocate(dat%tas_sigma)

    return 

end subroutine dataset_dealloc

subroutine reshape_to_4D(var4D,var3D,month0,mv)

    implicit none 

    real(wp), intent(OUT) :: var4D(:,:,:,:) 
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
    nyr = size(var4D,4)

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

