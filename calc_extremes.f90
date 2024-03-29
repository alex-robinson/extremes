!! TO COMPILE : 
!! gfortran -o calc_extremes.x -I/opt/local/include libs/ncio.f90 libs/index.f90 calc_extremes.f90 -L/opt/local/lib -lnetcdff -lnetcdf

program calc_extremes 

    use ncio 
    use index 

    use ieee_arithmetic  ! for nan-checks 

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
        
        ! Dimension sizes
        integer :: nx, ny, nm, ns, nyr, nt, nsig 

        ! Dimensions
        real(wp), allocatable :: lon(:)
        real(wp), allocatable :: lat(:)
        real(wp), allocatable :: month(:)
        real(wp), allocatable :: seas(:) 
        real(wp), allocatable :: year(:) 
        real(wp), allocatable :: sigma(:) 

        ! Grid variables 
        real(wp), allocatable :: cell_wt(:,:)
        real(wp), allocatable :: f_land(:,:)
        real(wp), allocatable :: wt70land(:,:)
        
        ! Variables
        real(wp), allocatable :: tas3D(:,:,:)
        real(wp), allocatable :: tas(:,:,:,:)
        real(wp), allocatable :: tas_ref(:,:,:)
        real(wp), allocatable :: tas_sm(:,:,:,:)
        real(wp), allocatable :: tas_detrnd(:,:,:,:)
        real(wp), allocatable :: tas_sd(:,:,:)
        real(wp), allocatable :: tas_sigma(:,:,:,:)

        real(wp), allocatable :: tas_rec_hi(:,:,:,:)
        real(wp), allocatable :: tas_rec_hi_lin(:,:,:,:)
        real(wp), allocatable :: tas_rec_hi_1880(:,:,:,:)
        real(wp), allocatable :: tas_lin_m(:,:,:) 
        real(wp), allocatable :: tas_lin_r(:,:,:)

        real(wp), allocatable :: frac_sigma(:,:,:) 
        
    end type 

    type(dataset_class) :: dat 

    character(len=512) :: filename_in 
    character(len=512) :: filename_out_0
    character(len=512) :: filename_out_1 
    character(len=512) :: filename_out_2

    logical  :: load_tas 
    logical  :: load_stats_1 

    integer  :: year0 
    integer  :: year1

    ! Test time series calculations 
    ! call test_timeseries("test.nc",n=1000,mu=0.0_wp,sigma=2.0_wp,alpha=0.1_wp)
    ! stop 

    ! filename_in    = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1.nc"
    ! filename_out_0 = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats_0.nc"
    ! filename_out_1 = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats_1.nc"
    ! filename_out_2 = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats_2.nc"

    filename_in    = "data/BerkeleyEarth/2021-02_BEST/Complete_TAVG_LatLong1.nc"
    filename_out_0 = "data/BerkeleyEarth/2021-02_BEST/Complete_TAVG_LatLong1_stats_0.nc"
    filename_out_1 = "data/BerkeleyEarth/2021-02_BEST/Complete_TAVG_LatLong1_stats_1.nc"
    filename_out_2 = "data/BerkeleyEarth/2021-02_BEST/Complete_TAVG_LatLong1_stats_2.nc"

    ! Options 
    load_tas        = .TRUE.
    load_stats_1    = .FALSE. 

    year0           = 1750 
    year1           = 2020 

    ! Do not load original data if loading stats_1 
    if (load_stats_1) load_tas = .FALSE. 

    ! Load data including dimension info (nx,ny)
    call load_best(dat,filename_in,year0=year0,year1=year1,mv=mv,load_tas=load_tas)

if (load_stats_1) then 

    ! To do:
    !call stats_read_1(dat,filename_out_1)

    call nc_read(filename_out_0,"cell_wt",        dat%cell_wt)
    call nc_read(filename_out_0,"f_land",         dat%f_land)
    call nc_read(filename_out_0,"wt70land",       dat%wt70land)
    
    call nc_read(filename_out_0,"tas",            dat%tas,        missing_value=mv)
    call nc_read(filename_out_0,"tas_sm",         dat%tas_sm,     missing_value=mv)
    call nc_read(filename_out_0,"tas_detrnd",     dat%tas_detrnd, missing_value=mv)
    
    call nc_read(filename_out_1,"tas_lin_m",      dat%tas_lin_m,      missing_value=mv)
    call nc_read(filename_out_1,"tas_lin_r",      dat%tas_lin_r,      missing_value=mv)
    call nc_read(filename_out_1,"tas_sd",         dat%tas_sd,         missing_value=mv)
    call nc_read(filename_out_1,"tas_sigma",      dat%tas_sigma,      missing_value=mv)
    call nc_read(filename_out_1,"tas_rec_hi",     dat%tas_rec_hi,     missing_value=mv)
    ! call nc_read(filename_out_1,"tas_rec_hi_lin", dat%tas_rec_hi_lin, missing_value=mv)
    call nc_read(filename_out_1,"tas_rec_hi_1880",dat%tas_rec_hi_1880,missing_value=mv)

else 
    ! Calculate stats1 

    ! Perform stats calculations
    call stats_calc_1(dat,L=15.0_wp,mv=mv)
    
    ! === Write output ===

    ! Initialize file and write grid variables
    call write_extremes_init(filename_out_0,dat%lon,dat%lat,dat%sigma,year0=year0,year1=year1)
        
    call nc_write(filename_out_0,"cell_wt", dat%cell_wt, dim1="lon",dim2="lat")
    call nc_write(filename_out_0,"f_land",  dat%f_land,  dim1="lon",dim2="lat")
    call nc_write(filename_out_0,"wt70land",dat%wt70land,dim1="lon",dim2="lat")
    
    ! Write variables
    call nc_write(filename_out_0,"tas",       dat%tas,       dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out_0,"tas_ref",   dat%tas_ref,   dim1="lon",dim2="lat",dim3="month",missing_value=mv)
    call nc_write(filename_out_0,"tas_sm",    dat%tas_sm,    dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out_0,"tas_detrnd",dat%tas_detrnd,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    
    ! Initialize file and write grid variables
    call write_extremes_init(filename_out_1,dat%lon,dat%lat,dat%sigma,year0=year0,year1=year1)

    call nc_write(filename_out_1,"tas_lin_m", dat%tas_lin_m, dim1="lon",dim2="lat",dim3="month",missing_value=mv)
    call nc_write(filename_out_1,"tas_lin_r", dat%tas_lin_r, dim1="lon",dim2="lat",dim3="month",missing_value=mv)
    call nc_write(filename_out_1,"tas_sd",    dat%tas_sd,    dim1="lon",dim2="lat",dim3="month",missing_value=mv)
    call nc_write(filename_out_1,"tas_sigma", dat%tas_sigma, dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out_1,"tas_rec_hi",dat%tas_rec_hi,dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    ! call nc_write(filename_out_1,"tas_rec_hi_lin",dat%tas_rec_hi_lin, &
    !                 dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    call nc_write(filename_out_1,"tas_rec_hi_1880",dat%tas_rec_hi_1880, &
                    dim1="lon",dim2="lat",dim3="month",dim4="year",missing_value=mv)
    
end if 
    
    call stats_calc_2(dat,mv)

    ! Initialize file and write grid variables
    call write_extremes_init(filename_out_2,dat%lon,dat%lat,dat%sigma,year0=year0,year1=year1)
    call nc_write(filename_out_2,"frac_sigma",dat%frac_sigma,dim1="sigma",dim2="month",dim3="year",missing_value=mv)

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
    integer, allocatable :: idx_lin(:) 
    integer, allocatable :: idx_1880(:) 
    integer, allocatable :: idx_numeric(:) 

    integer  :: idx_mm(12,3)
    integer  :: m1, m2, m3 
    real(wp) :: lin_b 

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
    call which(dat%year .ge. 1971 .and. dat%year .le. 2020,idx_lin)
    call which(dat%year .ge. 1880 .and. dat%year .le. 2020,idx_1880)

    ! Initialize all variables with missing values 
    dat%tas_ref    = mv 
    dat%tas_sm     = mv 
    dat%tas_detrnd = mv 
    dat%tas_rec_hi = mv
    dat%tas_rec_hi_lin  = mv
    dat%tas_rec_hi_1880 = mv
    dat%tas_lin_m  = mv
    dat%tas_lin_r  = mv
    dat%tas_sd     = mv 
    dat%tas_sigma  = mv 

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
            
            ! Calculate records 
            call calc_series_records(dat%tas_rec_hi(i,j,m,:),dat%tas(i,j,m,:),low_records=.FALSE.,mv=mv)
            call calc_series_records(dat%tas_rec_hi_lin(i,j,m,:),dat%tas(i,j,m,:),low_records=.FALSE.,mv=mv,idx=idx_lin)
            call calc_series_records(dat%tas_rec_hi_1880(i,j,m,:),dat%tas(i,j,m,:),low_records=.FALSE.,mv=mv,idx=idx_1880)

            ! Calculate linear regression 
            call calc_linear_regression(dat%tas_lin_m(i,j,m),lin_b,dat%tas_lin_r(i,j,m), &
                                            dat%tas(i,j,m,idx_lin),dat%year(idx_lin),mv)

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

subroutine stats_calc_2(dat,mv)
    ! Given a 4D dataset (nx,ny,nm,nyr), perform 
    ! initial steps of statistics calculations. 

    implicit none 

    type(dataset_class), intent(INOUT) :: dat  
    real(wp), intent(IN) :: mv 

    ! Local variables 
    integer :: i, j, m, y 
    integer :: nx, ny, nm, nyr, n 

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
    
    write(*,*) "stats_calc_2..."

    call calc_frac_sigma_series(dat%frac_sigma,dat%tas_sigma,dat%wt70land,dat%sigma,mv)

    write(*,*) "stats_calc_2... done."

    return 

end subroutine stats_calc_2

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

subroutine load_best(dat,filename,year0,year1,mv,load_tas)

    implicit none 

    type(dataset_class), intent(OUT) :: dat 
    character(len=*), intent(IN) :: filename 
    integer,  intent(IN)  :: year0 
    integer,  intent(IN)  :: year1
    real(wp), intent(IN)  :: mv 
    logical,  intent(IN)  :: load_tas 

    ! Local variables  
    real(wp), allocatable :: tas3D(:,:,:)
    integer,  allocatable :: idx(:) 

    write(*,*) "load_best..."

    dat%nx = nc_size(filename_in,"longitude")
    dat%ny = nc_size(filename_in,"latitude")
    dat%nt = nc_size(filename_in,"time")

    allocate(dat%lon(dat%nx))
    allocate(dat%lat(dat%ny)) 
    allocate(tas3D(dat%nx,dat%ny,dat%nt))

    call nc_read(filename,"longitude",dat%lon)
    call nc_read(filename,"latitude", dat%lat)
    
    ! Allocate dataset variables to prepare for populating them.
    call dataset_alloc(dat,year0=year0,year1=year1)

    ! Calculate grid variables (cell_wt,mask)
    call calc_cell_weights(dat%cell_wt,dat%lon,dat%lat)

    call nc_read(filename,"land_mask",dat%f_land)

    call which(dat%lat < -60.0_wp .or. dat%lat > 70.0_wp,idx)
    dat%wt70land = dat%cell_wt*dat%f_land 
    dat%wt70land(:,idx) = 0.0 
    dat%wt70land = dat%wt70land / sum(dat%wt70land)

    ! If desired, read and reshape data from file 
    ! (otherwise, assume this will be done elsewhere)
    if (load_tas) then 
        call nc_read(filename,"temperature",tas3D,missing_value=mv)
        where(abs(tas3D) .gt. 1e10) tas3D = mv          ! For safety
        where(.not. ieee_is_finite(tas3D)) tas3D = mv 
        call reshape_to_4D(dat%tas,tas3D,month0=1,mv=mv)
    end if 

    write(*,*) "load_best... done."
    
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
    allocate(dat%cell_wt(dat%nx,dat%ny))
    allocate(dat%f_land(dat%nx,dat%ny))
    allocate(dat%wt70land(dat%nx,dat%ny))
    allocate(dat%tas(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_ref(dat%nx,dat%ny,dat%nm))
    allocate(dat%tas_sm(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_detrnd(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_sd(dat%nx,dat%ny,dat%nm))
    allocate(dat%tas_sigma(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_rec_hi(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_rec_hi_lin(dat%nx,dat%ny,dat%nm,dat%nyr))
    allocate(dat%tas_rec_hi_1880(dat%nx,dat%ny,dat%nm,dat%nyr))

    allocate(dat%tas_lin_m(dat%nx,dat%ny,dat%nm))
    allocate(dat%tas_lin_r(dat%nx,dat%ny,dat%nm))

    allocate(dat%frac_sigma(dat%nsig,dat%nm,dat%nyr))

    return 

end subroutine dataset_alloc

subroutine dataset_dealloc(dat)

    implicit none 

    type(dataset_class), intent(INOUT) :: dat 

    ! deallocate variables 
    if (allocated(dat%cell_wt))         deallocate(dat%cell_wt)
    if (allocated(dat%f_land))          deallocate(dat%f_land)
    if (allocated(dat%wt70land))        deallocate(dat%wt70land)
    if (allocated(dat%tas))             deallocate(dat%tas)
    if (allocated(dat%tas_ref))         deallocate(dat%tas_ref)
    if (allocated(dat%tas_sm))          deallocate(dat%tas_sm)
    if (allocated(dat%tas_detrnd))      deallocate(dat%tas_detrnd)
    if (allocated(dat%tas_sd))          deallocate(dat%tas_sd)
    if (allocated(dat%tas_sigma))       deallocate(dat%tas_sigma)
    if (allocated(dat%tas_rec_hi))      deallocate(dat%tas_rec_hi)
    if (allocated(dat%tas_rec_hi_lin))  deallocate(dat%tas_rec_hi_lin)
    if (allocated(dat%tas_rec_hi_1880)) deallocate(dat%tas_rec_hi_1880)
    
    if (allocated(dat%frac_sigma))      deallocate(dat%frac_sigma)

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
    real(wp) :: m, b, r 

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

    ! Test linear regression 
    call calc_linear_regression(m,b,r,y,x,mv)

    write(*, "(a, i0)") "sample size = ", n
    write(*, "(a, f10.3)") "Mean      : ", mean
    write(*, "(a, f10.3)") "Stddev    : ", stdev  
    write(*,*) "====="
    write(*, "(a, f10.3)") "Slope     : ", m 
    write(*, "(a, f10.3)") "Intercept : ", b 
    write(*, "(a, f10.3)") "r-value   : ", r 

    call nc_create(filename)
    call nc_write_dim(filename,"x",x)
    call nc_write(filename,"y",y,dim1="x")
    call nc_write(filename,"ysm",ysm,dim1="x")

    return 

end subroutine test_timeseries


subroutine calc_frac_sigma_series(frac_sigma,tas_sigma,wts,sigma,mv)
    ! Get spatial fraction of weighted area that is >= sigma-thresholds
    ! defined in the vector `sigma`. This variable is calculated
    ! for each month and each year. 

    implicit none 

    real(wp), intent(OUT) :: frac_sigma(:,:,:)      ! [nsig,nm,nyr]
    real(wp), intent(IN)  :: tas_sigma(:,:,:,:)     ! [nx,ny,nm,nyr]
    real(wp), intent(IN)  :: wts(:,:)               ! [nx,ny]
    real(wp), intent(IN)  :: sigma(:)
    real(wp), intent(IN)  :: mv 

    ! Local variables 
    integer :: y, m, q, nyr, nm, nsig, nx, ny, n
    integer,  allocatable :: idx(:)  
    real(wp), allocatable :: tas_sigma_now(:) 
    real(wp), allocatable :: wts_now(:) 
    real(wp) :: wt_tot 

    nx   = size(tas_sigma,1)
    ny   = size(tas_sigma,2)

    nsig = size(frac_sigma,1) 
    nm   = size(frac_sigma,2)
    nyr  = size(frac_sigma,3)
    
    allocate(tas_sigma_now(nx*ny)) 
    allocate(wts_now(nx*ny)) 

    ! Initially set missing values everywhere
    frac_sigma = mv 

    write(*,"(a)",advance="no") "calc_frac_sigma_series..."

    do y = 1, nyr 

        do m = 1, nm 

            ! Get vector variable and weights 
            tas_sigma_now = reshape(tas_sigma(:,:,m,y),shape=[nx*ny])
            wts_now       = reshape(wts,shape=[nx*ny])

            where(tas_sigma_now .eq. mv) wts_now = 0.0_wp 

            wt_tot = sum(wts_now) 

            if (wt_tot .gt. 0.0_wp) then 

                do q = 1, nsig 

                    call which( tas_sigma_now .ne. mv       .and. &
                                tas_sigma_now .ge. sigma(q) .and. &
                                wts_now .gt. 0.0_wp, idx, n )

                    if (n .gt. 0) then 
                        frac_sigma(q,m,y) = sum(wts_now(idx)) / wt_tot 
                    else 
                        frac_sigma(q,m,y) = 0.0_wp 
                    end if 

                end do 

            end if 

        end do 
    end do 

    write(*,*) "done."

    return 

end subroutine calc_frac_sigma_series

subroutine calc_series_records(yrec,y,low_records,mv,idx)

    implicit none 

    real(wp), intent(OUT) :: yrec(:) 
    real(wp), intent(IN)  :: y(:) 
    logical,  intent(IN)  :: low_records 
    real(wp), intent(IN)  :: mv
    integer,  intent(IN), optional :: idx(:)

    ! Local variables 
    integer :: i, k, nt, n 
    integer,  allocatable :: idx_now(:) 
    real(wp), allocatable :: y_now(:) 
    real(wp) :: ylim 

    nt = size(y) 
    allocate(y_now(nt)) 
    y_now = y 

    if (present(idx)) then 
        ! Remove values that should not be treated 
        do k = 1, nt 
            if (.not. any(idx .eq. k)) then 
                y_now(k) = mv 
            end if 
        end do 
    end if 

    ! Set yrec initially to zero everywhere 
    yrec = 0

    ! Get indices of available values 
    call which(y_now .ne. mv, idx_now, n)

    if (n .gt. 0) then 
        ! Data available, proceed with record counting 

        ! Initial value is always a record
        k = idx_now(1) 
        yrec(k) = 1
        ylim    = y_now(k)

        if (low_records) then 
            ! Count low records 

            do i = 2, n 
                k = idx_now(i) 
                if (y_now(k) .lt. ylim) then 
                    yrec(k) = 1 
                    ylim    = y_now(k) 
                end if 
            end do 

        else 
            ! Count high records 

            do i = 2, n
                k = idx_now(i)  
                if (y_now(k) .gt. ylim) then 
                    yrec(k) = 1 
                    ylim    = y_now(k) 
                end if 
            end do 

        end if

    end if 

    return 

end subroutine calc_series_records

subroutine calc_linear_regression(m,b,r,y,x,mv)
    ! Calculate linear regression of the form
    ! y = mx + b, with correlation coefficient r 
    ! using least squares regression algorithm: 
    ! 1. compute sum of x
    ! 2. compute sum of x**2
    ! 3. compute sum of x * y
    ! 4. compute sum of y
    ! 5. compute sum of y**2
            
    implicit none 

    real(wp), intent(OUT) :: m 
    real(wp), intent(OUT) :: b 
    real(wp), intent(OUT) :: r 
    real(wp), intent(IN)  :: x(:) 
    real(wp), intent(IN)  :: y(:) 
    real(wp), intent(IN)  :: mv
    
    ! Local variables 
    integer :: i, k, n 
    integer, allocatable :: idx(:) 
    real(wp) :: npts
    real(wp) :: sumx, sumx2, sumxy, sumy, sumy2 

    call which(x .ne. mv .and. y .ne. mv,idx,n)

    if (n .ge. 2) then 
        ! Points are available, perform linear regression 

        sumx  = 0.0_wp 
        sumx2 = 0.0_wp 
        sumxy = 0.0_wp 
        sumy  = 0.0_wp 
        sumy2 = 0.0_wp 

        do i = 1, n
            k = idx(i)                                              
            sumx  = sumx  + x(k)
            sumx2 = sumx2 + x(k) * x(k)
            sumxy = sumxy + x(k) * y(k)
            sumy  = sumy  + y(k)                                                              
            sumy2 = sumy2 + y(k) * y(k)
        end do

        ! Specify number of points in real format for easy division
        npts = real(n,wp)

        ! Calculate the slope, the y-intercept and the correlation coefficient r
        ! 1. compute slope
        ! 2. compute y-intercept
        ! 3. compute correlation coefficient
        m = (npts * sumxy  -  sumx * sumy)  / (npts * sumx2 - sumx**2)
        b = (sumy * sumx2  -  sumx * sumxy) / (npts * sumx2 - sumx**2)
        r = (sumxy - sumx * sumy / npts) /  &
                     sqrt((sumx2 - sumx**2/npts) * (sumy2 - sumy**2/npts))

    else 
        ! Assign missing values to output

        m = mv 
        b = mv 
        r = mv

    end if

    return 

end subroutine calc_linear_regression

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

subroutine calc_cell_weights(cell_wt,lon,lat)

    implicit none 

    real(wp), intent(OUT) :: cell_wt(:,:) 
    real(wp), intent(IN)  :: lon(:) 
    real(wp), intent(IN)  :: lat(:) 
    
    ! Local variables 
    integer :: i, j, nx, ny 

    nx = size(cell_wt,1)
    ny = size(cell_wt,2) 

    do j = 1, ny 
        cell_wt(:,j) = cos(lat(j)*degrees_to_radians)
    end do 

    return 

end subroutine calc_cell_weights

subroutine write_extremes_init(filename,lon,lat,sigma,year0,year1)

    implicit none 

    character(len=*), intent(IN) :: filename 
    real(wp), intent(IN) :: lon(:)
    real(wp), intent(IN) :: lat(:)
    real(wp), intent(IN) :: sigma(:)
    integer,  intent(IN) :: year0 
    integer,  intent(IN) :: year1 
    
    call nc_create(filename)
    call nc_write_dim(filename,"lon",lon,units="degrees_east", &
                        long_name="Longitude")
    call nc_write_dim(filename,"lat",lat,units="degrees_north", &
                        long_name="Latitude")
    
    call nc_write_dim(filename,"sigma",sigma,units="", &
                        long_name="sigma-level")

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

