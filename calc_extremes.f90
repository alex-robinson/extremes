!! TO COMPILE : 
!! gfortran -o calc_extremes.x -I/opt/local/include libs/ncio.f90 calc_extremes.f90 -L/opt/local/lib -lnetcdff -lnetcdf

program calc_extremes 

    use ncio 

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

    integer, parameter :: ntest = 1000
    real(wp) :: var(ntest)
    real(wp) :: stdev 

    integer :: i 

    filename_in  = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1.nc"
    filename_out = "data/BerkeleyEarth/2020-08_BEST/Land_and_Ocean_LatLong1_stats.nc"


    ! Test standard deviation 
    call random_seed
    do i = 1, ntest 
        call random_number(var(i))
    end do 
    call calc_stdev(stdev,var,mv)
    stop 

    ! Load data 
    call load_best(dat,filename_in,year0=1850,year1=2020,mv=mv)

    ! Write output
    call write_extremes_init(filename_out,dat%lon,dat%lat,year0=1850,year1=2020)
    call nc_write(filename_out,"tas",dat%tas,dim1="lon",dim2="lat",dim3="month",dim4="year")

contains 

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


subroutine calc_stdev(stdev,var,mv)

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

        write(*, "(a, i0)") "sample size = ", n
        write(*, "(a, f17.15)") "Mean :   ", mean
        write(*, "(a, f17.15)") "Stddev : ", stdev  

    else 
        ! Set standard deviation equal to missing value 

        stdev = mv 

    end if 

    return 

end subroutine calc_stdev


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

