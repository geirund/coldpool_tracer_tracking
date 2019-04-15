MODULE irt_parameters
implicit none

INTEGER               :: dsize_x 
INTEGER               :: dsize_y 
REAL                  :: res  ! resolution in m
REAL                  :: dt  ! time step in sec
CHARACTER(LEN=100)    :: odir

LOGICAL, PARAMETER    :: lperiodic = .TRUE.

INTEGER, PARAMETER    :: n_fields =1 ! nfield   ! number of additional averaging fields


! bins of coarse velocity field
INTEGER, PARAMETER    :: time_steps = 267    ! total number of timesteps
INTEGER, PARAMETER    :: nt_bins = 10
INTEGER, PARAMETER    :: nx_bins = 2
INTEGER, PARAMETER    :: ny_bins = 2

REAL                  :: cutoff             ! for rain intensity
REAL                  :: minimum_size        ! events smaller than that will be sorted out
! for cp_tracking: We are interested in deep convective cps. Thus precipi event
! causing cp should be larger than RB scale (2km) ~ 10x10 gps 
REAL, PARAMETER       :: threshold_ratio=1      ! Choose value between 0.0 and 1.0

REAL, PARAMETER       :: max_velocity = 10.   ! adjust acordingly
                                              ! velocities>max_velocity will be ignored to remove outliers
! define a minimal number of cells required for a coarse grained coordinate to be evaluated 
! if there are less, missing value will be assigned to that coarse cell
REAL, PARAMETER       :: edge_fraction = .01  ! the estimated fraction of the domain located on the edge of cold pools
INTEGER, PARAMETER    :: min_cells = 20 !256 !4 ! for clouds !256 for rain
INTEGER, PARAMETER    :: tracer_steps = 20
INTEGER, PARAMETER    :: max_tfields = 20     ! maximum number of distinct tracer fields to keep track of

INTEGER, PARAMETER    :: max_no_of_cells=10000  ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_no_of_tracks=10000    ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_length_of_track=267  ! buffer size, increase if necessary

REAL, PARAMETER       :: miss=-9999.           ! value<miss ==> missing_value
INTEGER, PARAMETER    :: max_tracers=500000     ! max number of tracers to track


END MODULE irt_parameters
