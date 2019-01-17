MODULE cp_parameters
implicit none

INTEGER               :: dsize_x, dsize_y  ! domain size
INTEGER               :: max_age           ! max # of timesteps for tracking one CP
INTEGER               :: nTrSet            ! # of tracer set per timestep 
INTEGER               :: max_tracer_CP     ! produkt of nTrSet and ntset 
INTEGER               :: ntset             ! number of timesteps at which tracer are set after precip start
REAL                  :: res               ! resolution in m
REAL                  :: dt                ! time step in sec
CHARACTER(LEN=200)    :: odir
logical               :: rad               ! If True tracer follow radial vel
character(LEN=3)      :: lformat
END MODULE cp_parameters
