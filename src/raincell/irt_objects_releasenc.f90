! identifies objects and establishes links
! needs input data as SRV file
! Compile: ifort -no-wrap-margin -o irt_objects_release.x irt_objects_release.f90

PROGRAM irt_objects
USE netcdf

USE irt_parameters, ONLY: dsize_x, dsize_y, lperiodic, &
    n_fields, time_steps, nt_bins, nx_bins, ny_bins, &
    max_no_of_cells, miss, cutoff ,minimum_size, dt, res, odir 

IMPLICIT NONE

INTEGER              :: domsize_x, domsize_y
INTEGER              :: nt,nx,ny,nz,imain
INTEGER              :: ii, ij,ix,iy,idx,idy
REAL, ALLOCATABLE    :: input_field(:,:,:)
REAL, ALLOCATABLE    :: overlay_field(:,:)
INTEGER, ALLOCATABLE :: event_number(:,:,:)
LOGICAL, ALLOCATABLE :: occupied(:,:)
INTEGER              :: i,j,a,b,c,t,it,fileid
INTEGER              :: counter_actual, counter_previous, counter_mixed
INTEGER              :: counter_total_actual, counter_total_previous
LOGICAL              :: delete_cell
REAL                 :: totarea_sum

! for cell statistics
REAL                 :: totarea(max_no_of_cells,2)
REAL                 :: field_mean(max_no_of_cells,2,n_fields+1)
REAL                 :: field_min(max_no_of_cells,2,n_fields+1)
REAL                 :: field_max(max_no_of_cells,2,n_fields+1)
REAL                 :: center_of_mass_x(max_no_of_cells,2),center_of_mass_y(max_no_of_cells,2)
INTEGER              :: first_point_x(max_no_of_cells,2),first_point_y(max_no_of_cells,2)
INTEGER              :: xfirst(max_no_of_cells,2), xlast(max_no_of_cells,2)
INTEGER              :: yfirst(max_no_of_cells,2), ylast(max_no_of_cells,2)

! variables for the overlay field
INTEGER              :: xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed

! for link statistics
INTEGER              :: largest_forward_link(max_no_of_cells)
REAL                 :: largest_forward_link_size(max_no_of_cells)
INTEGER              :: second_largest_forward_link(max_no_of_cells)
REAL                 :: second_largest_forward_link_size(max_no_of_cells)
INTEGER              :: largest_backward_link(max_no_of_cells)
REAL                 :: largest_backward_link_size(max_no_of_cells)
INTEGER              :: second_largest_backward_link(max_no_of_cells)
REAL                 :: second_largest_backward_link_size(max_no_of_cells)

! diagnosed cell velocities
REAL                 :: velocity_x(max_no_of_cells)
REAL                 :: velocity_y(max_no_of_cells)

! dummy variables for velocity diagnosis
REAL                 :: area_weight,dx,dy

INTEGER              :: srv_header_input(8)
INTEGER              :: srv_header_coarse(8)

! file names
CHARACTER (len=290)   :: input_filename(n_fields+1)
CHARACTER (len=190)   :: output_filename
CHARACTER (len=190)   :: mask_filename
CHARACTER (len=190)   :: coarsevel_filename

CHARACTER (len=1)    :: iteration_str
INTEGER              :: headlen

! time handling
INTEGER              :: n_actual,n_previous,n_previous_it
LOGICAL              :: previous_exists
INTEGER              :: previousdate, previoustime, previoustimestep

! which iteration?
INTEGER              :: iteration = 1   ! 1: first iteration
                                        ! 2: second iteration

! for advection velocity field
REAL                 :: coarse_vel_x(nx_bins,ny_bins,nt_bins)
REAL                 :: coarse_vel_y(nx_bins,ny_bins,nt_bins)
REAL                 :: coarse_dummy(nx_bins,ny_bins)
REAL                 :: vx,vy
INTEGER              :: ierr

! GET NAMELIST INPUT
namelist /INPUTgeneral/ dsize_x, dsize_y, dt, res
namelist /INPUTrain/  cutoff,minimum_size
namelist /INPUTIO/ odir
open(10,file='job/namelist.dat')
read(10,nml=INPUTgeneral)
read(10,nml=INPUTIO)
read(10,nml=INPUTrain)

write(*,*) odir
write(*,*) dsize_x, dsize_y, dt, res
write(*,*) cutoff,minimum_size

DO i=0,n_fields
  WRITE(input_filename(i+1),"(A33,I2.2,A4)") "input/raincell/irt_objects_input_",i,".srv"
ENDDO
WRITE(*,*) input_filename(1)
output_filename    =trim(odir)// "output/raincell/irt_objects_output.txt"
mask_filename      = trim(odir)// "output/raincell/irt_objects_mask.srv"
coarsevel_filename = trim(odir)// "output/raincell/irt_advection_field.srv"
write(*,*) "hier 1"
! extracting the argument from the command line
CALL getarg(1,iteration_str)

! checking if this is the first or a higher iteration
IF (iteration_str .EQ. "1") THEN
   iteration = 1
ELSEIF (iteration_str .EQ. "2") THEN
   iteration = 2
ELSE
  WRITE (*,*) "ERROR: 1 for first iteration, 2 for subsequent iterations"
  STOP
ENDIF
write(*,*) "hier 2"

! open the input filenames
DO fileid=1,n_fields+1
  OPEN(fileid,FILE=trim(odir)//trim(input_filename(fileid)),FORM='unformatted', ACTION='read')
ENDDO

! open the output filenames
OPEN(20,FILE=trim(output_filename),FORM='formatted', ACTION='write')
OPEN(30,FILE=trim(mask_filename),FORM='unformatted', ACTION='write')
write(*,*) "hier 3"
! if not first iteration, read in the coarse velocity field
IF (iteration == 2) THEN
  OPEN (50,FILE=trim(coarsevel_filename),FORM='unformatted', STATUS='old', ACTION='read')
  DO i=1,nt_bins
    READ (50, END=200) srv_header_coarse
    READ (50) coarse_vel_x(:,:,i)
    READ (50) srv_header_coarse
    READ (50) coarse_vel_y(:,:,i)
    READ (50) srv_header_coarse
    READ (50) coarse_dummy(:,:)
  ENDDO
  CLOSE (50)
ENDIF

n_actual     = 1
n_previous   = 2
write(*,*) "hier 4"
previous_exists = .FALSE.  ! because it's the first time step
counter_actual = 0
counter_total_actual = 0
counter_total_previous = 0
previoustimestep=-1

! If periodic boundary conditions are switched off, a 1-gridbix thick 
! frame of missing values will be laid around the field, and the domainsize
! has to be increased by 2 in both dimensions
IF (lperiodic) THEN
  domsize_x = dsize_x
  domsize_y = dsize_y
ELSE
  domsize_x = dsize_x+2
  domsize_y = dsize_y+2
ENDIF

! allocating fields for the different quantities
ALLOCATE(input_field(domsize_x,domsize_y,n_fields+1))
ALLOCATE(overlay_field(domsize_x,domsize_y))
ALLOCATE(event_number(domsize_x,domsize_y,4))
ALLOCATE(occupied(domsize_x,domsize_y))
write(*,*) "hier 5"
! beginning of main loop
 CALL get_dim(nt,nx,ny,nz,trim(odir) //'input/raincell/irt_objects_input_00.nc') 
write(*,*) "hier 6", 'nt',nt,'nx',nx,'ny',ny,'nz',nz
DO imain=1,nt,1 

previousdate=srv_header_input(3)
previoustime=srv_header_input(4)
previoustimestep=previoustimestep+1
! reading the standard input fields
DO fileid=1,n_fields+1
  IF (lperiodic) THEN
        CALL read_nc_2D (input_field(:,:,1),'TOT_PREC',trim(odir) //'input/raincell/irt_objects_input_00.nc',previoustimestep+1)
        WRITE(*,*) imain, input_field(40,40,1),nt,"testnc"
        CALL WRITEMAP(input_field(:,:,1),imain,"testnc")
  ELSE
    input_field(:,:,fileid) = miss-1.
    READ (fileid) input_field(2:domsize_x-1,2:domsize_y-1,fileid)
  ENDIF
ENDDO

occupied(:,:)=.FALSE.
event_number(:,:,n_actual)=0
counter_previous = MAX(1,counter_actual)
counter_actual=1

!WRITE(*,*) "domsize: ",domsize_x,domsize_y
! identification of patches
DO iy=1, domsize_y
  DO ix=1, domsize_x
    IF (input_field(ix,iy,1) .GE. cutoff .AND. .NOT. occupied(ix,iy)) THEN
      ii = ix
      ij = iy
      totarea(counter_actual,n_actual)=0
      field_mean(counter_actual,n_actual,:)=0
      field_min(counter_actual,n_actual,:)=1E+9
      field_max(counter_actual,n_actual,:)=0
      center_of_mass_x(counter_actual,n_actual)=0
      center_of_mass_y(counter_actual,n_actual)=0
      xfirst(counter_actual,n_actual)=ii
      xlast(counter_actual,n_actual)=ii
      yfirst(counter_actual,n_actual)=ij
      ylast(counter_actual,n_actual)=ij
      delete_cell = .FALSE.
      !WRITE(*,*) "test",ii,ij
      ! calling the recursive method for rain area identification
      CALL area(miss, cutoff, ii, ij, domsize_x,domsize_y,n_fields, input_field, &
         occupied, counter_actual,event_number(:,:,n_actual),totarea(counter_actual,n_actual), &
         field_mean(counter_actual,n_actual,:),field_min(counter_actual,n_actual,:),field_max(counter_actual,n_actual,:), &
         center_of_mass_x(counter_actual,n_actual),center_of_mass_y(counter_actual,n_actual), &
         delete_cell,xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
         yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))

      IF (delete_cell) THEN
        ! delete this cell by overwriting it with -1 in event_number
        CALL set_event_number_to_value(domsize_x,domsize_y,-1,event_number(:,:,n_actual),counter_actual, &
                                       xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
				       yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))
      ELSEIF (counter_actual .GT. max_no_of_cells) THEN
        WRITE(*,*) "ERROR: number of cells >",max_no_of_cells
        STOP
      ELSEIF (totarea(counter_actual,n_actual) .GT. minimum_size) THEN
        center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)/ &
                                                  field_mean(counter_actual,n_actual,1)
        center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)/ &
	                                          field_mean(counter_actual,n_actual,1)
	field_mean(counter_actual,n_actual,:)=field_mean(counter_actual,n_actual,:)/&
	                                      totarea(counter_actual,n_actual)

	! take care for periodic boundary conditions
	IF (center_of_mass_x(counter_actual,n_actual) .GE. domsize_x+1) THEN
	  center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)-domsize_x
	ENDIF
	IF (center_of_mass_x(counter_actual,n_actual) .LE. 1) THEN
	  center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)+domsize_x
	ENDIF
	IF (center_of_mass_y(counter_actual,n_actual) .GE. domsize_y+1) THEN
	  center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)-domsize_y
	ENDIF
	IF (center_of_mass_y(counter_actual,n_actual) .LE. 1) THEN
	  center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)+domsize_y
	ENDIF

	first_point_x(counter_actual,n_actual)=ix
	first_point_y(counter_actual,n_actual)=iy
        counter_actual=counter_actual+1  ! the event counter is increased by one
      ELSE   ! if cell is too small, delete it

	CALL set_event_number_to_value(domsize_x,domsize_y,-1,event_number(:,:,n_actual),counter_actual, &
	                               xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
				       yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))
      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (lperiodic) THEN
  CALL write_srv(REAL(event_number(:,:,n_previous)),previousdate,previoustime,30,previoustimestep)
ELSE
  CALL write_srv(REAL(event_number(2:domsize_x-1,2:domsize_y-1,n_previous)),&
                 previousdate,previoustime,30,previoustimestep)
ENDIF

counter_total_previous = counter_total_actual
counter_total_actual = counter_total_actual+counter_previous-1

! do we have to read in the next step also?
IF (.NOT. previous_exists) THEN
  previous_exists = .TRUE.
  ! flip time indices and advance 5 minutes
  n_actual   = 3-n_actual
  n_previous = 3-n_previous
  !WRITE (*,*) "JUMP"
  CYCLE
ENDIF

! perform second iteration on "previous" field
IF (iteration .EQ. 2) THEN
  event_number(:,:,4)=0
  DO iy=1,domsize_y
    DO ix=1,domsize_x
      i = event_number(ix,iy,n_previous)
      IF (i .LT. 1) CYCLE
      ii  = NINT(center_of_mass_x(i,n_previous))
      ij  = NINT(center_of_mass_y(i,n_previous))
      CALL time_interpolation(domsize_x,domsize_y,time_steps,nx_bins,ny_bins,nt_bins,lperiodic, &
                              coarse_vel_x,coarse_vel_y,ii,ij,MAX(1,previoustimestep),vx,vy)
      idx = NINT(vx)
      idy = NINT(vy)
      IF (idx .LE. miss+1 .OR. idy .LE. miss+1) CYCLE
      IF (ix+idx .LT. 1 .OR. ix+idx .GT. domsize_x .OR. iy+idy .LT. 1 .OR. iy+idy .GT. domsize_y) CYCLE
      event_number(ix+idx,iy+idy,4) = i
      first_point_x(i,n_previous) = ix+idx
      first_point_y(i,n_previous) = iy+idy
    ENDDO
  ENDDO
  n_previous_it = 4
  !CALL write_srv(domsize_x,domsize_y,REAL(event_number(:,:,4)),previousdate,previoustime+30,30)
ELSE
  n_previous_it = n_previous
ENDIF

! calculate events for the mixed field
overlay_field = 0
DO iy=1,domsize_y
  DO ix=1,domsize_x
    IF (input_field(ix,iy,1) .LE. miss) THEN
      overlay_field(ix,iy) = miss-1
    ELSE IF(event_number(ix,iy,n_actual) .NE. 0 .OR. event_number(ix,iy,n_previous_it) .NE. 0) THEN
      overlay_field(ix,iy) = 1.
    ENDIF
  ENDDO
ENDDO

occupied(:,:)=.FALSE.

counter_mixed=1

!identification of patches for the overlaid field
DO iy=1, domsize_y
  DO ix=1, domsize_x
    IF (overlay_field(ix,iy) .GE. 0.5 .AND. .NOT. occupied(ix,iy)) THEN
      ii = ix
      ij = iy
      xfirst_mixed=ii
      xlast_mixed=ii
      yfirst_mixed=ij
      ylast_mixed=ij
      delete_cell = .FALSE.
      CALL overlay_area(miss, ii, ij, domsize_x,domsize_y, overlay_field, &
         occupied,counter_mixed,event_number(:,:,3),  &
	 delete_cell,xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed)
      IF (delete_cell) THEN
        CALL set_event_number_to_value(domsize_x,domsize_y,0,event_number(:,:,3),counter_mixed, &
	                               xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed)
       ELSEIF (counter_mixed .LE. max_no_of_cells) THEN
        counter_mixed=counter_mixed+1
      ELSE
        WRITE(*,*) "ERROR: number of cells >",max_no_of_cells
        STOP
      ENDIF
    ENDIF
  ENDDO
ENDDO

!WRITE (*,*) 'timestep:',previoustimestep,', time:',srv_header_input(4)
! identify foreward links and velocity
largest_forward_link=0
largest_forward_link_size=0
second_largest_forward_link=0
second_largest_forward_link_size=0
DO i=1,counter_previous-1
  a = event_number(first_point_x(i,n_previous),first_point_y(i,n_previous),3)
  IF (a==0) THEN 
    largest_forward_link(i) = -1    ! track interrupted by missing values
    CYCLE
  ENDIF
  velocity_x(i) = 0.
  velocity_y(i) = 0.
  area_weight = 0.
  DO j=1,counter_actual-1
    b = event_number(first_point_x(j,n_actual),first_point_y(j,n_actual),3)
    IF (b==0) CYCLE
    IF (a==b) THEN
      dx = center_of_mass_x(j,n_actual)-center_of_mass_x(i,n_previous)
      dy = center_of_mass_y(j,n_actual)-center_of_mass_y(i,n_previous)
      ! take care for cyclic boundary conditions:
      IF (dx .GT. domsize_x/2)  dx = dx-domsize_x
      IF (dx .LT. -domsize_x/2) dx = dx+domsize_x
      IF (dy .GT. domsize_y/2)  dy = dy-domsize_y
      IF (dy .LT. -domsize_y/2) dy = dy+domsize_y
      velocity_x(i) = velocity_x(i)+dx*totarea(j,n_actual)
      velocity_y(i) = velocity_y(i)+dy*totarea(j,n_actual)
      area_weight = area_weight+totarea(j,n_actual)
      IF (totarea(j,n_actual) .GT. largest_forward_link_size(i)) THEN
        second_largest_forward_link(i) = largest_forward_link(i)
	second_largest_forward_link_size(i) = largest_forward_link_size(i)
	largest_forward_link(i) = counter_total_actual+j
	largest_forward_link_size(i) = totarea(j,n_actual)
      ELSEIF (totarea(j,n_actual) .GT. second_largest_forward_link_size(i)) THEN
	second_largest_forward_link(i) = counter_total_actual+j
	second_largest_forward_link_size(i) = totarea(j,n_actual)
      ENDIF
    ENDIF
  ENDDO
  ! determine cell velocity if possible, otherwise assign missing value
  IF (area_weight .GE. minimum_size) THEN
    velocity_x(i) = velocity_x(i)/area_weight
    velocity_y(i) = velocity_y(i)/area_weight
  ELSE
    velocity_x(i) = miss
    velocity_y(i) = miss
  ENDIF
ENDDO

! write sizes
IF (counter_previous==1) THEN
  WRITE(20,*) previoustimestep,0,counter_total_previous+i,totarea(1,n_previous), &
	      field_mean(1,n_previous,:),field_min(1,n_previous,:),field_max(1,n_previous,:), &
	      xfirst(1,n_previous),xlast(1,n_previous),yfirst(1,n_previous),ylast(1,n_previous), &
	      center_of_mass_x(1,n_previous),center_of_mass_y(1,n_previous), &
	      velocity_x(1),velocity_y(1), &
	      largest_forward_link(1),largest_forward_link_size(1), &
	      second_largest_forward_link(1),second_largest_forward_link_size(1), &
	      largest_backward_link(1),largest_backward_link_size(1), &
	      second_largest_backward_link(1),second_largest_backward_link_size(1)
ENDIF
DO i=1,counter_previous-1
  WRITE(20,*) previoustimestep,i,counter_total_previous+i,totarea(i,n_previous), &
	      field_mean(i,n_previous,:),field_min(i,n_previous,:),field_max(i,n_previous,:), &
	      xfirst(i,n_previous),xlast(i,n_previous),yfirst(i,n_previous),ylast(i,n_previous), &
	      center_of_mass_x(i,n_previous),center_of_mass_y(i,n_previous), &
	      velocity_x(i),velocity_y(i), &
	      largest_forward_link(i),largest_forward_link_size(i), &
	      second_largest_forward_link(i),second_largest_forward_link_size(i), &
	      largest_backward_link(i),largest_backward_link_size(i), &
	      second_largest_backward_link(i),second_largest_backward_link_size(i)
ENDDO

! identify backward links
largest_backward_link=0
largest_backward_link_size=0
second_largest_backward_link=0
second_largest_backward_link_size=0
DO i=1,counter_actual-1
  a = event_number(first_point_x(i,n_actual),first_point_y(i,n_actual),3)
  IF (a==0) THEN
    largest_backward_link(i) = -1    ! track interrupted by missing values
    CYCLE
  ENDIF
  DO j=1,counter_previous-1
    b = event_number(first_point_x(j,n_previous),first_point_y(j,n_previous),3)
    IF (b==0) CYCLE
    IF (a==b) THEN
      IF (totarea(j,n_previous) .GT. largest_backward_link_size(i)) THEN
        second_largest_backward_link(i) = largest_backward_link(i)
	second_largest_backward_link_size(i) = largest_backward_link_size(i)
	largest_backward_link(i) = counter_total_previous+j
	largest_backward_link_size(i) = totarea(j,n_previous)
      ELSEIF (totarea(j,n_previous) .GT. second_largest_backward_link_size(i)) THEN
	second_largest_backward_link(i) = counter_total_previous+j
	second_largest_backward_link_size(i) = totarea(j,n_previous)
      ENDIF
    ENDIF
  ENDDO
ENDDO

! flip time indices
n_actual   = 3-n_actual
n_previous = 3-n_previous

! end main loop
ENDDO

200 CONTINUE

! close input/output files
CLOSE(20)

DO fileid=1,n_fields+1
  CLOSE(fileid)
ENDDO


CONTAINS

! -----------------------------------------------------------------------
!subroutines for reading and writing netcdf
! -----------------------------------------------------------------------

!! -----------------------------------------------------------------------
!! input 2D nc data
!! -----------------------------------------------------------------------
  SUBROUTINE read_nc_2D (poutput,varname, filename,ctime)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  real, dimension(:,:), intent(inout) :: poutput
  integer :: ncId, rhVarId, nx, ny, nz, nt, ctime
  integer, dimension(nf90_max_var_dims) :: dimIDs
!  real, allocatable, dimension(:,:,:) ::  zvar
    CALL check(nf90_open(filename, nf90_NoWrite, ncid))
    CALL check(nf90_inq_varid(ncid,varname, rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
! for reading the divergence file
   ! CALL check(nf90_inquire_dimension(ncid, dimIDs(4), len = nt))
   ! CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
   ! CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
   ! CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nz))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
!   CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nz))
write(*,*) ctime
  !  CALL check(nf90_get_var(ncid, rhVarId, poutput, start =(/1,1,1,ctime/),&
  !                                                  count= (/nx, ny,nz, 1/)))
    CALL check(nf90_get_var(ncid, rhVarId, poutput, start =(/1,1,ctime/),&
                                                    count= (/ny,nx,1/)))

    CALL check(nf90_close(ncid))
write(*,*) MAXVAL(poutput)
!
  END SUBROUTINE read_nc_2D
!12.Juni 2019>
! zum test
  SUBROUTINE WRITEMAP(map,t,FILE_NAME2)

  real, intent(in) :: map(350,350)
  integer, intent(in) :: t
  integer :: ncid, x_dimid, y_dimid, rec_dimid, dimids(3),count(3),start(3), varid
  character (len = *), intent(in):: FILE_NAME2
  character (len = 200):: FILE_NAME

   FILE_NAME =trim(FILE_NAME2)//trim(str(t))//'.nc' !'CPout'//trim(str(t)//'.nc'
    call check( nf90_create(trim(FILE_NAME), NF90_CLOBBER, ncid) )
    call check( nf90_def_dim(ncid, "x", dsize_x, x_dimid) )
    call check( nf90_def_dim(ncid, "y", dsize_y, y_dimid) )
    call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, rec_dimid) )
    dimids =  (/ y_dimid, x_dimid,rec_dimid /)
    call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
    call check( nf90_enddef(ncid) )
    count = (/ dsize_x, dsize_y, 1 /)
    start = (/ 1, 1, 1 /)
    call check( nf90_put_var(ncid, varid, map, start=start, count=count) )
    call check( nf90_close(ncid) )
  END SUBROUTINE WRITEMAP

! -----------------------------------------------------------------------
! function to shorten characters for output names
! -----------------------------------------------------------------------
character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    str = adjustl(str)
    write (str, *) k
    str = adjustl(str)

end function str

! end
!---------------------------------------------------
! CHECK NETCD DATA
!--------------------------------------------------
  SUBROUTINE check(status)

    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  END SUBROUTINE check

  SUBROUTINE get_dim(nt,nx,ny,nz,filename)

  IMPLICIT NONE
  character(*), intent(in) :: filename
  integer, intent(out)     :: nt, nx, ny,nz
  integer                  :: ncId, rhVarId
  integer, dimension(nf90_max_var_dims) :: dimIDs
    CALL check(nf90_open(filename, nf90_NoWrite, ncid))
    CALL check(nf90_inq_varid(ncid,"TOT_PREC", rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nt))
    nz =0
    !CALL check(nf90_inquire_dimension(ncid, dimIDs(4), len = nt))
    !CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    !CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    !CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nz))
    !CALL check(nf90_close(ncid))
END SUBROUTINE get_dim

END PROGRAM irt_objects


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE write_srv(field,date,time,file_id,timestep)
USE irt_parameters, ONLY: dsize_x, dsize_y
  IMPLICIT NONE
  REAL, INTENT(IN)      :: field(dsize_x,dsize_y)
  INTEGER, INTENT(IN)   :: date,time,timestep
  INTEGER, INTENT(IN)   :: file_id
  INTEGER               :: srv_header(8)

  srv_header(1) = 1	      ! Code
  srv_header(2) = 1	      ! Level
  srv_header(3) = date        ! Datum
  srv_header(4) = time        ! Zeitinkrement
  srv_header(5) = dsize_x
  srv_header(6) = dsize_y
  srv_header(7) = timestep !timestep
  srv_header(8) = 0

  WRITE (file_id) srv_header
  WRITE (file_id) field
  
  RETURN

END SUBROUTINE write_srv


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE set_event_number_to_value(domsize_x,domsize_y,value,event_number,nevent, &
                                     xfirst,xlast,yfirst,ylast)

  IMPLICIT NONE
  INTEGER, INTENT(IN)	    :: domsize_x,domsize_y
  INTEGER, INTENT(IN)       :: value
  INTEGER, INTENT(INOUT)    :: event_number(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: nevent
  INTEGER, INTENT(IN)       :: xfirst,xlast,yfirst,ylast
  INTEGER                   :: ix,iy,ix_mod,iy_mod
  DO iy=yfirst,ylast
    DO ix=xfirst,xlast
      ix_mod = MOD(ix-1+domsize_x,domsize_x)+1
      iy_mod = MOD(iy-1+domsize_y,domsize_y)+1
      IF (event_number(ix_mod,iy_mod)==nevent) event_number(ix_mod,iy_mod)=value
    ENDDO
  ENDDO
  
  RETURN

END SUBROUTINE set_event_number_to_value

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE decrease_resolution(domsize_x,domsize_y,miss,field)

  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: domsize_x,domsize_y
  REAL, INTENT(INOUT)   :: field(domsize_x,domsize_y)
  REAL                  :: value_mean
  REAL, INTENT(IN)      :: miss            ! value<miss ==> missing_value
  INTEGER               :: ix,iy

  ! build 2x2 km grid out of 1x1 km grid
  DO ix=1,domsize_x-1,2
    DO iy=1,domsize_y-1,2
      IF (field(ix,iy) .LT. miss .OR. field(ix,iy+1) .LT. miss .OR. field(ix+1,iy) .LT. miss .OR. field(ix+1,iy+1) .LT. miss) THEN
	value_mean = miss-1
      ELSE
        value_mean = NINT((field(ix,iy)+field(ix,iy+1)+field(ix+1,iy)+field(ix+1,iy+1))*0.25)
      ENDIF
      field(ix  ,iy  ) = value_mean
      field(ix  ,iy+1) = value_mean
      field(ix+1,iy  ) = value_mean
      field(ix+1,iy+1) = value_mean
    ENDDO
  ENDDO

END SUBROUTINE decrease_resolution

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


RECURSIVE SUBROUTINE area(miss, cutoff, ii,ij,domsize_x,domsize_y,n_fields,input_field, &
                  occupied,nevent,event_number,totarea,field_mean,field_min,field_max, &
		  COMx,COMy,delete_cell,xfirst,xlast,yfirst,ylast)
  IMPLICIT NONE
  INTEGER, INTENT(IN)	 :: ii, ij
  INTEGER, INTENT(IN)	 :: domsize_x,domsize_y
  REAL, INTENT(IN)	 :: miss
  INTEGER, INTENT(IN)	 :: n_fields
  REAL, INTENT(IN)       :: input_field(domsize_x,domsize_y,n_fields+1)
  INTEGER,INTENT(INOUT)  :: event_number(domsize_x,domsize_y)
  LOGICAL, INTENT(INOUT) :: occupied(domsize_x,domsize_y)
  INTEGER, INTENT(INOUT) :: nevent
  REAL, INTENT(INOUT)    :: field_mean(n_fields+1),field_min(n_fields+1),field_max(n_fields+1)
  REAL, INTENT(INOUT)    :: totarea
  REAL, INTENT(INOUT)	 :: COMx,COMy
  INTEGER, INTENT(INOUT) :: xfirst,xlast,yfirst,ylast
  INTEGER		 :: i, ii_mod, ij_mod,fieldid
  INTEGER		 :: icell(4), jcell(4)
  REAL, INTENT(IN)       :: cutoff
  LOGICAL, INTENT(INOUT) :: delete_cell
  
  ! Indices for all 4 flow directions
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  
  ii_mod = MOD(ii-1+domsize_x,domsize_x)+1
  ij_mod = MOD(ij-1+domsize_y,domsize_y)+1
  
  IF (.NOT. occupied(ii_mod,ij_mod)) THEN
     ! take care for periodic boundary conditions:
     IF (ii .LT. xfirst) xfirst=ii
     IF (ii .GT. xlast) xlast=ii
     IF (ij .LT. yfirst) yfirst=ij
     IF (ij .GT. ylast) ylast=ij
     occupied(ii_mod,ij_mod) = .TRUE.
     IF (input_field(ii_mod,ij_mod,1) .GE. cutoff) THEN
	! center of mass is now weighted by intensity!!!
        COMx = COMx + ii*input_field(ii_mod,ij_mod,1)
        COMy = COMy + ij*input_field(ii_mod,ij_mod,1)
        event_number(ii_mod,ij_mod) = nevent
	DO fieldid=1,n_fields+1
	  field_mean(fieldid) = field_mean(fieldid) + input_field(ii_mod,ij_mod,fieldid)
	  field_min(fieldid)  = MIN(field_min(fieldid),input_field(ii_mod,ij_mod,fieldid))
	  field_max(fieldid)  = MAX(field_max(fieldid),input_field(ii_mod,ij_mod,fieldid))
	ENDDO
        totarea = totarea + 1
        DO i=1,4
           CALL area(miss, cutoff, icell(i),jcell(i),domsize_x,domsize_y,n_fields,input_field,&
	             occupied,nevent,event_number,totarea,field_mean,field_min,field_max, &
		     COMx,COMy,delete_cell,xfirst,xlast,yfirst,ylast)
        ENDDO
     ELSEIF (input_field(ii_mod,ij_mod,1) .LE. miss) THEN
        ! if cell touches missing value, this cell will be deleted
        delete_cell = .TRUE.
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE area

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RECURSIVE SUBROUTINE overlay_area(miss, ii,ij,domsize_x,domsize_y,input_field,occupied,nevent, &
                                  event_number,delete_cell,xfirst,xlast,yfirst,ylast)
  IMPLICIT NONE
  INTEGER, INTENT(IN)	 :: ii, ij
  INTEGER, INTENT(IN)	 :: domsize_x,domsize_y
  REAL, INTENT(IN)	 :: miss
  REAL, INTENT(IN)       :: input_field(domsize_x,domsize_y)
  INTEGER,INTENT(INOUT)  :: event_number(domsize_x,domsize_y)
  LOGICAL, INTENT(INOUT) :: occupied(domsize_x,domsize_y)
  INTEGER, INTENT(INOUT) :: nevent
  INTEGER, INTENT(INOUT) :: xfirst,xlast,yfirst,ylast
  INTEGER		 :: i, ii_mod, ij_mod
  INTEGER		 :: icell(4), jcell(4)
  LOGICAL, INTENT(INOUT) :: delete_cell
  
  ! Indices for all 4 flow directions
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  
  ii_mod = MOD(ii-1+domsize_x,domsize_x)+1
  ij_mod = MOD(ij-1+domsize_y,domsize_y)+1
  
  IF (.NOT. occupied(ii_mod,ij_mod)) THEN
     ! take care for periodic boundary conditions:
     IF (ii .LT. xfirst) xfirst=ii
     IF (ii .GT. xlast) xlast=ii
     IF (ij .LT. yfirst) yfirst=ij
     IF (ij .GT. ylast) ylast=ij
     occupied(ii_mod,ij_mod) = .TRUE.
     IF (input_field(ii_mod,ij_mod) .GE. 0.5) THEN
	! center of mass is now weighted by intensity!!!
        event_number(ii_mod,ij_mod) = nevent
	DO i=1,4
           CALL overlay_area(miss, icell(i),jcell(i),domsize_x,domsize_y,input_field,occupied,nevent, &
	             event_number,delete_cell,xfirst,xlast,yfirst,ylast)
        ENDDO
     ELSEIF (input_field(ii_mod,ij_mod) .LE. miss) THEN
        ! if cell touches missing value, this cell will be deleted
        delete_cell = .TRUE.
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE overlay_area

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! linear time interpolation of velocity field

SUBROUTINE time_interpolation(domsize_x,domsize_y,time_steps,nx_bins,ny_bins,nt_bins,lperiodic, &
                              coarse_vel_x,coarse_vel_y,ix,iy,it,vel_x,vel_y)

INTEGER, INTENT(IN)  :: domsize_x,domsize_y
INTEGER, INTENT(IN)  :: time_steps
INTEGER, INTENT(IN)  :: nx_bins
INTEGER, INTENT(IN)  :: ny_bins
INTEGER, INTENT(IN)  :: nt_bins
INTEGER, INTENT(IN)  :: ix,iy,it
LOGICAL, INTENT(IN)  :: lperiodic        ! periodic boundaries?

REAL, INTENT(IN)     :: coarse_vel_x(nx_bins,ny_bins,nt_bins)
REAL, INTENT(IN)     :: coarse_vel_y(nx_bins,ny_bins,nt_bins)

REAL, INTENT(OUT)    :: vel_x,vel_y

REAL                 :: vel_x_prev,vel_y_prev
REAL                 :: vel_x_next,vel_y_next
REAL                 :: weight
INTEGER              :: nt

nt = MIN(FLOOR((REAL(it-1)/REAL(time_steps))*nt_bins)+1,nt_bins)
weight = (REAL(it-1)/REAL(time_steps))*nt_bins
weight = weight-FLOOR(weight)

IF ((nt .EQ. 1 .AND. weight .LT. 0.5) .OR. (nt .EQ. nt_bins .AND. weight .GT. 0.5)) THEN
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_x(:,:,nt),ix,iy,vel_x)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_y(:,:,nt),ix,iy,vel_y)
ELSEIF (weight .LT. 0.5) THEN
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_x(:,:,nt-1),ix,iy,vel_x_prev)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_y(:,:,nt-1),ix,iy,vel_y_prev)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_x(:,:,nt),ix,iy,vel_x_next)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_y(:,:,nt),ix,iy,vel_y_next)
  vel_x = (0.5-weight)*vel_x_prev + (weight+0.5)*vel_x_next
  vel_y = (0.5-weight)*vel_y_prev + (weight+0.5)*vel_y_next
ELSE
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_x(:,:,nt),ix,iy,vel_x_prev)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_y(:,:,nt),ix,iy,vel_y_prev)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_x(:,:,nt-1),ix,iy,vel_x_next)
  CALL bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_vel_y(:,:,nt-1),ix,iy,vel_y_next)
  vel_x = (1.5-weight)*vel_x_prev + (weight-0.5)*vel_x_next
  vel_y = (1.5-weight)*vel_y_prev + (weight-0.5)*vel_y_next
ENDIF

END SUBROUTINE time_interpolation


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! bilinear spatial interpolation of coarse velocity field

SUBROUTINE bilinear_interpolation(domsize_x,domsize_y,nx_bins,ny_bins,lperiodic,coarse_field,ix,iy,out_value)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: domsize_x,domsize_y
INTEGER, INTENT(IN)  :: nx_bins
INTEGER, INTENT(IN)  :: ny_bins
INTEGER, INTENT(IN)  :: ix,iy
LOGICAL, INTENT(IN)  :: lperiodic        ! periodic boundaries?

REAL, INTENT(IN)     :: coarse_field(nx_bins,ny_bins)

REAL, INTENT(OUT)    :: out_value

REAL                 :: halfbin_x,halfbin_y
REAL                 :: weight_x, weight_y
INTEGER              :: nx,ny
INTEGER              :: nx_west, nx_east, ny_south, ny_north

halfbin_x = NINT(REAL(domsize_x)/REAL(nx_bins)/2.)
halfbin_y = NINT(REAL(domsize_y)/REAL(ny_bins)/2.)

nx = FLOOR((REAL(ix+halfbin_x)/REAL(domsize_x))*nx_bins)
weight_x = (REAL(ix+halfbin_x)/REAL(domsize_x))*nx_bins
weight_x = weight_x-FLOOR(weight_x)

ny = FLOOR((REAL(iy+halfbin_y)/REAL(domsize_y))*ny_bins)
weight_y = (REAL(iy+halfbin_y)/REAL(domsize_y))*ny_bins
weight_y = weight_y-FLOOR(weight_y)

!catch
IF (nx<0 .OR. nx>nx_bins .OR. ny<0 .OR. ny>ny_bins) THEN
  WRITE (*,*) "ERROR in bilinear_interpolation: indices out of range."
  STOP
ENDIF

IF (lperiodic) THEN ! periodix boundaries

  IF (nx==0) THEN
    nx_west = nx_bins
    nx_east = nx+1
  ELSEIF (nx==nx_bins) THEN
    nx_west = nx
    nx_east = 1
  ELSE
    nx_west = nx
    nx_east = nx+1
  ENDIF

  IF (ny==0) THEN
    ny_south = ny_bins
    ny_north = ny+1
  ELSEIF (ny==ny_bins) THEN
    ny_south = ny
    ny_north = 1
  ELSE
    ny_south = ny
    ny_north = ny+1
  ENDIF

ELSE ! no periodic boundaries

  IF (nx==0) THEN
    nx_west = nx+1
    nx_east = nx+1
  ELSEIF (nx==nx_bins) THEN
    nx_west = nx
    nx_east = nx
  ELSE
    nx_west = nx
    nx_east = nx+1
  ENDIF

  IF (ny==0) THEN
    ny_south = ny+1
    ny_north = ny+1
  ELSEIF (ny==ny_bins) THEN
    ny_south = ny
    ny_north = ny
  ELSE
    ny_south = ny
    ny_north = ny+1
  ENDIF

ENDIF

out_value = coarse_field(nx_west,ny_south)*(1.-weight_x)*(1.-weight_y) + &
            coarse_field(nx_west,ny_north)*(1.-weight_x)*weight_y + &
            coarse_field(nx_east,ny_south)*weight_x*(1.-weight_y) + &
            coarse_field(nx_east,ny_north)*weight_x*weight_y

END SUBROUTINE


