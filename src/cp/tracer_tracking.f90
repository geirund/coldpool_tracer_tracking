! 2018 Sep 08 O Henneb
! read rain tracking output txt for COGs (loop trough every line)
! read 2D velocity field (time-stepwise)
! read rain track cells for precipitation boundaries (timestepwise)
! SUBROUTINE neigh idetifies precipitation boundaries
! SUBROUTINE set_tracer sets tracer at the boundaries
! SUBROUTINE ....
! SUBROUTINE ... 
! 
! 2018 Oct:
! change initial tracer placement from precip boundary to circles
! calculate effective radius r = sqrt(A)/pi 
! Problem: tracers which are inside the precip area never chatch up
! alternatives: 
!    set tracer at outremost circle
!    distribute tracer equally around outline
!      * sort boundary points by angle 
!      * calculate distance from one gp to another and sum it up ->
!        circulferance    
!      * divide by circum by no of desired tracer to get distance 
!
! 2018 Nov:
! change order of CALL routine: 
! output before update of tracer -> output of windvector match the position of
! tracer at output now, velocities before subtimestepping
!
! 2018 Dec:
! include option to track tracer base don radial velocity instead of full vector
! output of radial and tangential velocity also included
!
! 2018 Dec:
! allow precipitation tracking starting with smaller precipitation area, but set
! tracer not untill a larger size is reached which is given in jobfile  
! number of precipitation cells is given as input now. But can be done nicer
!
! 2018 Dec
! set tracer at multiple times
!
! OUTPUT: 
! traced(CP-ID,int tracer ID, property) 
! traced(:,:, 1)   x pos of tracer 
! traced(:,:, 2)   y pos of tracer
! traced(:,:, 3)   x gp 
! traced(:,:, 4)   y gp
! traced(:,:, 5)   distance to COG
! traced(:,:, 6)   timestep 
! traced(:,:, 7)   CP age defined by onset of precip ??? check this 
! traced(:,:, 8)   angle to COG
! traced(:,:, 9)   CP ID
! traced(:,:,10)   gloabl tracer number 
! traced(:,:,11)   active 0-1 
! traced(:,:,12)   ID merging 
! traced(:,:,13)   u 
! traced(:,:,14)   v
! traced(:,:,15)   distance from tracer to COG in x direction 
! traced(:,:,16)   merger 
! traced(:,:,17)   distance from tracer to COG in x direction
! traced(:,:,18)   timesteps set after precip start 
! traced(:,:,19)   vrad
! traced(:,:,20)   vtan
PROGRAM cp_tracking

USE netcdf
USE cp_parameters !, ONLY: dsize_x, dsize_y, max_tracer_CP, dt ! &
!                          resolution, dt ! resolution in m! dt in sec

IMPLICIT NONE
REAL, ALLOCATABLE    :: input_field(:,:)      ! eg precip field to find COGs
REAL, ALLOCATABLE    :: vel(:,:,:)               ! velocity field
REAL, ALLOCATABLE    :: nneighb(:,:)             ! identifier for cell boundaries
REAL, ALLOCATABLE    :: COMx(:), COMy(:)         ! store COG
REAL, ALLOCATABLE    :: rmax(:)                  ! area of precip
!INTEGER, ALLOCATABLE :: IDstart(:)               ! first timestep of precip 
REAL, ALLOCATABLE    :: track_numbers(:,:)       ! ID for precip and coldpoolsobjects
INTEGER, ALLOCATABLE :: already_tracked(:)       ! memory of cell counter
REAL, ALLOCATABLE    :: traced(:,:,:)            ! traced information, CP ID x internal tarcer ID x properties
REAL, ALLOCATABLE    :: traced_dummy(:,:)        ! dummy to sort tracer by angle
REAL, ALLOCATABLE    :: traced_prev(:,:,:)       ! traced information from previous timestep 
INTEGER,ALLOCATABLE  :: cpio(:,:)                ! to identify merger (dim1) and start time (dim2) splitting events have start time 0 to avoid them
INTEGER              :: ID                       ! local ID from rain track
INTEGER              :: tts                      ! timestep
!REAL                 :: areain                   ! size of precip cell
INTEGER              :: i                        ! running index
INTEGER              :: ierr                     ! error index for reading files
REAL                 :: xCOG, yCOG               ! buffer variable for read COGs
INTEGER              :: onset,begin              ! first cell reaches treshold, beginning of tracking
INTEGER              :: max_no_of_cells          ! maximum number if CPs (can be retrieved from the number of rain cells
!INTEGER              :: max_tracer_CP            ! max no of tracers per CP
INTEGER              :: max_tracers              ! max no of commulative tracers 
INTEGER,ALLOCATABLE  :: tracpo(:,:)              ! keeps track of first two indices in traced (CP ID, internal tracer ID) for every tracer 
INTEGER              :: srv_header_input(8)
INTEGER              :: timestep
!INTEGER              :: ntracer
INTEGER              :: count_tracer             ! counts the internal tracer in an individual CP

!INITIALIZE some values
namelist /INPUTgeneral/ dsize_x, dsize_y, dt, res, lformat 
namelist /INPUTtracer/ max_tracer_CP, nTrSet, ntset, max_age, rad
namelist /INPUTIO/ odir
open(100,file='job/namelist.dat')
read(100,nml=INPUTgeneral)
read(100,nml=INPUTIO)
read(100,nml=INPUTtracer)


!get number of tracked precip events
!CALL getarg(max_no_of_cells)
OPEN(1111,FILE='na.txt')
READ(1111,*) max_no_of_cells
!max_tracer_CP = nTrSet * ntset
max_tracers = max_no_of_cells*max_tracer_CP
write(*,*) "max_tracers", max_tracers, "max_no_of_cells", max_no_of_cells, "max_tracer_CP", max_tracer_CP 
! allocate fields
ALLOCATE(cpio(max_no_of_cells,3))
ALLOCATE(traced(max_no_of_cells,max_tracer_CP,20))
ALLOCATE(traced_dummy(max_tracer_CP,20))
ALLOCATE(traced_prev(max_no_of_cells,max_tracer_CP,20))
!ALLOCATE(IDstart(max_no_of_cells))
ALLOCATE(COMx(max_no_of_cells))
ALLOCATE(COMy(max_no_of_cells))
ALLOCATE(rmax(max_no_of_cells))
ALLOCATE(vel(dsize_x,dsize_y,2))
ALLOCATE(track_numbers(dsize_x,dsize_y))
ALLOCATE(nneighb(dsize_x,dsize_y))
ALLOCATE(already_tracked(max_no_of_cells))
ALLOCATE(input_field(dsize_x,dsize_y))
ALLOCATE(tracpo(2,max_tracers))

 i =1
 OPEN(1,FILE=trim(odir) // '/input/cp/mergingCPs.txt',FORM='formatted',ACTION='read',IOSTAT=ierr) 
 IF ( ierr == 0) then
   DO i =1, max_no_of_cells,1
write(*,*) i
     READ(1,*,END=400) cpio(i,1), cpio(i,2), cpio(i,3)
write(*,*) cpio(i,1), cpio(i,2), cpio(i,3)
   END DO
 ELSE 
     write(*,*) 'Beim Oeffnen der Datei ist ein Fehler Nr.', ierr,' aufgetreten'
 END IF
400 CONTINUE
write(*,*) 'read data'
OPEN(2,FILE=trim(odir) // '/output/raincell/irt_tracks_mask.srv',    FORM='unformatted', ACTION='read')
OPEN(40,FILE=trim(odir) // '/output/cp/coldpool_tracer_out.txt',FORM='formatted', ACTION='write')
if (trim(lformat) == 'srv') then
  OPEN(4,FILE=trim(odir) // '/input/cp/input_u.srv',FORM='unformatted',ACTION='read')
  OPEN(5,FILE=trim(odir) // '/input/cp/input_v.srv',FORM='unformatted',ACTION='read')
end if

write(*,*) i, 'rain cells found'

!INITIALIZE some values
already_tracked(:) = 0
traced(:,:,:) = 0.
traced_prev(:,:,:) = 0.
count_tracer = 1 ! counts individual pixels !OCH was ist mit pixeln gemeint? der
!IDstart(:) = 0
timestep=0 

! read when first precip is tracked
 OPEN(3,FILE=trim(odir) // '/input/cp/tracks_body_time_sorted.txt',FORM='formatted',ACTION='read')
 153 FORMAT (7X,I5,9X,I3,94X,F10.4,7X,F10.4)
 onset = minval(cpio(:,3),1)  !first tracers will be set when first raincell is larger than agiven trashold
 write(*,*) 'tracer start at' , onset
 ! read first line of precip tracking, to get first timestep of precip
 READ(3,153,END=200)  ID, tts, xCOG, yCOG
 begin = tts

 DO !start main loop
 write(*,*) 'at timestep', timestep
 IF (timestep .ge. maxval(cpio(:,3),1)) then
 GOTO 5000
 END IF 
! only for testing
   timestep=timestep+1 ! OCH: starts with 0 ?
   !if (IDstart(ID) .eq. 0) IDstart(ID) = cpio(ID,3) !new CP 
   track_numbers(:,:)=0 
! read velocity files
   if (trim(lformat) == 'srv') then
     READ (4,END=200) srv_header_input
     READ (4) vel(:,:,1)
     READ (5,END=200) srv_header_input
     READ (5) vel(:,:,2)
   else if (trim(lformat) == 'nc') then
     CALL read_nc_2D (vel(:,:,1),'U',trim(odir) // '/input/cp/input_u.nc',timestep)
     CALL read_nc_2D (vel(:,:,2),'V',trim(odir) // '/input/cp/input_v.nc',timestep)
     CALL read_nc_2D (track_numbers,'trackID',trim(odir) // '/output/raincell/irt_tracks_mask_fulltime.nc',timestep)
   end if
   ! read velocity files until start of tracking is reached, if so, also read
   ! tracking
   IF (timestep .GE. begin) THEN !max(onset,time_step_event)) THEN
     ! store the initially read COG at first timestep 
     ! or previously read 
     COMx(ID) = xCOG
     COMy(ID) = yCOG
     !area(ID) = areain
     DO WHILE (tts .eq. timestep ) !as long as CPs at the same tiemstep are read
       READ(3,153,END=200)  ID, tts, xCOG,yCOG ! read next line 
       IF (tts .eq. timestep ) THEN ! and store COG if still at the same timestep
         COMx(ID) = xCOG
         COMy(ID) = yCOG
         !area(ID) = areain
       END IF
       ! now the next timestep is already read
       ! dont overwrite COG (will be stored in xCOG until the next timestep)
     END DO ! all COGs for this timestep are read
     ! reading velocity field and passive tracer
     ! reading the track input files
   if (trim(lformat) == 'srv') then
     READ(2,END=200) srv_header_input
     READ(2) track_numbers(:,:)
   end if
     IF (timestep .GE. onset) THEN !not neccesarry to call any of the routines before onset (=first track has reached minimum size for tracer tracking)  
       ! identification of edges, and set tracer
       nneighb(:,:)=0 
       !CALL neigh(track_numbers,nneighb)
       !CALL set_tracer(nneighb, track_numbers,max_no_of_cells, timestep,traced,&
       !            count_tracer,already_tracked,tracpo,max_tracers,cpio(:,1:2))
       
       ! alterbative to set tracers at edges of precip, set tracer at circle
       ! around precip
       CALL maxcell(track_numbers,COMx,COMy,rmax,max_no_of_cells)
       CALL initCircle(max_no_of_cells,timestep,traced, COMx,COMy,&
                       rmax,cpio,max_tracers,tracpo,count_tracer)
       CALL velocity_interpol(vel(:,:,1),vel(:,:,2),timestep,traced, count_tracer, &
                          max_no_of_cells,tracpo,max_tracers)
       CALL geometry(traced,COMx,COMy,already_tracked,max_no_of_cells)
       CALL radvel(traced,already_tracked,max_no_of_cells)
  
  !     DO i =1,max_no_of_cells ! loop trough all cps with tracer
  !       IF (already_tracked(i) .gt. 1 ) THEN
  !       ! reset dummy first
  !         traced_dummy(:,:) = 0.
  !         traced_dummy = traced(i,:,:)
  !         CALL sort(traced_dummy(1:already_tracked(i),:),already_tracked(i))
  !!         if (traced(i,1,11)  .eq. 0) then !stop tracer only if precip has stoped
  !           !CALL oneside(traced_dummy(1:already_tracked(i),8),traced(i,:,:),already_tracked(i))
  !!         end if
  !       END IF
  !     END DO
       !CALL time_dev(traced,traced_prev,max_no_of_cells,max_tracer_CP,count_tracer,tracpo,max_tracers)
       !traced_prev = traced
       CALL write_output(traced,max_tracers,count_tracer,timestep,tracpo,&
                         max_no_of_cells,COMx,COMy)
  !     if rad then
  !       CALL radial_update(timestep,traced,count_tracer,max_no_of_cells,tracpo,max_tracers) 
  !     else
         CALL update_tracer(vel(:,:,1),vel(:,:,2),timestep,traced, count_tracer, &
                          max_no_of_cells,tracpo,max_tracers)
  !     end if
     END IF !if onset is reached
   END IF ! if begin is reached  
 END DO
 200 CONTINUE
 WRITE(*,*) "finished main loop" 
 5000 CONTINUE
 

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
    CALL check(nf90_inquire_dimension(ncid, dimIDs(4), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nz))
write(*,*) ctime
    CALL check(nf90_get_var(ncid, rhVarId, poutput(1:nx,1:ny), start =(/1,1,1,ctime/),&
                                                    count= (/nx, ny,1, 1/)))
    CALL check(nf90_close(ncid))
!write(*,*) poutput(1:3,1:3)
!
  END SUBROUTINE read_nc_2D

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

!--------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------------------------------------
! this routine is replaced by initCircle 
!--------------------------------------------------
SUBROUTINE neigh(track_numbers,nneighb)

  USE cp_parameters, ONLY : dsize_x, dsize_y
  IMPLICIT NONE
  INTEGER                :: i,j, imodm, imodp, jmodm,jmodp
  REAL, INTENT(IN)       :: track_numbers(dsize_x,dsize_y)
  REAL, INTENT(INOUT)    :: nneighb(dsize_x,dsize_y)
  nneighb(:,:) = 0
  DO i =2,dsize_x-1
    DO j =2,dsize_y-1
      IF (track_numbers(i,j) .gt. 0) THEN !check only gps with rain cell 
        imodp = mod(i+dsize_x,dsize_x)+1
        imodm = mod(i-2+dsize_x,dsize_x)+1
        jmodp = mod(j+dsize_y,dsize_y)+1
        jmodm = mod(j-2+dsize_y,dsize_y)+1
        IF (track_numbers(i,j) .ne. track_numbers(imodp,j)) nneighb((/i/),j)=1
        IF (track_numbers(i,j) .ne. track_numbers(imodm,j)) nneighb((/i/),j)=1
        IF (track_numbers(i,j) .ne. track_numbers(i,jmodp)) nneighb(i,(/j/))=1
        IF (track_numbers(i,j) .ne. track_numbers(i,jmodm)) nneighb(i,(/j/))=1
  
        ! diagonal neigh
        IF (track_numbers(i,j) .ne. track_numbers(imodp,jmodp)) nneighb(i,j) =1
        IF (track_numbers(i,j) .ne. track_numbers(imodm,jmodp)) nneighb(i,j) =1
        IF (track_numbers(i,j) .ne. track_numbers(imodp,jmodm)) nneighb(i,j) =1
        IF (track_numbers(i,j) .ne. track_numbers(imodm,jmodm)) nneighb(i,j) =1
      END IF
    END DO
  END DO
END SUBROUTINE neigh

SUBROUTINE maxcell(track_numbers,COMx,COMy,rmax,max_no_of_cells)
  USE cp_parameters, ONLY : dsize_x, dsize_y
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: max_no_of_cells
  REAL, INTENT(OUT)      :: rmax(max_no_of_cells)
  INTEGER                :: i,j
  REAL, INTENT(IN)       :: track_numbers(dsize_x,dsize_y)
  REAL, INTENT(IN)       :: COMx(max_no_of_cells), COMy(max_no_of_cells)
  REAL                   :: rt, dx, dy 
  rmax(:) = 0
  
  DO i =2,dsize_x-1
    DO j =2,dsize_y-1
      IF (track_numbers(i,j) .gt. 0) THEN
        !write(*,*) track_numbers(i,j)
        dx = mod(i +(dsize_x/2. - COMx(int(track_numbers(i,j))))+dsize_x-1,float(dsize_x))+1  - dsize_x/2.
        dy = mod(j +(dsize_y/2. - COMy(int(track_numbers(i,j))))+dsize_y-1,float(dsize_y))+1  - dsize_y/2.
        rt = sqrt(dx**2 + dy**2) 
        rmax(int(track_numbers(i,j))) =max(rmax(int(track_numbers(i,j))),rt)       
      END IF
    END DO
  END DO
END SUBROUTINE maxcell
!--------------------------------------------------------------------------------------
! Calculate  circle around COG dependent on size for initial tracer placement
! replaces routine neighbours
!--------------------------------------------------------------------------------------
!SUBROUTINE initCircle(max_no_of_cells,ts,IDstart,traced, COMx,COMy, &
!                      rmax,cpio,max_tracers,tracpo,count_tracer)
SUBROUTINE initCircle(max_no_of_cells,ts,traced, COMx,COMy, &
                      rmax,cpio,max_tracers,tracpo,count_tracer)
   USE cp_parameters, ONLY : dsize_x, dsize_y,max_tracer_CP, nTrSet, ntset
  
  INTEGER, INTENT(IN)       :: max_no_of_cells, max_tracers, ts
!  INTEGER, INTENT(IN)       :: IDstart(max_no_of_cells)
  REAL, INTENT(IN)          :: rmax(max_no_of_cells),&
                               COMx(max_no_of_cells), COMy(max_no_of_cells)
  INTEGER, INTENT(IN)       :: cpio(max_no_of_cells,3)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER, INTENT(INOUT)    :: tracpo(2,max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: pi, inc !, reff 
  INTEGER                   :: i, j, pj, it
  pi = 2.*asin(1.)
  inc = 2*pi/nTrSet !max_tracer_CP
  do i = 1,max_no_of_cells,1 ! loop trough precip cells

   if (cpio(i,2) .ne. 0 ) then ! avoid splitting events 
!   if (IDstart(i) == ts) then ! set new circle when precip begins
   if (ts .ge. cpio(i,3) .and. ts .le. cpio(i,3)+ntset) then ! set new circle when precip begins
     !reff = sqrt(area(i))/pi 
      do pj = 1,nTrSet,1 !max_tracer_CP,1 ! loop trough tracer number

        j =pj+(ts-cpio(i,3)) * nTrSet
        tracpo(:,count_tracer)=(/i,j/)
        count_tracer           = count_tracer + 1
        traced(i,j, 1) = MOD(COMx(i) + rmax(i)*cos(pj*inc)-1.,float(dsize_x))+1.
        traced(i,j, 2) = MOD(COMy(i) + rmax(i)*sin(pj*inc)-1.,float(dsize_y))+1.
        traced(i,j, 3) = nint(traced(i,j, 1))
        traced(i,j, 4) = nint(traced(i,j, 2))
        traced(i,j, 5) = rmax(i)
        traced(i,j, 6) = ts ! current time
        if (cpio(i,1) .ne. cpio(i,2)) then
          traced(i,j,16) = 0 
          traced(i,j,7) = ts-cpio(i,3)  !traced(cpio(i,2),j,7)  !set to age that the larger event of the merger has 
        else
          traced(i,j, 7) = ts-cpio(i,3)   ! age 
          traced(i,j, 16) = 1
          traced(i,j,12) = cpio(i,2)
        end if
        traced(i,j, 8) = pj*inc 
        traced(i,j, 9) = i 
        traced(i,j,10) = count_tracer 
        traced(i,j,11) = 1   ! active tracer 
        traced(i,j,12) = cpio(i,2)  ! start
        traced(i,j,18) =  ts-cpio(i,3)
        !traced(i,j,13) = 0   ! u 
        !traced(i,j,14) = 0   ! v
      end do
   end if
   end if
  end do
!  do i = 1,max_no_of_cells,1
!    if (cpio(i,1) .ne. cpio(i,2)) then
!      traced(i,:,12) = cpio(i,2)
!    end if
!  end do
END SUBROUTINE 


!--------------------------------------------------------------------------------------
! SET TRACER at their initial position acording to precipitation outlines
!-----------------------------------------------------------------------------------
SUBROUTINE set_tracer(nneighb,track_numbers, max_no_of_cells,&
                      timestep,traced,count_tracer,&
                      already_tracked,tracpo,max_tracers,cpio)
USE cp_parameters, ONLY : dsize_x, dsize_y,max_tracer_CP
  REAL, INTENT(INOUT)       :: nneighb(dsize_x,dsize_y)
  REAL, INTENT(IN)          :: track_numbers(dsize_x,dsize_y)
  INTEGER, INTENT(IN)       :: cpio(1700,2)
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells, max_tracers
!  INTEGER, INTENT(IN)       :: max_tracer_CP
  INTEGER, INTENT(INOUT)    :: tracpo(2,max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER                   :: ix, iy
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  DO iy=1, dsize_y
    DO ix=1, dsize_x
      IF (nneighb(ix,iy) .EQ. 1 ) THEN ! precip edge found
        IF (already_tracked(INT(track_numbers(ix,iy))) .LT. max_tracer_CP) THEN ! upper bound of tracer no
          IF (traced(INT(track_numbers(ix,iy)),1,7) .lt.1 )  THEN ! set tracer only when CP has not been tracked yet
            if  (cpio(INT(track_numbers(ix,iy)),2) .ne. 0 ) then ! avoid splitting events 
              ! ---- set tracer at gp center ----
              ! increase internal tracer counter for individual CP
              already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1
              ! set initial position of tracer
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1) = ix
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2) = iy
              ! keep track of indices for CP and interal tracer
              tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
              count_tracer           = count_tracer + 1  ! global tracer number
              ! ---- set tracer in between ----
              already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1)=ix + 1./3.
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2)=iy
              tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
              count_tracer           = count_tracer + 1

              already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1)=ix- 1./3.
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2)=iy
              tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
              count_tracer           = count_tracer + 1

              already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1)=ix
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2)=iy+ 1./3.
              tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
              count_tracer           = count_tracer + 1

              already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1)=ix
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2)=iy- 1./3.
              tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
              count_tracer           = count_tracer + 1

              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))-4: &
                                               already_tracked(INT(track_numbers(ix,iy))),6) = timestep     ! timestep when tracking begins
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))-4: &
                                               already_tracked(INT(track_numbers(ix,iy))),7) = 0   ! tracer age ???timestep 
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))-4: &
                                               already_tracked(INT(track_numbers(ix,iy))),9)=track_numbers(ix,iy) ! TRACK ID 
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))-4: &
                                               already_tracked(INT(track_numbers(ix,iy))),12)= cpio(INT(track_numbers(ix,iy)),2)
              ! set all new tracer active
              traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))-4: &
                                               already_tracked(INT(track_numbers(ix,iy))),11)= 1
            end if ! avoid splitting events
          END IF ! set tracer only at beginning
        ELSE 
          write(*,*) "maximum number of tracers per CP was exceeded, may affect the identification of CP"
          write(*,*) "CP ID:", INT(track_numbers(ix,iy)),"tracers:", already_tracked(INT(track_numbers(ix,iy)))
        END IF ! end upper bound of tracer 
      ENDIF ! end precip edge .eq. 1
    ENDDO ! end loop trough x valus
  ENDDO ! end loop trough y values

END SUBROUTINE set_tracer

!--------------------------------------------------------------------------------------
! UPDATE TRACER along the horizontal wind field
!-----------------------------------------------------------------------------------
SUBROUTINE velocity_interpol(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers)
USE cp_parameters, ONLY : dsize_x, dsize_y, res, dt, max_tracer_CP

  INTEGER, INTENT(IN)       :: timestep,max_no_of_cells, max_tracers
  REAL, INTENT(IN)          :: velx(dsize_x,dsize_y), &
                               vely(dsize_x,dsize_y)
  INTEGER, INTENT(IN)       :: count_tracer !,max_tracer_CP
  !INTEGER, INTENT(IN)       :: cp_field(dsize_x,dsize_y)
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: ix_new, iy_new, vx_intp, vy_intp
  REAL                      :: ix, iy, ixt, iyt    !t are half level (interface)
!  INTEGER                   :: ix_round, iy_round, &
  INTEGER                   :: ix_new_round, iy_new_round, start_time
  INTEGER                   :: ix_l, iy_l, ix_r, iy_r
  REAL                      :: wgt_x, wgt_y, wgt_xt, wgt_yt
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  INTEGER                   :: sub_dt   ! subtimstepping

 sub_dt =60 ! update evry min
   
  it = 1 ! counter trough traced
  DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until here
    ! determining the first time step of the event
     start_time=INT(traced(tracpo(1,it),tracpo(2,it),6))
     ! determining how many timesteps have passed since then
     tracer_ts =timestep-start_time+1

     ! getting the previous positions
     ix=traced(tracpo(1,it),tracpo(2,it),1)
     iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints

     !for u-Wind defined on xt (half level in x direction)
     ixt = ix-0.5

     wgt_xt  = MOD(ixt,1.)
     wgt_y  = MOD(iy,1.)

     ix_l = MOD(INT(ixt-wgt_xt)-1+dsize_x,dsize_x)+1
     iy_l = MOD(INT(iy-wgt_y)-1+dsize_y,dsize_y)+1
     ix_r = MOD(ix_l +dsize_x,dsize_x)+1
     iy_r = MOD(iy_l +dsize_y,dsize_y)+1

     vx_intp = velx(ix_l,iy_l)*(1-wgt_xt)*(1-wgt_y) &
             + velx(ix_l,iy_r)*(1-wgt_xt)*(  wgt_y) &
             + velx(ix_r,iy_l)*(  wgt_xt)*(1-wgt_y) &
             + velx(ix_r,iy_r)*(  wgt_xt)*(  wgt_y)

     !for v-Wind defined on xt (half level in x direction)
     iyt = iy-0.5

     wgt_x  = MOD(ix,1.)
     wgt_yt  = MOD(iyt,1.)

     ix_l = MOD(INT(ix-wgt_x)-1+dsize_x,dsize_x)+1
     iy_l = MOD(INT(iyt-wgt_yt)-1+dsize_y,dsize_y)+1

     ix_r = MOD(ix_l +dsize_x,dsize_x)+1
     iy_r = MOD(iy_l +dsize_y,dsize_y)+1

     vy_intp = vely(ix_l,iy_l)*(1-wgt_x)*(1-wgt_yt) &
             + vely(ix_l,iy_r)*(1-wgt_x)*(  wgt_yt) &
             + vely(ix_r,iy_l)*(  wgt_x)*(1-wgt_yt) &
             + vely(ix_r,iy_r)*(  wgt_x)*(  wgt_yt)
     !save the velocities of the tracer
     traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
     traced(tracpo(1,it),tracpo(2,it),14) = vy_intp

     it = it+1
   END DO

  RETURN
END SUBROUTINE velocity_interpol 

SUBROUTINE update_tracer(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers)
!OCH TO DO: dt and resolution should be parameter read by USE from module
USE cp_parameters, ONLY : dsize_x, dsize_y, res, dt, max_tracer_CP

  INTEGER, INTENT(IN)       :: timestep,max_no_of_cells, max_tracers
  REAL, INTENT(IN)          :: velx(dsize_x,dsize_y), &
                               vely(dsize_x,dsize_y)
  INTEGER, INTENT(IN)       :: count_tracer !,max_tracer_CP
  !INTEGER, INTENT(IN)       :: cp_field(dsize_x,dsize_y)
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: ix_new, iy_new, vx_intp, vy_intp
  REAL                      :: ix, iy, ixt, iyt    !t are half level (interface)
!  INTEGER                   :: ix_round, iy_round, &
  INTEGER                   :: ix_new_round, iy_new_round, start_time
  INTEGER                   :: ix_l, iy_l, ix_r, iy_r
  REAL                      :: wgt_x, wgt_y, wgt_xt, wgt_yt
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  INTEGER                   :: sub_dt   ! subtimstepping
  sub_dt =60 ! update evry min
  it = 1 ! counter trough traced
  DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until here
    ! determining the first time step of the event
     start_time=INT(traced(tracpo(1,it),tracpo(2,it),6))
     ! determining how many timesteps have passed since then
     tracer_ts =timestep-start_time+1

     ! getting the previous positions
     ix=traced(tracpo(1,it),tracpo(2,it),1)
     iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints

     !for u-Wind defined on xt (half level in x direction)
     do i = 1,5
       ixt = ix-0.5

       wgt_xt  = MOD(ixt,1.)
       wgt_y  = MOD(iy,1.)

       ix_l = MOD(INT(ixt-wgt_xt)-1+dsize_x,dsize_x)+1
       iy_l = MOD(INT(iy-wgt_y)-1+dsize_y,dsize_y)+1 

       ix_r = MOD(ix_l +dsize_x,dsize_x)+1   
       iy_r = MOD(iy_l +dsize_y,dsize_y)+1 

       vx_intp = velx(ix_l,iy_l)*(1-wgt_xt)*(1-wgt_y) &
               + velx(ix_l,iy_r)*(1-wgt_xt)*(  wgt_y) &
               + velx(ix_r,iy_l)*(  wgt_xt)*(1-wgt_y) &
               + velx(ix_r,iy_r)*(  wgt_xt)*(  wgt_y)

       !for v-Wind defined on xt (half level in x direction)
       iyt = iy-0.5

       wgt_x  = MOD(ix,1.)
       wgt_yt  = MOD(iyt,1.)

       ix_l = MOD(INT(ix-wgt_x)-1+dsize_x,dsize_x)+1
       iy_l = MOD(INT(iyt-wgt_yt)-1+dsize_y,dsize_y)+1    

       ix_r = MOD(ix_l +dsize_x,dsize_x)+1   
       iy_r = MOD(iy_l +dsize_y,dsize_y)+1

       vy_intp = vely(ix_l,iy_l)*(1-wgt_x)*(1-wgt_yt) &
               + vely(ix_l,iy_r)*(1-wgt_x)*(  wgt_yt) &
               + vely(ix_r,iy_l)*(  wgt_x)*(1-wgt_yt) &
               + vely(ix_r,iy_r)*(  wgt_x)*(  wgt_yt)

       !save the velocities of the tracer
       traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
       traced(tracpo(1,it),tracpo(2,it),14) = vy_intp
       ! update to new location
!OCH TO DO 
!first? update to temporary new position after 1m
!second: interpolate again 
!repeat 5 times
       ix_new = MOD((ix + sub_dt*vx_intp/res)-1.+FLOAT(dsize_x),FLOAT(dsize_x))+1.
       iy_new = MOD((iy + sub_dt*vy_intp/res)-1.+FLOAT(dsize_y),FLOAT(dsize_y))+1.
       ! uebergabe
       ix = ix_new
       iy = iy_new
     end do
     !save the velocities of the tracer
     traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
     traced(tracpo(1,it),tracpo(2,it),14) = vy_intp
!     ix_new = MOD((ix + dt*vx_intp/res)-1.+FLOAT(dsize_x),FLOAT(dsize_x))+1.
!     iy_new = MOD((iy + dt*vy_intp/res)-1.+FLOAT(dsize_y),FLOAT(dsize_y))+1.
     ix_new_round = nint(ix_new) ! maybe mod required
     iy_new_round = nint(iy_new)

     ! save new values for next loop
     traced(tracpo(1,it),tracpo(2,it),1) = ix_new
     traced(tracpo(1,it),tracpo(2,it),2) = iy_new
     traced(tracpo(1,it),tracpo(2,it),3) = float(ix_new_round) 
     traced(tracpo(1,it),tracpo(2,it),4) = float(iy_new_round)

     traced(tracpo(1,it),tracpo(2,it),6) = start_time !timestep
     traced(tracpo(1,it),tracpo(2,it),7) = traced(tracpo(1,it),tracpo(2,it),7) +1 !tracer_ts  !age
     traced(tracpo(1,it),tracpo(2,it),10) = it
     it = it + 1
!make separate routine
!     if (cp_field(ix_new_round,ix_new_round) .eq. 0 & ! no or no other cp at this gp
!         .or. cp_field(ix_new_round,ix_new_round) .eq. traced(tracpo(1,it),tracpo(2,it),9) ) then
!       cp_field(ix_new_round,ix_new_round) = traced(tracpo(1,it),tracpo(2,it),9)
!     else if ((cp_field(ix_new_round,ix_new_round) .le. -2) then  !already two cps or more
!       cp_field(ix_new_round,ix_new_round) = cp_field(ix_new_round,ix_new_round) -1
!     else ! already one cp 
!       cp_field(ix_new_round,ix_new_round) = -2
!     end if
   ENDDO
  RETURN
END SUBROUTINE update_tracer

! --------------------------------------------------------------------------------------
! CALCULATE ANGLE AND RADIUS FROM COG OF EACH TRACER 
!-----------------------------------------------------------------------------------
SUBROUTINE geometry(traced,COMx,COMy,already_tracked,max_no_of_cells)
USE cp_parameters, ONLY : dsize_x, dsize_y, max_tracer_CP

  INTEGER, INTENT(IN)       :: max_no_of_cells
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: DELTAx(max_no_of_cells,max_tracer_CP), &
                               DELTAy(max_no_of_cells,max_tracer_CP), pi
  INTEGER                   :: dx(max_no_of_cells), dy(max_no_of_cells) ! shift everything to have CP in thecenter
  INTEGER                   :: i,j
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  REAL, INTENT(IN)          :: COMx(max_no_of_cells),COMy(max_no_of_cells)
   pi = 2.*asin(1.)
  ! move COG into middle of center first 
  dx = INT(dsize_x)/2 - INT(COMx)
  dy = INT(dsize_y)/2 - INT(COMy)
  
  DO i = 1,max_no_of_cells ! loop trough every cp
!   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already tracerfor the CP
    do j = 1,max_tracer_CP,1
    !DO j=1,already_tracked(i)
      ! shift to center
      traced(i,j,1) = mod(traced(i,j,1)+dx(i)-1+dsize_x,float(dsize_x))+1
      traced(i,j,2) = mod(traced(i,j,2)+dy(i)-1+dsize_y,float(dsize_y))+1

      DELTAx(i,j) = traced(i,j,1) -dsize_x/2. !(COMx(i)+dx(i))-1.
      DELTAy(i,j) = traced(i,j,2) -dsize_y/2. !(COMy(i)+dy(i))-1 !dsize_y/2.
      traced(i,j,15) = DELTAx(i,j)
      traced(i,j,17) = DELTAy(i,j)

!      DELTAx(i,j) = traced(i,j,1) -COMx(i)-1. !dsize_x/2.
!      !mod(traced(i,j,1) -COMx(i)-1.+dsize_x,dsize_x)+1.
!      DELTAy(i,j) = traced(i,j,2) -COMy(i)-1. !dsize_y/2.
!      !mod(traced(i,j,2) - COMy(i)-1.+dsize_y,dsize_y)+1.

!      IF      (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .gt. 0) THEN
!       traced(i,j,8) = atan(DELTAx(i,j)/DELTAy(i,j))
!      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 2nd
!       traced(i,j,8) = atan(DELTAx(i,j)/DELTAy(i,j)) + pi
!      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 3nd
!       traced(i,j,8) = atan(DELTAx(i,j)/DELTAy(i,j)) + pi
!      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .gt. 0) THEN  ! 4nd
!       traced(i,j,8) = atan(DELTAx(i,j)/DELTAy(i,j)) +2.* pi
!      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .eq. 0)THEN
!       traced(i,j,8) = pi/2.
!      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .eq. 0) THEN
!       traced(i,j,8) = (3./2.) * pi
!      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .gt. 0) THEN
!       traced(i,j,8) = 0.
!      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .lt. 0) THEN
!       traced(i,j,8) =  pi
!      END IF
      IF      (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,8) = atan(DELTAy(i,j)/DELTAx(i,j))

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .gt. 0) THEN  ! 2nd
       traced(i,j,8) = atan(abs(DELTAx(i,j))/DELTAy(i,j)) + 1./2.*pi

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 3nd
       traced(i,j,8) = atan(DELTAy(i,j)/DELTAx(i,j)) + pi

      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 4nd
       traced(i,j,8) = atan(abs(DELTAx(i,j))/abs(DELTAy(i,j))) +3./2.* pi

      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .eq. 0)THEN
       traced(i,j,8) = 0.

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .eq. 0) THEN
       traced(i,j,8) =  pi

      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,8) = 1./2. *pi
      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .lt. 0) THEN
       traced(i,j,8) =  3./2.* pi
      END IF
      traced(i,j,5) = sqrt(DELTAy(i,j)**2.+DELTAx(i,j)**2.)
      ! shift back 
      traced(i,j,1) = mod(traced(i,j,1)-dx(i)-1+float(dsize_x),float(dsize_x))+1
      traced(i,j,2) = mod(traced(i,j,2)-dy(i)-1+float(dsize_y),float(dsize_y))+1
     END DO
!    END IF
   END DO
END SUBROUTINE geometry

! --------------------------------------------------------------------------------------
! Calculate radial velocities 
!-----------------------------------------------------------------------------------
SUBROUTINE radvel(traced, already_tracked,max_no_of_cells)
USE cp_parameters, ONLY : max_tracer_CP
  INTEGER, INTENT(IN)       :: max_no_of_cells
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: pi
  INTEGER                   :: i,j
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)

  pi = 2.*asin(1.)
  DO i = 1,max_no_of_cells ! loop trough every cp
!   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already
!   tracerfor the CP
    do j = 1,max_tracer_CP,1
      traced(i,j,19)   = traced(i,j,13)*cos(traced(i,j,8)) + traced(i,j,14) * sin(traced(i,j,8))
      traced(i,j,20)   = traced(i,j,14)*cos(traced(i,j,8)) - traced(i,j,13) * sin(traced(i,j,8))
    END DO
!    END IF
  END DO
 
END SUBROUTINE radvel

! --------------------------------------------------------------------------------------
! Update tracer along radial direction
!-----------------------------------------------------------------------------------
SUBROUTINE radial_update(timestep,traced,count_tracer,max_no_of_cells,tracpo,max_tracers)
USE cp_parameters, ONLY : dsize_x, dsize_y, res, dt, max_tracer_CP
  INTEGER, INTENT(IN)       :: timestep,max_no_of_cells, max_tracers
  INTEGER, INTENT(IN)       :: count_tracer !,max_tracer_CP
  !INTEGER, INTENT(IN)       :: cp_field(dsize_x,dsize_y)
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: ix_new, iy_new, vx_intp, vy_intp
  REAL                      :: ix, iy, ixt, iyt    !t are half level (interface)
!  INTEGER                   :: ix_round, iy_round, &
  INTEGER                   :: ix_new_round, iy_new_round, start_time
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)

  it = 1 ! counter trough traced
  DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until here
    ! determining the first time step of the event
     start_time=INT(traced(tracpo(1,it),tracpo(2,it),6))
     ! determining how many timesteps have passed since then
     tracer_ts =timestep-start_time+1

     ! getting the previous positions
     ix=traced(tracpo(1,it),tracpo(2,it),1)
     iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints
     !vx_intp = vr * cos(alpha)
     vx_intp = traced(tracpo(1,it),tracpo(2,it),19) * cos(traced(tracpo(1,it),tracpo(2,it),8))
     vy_intp = traced(tracpo(1,it),tracpo(2,it),19) * sin(traced(tracpo(1,it),tracpo(2,it),8))
     ix_new = MOD((ix + dt*vx_intp/res)-1.+FLOAT(dsize_x),FLOAT(dsize_x))+1.
     iy_new = MOD((iy + dt*vy_intp/res)-1.+FLOAT(dsize_y),FLOAT(dsize_y))+1.

     ix_new_round = nint(ix_new) ! maybe mod required
     iy_new_round = nint(iy_new)
     ! save new values for next loop
     traced(tracpo(1,it),tracpo(2,it),1) = ix_new
     traced(tracpo(1,it),tracpo(2,it),2) = iy_new
     traced(tracpo(1,it),tracpo(2,it),3) = float(ix_new_round)
     traced(tracpo(1,it),tracpo(2,it),4) = float(iy_new_round)

     traced(tracpo(1,it),tracpo(2,it),6) = start_time !timestep
     traced(tracpo(1,it),tracpo(2,it),7) = traced(tracpo(1,it),tracpo(2,it),7) +1 !tracer_ts  !age
     traced(tracpo(1,it),tracpo(2,it),10) = it
     it = it + 1
   ENDDO
  RETURN
END SUBROUTINE radial_update
! --------------------------------------------------------------------------------------
! SORT traced 
!-----------------------------------------------------------------------------------
SUBROUTINE sort(traced_dummy,jc)

  INTEGER, INTENT(IN)      :: jc
  REAL, INTENT(INOUT)      :: traced_dummy(jc,20)
  INTEGER                  :: j,k
  REAL                     :: a(20)
 ! sort by angle for every cpID:
  ! initialize
    do j=1,jc !   count_tracer!  
      a=traced_dummy(j,:)  !save value at j 
      do k=j-1,1,-1 ! go backwarts from value below current
        ! if  val before current value is smaller than actual value
        if (traced_dummy(k,8)<= a(8)) goto 10
        ! set this value at the current position 
        traced_dummy(k+1,:)=traced_dummy(k,:)

        ! DELTA X abspeichern
      end do
      k=0
      10 traced_dummy(k+1,:)=a
    end do
! write output sorted by angle for every cold pool
 j=1
  ! updating previous tracers
  DO WHILE (j .LT. jc)

    WRITE(50,200) INT(traced_dummy(j,6)), INT(traced_dummy(j,9)), & !timestep,CP ID
                  traced_dummy(j,1), traced_dummy(j,2),traced_dummy(j,3), &!position
                  traced_dummy(j,4), traced_dummy(j,5), traced_dummy(j,10)!, &
    j = j+1
  end do
 200 FORMAT   (2(2X,I4) ,   3(2X,F11.5), 2X,F5.3, 2(2X,F11.3))

END SUBROUTINE sort

! --------------------------------------------------------------------------------------
! STOP CP when all tracer are on the same side of the COG  
!-----------------------------------------------------------------------------------
SUBROUTINE oneside(tdummy,traced,tracked_no)

USE  cp_parameters, ONLY :max_tracer_CP
  INTEGER, INTENT(IN)       :: tracked_no
  REAL, INTENT(INOUT)       :: traced(max_tracer_CP,20)
  REAL, INTENT(IN)          :: tdummy(tracked_no)
  INTEGER                   :: j
  REAL                      :: dif 
  REAL                      :: pi
   pi = 2.*asin(1.)

!   DO i = 1,max_no_of_cells ! loop trough every cp
     IF (tracked_no .gt. 0) THEN ! do only sth if there are alreadytracerfor the CP
       DO j=1,tracked_no-1
         dif = tdummy(j+1) - tdummy(j)
         if (dif .lt. 0) then 
           write(*,*) "sorted not correctly"
         else if(dif .gt. 2./3.*pi) Then
            !write(*,*) "no tracers on both sides of CP, stop CP"
            !write(*,*) tdummy(j), tdummy(j+1), dif
            traced(:,11) = 0
         end if
       END DO
       dif = 2.*pi -tdummy(tracked_no) + tdummy(1)
       if (traced(1,11) .eq. 1) then
       if (dif .gt. 2./3.* pi ) then
         !write(*,*) "no tracers on both sides of CP, stop CP"
         !write(*,*) tdummy(1), tdummy(tracked_no)
         traced(:,11) = 0
       end if
       end if
       
     END IF
!   END DO

END SUBROUTINE oneside

! --------------------------------------------------------------------------------------
! STOP tarcer, when they accelerate after decelerating previously 
!-----------------------------------------------------------------------------------
SUBROUTINE time_dev(traced,traced_prev,max_no_of_cells,count_tracer,tracpo,max_tracers)
USE  cp_parameters, ONLY :max_tracer_CP
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells, max_tracer_CP,20)
  REAL, INTENT(IN)          :: traced_prev(max_no_of_cells, max_tracer_CP,20)
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  INTEGER, INTENT(IN)       :: max_no_of_cells, max_tracers,count_tracer !,max_tracer_CP
  REAL                      :: v0, v1, dv
  INTEGER                   :: it 
  it =1
    DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until
     v0=sqrt(traced_prev(tracpo(1,it),tracpo(2,it),14)**2 + &
             traced_prev(tracpo(1,it),tracpo(2,it),14)**2) 
     v1=sqrt(traced(tracpo(1,it),tracpo(2,it),14)**2 &
            +traced(tracpo(1,it),tracpo(2,it),14)**2)
     dv = v0-v1
     if (dv .gt. 1.) THEN !decelerate
      traced(tracpo(1,it),tracpo(2,it),11) = 1
     !if accelerated and decelerated already before 
     else if (dv .lt. 1. .and. traced(tracpo(1,it),tracpo(2,it),11) .eq. 1) THEN
       traced(tracpo(1,it),tracpo(2,it),11) = 0  ! set tracer inactive 
     end if
     it = it+1
   END DO 
END SUBROUTINE

SUBROUTINE write_output(traced,max_tracers,count_tracer,timestep,tracpo,&
                       max_no_of_cells,COMx,COMy)
USE  cp_parameters, ONLY :max_tracer_CP, max_age
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers, timestep, &
                               max_no_of_cells !, max_tracer_CP 
  INTEGER                   :: it
  INTEGER,INTENT(IN)        :: tracpo(2,max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells, max_tracer_CP,20)
  REAL, INTENT(IN)          :: COMx(max_no_of_cells), COMy(max_no_of_cells)
  it=1

  150 FORMAT    (2(4X,I4),    & !timestep, age
                2X,I6, 4X,I4, & !tracer and CP ID
                2(2X,F10.5),  & !pos 1-2
                2(4X,I4),     & !rounded
                2(2X,F7.3),   & !distance and angle
                2(2X,F8.2),   & !velocity                
                2(2X,F8.2),   & !x dist
                2(2X,F8.2),   & !COG
                1(2X,I1),     & ! merger
                2(2X,I4))       ! precip ID
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)
    IF (traced(tracpo(1,it),tracpo(2,it),11)  .eq. 1.) THEN  !trace only if tracer is active
          IF(INT(traced(tracpo(1,it),tracpo(2,it),7)) .le. max_age)then   ! output only up to  3hours
            WRITE(40,150) INT(timestep),INT(traced(tracpo(1,it),tracpo(2,it),7)),& !timestep, age
                        it,INT(traced(tracpo(1,it),tracpo(2,it),12)),& !tracer and CP ID
                        traced(tracpo(1,it),tracpo(2,it),1),traced(tracpo(1,it),tracpo(2,it),2),& !pos 1-2
                        INT(traced(tracpo(1,it),tracpo(2,it),3)),INT(traced(tracpo(1,it),tracpo(2,it),4)),& !rounded
                        traced(tracpo(1,it),tracpo(2,it),5),traced(tracpo(1,it),tracpo(2,it),8),& !distance and angle
                        traced(tracpo(1,it),tracpo(2,it),13),traced(tracpo(1,it),tracpo(2,it),14),& ! u, v Wind component
                        traced(tracpo(1,it),tracpo(2,it),15),traced(tracpo(1,it),tracpo(2,it),17),& ! x and y distance  
                        COMx(tracpo(1,it)), COMy(tracpo(1,it))                                   ,& ! center
                        INT(traced(tracpo(1,it),tracpo(2,it),16)), INT(traced(tracpo(1,it),tracpo(2,it),9)),&   ! merger, prec ID
                        INT(traced(tracpo(1,it),tracpo(2,it),18))   !timesteps set after precip start
    END IF
          END IF
    it = it+1
  END DO

END SUBROUTINE write_output

END PROGRAM cp_tracking
 
