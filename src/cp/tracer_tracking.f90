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
INTEGER, ALLOCATABLE :: map(:,:,:)
INTEGER,ALLOCATABLE  :: active_CP(:)
INTEGER,ALLOCATABLE  :: CPsize(:,:)
INTEGER, ALLOCATABLE :: prec_active(:,:)
INTEGER, ALLOCATABLE :: CPsets(:,:,:)
INTEGER, ALLOCATABLE :: counter(:), grdpnts(:,:)
integer :: j
!INITIALIZE some values
namelist /INPUTgeneral/ dsize_x, dsize_y, dt, res
namelist /INPUTtracer/ max_tracer_CP, nTrSet, ntset, max_age, rad, lformat
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
ALLOCATE(CPsize(max_no_of_cells,max_age+1))
ALLOCATE(prec_active(max_no_of_cells,2))
ALLOCATE(rmax(max_no_of_cells))
ALLOCATE(active_CP(max_no_of_cells))
ALLOCATE(vel(dsize_x,dsize_y,2))
ALLOCATE(track_numbers(dsize_x,dsize_y))
ALLOCATE(nneighb(dsize_x,dsize_y))
ALLOCATE(already_tracked(max_no_of_cells))
ALLOCATE(input_field(dsize_x,dsize_y))
ALLOCATE(tracpo(2,max_tracers))
ALLOCATE(map(dsize_x,dsize_y,2))
ALLOCATE(grdpnts(2,dsize_x*dsize_y))
ALLOCATE(counter(max_no_of_cells))
ALLOCATE(CPsets(max_no_of_cells,nTrSet,2))
! initialize arrays
map(:,:,:) = 0
active_cp(:) = 0
write(*,*) 'start random fnc'
CALL randomgrid(grdpnts)
write(*,*) 'finished random fnc'
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
OPEN(40,FILE=trim(odir) // '/output/cp/coldpool_tracer_out_all.txt',FORM='formatted', ACTION='write')
OPEN(41,FILE=trim(odir) //'/output/cp/coldpool_tracer_out.txt',FORM='formatted', ACTION='write')

if (trim(lformat) == 'srv') then
  OPEN(4,FILE=trim(odir) // '/input/cp/input_u.srv',FORM='unformatted',ACTION='read')
  OPEN(5,FILE=trim(odir) // '/input/cp/input_v.srv',FORM='unformatted',ACTION='read')
end if

write(*,*) i, 'rain cells found'

!INITIALIZE some values
already_tracked(:) = 0
traced(:,:,:) = 0.
traced_prev(:,:,:) = 0.
count_tracer = 0 ! counts individual pixels !OCH was ist mit pixeln gemeint? der
!IDstart(:) = 0
timestep=0 
prec_active(:,2) = 1
prec_active(:,1) = 0
! read when first precip is tracked
 OPEN(3,FILE=trim(odir) // '/input/cp/tracks_body_time_sorted.txt',FORM='formatted',ACTION='read')
 153 FORMAT (7X,I5,9X,I3,94X,F10.4,7X,F10.4)
 onset = minval(cpio(:,3),1)  !first tracers will be set when first raincell is larger than agiven trashold
 write(*,*) 'tracer start at' , onset
 ! read first line of precip tracking, to get first timestep of precip
 READ(3,153,END=200)  ID, tts, xCOG, yCOG
 prec_active(ID,1) = 1
 begin = tts

 DO !start main loop
! runs until all timesteps are read, or uncomment lines below
! IF (timestep .ge. 150) THEN !maxval(cpio(:,3),1)) then
! GOTO 5000
! END IF 
! only for testing
   timestep=timestep+1
  write(*,*) 'at timestep', timestep
 
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
       prec_active(ID,1) = 1
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

! start when first precip appears
     IF (timestep .GE. onset) THEN 
       ! identification of edges, and set tracer
       nneighb(:,:)=0
       CALL neigh2(track_numbers,max_no_of_cells,grdpnts,CPsets,counter) 
      ! this routine might be rewriten that it can be called in set_tracer
       CALL set_tracer2(max_no_of_cells,cpio,traced,count_tracer,counter,CPsets,&
                        timestep,active_CP,tracpo,max_tracers,already_tracked)
       !CALL neigh(track_numbers,nneighb)
       !CALL set_tracer(nneighb, track_numbers,max_no_of_cells, timestep,traced,&
       !            count_tracer,already_tracked,tracpo,max_tracers,cpio(:,1:2))
       ! alternative to set tracers at edges of precip, set tracer at circle
       ! around precip
       !CALL maxcell(track_numbers,COMx,COMy,rmax,max_no_of_cells)
       !CALL initCircle(max_no_of_cells,timestep,traced, COMx,COMy,&
       !                rmax,cpio,max_tracers,tracpo,count_tracer,active_CP,already_tracked)
       CALL velocity_interpol(vel(:,:,1),vel(:,:,2),timestep,traced, count_tracer, &
                          max_no_of_cells,tracpo,max_tracers)
       CALL geometry(traced,COMx,COMy,already_tracked,max_no_of_cells)
       CALL radvel(traced,already_tracked,max_no_of_cells)
       ! new circular tracer should imeditly be updated untill they match with
       ! previous tracer, sont have neg rad vel
        CALL update_newtracer(vel(:,:,1),vel(:,:,2),timestep,traced, count_tracer,&
                          max_no_of_cells,tracpo,max_tracers)
  
       DO i =1,max_no_of_cells ! loop trough all cps with tracer
         IF (traced(i,1,11) .eq. 1) then !already_tracked(i) .gt. 1 ) THEN
!  !       ! reset dummy first
           traced_dummy(:,:) = 0.
           traced_dummy = traced(i,:,:)
           CALL sort(traced_dummy(1:already_tracked(i),:),already_tracked(i))
  !!         if (traced(i,1,11)  .eq. 0) then !stop tracer only if precip has stoped
  !           !CALL oneside(traced_dummy(1:already_tracked(i),8),traced(i,:,:),already_tracked(i))
  !!         end if
           IF (traced(i,1,7) .lt. max_age) then
             CALL setongrid(traced_dummy,map(:,:,1),map(:,:,2),i,COMx(i),COMy(i),CPsize(i,int(traced(i,1,7)+1)))
!             IF (prec_active(ID,1) .eq. 0 .and. traced(i,1,7) .gt. 2) then ! precipitation has terminates 
!             IF ( traced(i,1,7) .gt. 2) then ! precipitation has terminates 
!               IF (CPsize(i,int(traced(i,1,7)+1)) .lt. 0.8*maxval(CPsize(i,:)))  then 
!                 traced(i,:,11) = 0  ! set cold pool inactive
!                 write(*,*) 'set ', i,'inactive'
!               END IF
!             END IF 
           END IF
         !ELSE
         !  write(*,*) i ,'is inactive' 
         END IF
       END DO
       CAll WRITEMAP(map(:,:,1),timestep,trim(odir) // '/output/cp/object')
       CAll WRITEMAP(map(:,:,2),timestep,trim(odir) // '/output/cp/tracer')

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
 map(:,:,:) = 0
 prec_active(:,2) = prec_active(:,1)
 prec_active(:,1) = 0 
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
    CALL check(nf90_inq_varid(ncid,varname, rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(4), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nz))
    CALL check(nf90_get_var(ncid, rhVarId, poutput(1:nx,1:ny), start =(/1,1,1,ctime/),&
                                                    count= (/nx, ny,1, 1/)))
    CALL check(nf90_close(ncid))
!write(*,*) poutput(1:3,1:3)
!
  END SUBROUTINE read_nc_2D

  SUBROUTINE WRITEMAP(map,t,FILE_NAME2)
  USE cp_parameters, ONLY : dsize_x, dsize_y
  
  integer, intent(in) :: map(dsize_x, dsize_y)
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

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!---------------------------------------------------
! CHECK NETCD DATA
!--------------------------------------------------
  SUBROUTINE check(status)

    integer, intent (in) :: status

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

SUBROUTINE randomgrid(grdpnts)
  USE cp_parameters, ONLY : dsize_x, dsize_y
  INTEGER, INTENT(OUT) :: grdpnts(2,dsize_x * dsize_y)
  REAL :: randy 
  INTEGER :: pos, i, j ,c,k  
  INTEGER :: numbers(dsize_x*dsize_y), &
             numbers2(dsize_x*dsize_y)

  grdpnts(:,:) = 0
  c = dsize_x*dsize_y
  numbers = (/(k, k = 1,c )/)

  DO k = 1, dsize_x*dsize_y
   CALL RANDOM_NUMBER(randy)
       pos = ceiling(randy*c)
       if (pos .eq. c ) then
         numbers2 = numbers
       else
         numbers2(1:pos-1) = numbers(1:pos-1)
         numbers2(pos:c-1) = numbers(pos+1:c)
         numbers2(c) = numbers(pos)
       end if
       numbers = numbers2
       c = c-1
  END DO
  c = 0
  DO i =1,dsize_x
    DO j =1,dsize_y
       c = c +1
       pos = numbers(c) 
       grdpnts(:,pos) = (/i,j/)
    END DO
  END DO
END SUBROUTINE randomgrid

SUBROUTINE neigh2(track_numbers,max_no_of_cells,grdpnts,CPsets, counter)
  USE cp_parameters, ONLY : dsize_x, dsize_y, nTrSet
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: max_no_of_cells
  INTEGER, INTENT(OUT)   :: CPsets(max_no_of_cells,nTrSet,2)
  INTEGER, INTENT(OUT)   :: counter(max_no_of_cells)
  INTEGER                :: i,j,ii, imodm, imodp, jmodm,jmodp
  INTEGER, INTENT(IN)    :: grdpnts(2,dsize_x*dsize_y)
  REAL, INTENT(IN)       :: track_numbers(dsize_x,dsize_y)
  CPsets(:,:,:) = 0
  counter(:) = 0
  DO ii = 1,dsize_x*dsize_y
   i = grdpnts(1,ii)
   j = grdpnts(2,ii)
      IF (track_numbers(i,j) .gt. 0) THEN !check only gps with rain cell 
        imodp = mod(i+dsize_x,dsize_x)+1
        imodm = mod(i-2+dsize_x,dsize_x)+1
        jmodp = mod(j+dsize_y,dsize_y)+1
        jmodm = mod(j-2+dsize_y,dsize_y)+1
        IF((track_numbers(i,j) .ne. track_numbers(imodp,j)) .or. & 
           (track_numbers(i,j) .ne. track_numbers(imodm,j)) .or. &
           (track_numbers(i,j) .ne. track_numbers(i,jmodp)) .or. & 
           (track_numbers(i,j) .ne. track_numbers(i,jmodm)) .or. &
        ! diagonal neigh
           (track_numbers(i,j) .ne. track_numbers(imodp,jmodp)) .or. &
           (track_numbers(i,j) .ne. track_numbers(imodm,jmodp)) .or. &
           (track_numbers(i,j) .ne. track_numbers(imodp,jmodm)) .or. &
           (track_numbers(i,j) .ne. track_numbers(imodm,jmodm))) THEN 
             counter(int(track_numbers(i,j))) = counter(int(track_numbers(i,j))) +1
             if (counter(int(track_numbers(i,j))) .le. nTrSet) then 
               CPsets(int(track_numbers(i,j)),counter(int(track_numbers(i,j))),:) = (/i,j/)
             end if
         END IF
       END IF
  END DO
  write(*,*) 'counter 72',counter(72)
END SUBROUTINE neigh2

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
                      rmax,cpio,max_tracers,tracpo,count_tracer,active_CP,already_tracked)
   USE cp_parameters, ONLY : dsize_x, dsize_y,max_tracer_CP, nTrSet, ntset
  
  INTEGER, INTENT(IN)       :: max_no_of_cells, max_tracers, ts
!  INTEGER, INTENT(IN)       :: IDstart(max_no_of_cells)
  REAL, INTENT(IN)          :: rmax(max_no_of_cells),&
                               COMx(max_no_of_cells), COMy(max_no_of_cells)
  INTEGER, INTENT(INOUT) ::already_tracked(max_no_of_cells)
  INTEGER, INTENT(IN)       :: cpio(max_no_of_cells,3)
  INTEGER, INTENT(INOUT)    :: count_tracer, active_CP(max_no_of_cells)
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
     active_CP(i) = 1
     already_tracked(i) =(ts-cpio(i,3)+1) * nTrSet             

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
SUBROUTINE set_tracer2(max_no_of_cells,cpio,traced,count_tracer,counter,&
                       CPsets,ts,active_CP,tracpo,max_tracers,already_tracked)
USE cp_parameters, ONLY: nTrSet, max_tracer_CP

INTEGER, INTENT(IN) :: max_no_of_cells, ts,max_tracers
INTEGER, INTENT(IN) :: CPsets(max_no_of_cells,nTrSet,2)
INTEGER, INTENT(IN) :: counter(max_no_of_cells)
INTEGER, INTENT(INOUT) :: already_tracked(max_no_of_cells)
REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
INTEGER, INTENT(INOUT)    :: count_tracer, active_CP(max_no_of_cells)
INTEGER, INTENT(IN)       :: cpio(max_no_of_cells,3)
INTEGER, INTENT(INOUT)    :: tracpo(2,max_tracers)
REAL :: randy
INTEGER :: setn,plus,i,j,ii,cc
  DO i = 1,max_no_of_cells,1 
  if (i .eq. 72) write(*,*) 'for 72 at timestep',ts, cpio(72,:)
    IF (cpio(i,2) .ne. 0 ) then  ! no splittig event
      setn = ts-cpio(i,3) 
      if (i .eq. 72) write(*,*) 'in first loop ', setn, ntset
      plus = setn*nTrSet 
      cc = min(counter(i),nTrSet)
      IF (setn .ge. 0 .and. setn .le. ntset) then  ! in timerange of setting 
        already_tracked(i) = nTrSet * (setn+1) ! number of tracer which will have been set after this routine
        if (i .eq. 72) write(*,*) 'should be here', nTrSet, counter(i)
        IF (counter(i) .gt. 0) then
          active_CP(i) = 1
          ! set tracer at center of all outer gridpints
          count_tracer = count_tracer +cc
          traced(i,1+plus:cc+plus,1) = CPsets(i,1:cc,1) ! xpos
          traced(i,1+plus:cc+plus,2) = CPsets(i,1:cc,2) ! ypos
          traced(i,1+plus:cc+plus,3) = CPsets(i,1:cc,1) !xpos
          traced(i,1+plus:cc+plus,4) = CPsets(i,1:cc,2) !ypos
          IF (counter(i) .lt. nTrSet) then
            DO ii=cc+1,nTrSet
              count_tracer = count_tracer+1
              CALL RANDOM_NUMBER(randy)
              traced(i,ii+plus,1) = CPsets(i,mod(ii,cc)+1,1) + randy
              CALL RANDOM_NUMBER(randy)
              traced(i,ii+plus,2) = CPsets(i,mod(ii,cc)+1,2) + randy
              traced(i,ii+plus,3) = CPsets(i,mod(ii,cc)+1,1)
              traced(i,ii+plus,4) = CPsets(i,mod(ii,cc)+1,2)
            END DO
          END IF 
          traced(i,1+plus:nTrSet+plus,12)=cpio(i,2)  ! start
        ELSE   ! if precip does not exists, set tracer inbetween
          write(*,*) 'precip has stoped'
          count_tracer = count_tracer + nTrSet
!!!!! 
write(*,*) i, 'CP', 1+plus-1, setn, ts, cpio(i,3)
          traced_dummy = traced(i,1+(plus-1):nTrSet+(plus-1),:)
 !         CALL sort(traced_dummy(:,:),already_tracked(i)-plus)
          ! muss noch gemacht werten
          traced(i,1+plus:plus+nTrSet,1) = traced(i,1:nTrSet,1) ! xpos
          traced(i,1+plus:plus+nTrSet,2) = traced(i,1:nTrSet,2) ! xpos
          traced(i,1+plus:plus+nTrSet,3) = traced(i,1:nTrSet,3) ! xpos
          traced(i,1+plus:plus+nTrSet,4) = traced(i,1:nTrSet,4) ! xpos
          traced(i,1+plus:plus+nTrSet,12)= traced(i,1:nTrSet,12)! xpos
        END IF
        traced(i,1+plus:plus+nTrSet,11) = 1
        traced(i,1+plus:plus+nTrSet,6) = ts
        traced(i,1+plus:plus+nTrSet,(/7,18/)) =setn ! age, setn
        traced(i,1+plus:nTrSet+plus,9) = i
        traced(i,1+plus:plus+nTrSet,10) = (/(j, j =count_tracer+1-nTrSet,count_tracer )/)
        tracpo(1,count_tracer+1-nTrSet:count_tracer) = i
        tracpo(2,count_tracer+1-nTrSet:count_tracer) = (/(j, j =1+plus,plus+nTrSet )/)
      END IF !timestep within the timerange when tracers shell be set for thisCP
    END IF  ! precip event is no splitting
  END DO
END SUBROUTINE set_tracer2

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
     !write(*,*) it,count_tracer
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


SUBROUTINE update_newtracer(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers)
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
   IF (traced(tracpo(1,it),tracpo(2,it),18) .eq.traced(tracpo(1,it),tracpo(2,it),7) )THEN
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
     ix_new_round = mod(nint(ix_new)-1+dsize_x,dsize_x)+1 ! maybe mod required
     iy_new_round = mod(nint(iy_new)-1+dsize_y,dsize_y)+1 

     ! save new values for next loop
     traced(tracpo(1,it),tracpo(2,it),1) = ix_new
     traced(tracpo(1,it),tracpo(2,it),2) = iy_new
     traced(tracpo(1,it),tracpo(2,it),3) = float(ix_new_round) 
     traced(tracpo(1,it),tracpo(2,it),4) = float(iy_new_round)
    END IF
    it = it+1
   ENDDO
  RETURN
END SUBROUTINE update_newtracer


SUBROUTINE update_tracer(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers)
!OCH TO DO: dt and resolution should be parameter read by USE from module
USE cp_parameters, ONLY : dsize_x, dsize_y, res, dt, max_tracer_CP, max_age

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
!   IF (traced(tracpo(1,it),tracpo(2,it),11) .eq. 1) then
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
     ix_new_round = mod(nint(ix_new)-1+dsize_x,dsize_x)+1 ! maybe mod required
     iy_new_round = mod(nint(iy_new)-1+dsize_y,dsize_y)+1 

     ! save new values for next loop
     traced(tracpo(1,it),tracpo(2,it),1) = ix_new
     traced(tracpo(1,it),tracpo(2,it),2) = iy_new
     traced(tracpo(1,it),tracpo(2,it),3) = float(ix_new_round) 
     traced(tracpo(1,it),tracpo(2,it),4) = float(iy_new_round)

     traced(tracpo(1,it),tracpo(2,it),6) = start_time !timestep
     traced(tracpo(1,it),tracpo(2,it),7) = traced(tracpo(1,it),tracpo(2,it),7) +1 !tracer_ts  !age
     if (traced(tracpo(1,it),tracpo(2,it),7) .gt. max_age) then 
      traced(i,:,11) = 0
     end if
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
!    end if
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
SUBROUTINE setongrid(traced_dummy,map,map2,CPID,cogx,cogy,counter2)
USE cp_parameters, ONLY : dsize_x, dsize_y, max_tracer_CP
  INTEGER,INTENT(INOUT) :: map(dsize_x, dsize_y),map2(dsize_x, dsize_y)
  INTEGER               :: maptr(dsize_x, dsize_y),maptr2(dsize_x, dsize_y),maptr3(dsize_x, dsize_y) 
  INTEGER, INTENT(IN)      ::  CPID
  REAL, INTENT(IN) :: traced_dummy(max_tracer_CP,20)
  REAL :: dum(max_tracer_CP,20), dumdum(max_tracer_CP,20)
  REAL, INTENT(IN) :: cogx, cogy
  INTEGER :: i,j,jsave,k, l
  INTEGER ::ii,ij,jj, deltax,deltay,dx, dy, ix, iy,oldx,c,oldy,nx,ny,x,y,xtemp,ytemp,xp,yp
  REAL :: pi,  s, centerx, centery
  INTEGER :: d, cc, mapmap(dsize_x, dsize_y), checker, counter
  INTEGER, INTENT(OUT) :: counter2
  REAL :: rsave
  REAL :: dist, distnew
  INTEGER ::nextracer, startc
  INTEGER :: savei
  INTEGER :: firstx, firsty

  write(*,*) 'CP', CPID
  checker = 0
  counter2 = 0
  counter=0
  pi = 2.*asin(1.)
  c =1
  cc = 1
  centerx = 0
  centery = 0
  mapmap = 0
!  ! to find a point inside the object
!  ! cog of precip is not always inside
  ! SHIFT MAP TO CENTER COG
  deltax = dsize_x/2-int(cogx)
  deltay = dsize_y/2-int(cogy)
  maptr(:,:) = 0
  maptr2(:,:) = 0
  maptr3(:,:) = 0

!  do ix=1,dsize_x
!   do iy=1,dsize_y
!    maptr(ix,iy) =map(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1) 
!    maptr2(ix,iy) =map2(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1)
!   end do
!  end do    

  rsave = 0.
  DO i=1,max_tracer_CP
     IF (traced_dummy(i,7) .gt. traced_dummy(i,18)) THEN
       dumdum(c,:) = traced_dummy(i,:)  
       if (traced_dummy(i,5) .gt. rsave) then
         startc = c
         rsave = traced_dummy(i,5)
       end if

       c= c+1
     END IF
  END DO 
  c = c- 1
  dum(:,:) = dumdum(:,:)
  i =1
  savei = 0
  if (c .gt. 1) then ! at first setting no tracer for dum will be found
!    dum(1:c-startc+1,:) = dumdum(startc:c,:)
!    dum(c-startc+2:c,:) = dumdum(1:startc-1,:)
!    dum(c-startc+2:,8) = dumdum(c-startc+2:,8)+2*pi

    oldx = mod(int(dum(c-1,3))-1+deltax+dsize_x,dsize_x)+1
    oldy = mod(int(dum(c-1,4))-1+deltay+dsize_y,dsize_y)+1
    firstx = oldx  ! gp to which gap is closed the first time
    firsty = oldy

    DO WHILE (i .lt. c+1) !max_tracer_CP
      if (i .lt. savei) then 
        ! close circle and get out 
        x = mod(int(dum(i,3))-1+deltax+dsize_x,dsize_x)+1
        y = mod(int(dum(i,4))-1+deltay+dsize_y,dsize_y)+1
!        ! fill  from the previous to the current,which has passed the "end of the
!        ! circle" 
!        CALL filltracer(oldx, oldy,x,y,maptr, maptr2,maptr3, CPID, &
!                     centerx, centery, counter, checker)
!        ! and fill from current to the first one
!        CALL filltracer(x,y,firstx, firsty, maptr, maptr2,maptr3, CPID, &
!                     centerx, centery, counter, checker)
        goto 1288
      end if
      savei = i
      !write(*,*)  'at tracer', i
      x = mod(int(dum(i,3))-1+deltax+dsize_x,dsize_x)+1
      y = mod(int(dum(i,4))-1+deltay+dsize_y,dsize_y)+1
!      IF (maptr(x,y) .eq. 0 .or. maptr(x,y) .eq. CPID) THEN
!        IF (int(dum(i,3)) .ne. 0 .and. int(dum(i,4)) .ne. 0) THEN
          IF (dum(i,7) .gt. dum(i,18)) THEN ! .or. dum(i,7) .lt. 2) THEN
            IF (maptr2(x,y) .ne. CPID) THEN
              centerx = (centerx * counter +x)/(counter+1)
              centery = (centery * counter +y)/(counter+1)
              counter = counter+1
            END IF

            if((abs(oldx-x) .lt. 1.1 ) .and. (abs(oldy-y) .le. 1.)) then
               maptr3(x,y) =i
               maptr2(x,y) =-2
               maptr(x,y) =CPID

               i = i +1
               jsave = 0             
            else
              jsave =1
              dist=sqrt(real(oldx-x)**2+real(oldy-y)**2)! max(abs(oldx-xtemp),abs(oldy-ytemp))  
              do j=1,int(c-i)/2
                ij = mod(i+j-1,c-1)+1
                xtemp = mod(int(dum(ij,3))-1+deltax+dsize_x,dsize_x)+1  
                ytemp = mod(int(dum(ij,4))-1+deltay+dsize_y,dsize_y)+1
                if((abs(oldx-xtemp) .lt. 1.1 ) .and. (abs(oldy-ytemp) .le. 1.)) then
                  i = ij
                  maptr3(xtemp,ytemp) = ij
                  maptr2(xtemp,ytemp) = -2
                  maptr(xtemp,ytemp) = CPID
                  GOTO 2912  ! END OF LOOP
                else
                  distnew = sqrt(real(oldx-xtemp)**2+real(oldy-ytemp)**2) !max(abs(oldx-xtemp),abs(oldy-ytemp)) !sqrt(real(oldx-xtemp)**2+real(oldy-ytemp)**2)
                  if (distnew .lt. dist) then
                     jsave = j
                     dist = distnew
                     x = xtemp !mod(int(dum(i,3))-1+deltax+dsize_x,dsize_x)+1
                     y = ytemp !mod(int(dum(i,4))-1+deltay+dsize_y,dsize_y)+1
                   end if 
                 end if
               end do
               ! set inbetwen points
               xtemp = oldx
               ytemp = oldy
              maptr3(x,y) = ij
              maptr2(x,y) = -2
              maptr(x,y) = CPID

              CALL filltracer(xtemp,ytemp,x,y,maptr, maptr2,maptr3, CPID, &
                      centerx, centery, counter, checker) 
             end if ! end neighb yes or no
           END IF       
!         END IF
!        END IF
        oldx = x
        oldy = y
        i = i+jsave
        2912 CONTINUE
      END DO
      CALL filltracer(x,y,firstx, firsty, maptr, maptr2,maptr3, CPID, &
                     centerx, centery, counter, checker)

      1288 CONTINUE
ii = centerx
ij = centery
!!  d = 1
  maptr3 = maptr2
  maptr2 = maptr

  d = 0

  CALL getobj(maptr3,CPID,maptr,ii,ij,d)

! if object exceed domain boundaries, all gp are filled and map needs to be set
! back to outlines
  if ( ALL( maptr==CPID ) ) then
    maptr(:,:) = maptr2
  end if

! not neccessary any more
!  if (ANY(maptr(1:2,:) .eq. CPID)) maptr(:,:) = maptr2
!  if (ANY(maptr(dsize_x-2:dsize_x,:) .eq. CPID)) maptr(:,:) = maptr2
!  if (ANY(maptr(:,1:2) .eq. CPID)) maptr(:,:) = maptr2
!  if (ANY(maptr(:,dsize_y-2:dsize_y) .eq. CPID)) maptr(:,:) = maptr2

! just to controll where centers are placed
  maptr2(centerx,centery) = -2
  maptr2(dsize_x/2,dsize_y/2) = -4

! transfer nback   ! why + deltax???
  do ix=1,dsize_x
   do iy=1,dsize_y
    if  ((maptr(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1) .ne.0) .or. &
       (maptr2(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1) .ne.0)) then
! CPID  &
!      .or. maptr(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1) .eq. -2) then
       !if (map(ix,iy) .eq. 0) then
         counter2 = counter2+1
      ! end if

       map(ix,iy)=maptr(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1)
       map2(ix,iy)=maptr2(mod(ix+deltax-1+dsize_x,dsize_x)+1,mod(iy+deltay-1+dsize_y,dsize_y)+1)

    end if
   end do
  end do
 end if 

END SUBROUTINE setongrid

!------------------
!
!------------------
SUBROUTINE filltracer(startx, starty,x,y,maptr, maptr2,maptr3, CPID, &
                      centerx, centery, counter, checker)
USE cp_parameters, ONLY : dsize_x, dsize_y, max_tracer_CP

INTEGER, INTENT(IN) :: startx, starty,x, y, CPID
INTEGER, INTENT(INOUT) :: maptr(dsize_x, dsize_y), maptr2(dsize_x, dsize_y)
INTEGER, INTENT(INOUT) :: maptr3(dsize_x, dsize_y)
REAL, INTENT(INOUT) :: centerx, centery
INTEGER, INTENT(INOUT) :: counter, checker
INTEGER :: xtemp, ytemp, dx, dy
  xtemp = startx
  ytemp = starty
  do while (xtemp .ne. x .or. ytemp .ne. y)
    dx = xtemp-x
    dy = ytemp-y
    if (abs(dx) .gt. abs(dy)) then
      xtemp = xtemp - 1*((dx/abs(dx)))
      ytemp = ytemp -(dy/dx)*(dx/abs(dx))
      IF (maptr2(xtemp,ytemp) .ne. CPID) THEN
        centerx = (centerx * counter +xtemp)/(counter+1)
        centery = (centery * counter +ytemp)/(counter+1)
        counter = counter +1
        checker = 1
      END IF
      maptr(xtemp,ytemp) = CPID
      maptr2(xtemp,ytemp) = -2
      if (maptr3(xtemp,ytemp) .eq. 0) then
        maptr3(xtemp,ytemp) = -2
      end if
      checker = 1
    else
      ytemp = ytemp - 1*((dy/abs(dy)))
      xtemp = xtemp - (dx/dy)*(dy/abs(dy))
      IF (maptr2(xtemp,ytemp) .ne. CPID) THEN
        centerx = (centerx * counter +xtemp)/(counter+1)
        centery = (centery * counter +ytemp)/(counter+1)
        counter = counter +1
        checker = 1
      END IF
      maptr(xtemp,ytemp) = CPID
      maptr2(xtemp,ytemp) = -2
      if (maptr3(xtemp,ytemp) .eq. 0 )then
        maptr3(xtemp,ytemp) = -2
      end if
      checker = 1
    end if
  end do
END SUBROUTINE filltracer

!----------------------------------------------
! identifies all gp within an outline starting 
! at a given center
!----------------------------------------------
RECURSIVE SUBROUTINE getobj(objectin,CPID,objectout,ii,ij,c)
USE cp_parameters, ONLY : dsize_x, dsize_y, max_tracer_CP
  INTEGER, INTENT(IN) :: objectin(dsize_x, dsize_y)
  INTEGER, INTENT(IN) :: CPID
  INTEGER, INTENT(IN) :: ii, ij
  INTEGER :: icell(4), jcell(4)  ! looking in x and y direction for object points
  INTEGER, INTENT(INOUT) :: objectout(dsize_x, dsize_y)
  INTEGER, INTENT(INOUT) :: c
  INTEGER :: z
  INTEGER :: ii_mod, ij_mod

  ! gridpoint shifts
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  ii_mod = mod(ii-1+dsize_x,dsize_x) + 1
  ij_mod = mod(ij-1+dsize_y,dsize_y)+ 1

! avoid to loop unneccessarily often when boundaries are hit, what 
! happens when center is outside of the outline

  if (ii .lt. 2) goto 200
  if (ij .lt. 2) goto 200
  if (ii .gt. dsize_x-2) goto 200
  if (ij .gt. dsize_y-2) goto 200

 
  IF (objectout(ii_mod,ij_mod) .ne. CPID) then
    objectout(ii_mod,ij_mod) = CPID
    IF (objectin(ii_mod,ij_mod) .ne. -2) then
      DO z = 1,4
!        if (c .gt. 10000)  goto 200
        c = c+1 

        CALL getobj(objectin,CPID,objectout,icell(z),jcell(z),c)
      END DO
    ELSE

! if boundary of cell was found, must go back to previous gp not to jump over
! the boundaries when proceed searching from a bundoury point
      IF (z .eq. 1) icell(1) = icell(1)-1
      IF (z .eq. 2)  jcell(2) = jcell(2) -1
      IF (z .eq. 3)  icell(3) = icell(3) +1
      IF (z .eq. 4)  jcell(4) = jcell(4) + 1
    END IF
  END IF
goto 100
! if domain boundary was hit, set all gps to ID, what terminates the routine
200 CONTINUE
c = 0 
objectout(:,:) = CPID
100 CONTINUE
END SUBROUTINE getobj
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
      !IF(INT(traced(tracpo(1,it),tracpo(2,it),7)) .le. max_age)then   ! output only up to  3hours
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

        IF (INT(traced(tracpo(1,it),tracpo(2,it),7)) .gt. INT(traced(tracpo(1,it),tracpo(2,it),18)) &
           .or. INT(traced(tracpo(1,it),tracpo(2,it),7)) .eq. 0 ) THEN
          WRITE(41,150) INT(timestep),INT(traced(tracpo(1,it),tracpo(2,it),7)),&!timestep, age
                        it,INT(traced(tracpo(1,it),tracpo(2,it),12)),& !tracerand CP ID
                        traced(tracpo(1,it),tracpo(2,it),1),traced(tracpo(1,it),tracpo(2,it),2),&!pos 1-2
                        INT(traced(tracpo(1,it),tracpo(2,it),3)),INT(traced(tracpo(1,it),tracpo(2,it),4)),&!rounded
                        traced(tracpo(1,it),tracpo(2,it),5),traced(tracpo(1,it),tracpo(2,it),8),&!distance and angle
                        traced(tracpo(1,it),tracpo(2,it),13),traced(tracpo(1,it),tracpo(2,it),14),&! u, v Wind component
                        traced(tracpo(1,it),tracpo(2,it),15),traced(tracpo(1,it),tracpo(2,it),17),&! x and y distance  
                        COMx(tracpo(1,it)), COMy(tracpo(1,it)) ,& ! center
                        INT(traced(tracpo(1,it),tracpo(2,it),16)),INT(traced(tracpo(1,it),tracpo(2,it),9)),&   ! merger, prec ID
                        INT(traced(tracpo(1,it),tracpo(2,it),18))   !timestepsset after precip start

        END IF           
    !   END IF
    END IF
    it = it+1
  END DO

END SUBROUTINE write_output

END PROGRAM cp_tracking
 
