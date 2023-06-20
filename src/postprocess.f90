
subroutine write_data
  use class_general
  use class_parameters
  implicit none

  include "silo_f9x.inc"
  real, allocatable             :: x(:),y(:),z(:)
  integer, allocatable          :: nodelist(:)
  integer                       :: i,j,k
  integer                       :: dbfile, err, ierr, ndims, NSHAPETYPES
  parameter (NSHAPETYPES = 1) ! Only hexa meshes
  integer                       :: LNODELIST
  integer shapesize(NSHAPETYPES) /8/
  integer, dimension (1:NSHAPETYPES):: shapecounts
  character(len = 40)               :: name,name2
  character(len=6)                  :: fmt ! format descriptor
  character(len=5)                  :: x1 ! format descriptor
  logical                           :: exist
  real dtime
  integer optlistid
  
! Fields

! ========================== ALLOCATE AND FILL ARRAYS FOR SILO =====================

  k=NNODE
  allocate(x(1:NNODE))
  allocate(y(1:NNODE))
  allocate(z(1:NNODE))
  allocate(nodelist(1:NMESH*8))


  do i=1,NNODE
     x(i)=NODES_ARRAY(i)%x
     y(i)=NODES_ARRAY(i)%y
     z(i)=NODES_ARRAY(i)%z
  enddo



! Fill arrays
  do i=1,NMESH
    do j=1,8
        nodelist(8*(i-1)+j)=MESHES_ARRAY(i)%nodes(j)
    enddo
  enddo
! ====================== ARRAYS ARE READY TO BE PRINTED =============================

  dtime=STEP*2.*DT-DT
  if (STEP == 0) dtime=0.0
  print*,"==========================================================================="
  print*,"Printing data at time = ", real(dtime,4), "|     Time step = ", STEP
  print*,"==========================================================================="


! Naming the outputs
  fmt = '(I5.5)' ! an integer of width 5 with zeros at the left
  write (x1,fmt) STEP ! converting integer to string using a 'internal file'
  name='plots//step'//trim(x1)
  name2='plots/step'//trim(x1)
  shapecounts(1)=NMESH
  LNODELIST=NMESH*8
  ndims = 3
  name= trim(name)
  ierr = dbcreate(trim(name), LEN_TRIM(name), DB_CLOBBER, DB_LOCAL, &
                  "Navier Stokes Equation Solution", 31, DB_PDB, dbfile)
!  print*, ierr, dbfile
  if(dbfile.eq.-1) then
    write (6,*) 'Could not create Silo file!\n'
  endif

  dtime=real(dtime,4)
  err  = dbmkoptlist(2, optlistid)
  err = dbaddiopt(optlistid, DBOPT_CYCLE, STEP)
  if(kind(dtime) == 4) then
     err = dbadddopt(optlistid, DBOPT_TIME, dtime)
  else
     err = dbadddopt(optlistid, DBOPT_DTIME, dtime)
  endif
  err = dbsetdepwarn(0) 


  err = dbputzl(dbfile, "zonelist", 8, NMESH, ndims, nodelist,&
        LNODELIST, 1, shapesize, shapecounts, NSHAPETYPES, ierr)

  err = dbputum(dbfile, "mesh", 4, ndims, REAL(x,4), REAL(y,4), REAL(z,4), &
        "X", 1, "Y", 1, "Z", 1, DB_FLOAT, NNODE, NMESH, &
        "zonelist", 8, DB_F77NULL, 0, optlistid, ierr)

 
  err = dbputuv1(dbfile, "u", 1, "mesh", 4, REAL(meshes_uxyzt,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
  err = dbputuv1(dbfile, "v", 1, "mesh", 4, REAL(meshes_vxyzt,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
  err = dbputuv1(dbfile, "w", 1, "mesh", 4, REAL(meshes_wxyzt,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
  err = dbputuv1(dbfile, "p", 1, "mesh", 4, REAL(meshes_p,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
  err = dbputuv1(dbfile, "T", 1, "mesh", 4, REAL(meshes_Txyzt,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
  err = dbputuv1(dbfile, "D", 1, "mesh", 4, REAL(meshes_dila,4), NMESH, &
        DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)

  err = dbputdefvars(dbfile, "defvars", 7, 1, &
        "vel", 3, DB_VARTYPE_VECTOR, "{u,v,w}", 7, DB_F77NULL, ierr)


  err = dbclose(dbfile)

  inquire(file="plots.visit", exist=exist)
  if (exist) then
    open(12, file="plots.visit", status="old", position="append", action="write")
  else
    open(12, file="plots.visit", status="new", action="write")
  end if
  write(12, '(a)') trim(name2)
  close(12)

  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(nodelist)
end subroutine write_data

