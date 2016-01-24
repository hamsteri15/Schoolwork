! **********************************************************************


module CONSTANTS

  implicit none

  save
	!simulation control parameters and geometry specific info
	character :: CASE_NAME
	integer :: OUTPUT_N,n_room
	integer :: BOX_MOV,EXTRA,INCL_EXP !not used in RCM, read anyways...
	double precision :: IMPACT_FRIC,INFLOW_MOM,DROPPING_FLOW,PHI_CRIT,POS_CHANGE !not used in RCM, read anyways...

  double precision :: DT, T_MAX,RHO,GRAVITY,OUTPUT_INTERVAL,CFL_MAX
	double precision :: DIST_TIME,CROSS_SECTION
	double precision :: AMPLITUDE,FREQ
	double precision,allocatable,dimension(:) :: OUTPUT_TIMES

	!room info	
	double precision, allocatable :: x_room(:), y_room(:), z_room(:)
  double precision, allocatable :: l_room(:), b_room(:), h_room(:)
	double precision, allocatable :: INIT_MASS(:),DY(:)
	integer, allocatable :: JMAX(:)

	!opening info
	INTEGER :: n_opening
  INTEGER, ALLOCATABLE :: conn_open(:,:), type_open(:)
  DOUBLE PRECISION, ALLOCATABLE :: A_open(:), b_open(:), h_open(:)
  DOUBLE PRECISION, ALLOCATABLE :: r_open(:,:), norm_open(:,:), C_d(:)  

end module CONSTANTS

! **********************************************************************
