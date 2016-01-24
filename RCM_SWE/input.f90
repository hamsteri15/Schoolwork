MODULE INPUT
IMPLICIT NONE


CONTAINS

! **********************************************************************

subroutine READ_INPUTS()

  use CONSTANTS
	use WORKSPACE, ONLY: CUMULATIVE_MASS,TIME
	
	
  implicit none
	integer :: I

	!simulation info
	!OPEN(20,FILE='input_models/runparameters.txt') 
	OPEN(20,FILE='input/runparameters.txt') 
	READ(20,*), T_MAX
	READ(20,*), CFL_MAX
	READ(20,*), RHO
	READ(20,*), GRAVITY
	READ(20,*), BOX_MOV
	READ(20,*), AMPLITUDE
	READ(20,*), FREQ
	READ(20,*), EXTRA
	READ(20,*), POS_CHANGE
	READ(20,*), INFLOW_MOM
	READ(20,*), DROPPING_FLOW
	READ(20,*), IMPACT_FRIC
	READ(20,*), PHI_CRIT
	READ(20,*), OUTPUT_INTERVAL	
	READ(20,*), INCL_EXP
	READ(20,*), CASE_NAME
	CLOSE(20)

	!==!

	!room info
	OPEN(20,FILE="input/room_info.txt")
	!read number of rooms
	READ(20,*) n_room

	ALLOCATE( x_room(n_room),y_room(n_room),z_room(n_room) )
  ALLOCATE( l_room(n_room),b_room(n_room),h_room(n_room) )
	ALLOCATE( INIT_MASS(n_room),JMAX(n_room), DY(n_room) )

	 !read geometry specific info for the rooms
	 DO I = 1,n_room
      READ(20,*), x_room(I) ! room center x-coord in hull coordinates
      READ(20,*), y_room(I) ! room center y-coord in hull coordinates
      READ(20,*), z_room(I) ! room bottom z-coord in hull coordinates
      READ(20,*), l_room(I)
      READ(20,*), b_room(I)
      READ(20,*), h_room(I)
			READ(20,*), JMAX(I)
			DY(I) = b_room(I)/JMAX(I)
	ENDDO

	CLOSE(20)
 
	!==!

	OPEN(20,FILE='input/opening_info.txt') 
		 !get the amount of openings
		 READ(20,*) n_opening
		 write(*,'(/A,I3)') 'number of openings:', n_opening

     ALLOCATE( conn_open(n_opening,2), type_open(n_opening) )
     ALLOCATE( A_open(n_opening), b_open(n_opening), h_open(n_opening) )
     ALLOCATE( r_open(n_opening,3), norm_open(n_opening,3), C_d(n_opening) )

		 IF (n_opening.NE.0) THEN
          DO I = 1,2
             READ(20,*) 
          END DO
       END IF

			openings:DO I = 1,n_opening

          READ(20,*) conn_open(I,1:2), type_open(I), A_open(I), &
               b_open(I), h_open(I), r_open(I,:), norm_open(I,:), C_d(I)

        

          write(*,'(/3(A,I2,A,I2))') 'opening:', I, ', from:', conn_open(I,1), &
               ' to:', conn_open(I,2)
          write(*,'(5(A,G11.3))') 'type:', type_open(I), &
               ', Area:',A_open(I),'m^2, b:', b_open(I), 'm, h:', h_open(I), 'm, Cd:', C_d(I) 
          write(*,'(A,3G10.3)') 'location:',r_open(I,:)
          write(*,'(A,3G10.3)') 'normal:',norm_open(I,:)

					IF (type_open(I) .NE. 1) THEN
							PRINT *,"Do not use type 2 openings in RCM they are currently not applied, stopping program."
					ENDIF


          IF (n_room < conn_open(I,1) .OR. n_room < conn_open(I,2) ) THEN
             write(*,*) 'Entered conn_open to a non-existent room, stopping program'
             STOP
          END IF
       END DO openings



	CLOSE(20)

	!======!
 

	OUTPUT_N = nint(T_MAX/OUTPUT_INTERVAL)
	ALLOCATE (OUTPUT_TIMES(0:OUTPUT_N))
	TIME = 0.
  CUMULATIVE_MASS = 0.
 	

end subroutine READ_INPUTS


! **********************************************************************



subroutine INIT(JMAX,V_NOW,V_NEW,H_NOW,H_NEW,OUTPUT_TIMES,RHO,DY,n_room,&
		b_room,l_room,x_room,y_room,z_room,OUTPUT_INTERVAL,OUTPUT_N)
!initializes the water heights and velocities at the beginning of the simulation
!the driving source terms (A_Z and F_Y) are calculated in a separate routine
 
	USE ADDITIONAL, ONLY : function_surface_eq

  implicit none
	
	integer, intent(IN) :: n_room,OUTPUT_N
	double precision, intent(IN) :: RHO,OUTPUT_INTERVAL
	integer, intent(IN),dimension(1:n_room) :: JMAX
	double precision, intent(IN),dimension(1:n_room) :: b_room,l_room,x_room,y_room,z_room,DY
	!==!
	double precision, intent(OUT),dimension(0:OUTPUT_N) :: OUTPUT_TIMES
	double precision, dimension(1:n_room,0:MAXVAL(JMAX)+1),intent(OUT) :: V_NOW,V_NEW, H_NOW,H_NEW
	!==!
	integer :: I,J
	double precision :: H_INIT,DUMMY,dphi,phi,vol,mass,Y_COOR,temp
	double precision :: SURFACE(4)
	

	!open file and read the not used inputs to dummy variable
	OPEN(19,FILE="input/initstate.txt")
	DO I = 1,13
		READ(19,*) ,DUMMY
	ENDDO
	
	DO I = 1, n_room

		READ (19,*), dphi
		READ (19,*), phi
		READ (19,*), vol
		mass = RHO*vol 
		SURFACE = function_surface_eq(I , phi , mass)

		!temp = 0.
		DO J = 1, JMAX(I)

			!y-coordinates in the ship system
			Y_COOR = y_room(I)-(b_room(I)/2.)+(DY(I)/2.)+(DY(I)*J)

			!water heights in the _ROOM_ system, must be converted to positive values!
			H_NOW(I,J) = (-SURFACE(4)-(SURFACE(2)*Y_COOR))/(SURFACE(3))-z_room(I)		
		
			!take only positive water heights. NOTE! this causes a water loss when the water has
			!some initial angle. Has to be fixed later!!! use the variables mass and temp for comparison
			!of the input water mass and the one which is set to each room
			IF (H_NOW(I,J) .GT. 0.) THEN
				H_NOW(I,J) = 0.
			ENDIF
			H_NOW(I,J) = -H_NOW(I,J)

			!initial water velocity is set to zero
			V_NOW(I,J) = 0.
			
			!temp = temp+(RHO*(H_NOW(I,J)-z_room(I))*l_room(I)*DY(I))
		ENDDO		
		
	ENDDO
	CLOSE(19)
	H_NEW(:,:) = 0.; V_NEW(:,:) = 0.
		

	!initialize an array which consists the times at which the outputs are printed
	OUTPUT_TIMES(0) = 0.
	DO J = I,OUTPUT_N
		OUTPUT_TIMES(J) = OUTPUT_INTERVAL*J
	ENDDO
	

	
	
end subroutine INIT

! **********************************************************************

END MODULE INPUT
