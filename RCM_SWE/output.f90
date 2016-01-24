MODULE OUTPUT
  IMPLICIT NONE
	
CONTAINS

!********************************************************************************************
	subroutine OUTPUT_RESULTS()
	!Outputs the results. Called at every N but does not output at every time level.
		
		USE CONSTANTS
		USE ARRAYS
		USE WORKSPACE

		implicit none
		
		logical :: TIMETOWRITE		
		integer :: I
		double precision :: DIF

	
		IF (TIME .EQ. 0. ) THEN
			!open files and write headers
			CALL WRITE_HEADERS(CASE_NAME)
			
		ENDIF
		

		!check whether it is time to output
		TIMETOWRITE = .FALSE.
		DO I = 0,OUTPUT_N
			DIF = ABS(TIME-OUTPUT_TIMES(I))
			IF (DIF .LT. DT/2.) THEN
						TIMETOWRITE = .TRUE.
			ENDIF
		ENDDO

	
		!if time to output, output
		IF (TIMETOWRITE .EQV. .TRUE.) THEN

			CALL WRITE_MASS(JMAX,H_NOW,RHO,TIME,DY,n_room,l_room)
			CALL WRITE_GIF_FILE(JMAX,H_NOW,V_NOW,RHO,TIME,DY,PHI,n_room,y_room,b_room)
			CALL WRITE_FORCES(JMAX,H_NOW,A_Z,RHO,TIME,DY,n_room,l_room,b_room,y_room,conn_open,n_opening,A_open)

		ENDIF
		





			
			                                                                     		

	end subroutine OUTPUT_RESULTS




!********************************************************************************************
	subroutine WRITE_HEADERS(CASE_NAME)
	!opens the result files and writes headers to the first row
		
		

		implicit none
		character, intent(IN) :: CASE_NAME
		character(LEN=40) :: write2file		

		
			
		WRITE(write2file,'(A18)') 'output/massbal.out'
		
		OPEN(13,FILE=write2file)			!mass balance output

		WRITE(write2file,'(A21)') 'output/simulation.out'
		OPEN(14,FILE=write2file)			!.gif plotting output

		WRITE(write2file,'(A23)') 'output/waterheights.out'
		OPEN(15,FILE=write2file)			! water level at some cross section


		WRITE(write2file,'(A23)') 'output/distribution.out'
		OPEN(16,FILE=write2file)			! water level at some cross section

		WRITE(write2file,'(A17)') 'output/forces.out'
		OPEN(17,FILE=write2file)			!forces and moments exerted on the compartment


		WRITE(13,'(A,7X,A,10X,A,10X,A,10X,A,10X,A,10X,A)')    &
			  "#","T","ROOM 1","ROOM 2", "ROOM 3","ROOM 4","ROOM I"

		WRITE(14,'(A,7X,A,12X,A,12X,A,12X,A,12X,A)')    &
			  "#","X","H","V","T","PHI"


		WRITE(15,'(A,7X,A,12X,A,12X,A)')    &
			  "#","T","H","Y"

		WRITE(16,'(A,7X,A,12X,A,12X,A)')    &
			  "#","X","H","V"  


		WRITE(17,'(A,7X,A,12X,A,12X,A,12X,A)')    &
	    "#","T","M","F"


			                                                                     		

	end subroutine WRITE_HEADERS


!********************************************************************************************


	subroutine WRITE_MASS(JMAX,H_NOW,RHO,TIME,DY,n_room,l_room)
	!writes the total amount of mass in each room
		
		

		implicit none
		integer, intent(IN) :: n_room
		double precision, intent(IN) :: RHO,TIME
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		double precision, intent(IN),dimension(1:n_room) :: DY,l_room
		double precision,intent(IN),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: H_NOW
		!==!
		integer :: I,J
		double precision :: MASS
		double precision,allocatable,dimension(:) :: MASS_VECTOR

		ALLOCATE( MASS_VECTOR(1:n_room) )

		DO I = 1,n_room
			MASS = 0.
			DO J = 1, JMAX(I)
				MASS = MASS + H_NOW(I,J)*RHO*l_room(I)*DY(I)
			ENDDO		
			MASS_VECTOR(I) = MASS
		ENDDO

		WRITE(13,'(10E14.6)') TIME,MASS_VECTOR
		
		DEALLOCATE(MASS_VECTOR)


			                                                                     		

	end subroutine WRITE_MASS


!********************************************************************************************


	subroutine WRITE_GIF_FILE(JMAX,H_NOW,V_NOW,RHO,TIME,DY,PHI,n_room,y_room,b_room)
	!writes the water heights and velocities in a format which can
	!be used in the plot_simu.py script for creating a .gif of the 
	!simulation
		
		

		implicit none
		integer, intent(IN) :: n_room
		double precision, intent(IN) :: RHO,TIME,PHI
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		double precision, intent(IN),dimension(1:n_room) :: DY,y_room,b_room
		double precision,intent(IN),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: H_NOW,V_NOW
		!==!
		integer :: I,J
		double precision :: Y_COOR

		DO I = 1,n_room
			
			DO J = 1, JMAX(I)
				!y-coordinates in the ship system
				Y_COOR = y_room(I)-(b_room(I)/2.)+(DY(I)/2.)+(DY(I)*J)
				WRITE(14,'(10E14.6)') Y_COOR,H_NOW(I,J),V_NOW(I,J),TIME,PHI 
			ENDDO		
			
		ENDDO
	


			                                                                     		

	end subroutine WRITE_GIF_FILE


!********************************************************************************************


	subroutine WRITE_FORCES(JMAX,H_NOW,A_Z,RHO,TIME,DY,n_room,l_room,b_room,y_room,conn_open,n_opening,A_open)
	!writes the force and moment exerted on the compartment
	!the moment is calculated with respect to the origo of the ship coordinate system.
		
		USE ADDITIONAL , ONLY: open_side,open_height

		implicit none
		integer, intent(IN) :: n_room,n_opening
		double precision, intent(IN) :: RHO,TIME
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		integer, intent(IN),dimension(1:n_opening,2) :: conn_open	
		double precision, intent(IN),dimension(1:n_opening) :: A_open
		double precision, intent(IN),dimension(1:n_room) :: DY,l_room,b_room,y_room
		double precision,intent(IN),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: H_NOW,A_Z
		!==!
		integer :: I,J,SIDE
		double precision :: MOMENT,Y_COOR,FORCE_L,FORCE_R,FORCE,HEIGHT
	
		MOMENT = 0.
		FORCE = 0.
		DO I = 1,n_room

			!get the open side
			SIDE = open_side(i,conn_open,n_opening,y_room,n_room)

			!no open sides in room
			IF (SIDE .EQ. -1) THEN
				FORCE_L = 0.5*RHO*A_Z(I,1)*(H_NOW(I,1))**2	
				FORCE_R = 0.5*RHO*A_Z(I,JMAX(I))*H_NOW(I,JMAX(I))**2
			
			ELSE
				!left side of the room is open
				IF (SIDE .EQ. 1) THEN
					!get the opening height
					HEIGHT = open_height(i,conn_open,n_opening,n_room,A_open,l_room)
					FORCE_L = 0.5*RHO*A_Z(I,1)*(H_NOW(I,1)**2 - HEIGHT**2)
					FORCE_R = 0.5*RHO*A_Z(I,JMAX(I))*H_NOW(I,JMAX(I))**2
				ENDIF
				!right side of the room is open
				IF (SIDE .EQ. 2) THEN
					!get the opening height
					HEIGHT = open_height(i,conn_open,n_opening,n_room,A_open,l_room)
					FORCE_L = 0.5*RHO*A_Z(I,1)*(H_NOW(I,1))**2	
					FORCE_R = 0.5*RHO*A_Z(I,JMAX(I))*(H_NOW(I,JMAX(I))**2 - HEIGHT**2)
				ENDIF
			ENDIF
			
			!calculate the force 
			FORCE = FORCE + (FORCE_R-FORCE_L)*l_room(I)
			
			!calculate the moment
			DO J = 1, JMAX(I)
				Y_COOR = y_room(I)-(b_room(I)/2.)+(DY(I)/2.)+(DY(I)*J)
				MOMENT = MOMENT + (RHO*A_Z(I,J)*H_NOW(I,J)*l_room(I)*DY(I))*Y_COOR
			ENDDO

		ENDDO

		WRITE(17,'(10E14.6)') TIME,MOMENT,FORCE

			                                                                     		

	end subroutine WRITE_FORCES


!********************************************************************************************



END MODULE OUTPUT
