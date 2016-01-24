MODULE ADDITIONAL
  IMPLICIT NONE

CONTAINS

! **********************************************************************
  FUNCTION function_surface_eq(i , phi_i , m_i )

    ! room position and dimensions are given in module CONSTANTS
    USE CONSTANTS, ONLY : x_room, y_room, z_room, l_room, b_room, h_room,RHO
 

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: phi_i, m_i
    DOUBLE PRECISION :: b_w, h, vol_i, vol_tot, vol_a, h_w
    DOUBLE PRECISION :: n_surf(3), r_surf(3), d_surf
    DOUBLE PRECISION :: function_surface_eq(4)

    DOUBLE PRECISION :: vol_cond1, vol_cond2

    ! normal vector of the suface in ship coords, pointing above
    n_surf(1) =         0.0 ! x-direction, in 2D-model is zero
    n_surf(2) =  SIN(phi_i)
    n_surf(3) = -COS(phi_i)

    vol_i = m_i/RHO
    vol_tot = l_room(i)*b_room(i)*h_room(i) ! total volume of the room
    vol_a = vol_tot - vol_i ! air volume in the room i

    ! conditions for different surface positions
    vol_cond1 = 0.5* ABS( TAN(phi_i)*l_room(i)*b_room(i)**2 )
    vol_cond2 =  0.5*ABS(1.0/TAN(phi_i))*l_room(i)*h_room(i)**2
    IF (vol_cond2 > 0.5*vol_tot) THEN
       vol_cond2 = 0.5*vol_tot
    END IF
       


    ! bottom is partially covered by water
    bottom_covered:IF ( vol_i < vol_cond1 .AND. vol_a > vol_cond2 ) THEN

       ! breadth of the area in room covered by water
       b_w = SQRT(2.*vol_i/(ABS(TAN(phi_i))*l_room(i)) ) 

       ! maximum height of the water in room
       h_w = 2.*vol_i /(b_w*l_room(i) )

       ! bottom visible and water does not reach room top
       bottom_vis_topcond:IF ( h_w < h_room(i) ) THEN 
          ! surface goes through the point
          r_surf(1) = 0.0 + x_room(i)
          r_surf(2) = SIGN(1d0,phi_i)*( b_w - b_room(i)/2 ) + y_room(i)
          r_surf(3) = 0.0 + z_room(i)

          
       ELSE ! water reaches room top, but bottom stays partially covered

          ! average breadth of the water
          b_w = vol_i/(l_room(i)*h_room(i)) 
       
          ! surface goes through the point
          r_surf(1) = 0.0 + x_room(i)
          r_surf(2) = SIGN(1d0,phi_i)*(-0.5*b_room(i) + b_w) + y_room(i)
          r_surf(3) = -0.5*h_room(i)  + z_room(i)

       END IF bottom_vis_topcond


    ! bottom is fully covered by water, but water is not reaching the roof
    ELSE IF ( vol_a > vol_cond1 ) THEN 

       h = vol_i/(l_room(i)*b_room(i)) 
       
       ! surface goes through the point
       r_surf(1) = 0.0 + x_room(i)
       r_surf(2) = 0.0 + y_room(i)
       r_surf(3) = -h  + z_room(i)


    ELSE ! bottom fully covered and water reaches room top

       ! breadth of the airpocket area on room top
       b_w = SQRT(2.*vol_a/(ABS(TAN(phi_i))*l_room(i)) ) 

       ! surface goes through the point       
       r_surf(1) = 0.0 + x_room(i)
       r_surf(2) = -SIGN(1d0,phi_i)*( b_w - b_room(i)/2 ) + y_room(i)
       r_surf(3) = -h_room(i) + z_room(i)
       

    END IF bottom_covered

    ! surface equation in ship coordinates
    ! equation for the surface is ax + by + cz + d = 0,
    ! where {a,b,c} = n_surf and d = d_surf = -n_surf . r_surf 
    d_surf = -DOT_PRODUCT(n_surf,r_surf)

    function_surface_eq(1:3) = n_surf
    function_surface_eq(4)   = d_surf


  END FUNCTION function_surface_eq



! **********************************************************************

	FUNCTION randomnumber(N)

	!Creates a random point RPOINT = [0,1] to be used for
	!the rcm method in sampling. The algorithm uses the Van der Corput
	!sequence.
		
	
		implicit none
			
		integer,intent(IN) :: N	
		!==!
		integer ::  J,K1, K2, AJ, BIG_AJ,TEST
		double precision :: THETA,randomnumber
		                                                                                                                                            
		!coefficients for the sequence
		K1 = 5 
		K2 = 3 	
		                                                                     
		TEST = 100	!a_j, set to arbitrary value at first
		J = 0
		THETA = 0


		DO WHILE (TEST .GE. 1)
		
			TEST = N/(K1**J)
			AJ = MOD(TEST,K1)
			BIG_AJ = MOD(K2*AJ,K1)
			THETA = THETA + BIG_AJ*K1**(-REAL(J)-1)
			J = J+1
		
		ENDDO
		
		RANDOMNUMBER = THETA
		

	END FUNCTION randomnumber

! **********************************************************************

  FUNCTION pressure_on_room_side(i , y_loc , z_loc )
		!returns the pressure on either left or ride side of the room at water
		!height z_loc. If z_loc is over the water, 0 value is returned.
		!(y_loc,z_loc) must be in hull coordinates!
		 
    
    USE CONSTANTS, ONLY : x_room, y_room, z_room, l_room, b_room, h_room,RHO, DY,JMAX
		USE ARRAYS, ONLY: H_NOW,A_Z
 		

    IMPLICIT NONE

    integer, intent(IN) :: i
    double precision, intent(IN) :: y_loc,z_loc
    !==!
    double precision :: pressure_on_room_side
		double precision :: DELTA_Z,PRESSURE,Y_MIN,Y_MAX,Z_MIN,Z_MAX

		


		Y_MIN = y_room(i) - (b_room(i)/2.) - DY(i); Y_MAX = y_room(i) + (b_room(i)/2.) + DY(i)
		Z_MIN = z_room(i); Z_MAX = z_room(i) + (h_room(i))


		!IF (y_loc .GT. Y_MAX  .OR. y_loc .LT. Y_MIN &
			!.OR. z_loc .GT. Z_MAX .OR. z_loc .LT. Z_MIN) THEN
			!PRINT *,"attempting to get pressure values outside room, stopping program"
			!STOP
		!ENDIF


		DELTA_Z = abs(z_loc - z_room(i))

		IF (y_loc .LT. y_room(i)) THEN
			PRESSURE = (H_NOW(i,1)-DELTA_Z)*A_Z(i,1)*RHO
		ELSE
			PRESSURE = (H_NOW(i,JMAX(i))-DELTA_Z)*A_Z(JMAX(i),1)*RHO
		ENDIF

		IF (PRESSURE .LT. 0.) THEN
			PRESSURE = 0.
		ENDIF
	
		pressure_on_room_side = PRESSURE


    


	END FUNCTION pressure_on_room_side


! **********************************************************************

  FUNCTION open_side(i,conn_open,n_opening,y_room,n_room)
	!returns 1 if left side of the room i has an opening, 2 if the right side has opening
	!and -1 if no openings in the room
		 
 		

    IMPLICIT NONE

    integer, intent(IN) :: i,n_opening,n_room
		integer, intent(IN),dimension(1:n_opening,2) :: conn_open	
		double precision, intent(IN),dimension(1:n_room) :: y_room
    !==!
    integer :: open_side, N,NEIGHBOUR
		
		IF (n_opening .EQ. 0) THEN
				open_side = -1
				GO TO 100
		ENDIF

		
		DO N = 1,n_opening		

		  open_side = -1

			IF(conn_open(N,1) .EQ. i) THEN
				NEIGHBOUR = conn_open(N,2)
				IF (y_room(i) .LT. y_room(NEIGHBOUR)) THEN
						open_side = 2
				ELSE 
						open_side = 1		
				ENDIF 												
			ENDIF

		 IF(conn_open(N,2) .EQ. i) THEN
				NEIGHBOUR = conn_open(N,1)
				IF (y_room(i ) .LT. y_room(NEIGHBOUR)) THEN
						open_side = 2
				ELSE 
						open_side = 1		
				ENDIF
		 ENDIF

		ENDDO
		

    100 CONTINUE


	END FUNCTION open_side

! **********************************************************************

  FUNCTION open_height(i,conn_open,n_opening,n_room,A_open,l_room)
	!returns the height of the opening corresponding to room i
	!if no opening in room i, 0 is returned
		 
    IMPLICIT NONE

    integer, intent(IN) :: i,n_opening,n_room
		integer, intent(IN),dimension(1:n_opening,2) :: conn_open
		double precision, intent(IN),dimension(1:n_room) :: l_room	
		double precision, intent(IN),dimension(1:n_opening) :: A_open	
    !==!
		integer :: N
    double precision :: open_height,area,height
		

		
		DO N = 1,n_opening		
			height = 0.
			IF(conn_open(N,1) .EQ. i .OR. conn_open(N,2) .EQ. i) THEN
				area = A_open(N)
				height = area/l_room(i)											
			ENDIF
		ENDDO

		open_height = height


	END FUNCTION open_height

! **********************************************************************




END MODULE ADDITIONAL
