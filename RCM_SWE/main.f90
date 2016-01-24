!
!Program RCM:
!
!-Solves the one-dimensional shallow water equations using the Glimm's Random Choice Method
!
!-Allows the solution of multiple square tank domains which are excited via harmonic roll or sway motion.
!
!-Solves the possible water exchange between the domains using pressure difference and Bernoulli's equation 
!to calculate the flowrate from one tank to another
!
!-The Riemann solver in this code is written based on the solver presented in "Shock-Capturing methods for free-
!surface shallow water flow" by E.F. Toro
!
!Read the README.pdf file before setting the inputs and use tab width of two spaces when opening the code files!!
!
!
!
!	-=petteri.peltonen@aalto.fi=-
!========================================================================================================



PROGRAM RCM

  use WORKSPACE
  use ARRAYS
  use CONSTANTS
	use INPUT
	use ADDITIONAL
	use FLUID_SOLVER
	use OUTPUT

  implicit none
	
	INTEGER I,N
	logical :: exist

	!read input files
	CALL READ_INPUTS()
		
	!allocate memory for the arrays which are solved at each time level
	allocate (V_NOW(1:n_room,0:MAXVAL(JMAX)+1))
	allocate (V_NEW(1:n_room,0:MAXVAL(JMAX)+1))
	allocate (H_NOW(1:n_room,0:MAXVAL(JMAX)+1))
	allocate (H_NEW(1:n_room,0:MAXVAL(JMAX)+1))
	allocate (A_Z(1:n_room,0:MAXVAL(JMAX)+1))
	allocate (F_Y(1:n_room,0:MAXVAL(JMAX)+1))

	!initialize the room domains
	CALL INIT(JMAX,V_NOW,V_NEW,H_NOW,H_NEW,OUTPUT_TIMES,RHO,DY,n_room,&
		b_room,l_room,x_room,y_room,z_room,OUTPUT_INTERVAL,OUTPUT_N)

	WRITE (*,*)"-------------------------------------"
	WRITE (*,*)"STARTING SIMULATION"
	WRITE (*,*)"-------------------------------------"

	DO N=1,10000000
		
		TIME = TIME + DT

		!calculate source terms for a box movement case. If the motion of the box and water motion are coupled, another 
		!routine will have to be written

		CALL SOURCE_TERM(JMAX,V_NOW,H_NOW,F_Y,A_Z,GRAVITY,RHO,n_room,BOX_MOV,b_room,l_room,&
				x_room,y_room,z_room,DY,AMPLITUDE,FREQ,TIME,Q1,PHI)

		!calculate new time step
		CALL ADJUST_DT(N,JMAX,V_NOW,H_NOW,A_Z,DY,DT,TIME,CFL_MAX,n_room)
		
		!apply boundary conditions to ghost cells
		CALL BOUND(JMAX,V_NOW,H_NOW,A_Z,F_Y,DY,DT,GRAVITY,RHO,n_room,n_opening,conn_open,r_open,C_d,A_open,l_room)


		!solve the Riemann problems for all rooms
		CALL RIEMANN(JMAX,V_NOW,H_NOW,V_NEW,H_NEW,A_Z,DY,DT,n_room,N)
	
		!update the variables H_NOW,V_NOW to new time level
		!this is the operational splitting phase
		CALL UPDATE_VARIABLES(JMAX,V_NOW,H_NOW,H_NEW,V_NEW,F_Y,DT,n_room)

		!output called at every N but outputs only at certain time instances (see runparameters)
		CALL OUTPUT_RESULTS()

		!exit if simulation time is reached
		IF (TIME .GE. T_MAX) THEN
				EXIT
		ENDIF				

	ENDDO
	
 


END PROGRAM RCM

