MODULE FLUID_SOLVER
  IMPLICIT NONE
	
CONTAINS

	! **********************************************************************

	subroutine RIEMANN(JMAX,V_NOW,H_NOW,V_NEW,H_NEW,A_Z,DY,DT,n_room,N)
	!solves and samples the Riemann problems for all cells and all rooms.
	!calls multiple routines in the process.
		
		use ADDITIONAL , ONLY: randomnumber

		implicit none

		integer, intent(IN) :: n_room,N
		double precision, intent(IN) :: DT
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		double precision, intent(IN),dimension(1:n_room) :: DY
		double precision,intent(IN),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: V_NOW,H_NOW,A_Z
		!==!
		double precision,intent(OUT),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: V_NEW,H_NEW
		!==!
		integer :: I,J
		double precision :: RPOINT,HL,VL,AL,HR,VR,AR,Y,A,H,V

		!get a random variable for sampling
		!a single value is used for all cells in all rooms at time level N
		RPOINT = randomnumber(N)

		!loop over all rooms
		DO I = 1,n_room	

			!loop over all cells
			DO J = 1, JMAX(I)

						!the sample point is on the left side of the cell
						IF (RPOINT .LE. 0.5) THEN 

									!left initial values of the Riemann problem
									HL = H_NOW(I,J-1) 									!water height
				     			VL = V_NOW(I,J-1) 									!velocity 
				     			AL = SQRT(A_Z(I,J-1)*H_NOW(I,J-1))  !celerity

				          !right initial values of the Riemann problem                                    
				     			HR = H_NOW(I,J) 									
				     			VR = V_NOW(I,J)										
				     			AR = SQRT(A_Z(I,J)*H_NOW(I,J)) 			
								
									!the value point Y = (y,t) used in the sampling
									Y = RPOINT*(DY(I)/DT)

				 		!sample point is on the right side of the cell
						ELSE  
									!left initial values of the Riemann problem										
									HL = H_NOW(I,J) 
				     			VL = V_NOW(I,J) 
				     			AL = SQRT(A_Z(I,J)*H_NOW(I,J))  		
				                
									!right initial values of the Riemann problem                                                   
				     			HR = H_NOW(I,J+1) 
				     			VR = V_NOW(I,J+1)
				     			AR = SQRT(A_Z(I,J+1)*H_NOW(I,J+1))
		
									Y = (RPOINT-1.)*(DY(I)/DT)
				
						ENDIF
				
						A = A_Z(I,J)			!acceleration in cell J, explicit treatment
	
						!solve and sample the solution
						CALL SOLVE_CELL(HL,VL,AL,HR,VR,AR,Y,A,H,V)

						!print error message if NaN
						IF (H /= H .OR. V /= V) THEN
								PRINT *,"FLUID SOLVER CRASHES!"
								STOP
						ENDIF
				
						!set the values to temporary variables 'new'
						!the final values are obtained after the operational splitting
						H_NEW(I,J) = H
						V_NEW(I,J) = V

			ENDDO

		ENDDO	


	end subroutine RIEMANN


	! *********************************************************************


	subroutine SOLVE_CELL(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!determines if a dry region solver or a wet region solver must be called
	!and calls the appropriate solver
		
	
		implicit none
	
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V
		double precision :: DV_CRIT
	
		!this is the critical depth condition, which determines wheter
		!dry region is in the cell
		DV_CRIT = 2.0*(AL+AR) -VR + VL
	
		!call the correct solver based on the condition
		IF (DV_CRIT .LE. 0. .OR. HR .LE. 0. .OR. HL .LE. 0.) THEN
		 
			CALL DRY_CASES(HL,VL,AL,HR,VR,AR,Y,A,H,V)
			
		ELSE
		
			CALL WET_CASES(HL,VL,AL,HR,VR,AR,Y,A,H,V)
		
		ENDIF

	end subroutine SOLVE_CELL


	! **********************************************************************


	subroutine DRY_CASES(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!Riemann problem solution for dry cases, three possibilities.

		implicit none
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V


		IF (HL .LE. 0.) THEN
				!the dry region is on the left side				
		 		CALL DRY_LEFT(HL,VL,AL,HR,VR,AR,Y,A,H,V)

		ELSE 
				IF (HR .LE. 0.) THEN
				!the dry region is on the right side		
				CALL DRY_RIGHT(HL,VL,AL,HR,VR,AR,Y,A,H,V)								
		

			ELSE
				!the dry region is between left and right sides									
				CALL DRY_MIDDLE(HL,VL,AL,HR,VR,AR,Y,A,H,V)

			ENDIF
		ENDIF
	 	
	end subroutine DRY_CASES


	! **********************************************************************

	subroutine DRY_MIDDLE(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!when the middle section is dry, the solution consists of 
	!two rarefactions on the right and left side.


		implicit none
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V
		double precision :: SR,SR2,SL,SL2,A_TEMP 		!auxilary variables

		                                             
		!wave speeds corresponding to the left eigenvalue                                                       
		SL = VL - AL 			
		SL2 = VL + 2.0*AL 

		!wave speeds corresponding to the right eigenvalue
		SR2 = VR - 2.0*AR 
		SR = VR + AR 			
		
		
		                                                                     
		IF(Y.LE.SL) THEN                                                                        
			 !Sample point lies to the left of the left                    
		   !rarefaction                                                    
		                                                                     
		   H = HL 
		   V = VL 
		ENDIF 
		                                                                     
		IF(Y.GT.SL.AND.Y.LE.SL2) THEN 
			 !Sample point inside the left rarefaction 
		                                                               
		   V = (VL + 2.0*AL + 2.0*Y)/3.0 
		   A_TEMP = (VL + 2.0*AL - Y)/3.0 
		   H = A_TEMP*A_TEMP/A 
		ENDIF 
		                                                                    
		IF(Y.GT.SL2 .AND. Y.LE.SR2 ) THEN 
			 !Sample point in the dry region                                                                       
			                                                                    
		   H = 0. 
		   V = 0. 
		ENDIF 
		                                                                     
		IF(Y.GT.SR2 .AND. Y .LE. SR) THEN 
			 !Sample point inside the right rarefaction                                                                        
				           		                                                                      
		   V = ( VR - 2.0*AR + 2.0*Y)/3.0 
		   A_TEMP = (-VR + 2.0*AR + Y)/3.0 
		   H = A_TEMP*A_TEMP/A	 
		ENDIF 
		                                                                     
		IF(Y.GT.SR)THEN 
			 !Sample point to the left of the right rarefaction                                                                       
		                                                                                                                        
		   H = HR 
		   V = VR 
		ENDIF 
		                                                                     
		                                  
	

	
	end subroutine DRY_MIDDLE


	! **********************************************************************

	subroutine DRY_LEFT(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!if the left side is dry the solution consists of only single right 
	!rarefaction
		implicit none
	
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V
		double precision :: SR,SR2,A_TEMP 		!auxilary variables

		!wave speeds corresponding to the right eigenvalue
		SR = VR + AR
		SR2 = VR - 2.*AR

		IF (Y .LT. SR2) THEN
			!the sample point on the left of the rarefaction ! y < (u_R - 2a_R) !
			!this is the dry region
			H = HL
			V = VL
				 
		ELSE
			IF (Y .LE. SR) THEN
				!sample point inside the rarefaction  ! (u_R - 2a_R) < y < (u_R + a_R) !
				V = (VR - 2.*AR + 2.*Y)/3.
				A_TEMP = (-VR + 2.*AR + Y)/3.
				H = A_TEMP*A_TEMP/A
	

			ELSE
				!sample point on the right of the rarefaction ! y > (u_R + a_R) !
				H = HR
				V = VR
			ENDIF
		ENDIF
	

	
	
	end subroutine DRY_LEFT


	! **********************************************************************
	subroutine DRY_RIGHT(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!if the right side is dry the solution consists of only single left 
	!rarefaction

		implicit none
	
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V
		double precision :: SL,SL2,A_TEMP 		!auxilary variables
	
		!wave speeds corresponding to the left eigenvalue
		SL = VL - AL
		SL2 = VL + 2.*AL

		IF (Y .LE. SL) THEN
			!the sampling point is on the left of the rarefaction	! y < (u_L - a_L) !			
			H = HL
			V = VL
	
			ELSE		
				IF (Y .LE. SL2 ) THEN
					!the sampling point is inside the rarefaction ! (u_L - a_L) < y < (u_L + 2a_L) !		
					V = (VL + 2.*AL + 2.*Y)/3. 			
					A_TEMP = (VL + 2.*AL - Y)/3. 
					H = A_TEMP*A_TEMP/A

				ELSE

				!the sampling point is on the right of the rarefaction ! y > (u_L +2 a_L)
				!this is the dry region 
				H = HR 
				V = VR
			ENDIF

		ENDIF
		 	
	end subroutine DRY_RIGHT


	! **********************************************************************

	subroutine WET_CASES(HL,VL,AL,HR,VR,AR,Y,A,H,V)
	!Riemann problem solution for postivie water heights.
		
		!USE ADDITIONAL, ONLY:H_ZERO,F_AND_F_PRIME
		implicit none
		double precision,intent(IN) :: HL,VL,AL,HR,VR,AR,Y,A
		double precision,intent(OUT) :: H,V	

		!auxilary variables
		double precision :: H0,FL,FR
		double precision :: FL_PRIME,FR_PRIME,F,F_PRIME
		double precision :: H_STAR,DH,H_PREV,TOL,A_STAR,V_STAR
		double precision :: SL,QL,STL,A_TEMP,SHL 		
		double precision :: QR, SR, SHR, STR, REL
		double precision :: H_LOW,H_HIGH,DH_PREV
	
	
	
		integer :: I, ITER
		
		!get the initial value for Newton-Raphson
		CALL H_ZERO(HL,VL,AL,HR,VR,AR,Y,A,H0)

		!iteration tolerance
		TOL = EPSILON(H0)*1000.
	
		ITER = 1000

		H_LOW = -TOL
		H_HIGH = 5.*H0

		!calculate h* using the N-R iteration
		DO I = 1,100000

			!calculate the functions and their derivatives
			CALL F_AND_FPRIME(FL,FL_PRIME,H0,HL,AL,A)	!f_L & f'_L
			CALL F_AND_FPRIME(FR,FR_PRIME,H0,HR,AR,A) !f_R & f'_R
			F = FL + FR + (VR-VL)
			F_PRIME = FL_PRIME + FR_PRIME
	
			IF (F .GT. 0.) THEN
				H_HIGH = H0
			ELSE 
				H_LOW = H0
			ENDIF
		
			!Newton-Rapshon
			H0 = H0 - (F/F_PRIME)	

			!use bisection if iteration out of bounds, very slow, or starts to oscillate
			IF ((H0 - H_LOW)*(H0-H_HIGH) .GT. 0. .OR. &
					ABS(2*F) .GT. ABS(DH_PREV*F_PRIME) .OR. DH_PREV .EQ. DH) THEN
			
				H0 = 0.5*(H_HIGH+H_LOW)
				DH_PREV = DH
				DH = ABS(H0 - H_PREV)/(H0+H_PREV)/2

			ELSE
				!otherwise use Newton-Raphson
				DH_PREV = DH
				DH = ABS(H0 - H_PREV)/(H0+H_PREV)/2.
			ENDIF

			IF (I .EQ. ITER .AND. DH .GT. TOL*100) THEN
				!stop program if iteration is not converging at all
				PRINT *,"The N-R iteration in the fluid solver is not converging!"
				STOP
			ENDIF
	
			H_PREV = H0
			IF (DH .LT. TOL) EXIT
	
		ENDDO
	
		!after the iteration the mid region water height is stored in H0	

		!mid-region variables
		H_STAR = H0 !store the iterated value for h* 
		A_STAR = SQRT(H_STAR*A)  !calculate a*
		V_STAR = 0.5*(VL + VR) + 0.5*(FR - FL) !calculate v*  

	



		!at this point we have the solution for the star region
		!now we need to determine the solutions on the left
		!and right sides of the star region based on 
		!where the sample point Y lies.
	
		!==================================================!
		!									LEFT WAVE												 !
		!==================================================!
		IF (Y .LE. V_STAR) THEN
	
			IF (H_STAR .GE. HL) THEN 
			!the left wave is a shock

				QL = SQRT((H_STAR + HL)*H_STAR/(2.0*HL*HL))	!	eq. 5.26 
				SL = VL - AL*QL
	
					IF (Y .LE. SL) THEN 
						!the left wave is a shock and sample point on the left of the wave			
						H = HL
						V = VL	
							
					ELSE 							
						!the left is a shock but sample point in the star region		
						H = H_STAR
						V = V_STAR
					ENDIF


			ELSE 
			!the left wave is rarefaction wave	

					SHL = VL - AL					!eq. 5.27
					STL = V_STAR - A_STAR

					IF(Y.LE.SHL)THEN  
						!left wave is rarefaction and sample point on the left of the rarefaction                                                                                                                                   
						H = HL 
				 	 	V = VL 
		
					ELSE
				 				
						IF (Y .LE. STL)	THEN 
							 !left wave rarefaction and sample point inside the rarefaction
							 V = (VL + 2.0*AL + 2.0*Y)/3.0 
							 A_TEMP = (VL + 2.0*AL - Y)/3.0 
							 H = A_TEMP*A_TEMP/A 
							
						ELSE 									
							!left wave rarafaction but sample point inside star region
							H = H_STAR
							V = V_STAR
						ENDIF

					ENDIF

			ENDIF

	
		!==================================================!
		!									RIGHT WAVE											 !
		!==================================================!
		ELSE

			IF (H_STAR .GE. HR) THEN 	
			!the right wave is a shock

				QR = SQRT((H_STAR + HR)*H_STAR/(2.0*HR*HR)) !eq. 5.31
				SR = VR + AR*QR 
	
				IF (Y .GE. SR) THEN
					!right wave is a shock and sample point on the right of the wave 
					H = HR
					V = VR

				ELSE 
					!right wave is a shock but sample point in the star region
					H = H_STAR
					V = V_STAR
				ENDIF


			ELSE 
			!right wave is a rarefaction
	
				SHR = VR + AR 
				STR = V_STAR + A_STAR

				IF (Y .GE. SHR) THEN
					!right wave is a rarefaction and sample point on the right of the wave
					H = HR
					V = VR

				ELSE 

					IF(Y.GE.STR) THEN 
							!right wave is a rarefaction and sample point inside the rarefaction
							V = (VR  - 2.0*AR + 2.0*Y)/3.0 
				      A_TEMP = (-VR + 2.0*AR + Y)/3.0 
				      H = A_TEMP*A_TEMP/A 

					ELSE 
							!right wave is a rarefaction but sample point inside the star region
							H = H_STAR
							V = V_STAR

					ENDIF

				ENDIF

			ENDIF
			
		ENDIF
	 	
	end subroutine WET_CASES


	! **********************************************************************




		subroutine F_AND_FPRIME(F,FD,H,HK,AK,A) 
		!calculates the functions f_k(h*,h_k) and derivatives f_k(h*,h_k)'
		!the h* value comes from current iteration round 
		!of the N-R iteration algorithm this is called multiple times during the iteration   

						                           
			implicit none 
						                                                                                                                                    
			double precision ::  F,FD,H,HK,AK,A,A_STAR,GES
				                                                                  
			IF(H.LE.HK)THEN 
			!rarefaction wave                                                                                      
						                                                           
				A_STAR  = SQRT(A*H) 
				F  = 2.0*(A_STAR-AK) 
				FD = A/A_STAR 

			ELSE 
			!shock wave                                                                  
						                                                                                                 
				GES = SQRT(0.5*A*(H+HK)/(H*HK)) 
				F   = (H-HK)*GES 
				FD  = GES - 0.25*A*(H-HK)/(GES*H*H) 

			ENDIF 
				                                                                   
		end subroutine                  
		                      
	! **********************************************************************       

		subroutine H_ZERO(HL,VL,AL,HR,VR,AR,Y,A,H0)
		!initializes the Newton-Raphson iteration -> sets H* = H0
		
	
			implicit none
	
			double precision :: HL,VL,AL,HR,VR,AR,Y,A,HMIN,H0,GEL,GER
	
					                                                                   
			HMIN = MIN(HL,HR) 
				                                                                   
		! Use Two-Rarefaction (TRRS) solution as starting value             
				                                                                   
			H0 = (1.0/A)*(0.5*(AL+AR)-0.25*(VR-VL))**2 
				                                                                   
		 	IF(H0.LE.HMIN)THEN 
				                                                                   
		! Use Two-Rarefaction (TRRS) approximation as                    
		! starting value                                                 
				                                                                  
		
			ELSE 
				                                                                   
		! Use two-shock (TSRS) solution as starting value                
		! with H0 as computed from TRRS as estimate                      
				                                                                   
				                                                                   
			GEL = SQRT(0.5*A*(H0+HL)/(H0*HL)) 
			GER = SQRT(0.5*A*(H0+HR)/(H0*HR)) 
			H0  = (GEL*HL + GER*HR - (VR-VL))/(GEL + GER) 
				                                                                   
			ENDIF 


		!	PRINT *,H0


		end subroutine H_ZERO


	! **********************************************************************

	subroutine SOURCE_TERM(JMAX,V_NOW,H_NOW,F_Y,A_Z,GRAVITY,RHO,n_room,BOX_MOV,b_room,l_room,&
			x_room,y_room,z_room,DY,AMPLITUDE,FREQ,TIME,Q1,PHI)
	! calculates the source terms in a forced motion case. The source terms are calculated 
	!as in Dillingham 1981
		

		implicit none
	
		integer, intent(IN) :: n_room,BOX_MOV
		integer, intent(IN),dimension(1:n_room) :: JMAX
		double precision, intent(IN) :: RHO,GRAVITY,AMPLITUDE,FREQ,TIME
		double precision, intent(IN),dimension(1:n_room) :: b_room,l_room,x_room,y_room,z_room,DY
		double precision, dimension(1:n_room,0:MAXVAL(JMAX)+1),intent(IN) :: V_NOW,H_NOW
		!==!
		double precision, dimension(1:n_room,0:MAXVAL(JMAX)+1),intent(OUT) :: F_Y,A_Z
		double precision,intent(OUT) :: PHI,Q1
		!==!
		double precision :: ROLL_AMP,SWAY_AMP,PHI_PRIME,PHI_PRIME2,Q1_PRIME,Q1_PRIME2,Y_COOR
		integer :: I,J		
	


		IF (BOX_MOV .EQ. 2) THEN
			ROLL_AMP = 0.; SWAY_AMP = AMPLITUDE
			ELSE
	
				IF (BOX_MOV .EQ. 3) THEN	
					ROLL_AMP = AMPLITUDE; SWAY_AMP = 0.
				ELSE
					ROLL_AMP = 0.; SWAY_AMP = 0.
					IF (TIME .EQ. 0.) THEN
						PRINT *,"NO FORCED MOTION!"
					ENDIF

				ENDIF

		ENDIF
	
		!calculate the roll angle phi, roll speed phi', and roll acceleration phi''
		PHI = ROLL_AMP*COS(FREQ*TIME) 
		PHI_PRIME = -ROLL_AMP*FREQ*SIN(FREQ*TIME) 
		PHI_PRIME2 = -ROLL_AMP*FREQ*FREQ*COS(FREQ*TIME) 

		!calculate the location of the compartment, speed of the compartment and acceleration of the compartment 
		Q1 = SWAY_AMP*COS(FREQ*TIME)
		Q1_PRIME = -SWAY_AMP*FREQ*SIN(FREQ*TIME)
		Q1_PRIME2 = -SWAY_AMP*FREQ*FREQ*COS(FREQ*TIME)

	
		DO I = 1,n_room

			DO J = 1, JMAX(I)
				!y-coordinates in the ship system
				Y_COOR = y_room(I)-(b_room(I)/2.)+(DY(I)/2.)+(DY(I)*J)		
				!A_Z and F_Y 
				A_Z(I,J) = GRAVITY*COS(PHI) + PHI_PRIME2*(Y_COOR) - PHI_PRIME*PHI_PRIME*H_NOW(I,J) + 2.*PHI_PRIME*V_NOW(I,J) &
							 -Q1_PRIME2*SIN(PHI)
				F_Y(I,J) = -GRAVITY*SIN(PHI) + PHI_PRIME2*H_NOW(I,J) + PHI_PRIME*PHI_PRIME*(Y_COOR) - Q1_PRIME2*COS(PHI)
			
					
			ENDDO
		
		ENDDO                 
	
	end subroutine SOURCE_TERM


	! **********************************************************************

	subroutine ADJUST_DT(N,JMAX,V_NOW,H_NOW,A_Z,DY,DT,TIME,CFL_MAX,n_room)
	!adjusts time step according to the CFL < CFL_max input, called at every time
	!level
		
		implicit none
	
		integer, intent(IN) :: n_room,N
		double precision, intent(IN) :: TIME,CFL_MAX
		integer, intent(IN),dimension(1:n_room) :: JMAX
		double precision, intent(IN),dimension(1:n_room) :: DY
		double precision, dimension(1:n_room,0:MAXVAL(JMAX)+1),intent(IN) :: V_NOW,H_NOW,A_Z
		!==!
		double precision, intent(OUT) :: DT
		!==!
		integer :: I,J,ROOM_NUMBER
		double precision :: SN_MAX,TEMP,CFL
		real, dimension(3) :: OUT_ARRAY


		SN_MAX = 0.
	
		!find the largest wave speed
		DO I = 1, n_room

			DO J=1,JMAX(I)
				TEMP = ABS(V_NOW(I,J)+SQRT(H_NOW(I,J)*A_Z(I,J)))
				IF (TEMP .GE. SN_MAX) THEN
					SN_MAX = TEMP
					ROOM_NUMBER = I
				ENDIF	
			ENDDO
		
		ENDDO

		!use the largest wave speed to calcualte the smallest required DT
		DT = (CFL_MAX*DY(ROOM_NUMBER))/SN_MAX
		CFL = (DT*SN_MAX)/DY(ROOM_NUMBER)
		OUT_ARRAY(1) = TIME; OUT_ARRAY(2) = CFL; OUT_ARRAY(3) = DT

		!output TIME,CFL and DT for monitoring the progress
		IF (MOD(N,1000) .EQ. 0) THEN
			PRINT *,"TIME:", OUT_ARRAY(1),"CFL:", OUT_ARRAY(2),"DT:", OUT_ARRAY(3)
		ENDIF
	
	
	

	end subroutine ADJUST_DT


	! **********************************************************************

	subroutine BOUND(JMAX,V_NOW,H_NOW,A_Z,F_Y,DY,DT,GRAVITY,RHO,n_room,n_opening,conn_open,r_open,C_d,A_open,l_room)
	!updates ghost cells for boundary conditions, also takes care of the possible mass exchange 
	!if openings are given as inputs
		
		use ADDITIONAL, ONLY:pressure_on_room_side

		implicit none
		integer, intent(IN) :: n_room,n_opening
		double precision, intent(IN) :: DT,GRAVITY,RHO
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		integer, intent(IN),dimension(1:n_opening,2) :: conn_open		
		double precision, intent(IN),dimension(1:n_room) :: DY,l_room
		double precision, intent(IN),dimension(1:n_opening) :: C_d,A_open	
		double precision, intent(IN),dimension(1:n_opening,3) :: r_open	
		!==!
		double precision, dimension(1:n_room,0:MAXVAL(JMAX)+1),intent(OUT) :: V_NOW,H_NOW,A_Z,F_Y
		!==!
		integer :: I,J,FROM,TO
		double precision,  dimension(1:2) :: DH
		double precision :: VH,P1,P2,QH,DH_FROM,DH_TO,MOMFLUX


	

		DO I = 1,n_opening
			FROM = conn_open(I,1)
			TO  = conn_open(I,2)

			P1 = pressure_on_room_side(FROM , r_open(I,2) , r_open(I,3) )
			P2 = pressure_on_room_side(TO , r_open(I,2) , r_open(I,3) )

			!velocity from -> to
			VH = SIGN(1d0,(P1-P2))*SQRT((2./RHO)*ABS(P1-P2))

			!"volume" flowrate from -> to
			QH = C_d(I)*VH*(A_open(I)/l_room(FROM))
		
			!change of water height
			DH_FROM = (QH*DT)/DY(FROM) !ROOM 'FROM'
			DH_TO = (QH*DT)/DY(TO)   !ROOM 'TO'

			!change of momentum flux
			MOMFLUX = (rho*QH*VH*(A_open(I)/l_room(FROM)))

			!add the water height to the other room and subtract from the other
			!positive direction is from box 1 to box 2
			H_NOW(FROM,JMAX(FROM)) = -DH_FROM + H_NOW(FROM,JMAX(FROM)) 
			H_NOW(TO,1) = DH_TO + H_NOW(TO,1) 

			!add the momentum flux to the other room and subtract from the other
			F_Y(FROM,JMAX(FROM)) = F_Y(FROM,JMAX(FROM)) - MOMFLUX
			F_Y(TO,1) = F_Y(TO,1) + MOMFLUX
		
		
		ENDDO

		DO I = 1,n_room

			!closed left sides
			H_NOW(I,0) = H_NOW(I,1) 
			V_NOW(I,0) = -V_NOW(I,1)
			A_Z(I,0) = A_Z(I,1)
			F_Y(I,0) = 0.	
	
			!closed right sides
			H_NOW(I,JMAX(I)+1) = H_NOW(I,JMAX(I))
			V_NOW(I,JMAX(I)+1) = -V_NOW(I,JMAX(I))
			A_Z(I,JMAX(I)+1) = A_Z(I,JMAX(I))
			F_Y(I,JMAX(I)+1) = 0.

		ENDDO              
	 

	end subroutine BOUND


	! **********************************************************************
	subroutine UPDATE_VARIABLES(JMAX,V_NOW,H_NOW,H_NEW,V_NEW,F_Y,DT,n_room)
	!updates the variables to the next time level
	!uses the first order Euler scheme for operational splitting
	
		implicit none
		integer, intent(IN) :: n_room
		double precision, intent(IN) :: DT
		integer, intent(IN),dimension(1:n_room) :: JMAX	
		double precision,intent(IN),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: V_NEW,H_NEW,F_Y
		!==!
		double precision,intent(OUT),dimension(1:n_room,0:MAXVAL(JMAX)+1) :: V_NOW,H_NOW
		!==!
		integer :: I,J


		DO I = 1,n_room 

			DO J = 1,JMAX(I)

				H_NOW(I,J) = H_NEW(I,J) 
				V_NOW(I,J) = V_NEW(I,J) + ((DT)*F_Y(I,J))

			ENDDO     
		         
		ENDDO
	
	end subroutine UPDATE_VARIABLES

	! **********************************************************************



END MODULE FLUID_SOLVER

