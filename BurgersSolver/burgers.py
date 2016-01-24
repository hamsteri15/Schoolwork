#
#	#= FVM solver for 1D frictionless Burger's equation =#  	
#
# INTEM=1 --> UPWIND				
# INTEM=2 --> SECOND ORDER UPWIND		  
# INTEM=3 --> THIRD ORDER UPWIND		  
# INTEM=4 --> SECOND ORDER UPWIND BIASED	
# INTEM=5 --> CDS				
# INTEM=6 --> QUICK				
#
#POSITIVE INTEM VALUE IS USED FOR FLUX LIMITING	
#SEE THE AVAILABLE LIMITERS FROM THE FUNCTION limiter_func()			
#							petteri.peltonen@aalto.fi	
######################################################################################################



import numpy as np
import pylab as pl
import matplotlib
from sys import exit


def main():
	
	
	#ex. 1
	FVM(0.005,2,"albada2","LAX")
	FVM(0.005,2,"albada1","LAX")
	
	FVM(0.01,2,"albada2","LAX")
	FVM(0.01,2,"albada1","LAX")
	
	FVM(0.05,2,"albada2","LAX")
	FVM(0.05,2,"albada1","LAX")
	

	#ex. 2
	FVM(0.005,2,"vanleer","ROE")
	FVM(0.005,2,"monocent","ROE")
	
	FVM(0.01,2,"vanleer","ROE")
	FVM(0.01,2,"monocent","ROE")
	
	FVM(0.05,2,"vanleer","ROE")
	FVM(0.05,2,"monocent","ROE")

	#without limiters
	#FVM(0.005,-5,"vanleer","ROE")
	#FVM(0.005,-1,"monocent","ROE")
	
	#FVM(0.01,-5,"vanleer","ROE")
	#FVM(0.01,-1,"monocent","ROE")


#########################################################################################################




def FVM(timestep,intem,limiter,scheme):
	

	dx=0.1		#stepsize x
	nx=100
	L=dx*nx		#length			
	dt=timestep	#stepsize t
	mu=0		#viscosity
	time=4.		#simulation time
	
	
	S=1		#face areas, arbitrary value is used here
	V=S*dx		#cell volume
	
	#cell centre coordinates and time values for post processing
	xc=np.arange(dx/2,L,dx)
	t=np.arange(0,time+dt,dt)	

	
	nt=int(time/dt)
	
	
	#initialize the temperature "domain" U(x,t) and set the initial field 
	U=np.zeros((nx+4,nt+1))
	
	start=3+1./dx
	end=2+5./dx
	U[start:end,0]=2.
	
	#initialize flux vector F
	F=np.zeros(nx)
	
	n=0
	#time loop
	while (n<nt):

		U=bound(mu,U,dx,F,n,S,intem,dt)
		if (scheme=="ROE"):
			F=fluxROE(mu,U,dx,F,n,S,intem,limiter,dt)

		elif (scheme=="LAX"):
			F=fluxLAX(mu,U,dx,F,n,S,intem,limiter,dt)
		
		dU=dt/V*(F[:])
		
		U[2:-2,n+1]=U[2:-2,n]+dU[:]	

		n+=1
	

    #post processing comes here
	U1=U[2:-2,nt/2]
	U2=U[2:-2,-1]
	#np.savetxt('tulokset/{}_dt{}.out'.format(limiter,dt), np.c_[xc,U1,U2])
 	#np.savetxt('tulokset/upwind_dt{}.out'.format(dt), np.c_[xc,U1,U2])

	pl.plot(xc,U1)
	pl.plot(xc,U2)
	pl.show()

	
#calculates the net fluxes using the TVD Lax-Wendroff method.
#MUSCL schemes are _not_ used here
def fluxLAX(mu,U,dx,F,n,S,intem,limiter,dt):
	
	
	i=2
	j=0
	eps=1E-10
	while (j <= (len(F)-1)):

		diff0=U[i-1,n] - U[i-2,n]
		diff1=U[i,n]   - U[i-1,n]
		diff2=U[i+1,n]   - U[i,n]

		r_r =max(-1.E7,min(1E7,diff1/(diff2+eps)))
		r_l =max(-1.E7,min(1E7,diff0/(diff1+eps)))
				
		phi_r2,phi_l2  =  limiter_func(limiter,r_r,r_l)
	
		temp=(U[i,n]**2+U[i+1,n]**2)/2.
		u_avg=(U[i+1,n]+U[i,n])/2.
		#the velocities u_i+1/2 are calculated using central difference (u_avg)		
		f_p=0.5*(temp-abs(u_avg)*(1+phi_r2*((dt/dx)*abs(u_avg)-1))*(U[i+1,n]-U[i,n]))	


		temp=(U[i-1,n]**2+U[i,n]**2)/2.
		u_avg=(U[i-1,n]+U[i,n])/2.
		#the velocities u_i-1/2 are calculated using central difference (u_avg)	
		f_m=0.5*(temp-abs(u_avg)*(1+phi_l2*((dt/dx)*abs(u_avg)-1))*(U[i,n]-U[i-1,n]))
		
		F[j]=-S*(f_p-f_m)
		i=i+1
		j=j+1	
	
	return F



#calculates the net fluxes using the Roe scheme.
#Uses MUSCL schemes for velocity interpolation!
def fluxROE(mu,U,dx,F,n,S,intem,limiter,dt):

	i=2
	j=0
	eps=1E-10
	while (j <= (len(F)-1)):


		Ur_m,Ul_m=MUSCL(U,limiter,intem,i,n)
		Ur_p,Ul_p=MUSCL(U,limiter,intem,i+1,n)
		
		Am=0.5*(Ur_m+Ul_m)
		temp=(Ur_m**2+Ul_m**2)/2.
		f_m=0.5*(temp-abs(Am)*(Ur_m-Ul_m))
		
		
		Ap=0.5*(Ur_p+Ul_p)
		temp=(Ur_p**2+Ul_p**2)/2.
		f_p=0.5*(temp-abs(Ap)*(Ur_p-Ul_p))
		
		F[j]=-S*(f_p-f_m)
		i=i+1
		j=j+1	
	
	return F


def MUSCL(U,limiter,intem,i,n):
	
	
	inter=int(abs(intem))

	#second-order upwind
	if (inter==2):
		rk1=0.
		rk2=2.-rk1
	#third-order upwind
	if (inter==3):
		rk1=1.333
		rk2=2.-rk1
	#second-order upwind biased
	if (inter==4):
		rk1=1.
		rk2=2.-rk1
	#central difference
	if (inter==5):
		rk1=2.
		rk2=2.-rk1
	#QUICK SCHEME
	if (inter==6):
		rk1=1.5
		rk2=2.-rk1
	
	#first order upwind
	if (inter==1):
		rk1=0.
		rk2=0.      
	
      

      	r1      = 0.25*rk1
      	r2      = 0.25*rk2
      	eps     = 1.E-10

	#without limiters
	if (intem<=-2):
		diff0=U[i-1,n] - U[i-2,n]
		diff1=U[i,n]   - U[i-1,n]
		diff2=U[i+1,n] - U[i,n]
		Ur =U[i,n]   - (diff1*r1+diff2*r2)  
		Ul   =U[i-1,n] - (diff1*r1+diff0*r2)
		
	 
	#with limiter
	elif (inter>=2):
		diff0=U[i-1,n] - U[i-2,n]
		diff1=U[i,n]   - U[i-1,n]
		diff2=U[i+1,n]   - U[i,n]

		r_r =max(-1.E7,min(1E7,diff1/(diff2+eps)))
       		r_l =max(-1.E7,min(1E7,diff1/(diff0+eps)))	

		phi_r,phi_l  =  limiter_func(limiter,r_r,r_l)

		Ur   = U[i,n] - phi_r*(diff1*r1 + diff2*r2)
		Ul   = U[i-1,n] + phi_l*(diff1*r1+diff0*r2)


	#first order upwind
	else:
		Ur   = U[i,n]
		Ul   = U[i-1,n]
		phi_r=0
		phi_l=0
	
	return Ur,Ul

#Ghost cell update
#Two ghost cells on both ends.
def bound(mu,U,dx,F,n,S,intem,dt):
	

	
	U[1,n]=0. #U[-3,n]
	U[0,n]=0. #U[-4,n]

	U[-2,n]=0. #U[3,n]
	U[-1,n]=0. #U[2,n]


	return U



def limiter_func(limiter,r_right,r_left):

	
	if (limiter=="albada1"):
		phi_right=(r_right**2+r_right)/(r_right**2+1)
		phi_left=(r_left**2+r_left)/(r_left**2+1)


	if (limiter=="albada2"):
		phi_right=(2*r_right)/(1+r_right**2)
		phi_left=(2*r_left)/(1+r_left**2)
		
	
	if (limiter=="superbee"):

		phi_right=max(0,min(2*r_right,1),min(r_right,2))
		phi_left=max(0,min(2*r_left,1),min(r_left,2))


	if (limiter=="vanleer"):

		phi_right=(r_right+abs(r_right))/(1+r_right)
		phi_left=(r_left+abs(r_left))/(1+r_left)


	if (limiter=="monocent"):
	
		phi_right=max(0,min(2*r_right,0.5*(1+r_right),2))
		phi_left=max(0,min(2*r_left,0.5*(1+r_left),2))

	
	return phi_right,phi_left

		
main()
