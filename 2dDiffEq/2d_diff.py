import numpy as np
import pylab as pl
import matplotlib
#from sympy import Symbol
#from sympy.solvers import nsolve
#from sympy import sin, tan

def main():
	
	FVM(5.,5.)
	FVM(8.,2.)
	FVM(2.,8.)
	summation()

#########################################################################################################

def TDMAsolver(a, b, c, d):
	nf = len(a) # number of equations
	ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy the array
	
	for it in xrange(1, nf):
		mc = ac[it]/bc[it-1]
		bc[it] = bc[it] - mc*cc[it-1]
		dc[it] = dc[it] - mc*dc[it-1]
		 
	xc = ac
	xc[-1] = dc[-1]/bc[-1]
	 
	for il in xrange(nf-2, -1, -1):
		xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
	del bc, cc, dc # delete variables from memory
	return xc	


def FVM(kx,ky):
	Lx=1.0		#length 	x-dir
	Ly=2.0		#length 	y-dir

	rho=4000.	#density
	c0=1.25E3	#heat capacity
	mu=640.		#heat transfer coefficient
	kx=kx		#conductivity x
	ky=ky		#conductivity y
	q0=80.E3	#flux at the bottom of the plate

	alfa_x=kx/(rho*c0)
	alfa_y=ky/(rho*c0)
	
	Tinf=300.	#initial temperature
	T0=200.		#boundary temperature

	dx=0.02		#stepsize x
	dy=0.02		#stepsize y
	dt=3600		#stepsize t
	time=3600*400	#simulation time
					
	nx=(Lx/dx)	
	ny=(Ly/dy)
	nt=(time/dt)
	Sx=dy
	Sy=dx
	V=dx*dy
	#initialize the temperature "domain" T(x,y,t) and set the initial field 
	T=np.zeros((nx+2,ny+2,nt+1))
	T[1:-1,1:-1,0]=Tinf
	
	F1=np.zeros((nx,ny))
	F2=np.zeros((nx,ny))
	#time loop
	
	n=0
	while (n<nt):

		#update ghost cells
		T=boundary(T,nx,ny,dt,dx,dy,Sx,Sy,alfa_x,alfa_y,mu,V,q0,kx,ky,n,Tinf,T0)

		#calculate fluxes in both directions
		F1[:,:]=flux_x(alfa_x,T,dx,n,Sx,nx,ny)		
		F2[:,:]=flux_y(alfa_y,T,dy,n,Sy,nx,ny)

		#calculate the residual
		R=(F1[:,:]+F2[:,:])*(dt/(V))
		
		#calculate the implicit temperature change
		deltaT1=delta_Ty(Sy,Sx,alfa_y,dy,dt,R,nx,ny,V)	
		deltaT2=delta_Tx(Sy,Sx,alfa_x,dx,dy,dt,deltaT1,nx,ny,V,mu,kx)
		
		#march in time		
		T[1:-1,1:-1,n+1]=T[1:-1,1:-1,n]+deltaT2[:,:]
			
		n=n+1
	
	#generate grid for contour plotting
	Y, X = np.mgrid[slice(0, Ly, dy),
                slice(0, Lx, dx)]
	
	#plot contours and output results
	plotContour(X,Y,T[1:-1,1:-1,-1],nx,ny,Lx,Ly,kx,ky)
	#save(T,Lx,Ly,dx,dy,dt,time,kx,ky,nx,ny)
	
	
	
	
def flux_x(alfa_x,T,dx,time,Sx,nx,ny):
	
	F=np.zeros((nx,ny))
	n=time
	i=0
	j=0
	i2=1
	j2=1
	while (j<ny): 
		while (i < nx):
			F[i,j]=(alfa_x*Sx/dx)*((T[i2+1,j2,n]-T[i2,j2,n])\
			-(T[i2,j2,n]-T[i2-1,j2,n]))
			i=i+1		
			i2=i2+1
		
		i2=1
		i=0
		j2=j2+1
		j=j+1
	return F
	
def flux_y(alfa_y,T,dy,time,Sy,nx,ny):

	F=np.zeros((nx,ny))
	n=time
	i=0
	j=0
	i2=1
	j2=1
	while (i<nx): 
		while (j < ny):
			F[i,j]=(alfa_y*Sy/dy)*((T[i2,j2+1,n]-T[i2,j2,n])\
			-(T[i2,j2,n]-T[i2,j2-1,n]))
			j=j+1		
			j2=j2+1
		
		j2=1
		j=0
		i2=i2+1
		i=i+1

	
	return F

def delta_Tx(Sy,Sx,alfa_x,dx,dy,dt,F,nx,ny,V,mu,kx):

	a=np.zeros(nx)
	b=np.zeros(nx)
	c=np.zeros(nx)
	deltaT=np.zeros((nx,ny))

	a1=(-Sx*alfa_x/dx)*(dt/V)
	c1=(-Sx*alfa_x/dx)*(dt/V)
	
	a[:]=a1
	a[0]=0

	b[:]=1-(a1+c1)
	b[0]=1-(a1+c1)+a1
	b[-1]=1-(a1+c1)-c1
	
	c[:]=c1
	c[-1]=0

	j=0
	while (j<ny):
		deltaT[:,j]=TDMAsolver(a,b,c,F[:,j])
		j=j+1

	return deltaT



def delta_Ty(Sy,Sx,alfa_y,dy,dt,F,nx,ny,V):
	
	a=np.zeros(ny)
	b=np.zeros(ny)
	c=np.zeros(ny)
	deltaT=np.zeros((nx,ny))
	
	a1=(-Sy*alfa_y/dy)*(dt/V)
	c1=(-Sy*alfa_y/dy)*(dt/V)
	
	a[:]=a1
	a[0]=0

	b[:]=1-(a1+c1)
	b[0]=1-(a1+c1)+a1
	b[-1]=1-(a1+c1)-a1

	c[:]=c1
	c[-1]=0

	i=0
	while (i<nx):
		deltaT[i,:]=TDMAsolver(a,b,c,F[i,:])
		i=i+1

	return deltaT






def boundary(T,nx,ny,dt,dx,dy,Sx,Sy,alfa_x,alfa_y,mu,V,q0,kx,ky,n,Tinf,T0):

	#north
	T[1:-1,-1,n]=(1./3.)*(8*T0-6.*T[1:-1,-2,n]+T[1:-1,-3,n])
	
	#east 
	Tb=(-9*T[-2,1:-1,n]+T[-3,1:-1,n]+3.*Sx*Tinf*mu)/(3.*Sx*mu-8.)
	T[-1,1:-1,n]=(1./3.)*(8*Tb[:]-6*T[-2,1:-1,n]+T[-3,1:-1,n])	

	#south
	T[1:-1,0,n]=Sy*q0/ky+T[1:-1,1,n]

	#west 
	T[0,1:-1,n]=T[1,1:-1,n]

	return T


#calculates the analytic solution and outputs the results
def summation():
	
	l_a = 1.0
	l_b =2.0
	dx = 0.02
	dy = 0.02
	kx = 2.
	ky = 8.
	dtime= 3.6e3
	dt=3600.
	t_max = 3600*1000
	conv_lim = 1e-5
	u_beg = 2e2
	u_inf = 3e2
	rho = 4.00e3
	c_heat = 1.25e3
	mu = 0.64e3
	q_in = 8.0e4
	kappa = 5.0
	k_max = 8.0
	k_min = 2.0
	imax = int(l_a/dx)
	jmax = int(l_b/dy)
	
	nx=imax
	ny=jmax
	xc=np.arange(dx/2,l_a,dx)
	yc=np.arange(dy/2,l_b,dy)
	t=np.arange(0,t_max+dt,dt)


	x_pos = np.zeros(imax+2)
	y_pos = np.zeros(jmax+2)

	ro_cp = rho*c_heat
	beta = np.sqrt(ky/kx)
	muokx = mu/kx
	qioky = q_in/ky
	numiter = 100
	numeig = 50000
	tolerance = 1e-21

	par_1 = l_a*muokx
	par_2 = -l_b*qioky
	x_pos[1:imax+1] = 0.5*(dx/l_a)*np.array(range(1,2*imax,2))
	y_pos[1:jmax+1] = 0.5*(dy/l_b)*np.array(range(1,2*jmax,2))
	x_pos[0]=0.0
	x_pos[imax+1]=2*imax*0.5*(dx/l_a)
	y_pos[0]=0.0
	y_pos[jmax+1]=2*jmax*0.5*(dy/l_b)
	u_gap = u_beg - u_inf
	exact  = u_inf*np.ones((imax+2,jmax+2))



	for n2 in range(numeig,-1,-1):
	    temps = np.zeros((imax+2,jmax+2))
	    nc = 2*n2
	    for n1 in range(2,0,-1):
		n = nc+n1
		f0 = n*np.pi
		ff = f0 -np.pi
		
		for iter in range(1,numiter+1):
		    f1 = np.arctan(par_1/f0) + ff
		    if abs(f0/f1 -1.0) < tolerance:
		        break
		    f0 = f1
		f0 = par_1/f1
		
		by_pi = np.arctan(f0) + (n1 -1.0)*np.pi
		eigen = f1/beta/l_a
		coeff = 2.0*np.sin(by_pi)/(f1+f0/(1.0+f0*f0))
		f0 = -eigen*l_b
		f1 = par_2/f0
		ff = coeff/(1.0 + np.exp(2.0*f0))
		
		x_var_0= ff*np.cos(by_pi*x_pos)
		y_var_1=np.exp(f0*y_pos)
		y_var_2=np.exp(f0*(1.0-y_pos))
		
		
		for iter in range(0,jmax+2):
		    temps[0:imax+2, iter] = temps[0:imax+2,iter] + x_var_0[0:imax+2]*(u_gap*y_var_2[iter]*(1.0+y_var_1[iter]**2)+f1*y_var_1[iter]*(1.0-y_var_2[iter]**2))
	    exact = exact + temps        
		        
		    
	
	x_pos = x_pos*l_a
	y_pos = y_pos*l_b

	T=exact
	T1=T[1:-1,ny/2]
	T2=T[1:-1,ny/4]
	T3=T[nx/2,1:-1]
	T4=T[nx/4,1:-1]
	np.savetxt('analytic_k3_y1.out', (xc,T1))
	np.savetxt('analytic_k3_y05.out', (xc,T2))
	np.savetxt('analytic_k3_x025.out', (yc,T3))
	np.savetxt('analytic_k3_x05.out', (yc,T4))


	pl.plot(xc,T1)
	pl.show()
	pl.plot(xc,T2)
	pl.show()
	pl.plot(yc,T3)
	pl.show()
	pl.plot(yc,T4)
	pl.show()

def save(T,Lx,Ly,dx,dy,dt,time,kx,ky,nx,ny):

	#cell centre coordinates and time values for post processing
	xc=np.arange(dx/2,Lx,dx)
	yc=np.arange(dy/2,Ly,dy)
	t=np.arange(0,time+dt,dt)


	T1=T[1:-1,ny/2,-1]
	T2=T[1:-1,ny/4,-1]
	T3=T[nx/2,1:-1,-1]
	T4=T[nx/4,1:-1,-1]

	if (kx==5.):
		np.savetxt('k1_y1.out', (xc,T1))
		np.savetxt('k1_y05.out', (xc,T2))
		np.savetxt('k1_x025.out', (yc,T3))
		np.savetxt('k1_x05.out', (yc,T4))

	elif (kx==8.):
		np.savetxt('k2_y1.out', (xc,T1))
		np.savetxt('k2_y05.out', (xc,T2))
		np.savetxt('k2_x025.out', (yc,T3))
		np.savetxt('k2_x05.out', (yc,T4))

	elif (kx==2.):
		np.savetxt('k3_y1.out', (xc,T1))
		np.savetxt('k3_y05.out', (xc,T2))
		np.savetxt('k3_x025.out', (yc,T3))
		np.savetxt('k3_x05.out', (yc,T4))

	

def plotContour(X,Y,u,nx,ny,Lx,Ly,kx,ky):
		
	pl.pcolor(X, Y, np.transpose(u), cmap='RdBu')
	pl.colorbar()
	pl.title("kx = {}, ky = {}".format(kx,ky))
	pl.axis([X.min(), X.max(), Y.min(), Y.max()])
	pl.savefig("conotour{}.png".format(int(kx)))		
	pl.show()



main()
