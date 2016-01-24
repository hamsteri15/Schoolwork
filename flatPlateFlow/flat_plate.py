#																																			
#					
#            ##	       #		
#	    #  #       #	Finite difference solver for laminar boundary layer flow		
#	   #    #      #
#	  # #### #     #
#	 #	  #	   	
#	#	   #   #  Aalto University						
################################################################################################################################################
import numpy as np
import pylab as pl
import matplotlib
from scipy import meshgrid


 

def main():
	
	nu = 1e-5 		#kinematic viscosity
	U = 5 	 		#free stream velocity
	dU = 0			#Free stream acceleration, set 1 if you want that U(x)=U*(1-x/L) or 
				#2 if you want that U(x)=U*(1+x/L). Otherwise U(x)=U
 

	L = 0.6		        #length of the domain
	H = 0.05                #height of the domain
	nx =100             	#number of nodes in x-direction
	ny =100            	#number of nodes in y-direction

	impsi=1			#Implicit or explicit approach. See White ch. 4-7. 
				#1=implicit #2=explicit
				#If you use the explicit approach, pay more attention to 
				#step sizes and follow the stability criteria. If you have no idea what you are
				#doing, use the implicit approach!
			
	



	plot = True	        #True if you wish to see u-velocity profiles. The velocity
				#profiles are plotted from three different locations (x1,x2,x3)
	
	plot2 = False		#True if you wish to see velocity contours. For small step sizes
				#this may take a while!!
		
	output = True		#output raw data of the velocity profiles into profiles.dat file



	x1 = 3			#First _NODE_ for the u-velocity profile. 0 <= x1 <= nx-1 ! 0 is the 					#beginning
				#of the plate and nx-1 is the end of the
				#plate.						
	x2 = 5             	#Second _NODE_ for the u-velocity profile.		
	x3 = 7              	#Third _NODE_ for the u-velocity profile.
				#The coordinate of the node can be calculated from x_coor = xn * (L/nx)	

	

	if impsi==2:
		explicit(U,dU,nu,L,H,nx,ny,plot,plot2,impsi,output,x1,x2,x3)
	if impsi ==1:
		implicit(U,dU,nu,L,H,nx,ny,plot,plot2,impsi,output,x1,x2,x3)

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


def plotContour(X,Y,u,nx,ny,impsi,L,H):

	pl.pcolor(X, Y, np.transpose(u), cmap='RdBu')
	pl.colorbar()
	pl.title("u-velocity contour")
	pl.axis([X.min(), X.max(), Y.min(), Y.max()])		
	pl.show()
			

def outputData(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3):
	f = open('profiles.dat','w')
	
	U_plot=u[x1,0:ny]
	U_plot2=u[x2,0:ny]
	U_plot3=u[x3,0:ny]
	y=Y[:,0]

	f.write('#y coordinate	u-velocity\n')
	f.write('################### x1 #################\n')	
	i=0
	while (i<len(U_plot)):	
		f.write('%.5f	%.5f\n' % (y[i],U_plot[i]))
		i=i+1
	f.write('################### x2 ##################\n')
	i=0
	while (i<len(U_plot2)):	
		f.write('%.5f	%.5f\n' % (y[i],U_plot2[i]))
		i=i+1  
	f.write('################### x3 ##################\n')
	i=0
	while (i<len(U_plot3)):	
		f.write('%.5f	%.5f\n' % (y[i],U_plot3[i]))
		i=i+1 
	f.close() 

def plotUprofile(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3):
		U_plot=u[x1, 0:ny]
		U_plot2=u[x2, 0:ny]
		U_plot3=u[x3, 0:ny]
		y=Y[:,0]

		pl.plot(U_plot,y,color="black",label="x1")
		pl.plot(U_plot2,y,color="red",label="x2")
		pl.plot(U_plot3,y,color="blue",label="x3")
		pl.xlim(0,u.max()+1)
		pl.ylim(0,Y.max()/10)
		pl.xlabel("U")
		pl.ylabel("y")
		pl.title("u-velocity profiles")
		pl.legend(loc="upper left")
		
		pl.show()
		
	
def implicit(U,dU,nu,L,H,nx,ny,plot,plot2,impsi,output,x1,x2,x3):
	dx = L/(nx)
	dy = H/(ny)
	#create a grid for post processing
	Y, X = np.mgrid[slice(0, H, dy),
                slice(0, L, dx)]

	u=np.zeros((nx,ny)) 
	v=np.zeros((nx,ny))
	
	#0 is the first index in python --> ny-1 is the last index in y-direction	
	v[:,0]=0 
	v[0,:]=0
	u[:,0]=0 #no slip
	u[0,0:]=U #inlet 
		
	
	u[0:nx,ny-1]=U #free stream
			
	#vectors for the TDMA algorithm a_i*x_(i-1) + b_i*x_i + c_i*x_(i+1) = d_i		
	a=np.zeros(ny-2)
	b=np.zeros(ny-2)
	c=np.zeros(ny-2)
	d=np.zeros(ny-2)
	
	#free stream acceleration/deceleration
	U2=np.zeros(nx)	
	i=0
	if dU==1:
		while (i<nx-1):
			U2[i]=U*(1-(i*dx)/L)
			i=i+1
	elif dU==2:	
		while (i<nx-1):
			U2[i]=U*(1+(i*dx)/L)
			i=i+1
	
	i=0
	m=0
	n=1
	
	#loop x-nodes to the end and y-nodes from second to second to last 
	while (m<(nx-1)):
		while (n<=(ny-2)):
			
			#second node (first one after the known no-slip plate)  
			if n==1:
				alfa=(nu*dx)/(u[m,n]*dy*dy)		
				beta=(v[m,n]*dx)/(2*u[m,n]*dy)
				a[i]=0	
				b[i]=(1+2*alfa)
				c[i]=-alfa
				d[i]=u[m,n]-(beta*(u[m,n+1]-u[m,n-1]))+((U2[m+1]-U2[m])/(2*u[m,n]))
				
			#middle nodes
			if (n > 1 and n < (ny-2)):
				alfa=(nu*dx)/(u[m,n]*dy*dy)
				beta=(v[m,n]*dx)/(2*u[m,n]*dy)
				a[i]=-alfa
				b[i]=(1+2*alfa)
				c[i]=-alfa
				d[i]=u[m,n]-(beta*(u[m,n+1]-u[m,n-1])) +((U2[m+1]-U2[m])/(2*u[m,n]))
			
			#second to last node (last one before the known free stream node)
   			if (n==(ny-2)):
    			    	alfa=(nu*dx)/(u[m,n]*dy*dy)
			     	beta=(v[m,n]*dx)/(2*u[m,n]*dy)
			     	a[i]=-alfa	
		             	b[i]=(1+2*alfa)
			     	c[i]=0		
		             	d[i]=u[m,n]-(beta*(u[m,n+1]-u[m,n-1]))+(alfa*U)+((U2[m+1]-U2[m])/(2*u[m,n]))			
			
			i=i+1
			n=n+1
		
		#march in x-direction.			
		u[m+1,1:ny-1]=TDMAsolver(a,b,c,d)
		j=1
		while j<ny-1:
        		v[m+1,j] = v[m+1,j-1] - (dy/(2*dx))*(u[m+1,j]-u[m,j]+u[m+1,j-1]-u[m,j-1])
			j=j+1
		j=0
		m=m+1
		n=1
		i=0
	
	
	
	if (plot==True):
		plotUprofile(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3)
	if (plot2==True):
		plotContour(X,Y,u,nx,ny,impsi,L,H)

	if (output==True):	
		outputData(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3)


def explicit(U,dU,nu,L,H,nx,ny,plot,plot2,impsi,output,x1,x2,x3):
	dx = L/(nx)
	dy = H/(ny)
	#create a grid for post processing
	Y, X = np.mgrid[slice(0, H, dy),
                slice(0, L, dx)]

	u=np.zeros((nx,ny)) 
	v=np.zeros((nx,ny))
	
	#0 is the first index in python --> ny-1 is the last index in y-direction	
	v[:,0]=0 
	v[0,:]=0
	u[:,0]=0 #no slip	
	u[0,0:]=U #inlet 
	u[0:nx,ny-1]=U #free stream

	
	#free stream acceleration/deceleration
	U2=np.zeros(nx)
	i=0
	if dU==1:
		while (i<nx-1):
			U2[i]=U*(1-(i*dx)/L)
			i=i+1
	elif dU==2:	
		while (i<nx-1):
			U2[i]=U*(1+(i*dx)/L)
			i=i+1
	i=0
	
	n=1
	m=0
	
	#loop x-nodes to the end and y-nodes from second to second to last 
	while (m<(nx-1)):
		while (n<=(ny-2)):
			alfa=(nu*dx)/(u[m,n]*dy*dy)
			beta=(v[m,n]*dx)/(2*u[m,n]*dy)
			u[m+1,n]=((alfa-beta)*u[m,n+1])+((1-(2*alfa))*u[m,n])+((alfa+beta)*u[m,n-1])+((U2[m+1]-U2[m])/(2*u[m,n]))	
			v[m+1,n]=v[m+1,n-1]-(dy/(2*dx))*(u[m+1,n]-u[m,n]+u[m+1,n-1]-u[m,n-1])		
			
			n=n+1
		stabx=(u[m+1,1:].min()*dy**2)/(2*nu)
		staby=(2*nu)/(abs(v[m+1,1:].max()))
		
		print ("%.8f  	%.8f x-stability || y-stability\n" % (stabx, staby))
		#For some reason, there seems to be unstability at the beginning of the plate. Thats why
		#m>10
		if ((stabx<dx or staby<dy) and m>10):
			print "step size error"
			#return 0
		
		
		m=m+1
		n=1
		
	


		
	if (plot==True):		
		plotUprofile(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3)
	if (plot2==True):
		plotContour(X,Y,u,nx,ny,impsi,L,H)
	if (output==True):	
		outputData(X,Y,u,nx,ny,impsi,L,H,x1,x2,x3)

main()
