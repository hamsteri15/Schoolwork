import numpy as np
import pylab as pl
import sys
from shocktubecalc import sod #analytic solution


##############################
# LaxFriedrichs(self) == OK
# computeFluxes_F(self) == OK
##############################



class WENO5(object):

	
	def __init__(self, domain):
	
		
		self.domain = domain
		
		#two ghost cells required here
		self.fL = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) 
		self.fR = np.zeros((4,(domain.nx + 4) * (domain.ny + 4)))
		
		self.gL = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) 
		self.gR = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) 
		
		self.FLp05 = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) #fL_i+1/2
		self.FRm05 = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) #fR_i-1/2

		self.GLp05 = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) #gL_i+1/2
		self.GRm05 = np.zeros((4,(domain.nx + 4) * (domain.ny + 4))) #gR_i-1/2
		
		#residual at all cells
		self.dF = np.zeros((4,(domain.nx + 4) * (domain.ny + 4)))
		self.dG = np.zeros((4,(domain.nx + 4) * (domain.ny + 4)))
		
		self.R = np.zeros((4,(domain.nx + 4) * (domain.ny + 4)))
		
		
		
		#self.LaxFriedrichs()
		
		#print self.gL[0, self.domain.nx:self.domain.nx*2]
		
		
			
		
		
	def RK3step(self,dt):
	
			
			Utemp = np.array(self.domain.U)
			

			



			#1st step
			#self.LaxFriedrichs()
			#self.computedF()
			self.computeResidual()
			self.domain.updateConservative(Utemp-dt*self.R)
			self.domain.updatePrimitive()
			self.domain.updateFluxVector()
			self.domain.updateBCs()
			
			#print self.R[0,:].reshape(self.domain.nx + 4 , self.domain.ny + 4)

			
			#sys.exit()

			
			#2nd step
			#self.LaxFriedrichs()
			#self.computedF()
			self.computeResidual()
			self.domain.updateConservative(0.75*Utemp + 0.25*(self.domain.U - dt*self.R))
			self.domain.updatePrimitive()
			self.domain.updateFluxVector()
			self.domain.updateBCs()
		
			
			#3rd step
			#self.LaxFriedrichs()
			#self.computedF()
			self.computeResidual()
			self.domain.updateConservative( (Utemp + 2*(self.domain.U - dt*self.R))/3. )
			self.domain.updatePrimitive()
			self.domain.updateFluxVector()
			self.domain.updateBCs()

			
		

	def printArray(self, arr, nx, ny):

		for j in range(0,ny):
			for i in range(0, nx):
				print "{} ".format(arr[i + nx*j])

			print " "


	def computeResidual(self):

		self.LaxFriedrichs()
		self.computedF()


		nx = self.domain.nx + 4
		ny = self.domain.ny + 4
		for k in range(0,4):
			for jj in range(0, self.domain.ny):
				for ii in range(0, self.domain.nx):

					i = ii + 2
					j = jj + 2

					self.R[k, i + nx*j] = self.dF[k, i + nx*j] 
			
			
			
		
		
	def updateBC(self):
	
		"""
		Boundary conditions -> no flux from either side of the domain
		
		self.FRm05[:,1] = self.FRm05[:,2]
		self.FLp05[:,1] = self.FLp05[:,2]
		
		self.FRm05[:,-2] = self.FRm05[:,-3]
		self.FLp05[:,-2] = self.FLp05[:,-3]
		"""
		
		
		nx = self.domain.nx + 4
		ny = self.domain.ny + 4
		for k in range(0,4):
			for j in range(0, self.domain.ny):
			
				#WEST
				i = 1
				self.FRm05[k, i + nx * j] =  self.FRm05[k, (i + 1) + nx * j]
				self.FLp05[k, i + nx * j] =  self.FLp05[k, (i + 1) + nx * j]
				
		
				#EAST
				i = self.domain.nx
				self.FRm05[k, (i + 1) + nx * j] =  self.FRm05[k, (i - 1) + nx * j]
				self.FLp05[k, (i + 1) + nx * j] =  self.FLp05[k, (i - 1) + nx * j]
				
		
	def computedF(self):
		
		
		
		#fluxes at cell walls
		self.computeFluxes_F()
		
		
		self.updateBC()
		
		nx = self.domain.nx + 4
		ny = self.domain.ny + 4
	
		

		for k in range(0,4):
			for jj in range(0,self.domain.ny):
				for ii in range(0,self.domain.nx):
					
					i = ii + 2
					j = jj + 2
								 
					iii = i + 1
  				
					

					self.dF[k,i + nx*j] = (	self.FLp05[k,(i + 0) + nx*j] \
											+ self.FRm05[k,(i + 0)+ nx*j] \
											- self.FLp05[k,(i - 1) + nx*j] \
											- self.FRm05[k,(i - 1) + nx*j])/self.domain.dx 
			
		
					
		
	def computeFluxes_F(self):
	
		#flux-splitting for left and right states of all cells
		#self.LaxFriedrichs()
		
		nx = self.domain.nx + 4
		ny = self.domain.ny + 4
		for k in range(0, 4):
			for j in range(0, self.domain.ny):
				for i in range(0, self.domain.nx):
		
					#pass the ghost cells
					ii = i + 2 
					jj = j + 2

					"""
					Computation of the fL_i+1/2
					"""

					#Read only once to create the stencils to avoid cache misses
					fL_m2 = self.fL[k,(ii-2) + nx*jj]
					fL_m1 = self.fL[k,(ii-1) + nx*jj]
					fL_ = self.fL[k,ii + nx*jj]
					fL_p1 = self.fL[k,(ii+1) + nx*jj]
					fL_p2 = self.fL[k,(ii+2) + nx*jj]

					#Eno stencils
					f0 = (2*fL_m2 - 7*fL_m1 + 11*fL_)/6.
					f1 = (-fL_m1 + 5*fL_ + 2*fL_p1)/6.
					f2 = (2*fL_ + 5*fL_p1 - fL_p2)/6.
				
				
				
					#Smoothness indicators
					b0 = (13./12.) * (fL_m2 - 2*fL_m1 + fL_)**2 + (1./4.) * (fL_m2 - 4*fL_m1 + 3*fL_)**2 
					b1 = (13./12.) * (fL_m1 - 2*fL_ + fL_p1)**2 + (1./4.) * (fL_m1 - fL_p1)**2
					b2 = (13./12.) * (fL_ - 2*fL_p1 + fL_p2)**2 + (1./4.) * (3*fL_ - 4*fL_p1 + fL_p2)**2
				
				
					
					#Weights
					epsilon = 1E-6
				
					alpha_0 = 0.1/(epsilon + b0)**2
					alpha_1 = 0.6/(epsilon + b1)**2
					alpha_2 = 0.3/(epsilon + b2)**2
				
					alpha = alpha_0 + alpha_1 + alpha_2

					
				
					#Flux at fL_i+1/2
					self.FLp05[k, ii + nx*jj] = (alpha_0/alpha)*f0 + (alpha_1/alpha)*f1 + (alpha_2/alpha)*f2
			
				
				
			
					"""
					Computation of the fR_i-1/2
					"""
				
					#Read only once to create the stencils to avoid cache misses
					fR_m2 = self.fR[k,(ii-2) + nx*jj]
					fR_m1 = self.fR[k,(ii-1) + nx*jj]
					fR_ = self.fR[k,ii + nx*jj]
					fR_p1 = self.fR[k,(ii+1) + nx*jj]
					fR_p2 = self.fR[k,(ii+2) + nx*jj]
				
					#Eno stencils
					f0 = (-fR_m2 + 5*fR_m1 + 2*fR_)/6.
					f1 = (2*fR_m1 + 5*fR_ - fR_p1)/6.
					f2 = (11*fR_ - 7*fR_p1 + 2*fR_p2)/6.
				
				
					
				
					#Smoothness indicators
					b0 = (13./12.) * (fR_m2 - 2*fR_m1 + fR_)**2 + (1./4.) * (fR_m2 - 4*fR_m1 + 3*fR_)**2 
					b1 = (13./12.) * (fR_m1 - 2*fR_ + fR_p1)**2 + (1./4.) * (fR_m1 - fR_p1)**2
					b2 = (13./12.) * (fR_ - 2*fR_p1 + fR_p2)**2 + (1./4.) * (3*fR_ -4*fR_p1 + fR_p2)**2
				
						

					#Weights
					epsilon = 1E-6
				
					alpha_0 = 0.3/(epsilon + b0)**2
					alpha_1 = 0.6/(epsilon + b1)**2
					alpha_2 = 0.1/(epsilon + b2)**2
				
					alpha = alpha_0 + alpha_1 + alpha_2
		
		
					#Flux at fR_i-1/2
					self.FRm05[k, ii+nx*jj] = (alpha_0/alpha)*f0 + (alpha_1/alpha)*f1 + (alpha_2/alpha)*f2

					
		
	def LaxFriedrichs(self):
		"""
		Computes the Lax-Friedrichs flux split in all cell centers
		"""
		alpha = self.domain.getMaxEigenValue()
		
		
	
		nx = self.domain.nx + 4
		ny = self.domain.ny + 4
		
		

		for k in range(0,4):
			for jj in range(0,ny):
				
				for ii in range(0,nx - 1):
				
					i = ii
					j = jj
										

					self.fL[k, i + nx*j] = 0.5*(self.domain.F[k, i + nx*j] + alpha*self.domain.U[k, i + nx*j])
					self.fR[k, i + nx*j] = 0.5*(self.domain.F[k, (i + 1) + nx*j] - alpha*self.domain.U[k, (i + 1) + nx*j])



					
					#self.gL[k, i + nx*j] = 0.5*(self.domain.G[k, i + nx*j] + alpha*self.domain.U[k, i + nx*j])
					#self.gR[k, i + nx*j] = 0.5*(self.domain.G[k, i + nx*(j + 1)] - alpha*self.domain.U[k, i + nx*(j + 1)]) 




		
		
	
	

class Domain(object):


	"""
	Initialize with primitives, solve for conservative, update primitive, update flux vector based on new primitives
	"""

	def __init__(self,nx):
	
		self.L = 1.

		self.nx = nx
		self.ny = nx
		self.gamma = 1.4
		self.dx = self.L/nx
		self.dy = self.L/nx


				          
		#primitive variables
		self.rho, self.u, self.v, self.p, self.E = self.initializePrimitive(self.nx+4, self.ny+4, False)

		#conserved variables 4*(nx+4*ny+4) list
		self.U = self.initializeConservative()

		#flux vectors 3*(nx+4) list
		self.F, self.G = self.initializeFluxVector()
		
		
		
		
	def updateBCs(self):
		pass
		#self.U[:,0]=Utemp[:,0] #bc
		#self.U[:,-1]=Utemp[:,-1] #bc	
		
	def updatePrimitive(self):
		"""
		Updates the primitive variables
		assumes that U is up to date
		"""
		rho = self.U[0,:]
		rhoU = self.U[1,:]
		rhoV = self.U[2,:]
		rhoE	= self.U[3,:]
	
		self.rho = rho
		self.u = rhoU/rho
		self.v = rhoV/rho
		self.E = rhoE/rho
		self.p = (self.gamma-1.)*self.rho*(self.E - 0.5*(self.u**2 + self.v**2))
		
				
		
	def getMaxEigenValue(self):
		"""
		Gets the maximum sound speed at current state of the primitive variables
		q-c, q, q+c, where q = u + v
		"""
		c = np.sqrt(self.gamma*self.p/self.rho)
		

		q = self.u + self.v


		lam1 = np.amax(q + c)
		lam2 = np.amax(q)
		lam3 = np.amax(q - c)
		
		
		alpha = max(abs(lam1), abs(lam2), abs(lam3))
		
		
		return alpha
		
	def getFluxVector(self):
	
		return self.F
		
	def getConservative(self):
		
		return self.U		
	
	def updateConservative(self,new):
		self.U[:,:] = new
	
		
	def updateFluxVector(self):
		"""
		Updates the flux vector based on the current primitive variables
		"""
		
		self.F[0,:] = self.rho*self.u
		self.F[1,:] = self.rho*(self.u**2) + self.p
		self.F[2,:] = self.rho*self.u*self.v
		self.F[3,:] = self.u*(self.rho*self.E + self.p)
		
		self.G[0,:] = self.rho*self.v
		self.G[1,:] = self.rho*self.v*self.u + self.p
		self.G[2,:] = self.rho*self.v**2 + self.p
		self.G[3,:] = self.v*(self.rho*self.E + self.p)


		

		
	def initializeFluxVector(self):
		
		rho = self.U[0,:]
		rhoU = self.U[1,:]
		rhoV = self.U[2,:]
		rhoE	= self.U[3,:]



		u = rhoU/rho
		v = rhoU/rho
		E = rhoE/rho
		p = (self.gamma-1.)*rho*(E - 0.5*(u**2 + v**2))

		
		F = np.array([rho*u, rho*(u**2) + p, rho*u*v, u*(rho*E + p)]) 
		G = np.array([rho*v, rho*v*u + p, rho*v**2 + p, v*(rho*E + p)]) 

		
		return F,G	
		
            
	def initializeConservative(self):
    	
		
		U = np.array([self.rho, self.rho*self.u, self.rho*self.v, self.rho*self.E])   

		return U

	def initializePrimitive(self,nx,ny, plot=False):
		"""
		Initializes primitive variables for Sod's shock tube problem
		"""
		
		rho0 = np.zeros(nx*ny)
		u0 = np.zeros(nx*ny)
		p0 = np.zeros(nx*ny)
		v0 = np.zeros(nx*ny)

		x = np.linspace(0,self.L,nx); 
		y = np.linspace(0,self.L,nx)		


		self.x = x
		self.y = y

		for j in range(ny):
			for i in range(nx):
				
				if (x[i] <= self.L/2.):
			
					p0[i+nx*j] = 1.
					rho0[i+nx*j] = 1.
				else:
					p0[i+nx*j] = 0.1
					rho0[i+nx*j] = 0.125	
			
		if plot:		
		
			X,Y = np.meshgrid(x,y)
	
			rho_t = np.array(rho0).reshape(nx,ny)

			pl.pcolor(X, Y, rho_t, cmap='RdBu')
			pl.colorbar()
			pl.title("rho")
			pl.axis([X.min(), X.max(), Y.min(), Y.max()])
			pl.show()
		
			
			p_t = np.array(p0).reshape(nx,ny)

			pl.pcolor(X, Y, p_t, cmap='RdBu')
			pl.colorbar()
			pl.title("p")
			pl.axis([X.min(), X.max(), Y.min(), Y.max()])
			pl.show()
			



		
		#positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=0, gamma=self.gamma, npts=nx)
		
			
		
		
		E0 = p0/((self.gamma-1)*rho0) + 0.5*(u0**2 + v0**2) #total energy

		return rho0,u0,v0,p0,E0


	def plotRho1D(self, time):
        
		positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=time, gamma=self.gamma, npts=self.nx+4)

		rho_c = np.array(self.rho).reshape(self.nx + 4 , self.ny + 4)

		#print rho_c

		pl.plot(rho_c[3,:], label="WENO")
		pl.plot(values["rho"][0:-1], label="analytic")
          
	def plotRho2D(self,time):
	
		X,Y = np.meshgrid(self.x,self.y)
	
		rho_t = np.array(self.rho).reshape(self.nx+4,self.ny+4)

		pl.pcolor(X, Y, rho_t, cmap='RdBu')
		pl.colorbar()
		pl.title("rho")
		pl.axis([X.min(), X.max(), Y.min(), Y.max()])
		#pl.show()

		
		#positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=time, gamma=self.gamma, npts=self.nx+4)
		
		
		#rho_t = np.array(self.rho).reshape(nx,ny)

		#pl.plot(self.rho[5*self.nx:6*self.nx], label="WENO")
		#pl.plot(values["rho"][0:-1], label="analytic")
		          
    
	def plotAll(self,time):
	
		positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=time, gamma=self.gamma, npts=self.nx+4)
		
		
		pl.plot(self.rho, label="WENO-rho")
		pl.plot(values["rho"][0:-1], label="analytic-rho")          

		pl.plot(self.u, label="WENO-u")
		pl.plot(values["u"][0:-1], label="analytic-u")
		
		#pl.plot(self.p, label="WENO-E")
		#pl.plot(values["p"][0:-1], label="analytic-u")
		
		pl.legend(loc="best")

def main():

	domain = Domain(40)

	

	solver = WENO5(domain)
	
	

	dt = 0.001
	T = 0.08
	
	time = 0
	while (time<T):
		
		solver.RK3step(dt)
		time+=dt
		print time
	
	domain.plotRho1D(time)
	pl.show()
	


main()

