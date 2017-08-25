import numpy as np
import pylab as pl
import sys
from shocktubecalc import sod #analytic solution





class WENO5(object):

	
	def __init__(self, domain):
	
		
		self.domain = domain
		
		#two ghost cells required here
		self.fL = np.zeros((3,domain.nx + 4)) 
		self.fR = np.zeros((3,domain.nx + 4)) 
		
		self.FLp05 = np.zeros((3,domain.nx+4)) #fL_i+1/2
		self.FRm05 = np.zeros((3,domain.nx+4)) #fR_i-1/2
		
		#residual at all cells
		self.R = np.zeros((3,domain.nx+4))
		
		
		
		
		
		
		
		
			
		
		
	def RK3step(self,dt):
	
			
		Utemp = np.array(self.domain.U)
		
		

		#1st step
		self.computeResidual()
		self.domain.updateConservative(Utemp-dt*self.R)
		self.domain.updatePrimitive()
		self.domain.updateFluxVector()
		

		
		#2nd step
		self.computeResidual()
		self.domain.updateConservative(0.75*Utemp + 0.25*(self.domain.U - dt*self.R))
		self.domain.updatePrimitive()
		self.domain.updateFluxVector()
		
	
		
		#3rd step
		self.computeResidual()
		self.domain.updateConservative( (Utemp + 2*(self.domain.U - dt*self.R))/3. )
		self.domain.updatePrimitive()
		self.domain.updateFluxVector()
		
			
			
		
		
	def computeResidual(self):
		
		#flux-split
		self.LaxFriedrichs()
		
		#fluxes at cell walls
		self.computeFluxes()
	
	
		"""
		Boundary conditions -> no flux from either side of the domain
		"""
		self.FRm05[:,1] = self.FRm05[:,2]
		self.FLp05[:,1] = self.FLp05[:,2]
		
		self.FRm05[:,-2] = self.FRm05[:,-3]
		self.FLp05[:,-2] = self.FLp05[:,-3]
		
		
		
		for j in range(0,3):
			#range goes always from (a,b-1)
			for i in range(2,self.domain.nx+1):
					
								   				#i+1/2
				self.R[j,i] = ((self.FLp05[j,i] + self.FRm05[j,i]) - (self.FLp05[j,i-1] + self.FRm05[j,i-1]))/self.domain.dx 
			
		
		
		
		
	def computeFluxes(self):
	
		#flux-splitting for left and right states of all cells
		#self.LaxFriedrichs()
		

		#print self.fL[0,:]
		#sys.exit()

		for j in range(0,3):
			for i in range(0,self.domain.nx):
		
				ii = i + 2 #pass the ghost cells

				"""
				Computation of the fL_i+1/2
				"""

				#Read only once to create the stencils to avoid cache misses
				fL_m2 = self.fL[j,ii-2]
				fL_m1 = self.fL[j,ii-1]
				fL_ = self.fL[j,ii]
				fL_p1 = self.fL[j,ii+1]
				fL_p2 = self.fL[j,ii+2]

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
				self.FLp05[j,ii] = (alpha_0/alpha)*f0 + (alpha_1/alpha)*f1 + (alpha_2/alpha)*f2
			
				
				
			
				"""
				Computation of the fR_i-1/2
				"""
				
				#Read only once to create the stencils to avoid cache misses
				fR_m2 = self.fR[j,ii-2]
				fR_m1 = self.fR[j,ii-1]
				fR_ = self.fR[j,ii]
				fR_p1 = self.fR[j,ii+1]
				fR_p2 = self.fR[j,ii+2]
				
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
				self.FRm05[j,ii] = (alpha_0/alpha)*f0 + (alpha_1/alpha)*f1 + (alpha_2/alpha)*f2
				#print alpha_0

				#print self.FRm05[j,ii]

			

		
		
	def LaxFriedrichs(self):
		"""
		Computes the Lax-Friedrichs flux split in all cell centers
		"""
		alpha = self.domain.getMaxEigenValue()
		
		#print alpha
		for j in range(0,3):
			for i in range(0,self.domain.nx + 3):
				self.fL[j,i] = 0.5*(self.domain.F[j,i] + alpha*self.domain.U[j,i])
				self.fR[j,i] = 0.5*(self.domain.F[j,i+1] - alpha*self.domain.U[j,i+1])
				
				
				
		
				
				
	

class Domain(object):


	"""
	Initialize with primitives, solve for conservative, update primitive, update flux vector based on new primitives
	"""

	def __init__(self,nx):

		self.nx = nx
		self.gamma = 1.4
		self.dx = 1./nx
		
				          
		#primitive variables
		self.rho, self.u, self.p, self.E = self.initializePrimitive(nx+4)

		#conserved variables 3*(nx+4) list
		self.U = self.initializeConservative()

		#flux vector 3*(nx+4) list
		self.F = self.initializeFluxVector()
		
		
		
		
	def updateBCs(self, Utemp):
		#pass
		self.U[:,0]=Utemp[:,0] #bc
		self.U[:,1]=Utemp[:,1]
		self.U[:,-1]=Utemp[:,-1]
		self.U[:,-2]=Utemp[:,-2]
		
		
		
		#self.U[:,0] = self.U[:,-4]
		#self.U[:,1] = self.U[:,-3]
		
		
		
		#self.U[:,-1]= self.U[:,3]	
		#self.U[:,-2]= self.U[:,2]
		
	def updatePrimitive(self):
		"""
		Updates the primitive variables
		assumes that U is up to date
		"""
		rho = self.U[0,:]
		rhoU = self.U[1,:]
		rhoE	= self.U[2,:]
	
		self.rho = rho
		self.u = rhoU/rho
		self.E = rhoE/rho
		self.p = (self.gamma-1.)*self.rho*(self.E - 0.5*self.u*self.u)
				
		
	def getMaxEigenValue(self):
		"""
		Gets the maximum sound speed at current state of the primitive variables
		u-c, u, u+c
		"""
		c = np.sqrt(self.gamma*self.p/self.rho)
		
		lam1 = np.amax(self.u + c)
		lam2 = np.amax(self.u)
		lam3 = np.amax(self.u - c)
		
		
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
		self.F[2,:] = self.u*(self.rho*self.E + self.p)
		

		
	def initializeFluxVector(self):
		
		rho = self.U[0,:]
		rhoU = self.U[1,:]
		rhoE	= self.U[2,:]

		u = rhoU/rho
		E = rhoE/rho
		p = (self.gamma-1.)*rho*(E - 0.5*u*u)
		
		F = np.array([rho*u, rho*(u**2) + p, u*(rho*E + p)]) 
		
		return F	
		
            
	def initializeConservative(self):
    	
		E = self.p/((self.gamma-1)*self.rho) + 0.5*self.u**2 #total energy

		U = np.array([self.rho, self.rho*self.u, self.rho*E])   

		return U

	def initializePrimitive(self,nx):
		"""
		Initializes primitive variables for Sod's shock tube problem
		"""
		
		
		
		positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=0, gamma=self.gamma, npts=nx)
		
		
		rho0 = values["rho"]	
		u0 = values["u"]
		p0 = values["p"]
		
		rho0[-1] = rho0[-2]
		u0[-1] = u0[-2]
		p0[-1] = p0[-2]
		
		E0 = p0/((self.gamma-1)*rho0) + 0.5*u0**2 #total energy

		return rho0,u0,p0,E0        
          
          
	def plotRho(self,time):
	
		positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),geometry=(0., 1., 0.5), t=time, gamma=self.gamma, npts=self.nx+4)
		
		
		pl.plot(self.rho, label="WENO")
		pl.plot(values["rho"][0:-1], label="analytic")
		pl.ylim(0, 1.2)
		          
    
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
		#domain.plotRho(time)
		#pl.show()
	
	domain.plotRho(time)
	pl.show()
	


main()

