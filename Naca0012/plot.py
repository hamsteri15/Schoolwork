import numpy as np
import pylab as pl
import math

"""
CL VALUES
"""
cl_data = np.loadtxt("cl.out")
cl_measurements = np.loadtxt("cl_experimental.dat")

alfa_oma = cl_data[:,0]
cl_oma = cl_data[:,1]



alfa_measurements = cl_measurements[:,0]
cl_measurements = cl_measurements[:,1]

cl_white = []
alfa_white = np.arange(-5,15,0.1)
for alfas in alfa_white:
        cl_white.append(2*np.pi*np.sin(math.radians(alfas)))
       
cl_white = np.array(cl_white) 



pl.scatter(alfa_oma,cl_oma,label="simulated",color="black")
pl.plot(alfa_white,cl_white,label="$2\pi sin(\\alpha)$",color="black")
pl.scatter(alfa_measurements,cl_measurements,label="measured,NASA TM 4074",color="black",marker="x")

pl.legend(loc="best")
pl.xlabel("$\\alpha$",size=20)
pl.ylabel("$C_L$",size=20)
pl.xlim(-5,22)
pl.ylim(-0.5,2)
pl.grid(True)
pl.savefig("images/cl.png")
pl.show()


"""
CP VALUES  alfa = 0
"""

cp_data1 = np.loadtxt("cp_0.out")
cp_measured1 = np.loadtxt("cp_exp0.out")
x_oma = cp_data1[:,0]
cp_oma = cp_data1[:,1]

x_exp = cp_measured1[:,0]
cp_exp = cp_measured1[:,1]

pl.plot(x_oma,cp_oma,label="simulated",color="black")
pl.scatter(x_exp,cp_exp,label="measured, NASA R&M 3726",marker="x",color="black")


pl.ylim(1.,-0.6)
pl.xlim(0,1)
pl.xlabel("$x/c$",size=20)
pl.ylabel("$C_p$",size=20)
pl.grid(True)
pl.title("$\\alpha = 0$",size="16")
pl.legend(loc="best")
pl.savefig("images/cp_0.png")
pl.show()


"""
CP VALUES  alfa = 10
"""

cp_data1 = np.loadtxt("cp_10.out")
cp_measured1 = np.loadtxt("cp_exp10.out")
x_oma = cp_data1[:,0]
cp_oma = cp_data1[:,1]

x_exp = cp_measured1[:,0]
cp_exp = cp_measured1[:,1]

pl.plot(x_oma,cp_oma,label="simulated",color="black")
pl.scatter(x_exp,cp_exp,label="measured, NASA R&M 3726",marker="x",color="black")
pl.title("$\\alpha = 10$",size=16)


pl.ylim(2.,-6)
pl.xlim(0,1)
pl.xlabel("$x/c$",size=20)
pl.ylabel("$C_p$",size=20)
pl.grid(True)
pl.legend(loc="best")
pl.savefig("images/cp_10.png")
pl.show()


"""
CP VALUES  alfa = 15
"""

cp_data1 = np.loadtxt("cp_15.out")
cp_measured1 = np.loadtxt("cp_exp15.out")
x_oma = cp_data1[:,0]
cp_oma = cp_data1[:,1]

x_exp = cp_measured1[:,0]
cp_exp = cp_measured1[:,1]

pl.plot(x_oma,cp_oma,label="simulated",color="black")
pl.scatter(x_exp,cp_exp,label="measured, NASA R&M 3726",marker="x",color="black")
pl.title("$\\alpha = 15$",size=16)


pl.ylim(2.,-12)
pl.xlim(0,1)
pl.xlabel("$x/c$",size=20)
pl.ylabel("$C_p$",size=20)
pl.grid(True)
pl.legend(loc="best")
pl.savefig("images/cp_15.png")
pl.show()






















