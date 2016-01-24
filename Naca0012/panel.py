import math
import numpy as np
import pylab as pl
import scipy
import os

"""
USER NOTE!
WHEN RUNNING THE CODE FOR THE FIRST TIME, THE CODE SEARCHES FOR A cl.out FILE 
FROM THE SAME DIRECTORY AS THE CODE. CREATE THIS FILE MANUALLY BEFORE RUNNING THE CODE
"""



class Panel(object):

        def __init__(self,start_x,start_z,end_x,end_z,angle):

                """
                Rotate the points if alfa !=0
                """
                
                if (angle!=0.):
                    
                        theta = math.radians(-angle)
                        start_x = start_x*np.cos(theta)-start_z*np.sin(theta)
                        end_x   = end_x*np.cos(theta)-end_z*np.sin(theta)


                        start_z = start_x*np.sin(theta)+start_z*np.cos(theta)
                        end_z = end_x*np.sin(theta)+end_z*np.cos(theta)
                                              
                
                self.length = np.sqrt((end_x-start_x)**2 + (end_z - start_z)**2)
                self.xcoord = (start_x+end_x)/2.
                self.zcoord = (start_z+end_z)/2.

                self.start_x = start_x
                self.start_z = start_z
                
                self.end_x = end_x
                self.end_z = end_z

                self.sin_theta = (end_z - start_z)/self.length
                self.cos_theta = (end_x - start_x)/self.length


        
        def get_rij(self,neighbour):
                """
                returns the length of rij vector
                """

                x1 = neighbour.start_x
                x2 = self.xcoord

                z1 = neighbour.start_z
                z2 = self.zcoord

                diff_x = x2 - x1
                diff_z = z2 - z1
                
                return (np.sqrt(diff_x**2 + diff_z**2))                
                
              
       
        def get_rij_plus(self,neighbour):
                """
                returns the length of rij+1
                """ 

                x1 = neighbour.end_x
                x2 = self.xcoord

                z1 = neighbour.end_z
                z2 = self.zcoord

                diff_x = x2 - x1
                diff_z = z2 - z1

                return (np.sqrt(diff_x**2 + diff_z**2))                
                              

        def get_beta(self,neighbour):

                """
                returns the beta angle as radians
                """

                delta_y = (self.zcoord-neighbour.end_z)*(self.xcoord-neighbour.start_x)\
                        -(self.xcoord-neighbour.end_x)*(self.zcoord-neighbour.start_z)

                delta_x = (self.xcoord-neighbour.end_x)*(self.xcoord-neighbour.start_x)\
                        +(self.zcoord-neighbour.end_z)*(self.zcoord-neighbour.start_z)
               
                beta = math.atan2(delta_y, delta_x)
                 
                if (self == neighbour):
                        return math.pi

                else:
                        return beta

def main():
        
        U = 1.
        W = 0.
        panel_amount = 150
        last_digits = 12
        angle_list = np.arange(-5,16,1)

        
        
        for angle in angle_list:
                run_code(U,W,angle,panel_amount,last_digits)


def run_code(U,W,angle,panel_amount,last_digits):

        c = 1.
        alfa = 0.
        

        #panel_generator(last_digits_of_naca, chord length, number of panels, rotation angle)    
        panels = panel_generator(last_digits,c,panel_amount,angle)
        
        
        #forms the Aij matrix and calculates the values 1...N * 1...N
        mat_AIJ = matrix_Aij(panels)
        
        #calculates values A_i,n+1, A_n+1,j and A_n+1,n+1
        final_A = matrix_Aij_plus(panels,mat_AIJ)
              
        #calculates the right hand side of the system
        b = RHS(panels,U,W,alfa)
        
        #solves the system
        x = np.linalg.solve(final_A,b)

        #calculates the tangential velocities from the sigma and gamma values
        v=Velocities(panels,x,U,W,alfa)

        #output results
        post_process(v,panels,U,W,x,c,angle)


        


def post_process(v,panels,U,W,x,c,angle):

        cp = []
        xcoord = []
        circ = []
        
        for i in v:
                temp = 1-(i**2/(U**2 + W**2))           
                cp.append(temp)                       

        for i in panels:
                xcoord.append(i.xcoord)
        

        

        for i in panels:
                circ.append(i.length*x[-1])
                


        cp = np.array(cp)
        xcoord = np.array(xcoord)
        circ = np.array(circ)
        
        cl = (np.sum(circ))/(0.5*U*c)
        cl_theory = 2*np.pi*np.sin(math.radians(angle))
        
        
        print ("cl = {}, cl_t = {} angle = {}".format(cl,cl_theory,angle))

        
        
        """
        write the lift coefficients to a single file
        """
        if (angle == 0.):
                os.remove("cl.out")
                cl_file = open('cl.out', 'a')
                cl_file.write("#alfa    cl\n".format(angle,cl))
                cl_file.close()
    
        cl_file = open('cl.out', 'a')      
        cl_file.write("{0:5f} {1:5f}\n".format(angle,cl))  
        cl_file.close()


        """
        write the cp to separate files for alfa = 0,10 and 15
        """
        if (angle==0. or angle==10. or angle==15.):             
                np.savetxt('cp_{}.out'.format(angle), np.c_[xcoord,cp[1:]])
                

       
      

        


def matrix_Aij(panels):

        #Initialize the matrix
        A_ij = np.zeros((len(panels)+1,len(panels)+1))

        n=0; m=0;
        # i is now a panel object
        for i in panels:
                #loop over all neighbour panel objects                
                for j in panels:
                        
                        rij = i.get_rij(j)
                        rij_plus = i.get_rij_plus(j)
                        beta = i.get_beta(j)
                
                        #sin(0i - 0j)*ln(r_ij+/r_ij)
                        eka_termi = (j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)*   \
                                    math.log(rij_plus/rij)

                        #B*cos(0i-0j)
                        toka_termi = beta*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                                                 
                        aij = 1./(2*np.pi)*(eka_termi + toka_termi)
                     
                        A_ij[n,m] = aij
                        
                        m+=1
                
                m=0
                n+=1
                
        return A_ij
               

def matrix_Aij_plus(panels,A_ij):
        
        ###############################################################################################
                                                 #A_i,n+1
        #Calculates A_i,n+1
        n=0; m=0;

        for i in panels:               
                summa = 0.            
                for j in panels:
                        #calculates the sum of all neighbour panels                        
                        rij = i.get_rij(j)
                        rij_plus = i.get_rij_plus(j)
                        beta = i.get_beta(j)                        
                        eka_termi = math.log(rij_plus/rij)*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                        
                        toka_termi = beta*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)
                         
                        summa += (eka_termi - toka_termi)
                            

                A_ij[n,-1] = 1./(2.*np.pi)*(summa)
                n+=1
        

        ###############################################################################################
                                                 #A_n+1,j

        #Calculates the first summation (k=1) of A_n+1,j
        #k = 1
        m=0
        i = panels[0]        
        for j in panels:
                rij = i.get_rij(j)
                rij_plus = i.get_rij_plus(j)
                beta = i.get_beta(j)                        

                eka_termi = beta*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)

                toka_termi = math.log(rij_plus/rij)*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                
                temp = (eka_termi - toka_termi)
                A_ij[-1,m] = temp
                m+=1

        #Adds the second summation to the first summation (k=N)
        m=0
        i = panels[-1]       
        for j in panels:
                rij = i.get_rij(j)
                rij_plus = i.get_rij_plus(j)
                beta = i.get_beta(j)                        

                eka_termi = beta*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)

                toka_termi = math.log(rij_plus/rij)*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                
                temp = (eka_termi - toka_termi)
                A_ij[-1,m] += temp
                m+=1
                
        #Both now calculated, divide by 2*pi
        A_ij[-1,:]=A_ij[-1,:]/(2*np.pi)

        ###############################################################################################
                                                #A_n+1,n+1

        
        
        #k = 1
        #j=1...N
        summa1 = 0.      
        i = panels[0]        
        for j in panels:

                rij = i.get_rij(j)
                rij_plus = i.get_rij_plus(j)
                beta = i.get_beta(j)                        

                eka_termi = math.log(rij_plus/rij)*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)
              
                toka_termi = beta*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                         
                summa1 += (eka_termi + toka_termi)


        #k=N
        #j=1...N
        summa2 = 0.
        i = panels[-1]        
        for j in panels:

                rij = i.get_rij(j)
                rij_plus = i.get_rij_plus(j)
                beta = i.get_beta(j)                        

                eka_termi = math.log(rij_plus/rij)*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)
              
                toka_termi = beta*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                         
                summa2 += (eka_termi + toka_termi)
        

        summa = summa1+summa2

        A_np_np = summa/(2.*np.pi)
        
        A_ij[-1,-1] = A_np_np        
               


        return A_ij
                                
            
def RHS(panels,U,W,alfa):

        b = np.zeros(len(panels)+1)
        n = 0
        for i in panels:
                temp = i.sin_theta*np.cos(alfa)-i.cos_theta*np.sin(alfa)
                b[n] = np.sqrt((U**2 + W**2))*temp
                n+=1


        #b_N+1                
        #cos(theta_1 - alfa)
        temp1 = panels[0].cos_theta*np.cos(alfa) + panels[0].sin_theta*np.sin(alfa)

        #cos(theta_N - alfa)
        temp2 = panels[-1].cos_theta*np.cos(alfa) + panels[-1].sin_theta*np.sin(alfa)
        
        b[-1] = -np.sqrt(U**2 + W**2)*(temp1 + temp2)

        return b


              
def Velocities(panels,x,U,W,alfa):

        
        velocities = np.zeros(len(x))

        #n corrsepnds to i th index and m to j th.
        #new variables have to be used since i and j are now objects!
        n=0; m=0; #i,j

        
        for i in panels:
                #set both summations of i to zero
                summa1 = 0.
                summa2 = 0.
                m=0
                for j in panels:
                                      
                                 
                        rij = i.get_rij(j)
                        rij_plus = i.get_rij_plus(j)
                        beta = i.get_beta(j)                  

                        #calculates the two summations of i velocity        
                        
                        eka_termi_sigma = beta*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)
                        toka_termi_sigma = math.log(rij_plus/rij)*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)

                        #sigma(j)*(B*sin(0i - 0j) - ln(rij+/rij)*cos(0i-0j))
                        summa1 += x[m]*(eka_termi_sigma - toka_termi_sigma)


                        #ln(rij+/rij)*sin(0i-0j)+B*cos(0i-0j)
                        eka_termi_gamma = math.log(rij_plus/rij)*(j.cos_theta*i.sin_theta - i.cos_theta*j.sin_theta)
                        toka_termi_gamma = beta*(i.cos_theta*j.cos_theta + i.sin_theta*j.sin_theta)
                        summa2 += (eka_termi_gamma + toka_termi_gamma)                        

                        m+=1                                                  

                #temp = cos(0i-alfa)
                temp = panels[n].cos_theta*np.cos(alfa) + panels[n].sin_theta*np.sin(alfa)
                velocities[n] = np.sqrt(U**2 + W**2)*temp \
                                + 1./(2.*np.pi)*(summa1) \
                                + x[-1]/(2.*np.pi)*(summa2)

                                          
                                                       
                
                n+=1
                                
        
        return velocities
     
     

def panel_generator(last_digits, coord_l, num_points,angle):
        
                   
        t=(last_digits/100.)*coord_l
        c=coord_l
        x = np.arange(0.,c,c/num_points)
        
        x=np.append(x,c)
        
        #calculates the positive z-coordinates of naca00xx profile
        z_plus = 5.*t*c*(0.2969*np.sqrt(x/c)+(-0.1260)*(x/c)+(-0.3516)*(x/c)**2+(0.2843)*(x/c)**3+(-0.1015)*(x/c)**4)
        
        #fix the last index to have exactly zero value and set the lower side z-coordinates to be z_positive*(-1) 
        z_plus[-1] = 0.
        z_minus = -z_plus

       
        panel_list = []
        
              
        """
        Generates a list of panel objects. The panels are in exactly same order as in the lecture slides.
        Start of the fist panel is at x = coord length, z = 0
        
        """

        #pressure side panels (lower half)
        i=0
        while (i < len(z_minus)-1):
                
                start_x = x[len(x)-i-1]
                end_x = x[len(x)-i-2]
                start_z = z_minus[len(z_minus)-i-1]
                end_z = z_minus[len(z_minus)-i-2]                           
                panel_list.append(Panel(start_x,start_z,end_x,end_z,angle))
                i+=1

        #suction side panels (upper half)
        i=0
        while (i < len(z_plus)-2):

                start_x=x[i]
                end_x = x[i+1]

                start_z = z_plus[i]
                end_z = z_plus[i+1]

                panel_list.append(Panel(start_x,start_z,end_x,end_z,angle))
                
                i+=1

        """
        Closes the system = adds the final panel from the end of the last panel (at this point) 
        to the start of the first panel
        """
        panel_list.append(Panel(panel_list[-1].end_x,panel_list[-1].end_z,panel_list[0].start_x,panel_list[0].start_z,0.))



        """
        Rest of the lines are just for plotting

        """
        
        x_coordinaatit = []
        z_coordinaatit = []

        start_x = []
        end_x = []

        start_z = []
        end_z = []


        
        for i in panel_list:
                
                x_coordinaatit.append(i.xcoord)
                z_coordinaatit.append(i.zcoord)

                start_x.append(i.start_x)
                end_x.append(i.end_x)

                start_z.append(i.start_z)
                end_z.append(i.end_z)

        
        return panel_list










        
main()
