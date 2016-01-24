import unittest
from kartta import Kartta
from kartta import Seina
from yksilo import Yksilo
import os,sys,pygame,pygame.mixer,random
from pygame.locals import *


class Test(unittest.TestCase):


    def setUp(self):
        self.kartta = Kartta()

        
        
    
    def test_alusta_kartta(self): 
        """
        Testaa Kartta-luokan metodia alusta_kartta()
        """     

        seinien_maara = 0
   
        self.kartta.alusta_kartta(seinien_maara)
	
        seinien_maara2 = len(self.kartta.seinat)
	
        self.assertNotEqual(seinien_maara,seinien_maara2,"Seina-objekteja ei ole lisatty kartalle")
    


    def test_luo_valiseinat(self):
        """
        Testaa Kartta-luokan metodia luo_valiseinat()
        """  

        self.kartta.luo_valiseinat(0)
        temp = len(self.kartta.seinat)
        self.assertEqual(0,temp,"Virhe valiseinien luonnissa kun 0 valiseinaa")

        self.kartta.luo_valiseinat(1)
        temp2 = len(self.kartta.seinat)
        self.assertTrue(temp < temp2,"Virhe valiseinien luonnissa kun 1 valiseinaa")

        self.kartta.luo_valiseinat(2)
        temp3= len(self.kartta.seinat)
        self.assertTrue(temp2 < temp3,"Virhe valiseinien luonnissa kun 2 valiseinaa")

        self.kartta.luo_valiseinat(3)
        temp4= len(self.kartta.seinat)
        self.assertTrue(temp3 < temp4,"Virhe valiseinien luonnissa kun 3 valiseinaa")


    
    def test_piirra_yksilot(self):

        """
        Testaa Kartta-luokan metodia piirra_yksilot()
        Testit ovat huonoja ja aukaisevat ruudun kayttajan naytolle
        """  
        
        
        yksilo_lista = []
        
        for i in range(0,3):
                paikka_x = random.randint(10,390)
                paikka_y = random.randint(10,390)
                yksilo = Yksilo(paikka_x,paikka_y,5) 
                yksilo_lista.append(yksilo)


        white = 255,255,255
        self.kartta.ruutu.fill(white)
        self.kartta.piirra_yksilot(self.kartta.ruutu,yksilo_lista)

        pygame.display.flip()
        

    def test_paivita_kartta(self):
        """
        Testaa Kartta-luokan metodia paivita_kartta() paivittymista visuaalisesti.
        Arpoo yksiloille uusia paikkoja. Jalleen huono testi
        """  


        yksilo_lista = []
        
        for i in range(0,3):
                paikka_x = random.randint(10,390)
                paikka_y = random.randint(10,390)
                yksilo = Yksilo(paikka_x,paikka_y,5) 
                yksilo_lista.append(yksilo)
        


        while (self.kartta.paivita_kartta(yksilo_lista)):

                for yksilo in yksilo_lista:
                        yksilo.paikka[0] = random.randint(380,390)
                        yksilo.paikka[0] = random.randint(380,390)



                



        

    
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
