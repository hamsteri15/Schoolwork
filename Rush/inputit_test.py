import unittest
from inputit import Input
import os,sys,pygame,pygame.mixer,random
from pygame.locals import *

"""
"Testit" pitaisi menna lapi enteria hakkaamalla. Kaikki testit on huonoja
"""

class Test(unittest.TestCase):


    def setUp(self):
        self.input = Input()
        
        
        
    
    def test_pyyda_maara(self): 
        """
        Testaa pyyda maara metodia. Enteria painamalla pitaisi nousta ValueError
        """
        
        self.assertRaises(ValueError, lambda: self.input.pyyda_maara("paina enter",False))
    
    
    def test_laatikko(self):
        """
        Ajaa laatikko() metodin koodin. Ei siis testaa mitaan.
        """ 
        
        ruutu = pygame.display.set_mode((500,500))
        self.input.laatikko(ruutu,"paina enter",False)
    
        

    def test_pyyda_merkki(self):
        """
        Ajaa pyyda_merkki() metodin koodin. Ei siis testaa mitaan.
        """
        self.input.pyyda_merkki()
    

    def test_kysy_kysymys(self):
        """
        Ajaa kysy_kysymys() metodin koodin. 
        """
        ruutu = pygame.display.set_mode((500,500))
        self.input.kysy_kysymys(ruutu, "fdsa",False)
        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
