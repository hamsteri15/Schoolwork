import unittest
from pygame.math import Vector2 as Vector
from ruuhkasimulaatio import Simulaatio




class Test(unittest.TestCase):  


    def setUp(self):
        self.simulaatio = Simulaatio() 
        

    """
    Naita en osaa testata kuin graafisesti
    """
    
    def test_aloita_simulointi(self):      
        pass

    def test_kysy_yksilomaara(self):
        pass    

    def test_kysy_estemaara(self):
        pass


    def test_luo_yksilot(self):
        
        self.simulaatio.luo_yksilot(4)
        
        maara = len(self.simulaatio.kaytosmalli.yksilo_lista)
        self.assertEqual(4,maara,"yksiloiden lisaaminen ei onnistunut")
    

 

              



        

    
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
