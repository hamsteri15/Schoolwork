import unittest
from pygame.math import Vector2 as Vector
from yksilo import Yksilo



class Test(unittest.TestCase):


    def setUp(self):
        self.yksilo = Yksilo(100,100)
         
        
        
    
    def test_jarruta(self):
        
        """
        Tarkastaa etta voimavektori on pienempi jarrutuksen jalkeen.
        Jarrutus tarkoittaa siis ohjausvoimavektorin pituuden lyhentamista 
        """

        self.yksilo.voima = Vector(1,1)
        alkuperainen_pituus = self.yksilo.voima.length()
        
        self.yksilo.jarruta()
        jarrutus_pituus = self.yksilo.voima.length()

        self.assertTrue(alkuperainen_pituus > jarrutus_pituus, "Jarrutus ei toimi oikein")



    def test_katkaise_vektori(self):

        """
        Tarkastaa etta vektori lyhennetaan tarvittaessa ja etta sita ei lyhenneta kun ei tarvitse. 
        """


        vektori = Vector(3,3)
        alkuperainen_pituus = vektori.length()
        
        maksimi_pienempi = 0.5*alkuperainen_pituus
        maksimi_isompi = 2*alkuperainen_pituus

        
        self.yksilo.katkaise_vektori(vektori,maksimi_pienempi)        
        self.assertTrue(vektori.length() == maksimi_pienempi, "Vektorin katkaisu ei toimi")


        vektori = Vector(3,3)

        self.yksilo.katkaise_vektori(vektori,maksimi_isompi)        
        self.assertTrue(vektori.length() == alkuperainen_pituus, "Vektorin katkaisemattomuus ei toimi")




    def test_laske_uusi_paikka(self):
        """
        Testaa liikkeen positiivisiin ja negatiivisiin x,y suuntiin sekä että nollavektorilla ei tule virheitä
        """

        #positiivien x 
        self.yksilo.paikka = Vector(100,100)
        self.yksilo.voima = Vector(1,0)
        uusi_paikkavektori = self.yksilo.laske_uusi_paikka()
        self.assertEqual(uusi_paikkavektori,Vector(101,100), "oikealle liike ei toimi")
        self.yksilo.nopeus=Vector(0,0)

        #negatiivinen x
        self.yksilo.paikka = Vector(100,100)
        self.yksilo.voima = Vector(-1,0)
        uusi_paikkavektori = self.yksilo.laske_uusi_paikka()
        self.assertEqual(uusi_paikkavektori,Vector(99,100), "vasemmalle liike ei toimi")
        self.yksilo.nopeus=Vector(0,0)

        #positiivien y 
        self.yksilo.paikka = Vector(100,100)
        self.yksilo.voima = Vector(0,1)
        uusi_paikkavektori = self.yksilo.laske_uusi_paikka()
        self.assertEqual(uusi_paikkavektori,Vector(100,101), "ylospain liike ei toimi")
        self.yksilo.nopeus=Vector(0,0)

        #negatiivinen y 
        self.yksilo.paikka = Vector(100,100)
        self.yksilo.voima = Vector(0,-1)
        uusi_paikkavektori = self.yksilo.laske_uusi_paikka()
        self.assertEqual(uusi_paikkavektori,Vector(100,99), "alaspain liike ei toimi")
        self.yksilo.nopeus=Vector(0,0)
        
        #nollavoimavektori
        self.yksilo.paikka = Vector(100,100)
        self.yksilo.voima = Vector(0,0)
        uusi_paikkavektori = self.yksilo.laske_uusi_paikka()
        self.assertEqual(uusi_paikkavektori,Vector(100,100), "nollavektori ei toimi")
        self.yksilo.nopeus=Vector(0,0)
               
        


                  





    

                



        

    
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
