import unittest
import math
from pygame.math import Vector2 as Vector
from kaytosmalli import Kaytosmalli
from yksilo import Yksilo
from inputit import Input



class Test(unittest.TestCase):


    def setUp(self):
        self.kaytosmalli = Kaytosmalli()
        
        
        
    
    def test_lisaa_yksilo(self):
        """
        Testaa Kaytosmalli luokan metodia lisaa_yksilo()
        """

        eka = Yksilo(60,60)
        self.kaytosmalli.lisaa_yksilo(eka)
        self.assertEqual(len(self.kaytosmalli.yksilo_lista),1,"Yksilon lisaaminen ei toimi")


    
    def test_laske_voimat(self):
        """
        Testaa Kaytosmalli luokan metodia laske_voimat(). Tarkistaa vain etta voimat muuttuvat nollasta 
        poikkeavaksi.
        """
        eka = Yksilo(10,60)
        toka = Yksilo(70,60)
        
        self.kaytosmalli.lisaa_yksilo(eka)
        self.kaytosmalli.lisaa_yksilo(toka)

        self.kaytosmalli.laske_voimat()
        
        self.assertNotEqual(Vector(0,0),eka.voima, "Yksilon voima ei muutu vaikka pitaisi.")
        self.assertNotEqual(Vector(0,0),toka.voima, "Yksilon voima ei muutu vaikka pitaisi.")
        


    

    def test_laske_seeking(self):
        """
        Testaa Kaytosmalli luokan metodia laske_seeking()
        """
        eka = Yksilo(0,0)

        #ylospain maksimivoimalla       
        eka.paikka = self.kaytosmalli.KOHDE - Vector(0,-40)   
        ohjausvoima = self.kaytosmalli.laske_seeking(eka)
        ohjausvoima_tarkistus = eka.max_nopeus * (self.kaytosmalli.KOHDE - eka.paikka).normalize()                
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "liikkuminen maksiminopeudella ylospain ei toimi")
        
        #alaspain maksimivoimalla
        eka.paikka = self.kaytosmalli.KOHDE + Vector(0,-40)
        ohjausvoima = self.kaytosmalli.laske_seeking(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * (self.kaytosmalli.KOHDE - eka.paikka).normalize()              
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "liikkuminen maksiminopeudella alaspain ei toimi")

        #vinottain
        eka.paikka = self.kaytosmalli.KOHDE + Vector(-100,-40) 
        ohjausvoima = self.kaytosmalli.laske_seeking(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * (self.kaytosmalli.KOHDE - eka.paikka).normalize()              
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "liikkuminen vinottain ei toimi")
        

    def test_laske_seina_ja_seinamaetaisyys(self):

        """
        Testaa Kaytosmalli luokan metodeita laske_seina() ja seinamaetaisyys()
        """
        
        #seina yksilon oikealla puolella
        self.kaytosmalli.seinien_paikkavektorit.append(Vector(10,0))
        eka = Yksilo(0,0)
        ohjausvoima, temp = self.kaytosmalli.laske_seina(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * Vector(-1,0)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "seina yksilon oikealla puolella ei toimi")

        #seina yksilon vasemmalla puolella
        self.kaytosmalli.seinien_paikkavektorit.append(Vector(-9,0))
        eka = Yksilo(0,0)
        ohjausvoima, temp = self.kaytosmalli.laske_seina(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * Vector(1,0)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "seina yksilon vasemmalla puolella ei toimi")
        
        #seina yksilon ylapuolella
        self.kaytosmalli.seinien_paikkavektorit.append(Vector(0,-6))
        eka = Yksilo(0,0)
        ohjausvoima, temp = self.kaytosmalli.laske_seina(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * Vector(0,1)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "seina yksilon ylapuolella ei toimi")

        #seina yksilon alapuolella
        self.kaytosmalli.seinien_paikkavektorit.append(Vector(0,5))
        eka = Yksilo(0,0)
        ohjausvoima, temp = self.kaytosmalli.laske_seina(eka)
        ohjausvoima_tarkistus =  eka.max_nopeus * Vector(0,-1)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "seina yksilon alapuolella ei toimi")
        
        #======SEINAETAISYYDEN TESTIT======#

        #listassa useita seinia
        lyhin_etaisyys = self.kaytosmalli.seinamaetaisyys(eka)
        lyhin_etaisyys_tarkistus = Vector(0,5)
        self.assertEqual(lyhin_etaisyys,lyhin_etaisyys_tarkistus, "ei loyda lyhinta seinaa listasta")

        #tyhja lista
        del self.kaytosmalli.seinien_paikkavektorit[:]
        lyhin_etaisyys = self.kaytosmalli.seinamaetaisyys(eka)
        lyhin_etaisyys_tarkistus = Vector(1000,1000)
        self.assertEqual(lyhin_etaisyys,lyhin_etaisyys_tarkistus, "tyhjalla seinalistalla pitaisi palauttaa ylipitka vektori")
        


    def test_laske_ruuhkautumis(self):
        """
        Testaa Kaytosmalli luokan metodia laske_ruuhkautumis()
        """
        
        eka = Yksilo(60,60)
        toka = Yksilo(70,60)
        
        self.kaytosmalli.lisaa_yksilo(eka)
        self.kaytosmalli.lisaa_yksilo(toka)

        ohjausvoima = self.kaytosmalli.laske_ruuhkautumis(eka)
        ohjausvoima_tarkistus = Vector(-1,0)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "ruuhkautumisvoimaa ei laskettu oikein")

        ohjausvoima = self.kaytosmalli.laske_ruuhkautumis(toka)
        ohjausvoima_tarkistus = Vector(1,0)
        self.assertEqual(ohjausvoima_tarkistus,ohjausvoima, "ruuhkautumisvoimaa ei laskettu oikein")


    def test_jarruta_tarvittaessa(self):
        """
        Testaa Kaytosmalli luokan metodia jarruta_tarvittaessa()
        """
        
        eka = Yksilo(60,60)
        eka.voima = Vector(100,100)
        toka = Yksilo(70,60)
        
        self.kaytosmalli.lisaa_yksilo(eka)
        self.kaytosmalli.lisaa_yksilo(toka)
        self.kaytosmalli.jarruta_tarvittaessa(eka)
        self.assertEqual(Vector(0,0),eka.voima, "ei jarruta oikein")
        
        eka.voima = Vector(100,100)
        toka.paikka = Vector(85,60)
        self.kaytosmalli.jarruta_tarvittaessa(eka)
        self.assertEqual(Vector(100,100),eka.voima, "jarruttaa vaikka ei pitaisi")
        

    def test_paivita_paikat(self):
        """
        Testaa Kaytosmalli luokan metodia paivita_paikat()
        """

        eka = Yksilo(60,60)
        eka.voima = Vector(100,100)
        toka = Yksilo(70,60)

        self.kaytosmalli.lisaa_yksilo(eka)
        self.kaytosmalli.lisaa_yksilo(toka)
        self.kaytosmalli.paivita_paikat()

        self.assertEqual(toka.paikka, Vector(70,60), "ohausvoima 0, muuttaa silti paikkaa")        
        self.assertNotEqual(eka.paikka, Vector(60,60), "ohausvoima != 0, ei muuta silti paikkaa")
           
        eka.paikka = self.kaytosmalli.KOHDE + Vector(1,0)
        self.kaytosmalli.paivita_paikat()
        self.assertEqual(len(self.kaytosmalli.yksilo_lista), 1, "ei poista yksiloa listasta vaikka on lahella kohdetta") 

    def test_onko_edessa(self):
        """
        Testaa Kaytosmalli luokan metodia onko_edessa()
        """       

        eka = Yksilo(0,0)
        toka = Yksilo(10,0)
        
        self.kaytosmalli.lisaa_yksilo(eka)
        self.kaytosmalli.lisaa_yksilo(toka)
        
        #katsoo suoraan kohti
        eka.suunta = Vector(1,0)
        self.assertEqual(True,self.kaytosmalli.onko_edessa(eka,toka.paikka),"onko_edessa vaara arvo kun naapuri suoraan edessa")
        
        #katsoo vastakkaiseen suuntaan
        eka.suunta = Vector(-1,0)
        self.assertEqual(False,self.kaytosmalli.onko_edessa(eka,toka.paikka),"onko_edessa vaara arvo kun naapuri suoraan takana")

        #pitaisi olla viela katsealueella mutta ei suoraan edessa
        x = 2.
        eka.suunta = Vector(0,1)
        toka.paikka = Vector(x,10)
        self.assertEqual(True,self.kaytosmalli.onko_edessa(eka,toka.paikka),"onko_edessa vaara arvo kun naapuri viela katsealueella")

        #pitaisi olla viela katsealueella mutta ei suoraan edessa
        x = -2.
        eka.suunta = Vector(0,1)
        toka.paikka = Vector(x,10)
        self.assertEqual(True,self.kaytosmalli.onko_edessa(eka,toka.paikka),"onko_edessa vaara arvo kun naapuri viela katsealueella (negatiivinen)")

        #ei pitaisi olla enaa katsealueella mutta ei suoraan edessa

        x = math.tan(math.radians(Input.KATSEALUE/2.))*Input.MIN_JARRUTUSETAISYYS + 0.1
        eka.suunta = Vector(0,1)
        toka.paikka = Vector(x,Input.MIN_JARRUTUSETAISYYS)
        self.assertEqual(False,self.kaytosmalli.onko_edessa(eka,toka.paikka),"onko_edessa vaara arvo kun naapuri ei ole enaa katsealueella")

        
        
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
