import os,sys,math,pygame,pygame.mixer,numpy,random
from pygame.locals import *
from pygame.math import Vector2 as Vector
from inputit import Input

class Kaytosmalli(object):


        KOHDE = Input.KOHDE 
        MIN_SEINAETAISYYS = Input.MIN_SEINAETAISYYS 
        MIN_KOHDEETAISYYS = Input.MIN_KOHDEETAISYYS 
        MIN_VALIETAISYYS =  Input.MIN_VALIETAISYYS
        MIN_JARRUTUSETAISYYS = Input.MIN_JARRUTUSETAISYYS  
        KATSEALUE = Input.KATSEALUE 

        def __init__(self):
                                
            self.yksilo_lista = []
            self.seinien_paikkavektorit = []
            self.seinaneliot = []
            
                

        def lisaa_yksilo(self,yksilo):
            """
            Lisaa parametrin yksilo listaan self.yksilo_lista
            """
            self.yksilo_lista.append(yksilo)
        
          


        def laske_voimat(self):
            """
            Summaa ohausvoimat ja jarruttaa tarvittaessa
            """

            for yksilo in self.yksilo_lista:
                seeking_voima = self.laske_seeking(yksilo)
                seina_voima,C = self.laske_seina(yksilo) 
                ruuhkautumis_voima = self.laske_ruuhkautumis(yksilo)
                """
                Voimien painokertoimet. Metodi self.laske_ruuhkautumis() 
                palauttaa seinavoiman kertoimen
                """
                A = 1.0  #seeking
                B = 50.  #ruuhkautumis
                
                yksilo.voima = (A*seeking_voima + B*ruuhkautumis_voima + C*seina_voima)
            
                """
                Tarkastetaan pitaanko jarruttaa. Kaikille yksiloille ei tehda tarkastusta
                simuloinnin jumittumisen takia
                """
                if (random.randint(0, len(self.yksilo_lista)) > 1):
                        self.jarruta_tarvittaessa(yksilo)
                       


        def laske_seeking(self,yksilo):
            """
            Laskee ja palauttaa seeking-voiman parametrille yksilo
            """          
            suuntavektori = (self.KOHDE - yksilo.paikka).normalize()
            tavoitenopeus = suuntavektori * yksilo.max_nopeus
            ohjausvektori = tavoitenopeus + yksilo.nopeus       
            return ohjausvektori        
                

        def laske_seina(self,yksilo):
            """
            Laskee ja palauttaa seinavoiman parametrille yksilo.
            Jos yksilo ei ole seinan lahella palauttaa nollavektorin
            """
            ohjausvektori = Vector(0,0)          
            seinavektori = self.seinamaetaisyys(yksilo)
                                  
            if (seinavektori.length() < self.MIN_SEINAETAISYYS and \
                seinavektori.length() != 0):

                tavoitenopeus = -seinavektori.normalize()*yksilo.max_nopeus
                ohjausvektori = tavoitenopeus - yksilo.nopeus               
                C = (self.MIN_SEINAETAISYYS/seinavektori.length())**3
            else:
                C = 0.
            return ohjausvektori, C



        def laske_ruuhkautumis(self,yksilo):
            """
            Laskee ja palauttaa seeking voiman parametrille yksilo.
            """
            ohjausvektori = Vector(0,0)
            for naapuri in self.yksilo_lista:

                suunta = yksilo.paikka - naapuri.paikka
                etaisyys = suunta.length()

                if (etaisyys != 0 and etaisyys < self.MIN_VALIETAISYYS):
                        ohjausvektori += suunta/etaisyys
                    
            return ohjausvektori
            




        def seinamaetaisyys(self,yksilo):
            """
            Palauttaa vektorin, joka menee yksilosta lahimpaan seinaan
            """
            lyhin = Vector(1000,1000)

            for seina in self.seinien_paikkavektorit:

                etaisyys = seina - yksilo.paikka
                if (etaisyys.length() < lyhin.length()):
                    lyhin = etaisyys

            return lyhin
                    

        

        def alusta_seinien_paikkavektorit(self,Kartta):
            """
            Luo paikkavektorin jokaiselle seinaelementille seinavoimien laskentaa
            varten
            """
            for seina in Kartta.seinat:
                self.seinien_paikkavektorit.append(Vector(seina.keski_x,seina.keski_y))
                self.seinaneliot.append(seina)
      


        def jarruta_tarvittaessa(self,yksilo):
            """
            Jarruttaa jos joku naapuri on yksilon edessa
            """        
            for naapuri in self.yksilo_lista:
                
                if (self.onko_edessa(yksilo,naapuri.paikka)):                  
                    yksilo.jarruta()
                    break 
                                                     
                                      


        def onko_edessa(self,yksilo,paikkavektori):
            """
            Tarkistaa onko parametri paikkavektori yksilon edessa, 
            eli tarvitseeko jarruttaa.
            """        
            if (yksilo.paikka == paikkavektori):
                return False

            etaisyys = paikkavektori - yksilo.paikka
            kulma = yksilo.suunta.angle_to(etaisyys)
            if (etaisyys.length() < self.MIN_JARRUTUSETAISYYS and\
            abs(kulma)<self.KATSEALUE):
                return True

                    
            return False


        

        

        def paivita_paikat(self):
            """
            Paivittaa yksiloiden paikat ja jos yksilo on maalissa, poistaa sen listasta.
            Aikaisemmin tassa oli viela tarkistus mutta siita luovuttiin.
            """
            for yksilo in self.yksilo_lista:

                uusi_paikka = yksilo.laske_uusi_paikka()        
                yksilo.paikka = uusi_paikka

                if ((yksilo.paikka-self.KOHDE).length() < self.MIN_KOHDEETAISYYS):
                    self.yksilo_lista.remove(yksilo)
                      
                
            
            



        
                





        
