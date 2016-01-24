import os,sys,math,pygame,pygame.mixer,random
from pygame.locals import *
from kaytosmalli import Kaytosmalli
from yksilo import Yksilo
from kartta import Kartta
from inputit import Input


class Simulaatio(object):

        
    def __init__(self):     
                   

        self.input = Input()
        self.kaytosmalli = Kaytosmalli() 
        self.kartta = Kartta()                
        

    def aloita_simulointi(self):
        """
        Kysyy kayttajalta yksilo- ja seinamaarat. Alustaa kartan ja simuloi kunnes
        Kartta.paivita_kartta() palauttaa arvon False
        """        

        yksilomaara = self.kysy_yksilomaara()
        estemaara = self.kysy_estemaara()
        self.luo_yksilot(yksilomaara)
        self.kartta.alusta_kartta(estemaara)
        self.kaytosmalli.alusta_seinien_paikkavektorit(self.kartta)

        clock = pygame.time.Clock()
        
        
        while (self.kartta.paivita_kartta(self.kaytosmalli.yksilo_lista)):
        
            self.kaytosmalli.laske_voimat()                                   
            self.kaytosmalli.paivita_paikat()

            if (len(self.kaytosmalli.yksilo_lista)==0):
                    self.luo_yksilot(yksilomaara)
            
            
            clock.tick(Input.FPS)
                    
                    
    def kysy_yksilomaara(self):
        """
        Kysyy ja palauttaa kayttajalta yksilomaaran ja virheellisen syotteen jalkeen
        kysyy uudestaan virheilmoituksen kanssa
        """
        virheilmoitus = False
        while (True):
            
            try:
                yksilomaara = self.input.pyyda_maara("Anna yksilomaara (1 - 100)",virheilmoitus)
                if (yksilomaara<1 or yksilomaara>100):               
                    raise ValueError
                else:
                    break
            except ValueError:
                virheilmoitus = True
                
        return yksilomaara


    def kysy_estemaara(self):
        """
        Kysyy ja palauttaa kayttajalta valiseinamaaran ja virheellisen syotteen jalkeen
        kysyy uudestaan virheilmoituksen kanssa
        """

        virheilmoitus = False
        while (True):
            try:
                estemaara = self.input.pyyda_maara("Anna esteiden maara (0 - 3)",virheilmoitus)
                if (estemaara<0 or estemaara>3):
                    raise ValueError
                else:
                    break
            except ValueError:
                virheilmoitus = True
                      
    
        return estemaara
            

    def luo_yksilot(self,yksilomaara):
        """
        Luo parametrin yksilomaaran verran Yksiloita ja arpoo niille aloituspaikat
        ulkoseinien sisapuolelta
        """
        
        for i in range(0,yksilomaara):
            paikka_x = random.randint(10,Input.RUUDUN_LEVEYS-10)
            paikka_y = random.randint(45,Input.RUUDUN_KORKEUS-10)            
            yksilo = Yksilo(paikka_x,paikka_y)        
                     
            self.kaytosmalli.lisaa_yksilo(yksilo)



if __name__ == "__main__":
    simulaatio = Simulaatio()
    simulaatio.aloita_simulointi()

            

        
                

                
                
        



        






