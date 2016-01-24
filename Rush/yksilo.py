from pygame.math import Vector2 as Vector


class Yksilo(object):

    def __init__(self, x ,y , size=5, color = (0,0,255), width = 0):         
    
        self.x = x
        self.y = y
        self.size = size   
        self.color = color 
        self.width = width 


        self.massa = 1.
        
        self.paikka = Vector(x,y)
        
        self.nopeus = Vector(0,0)
        self.suunta = Vector(0,0)
        self.voima = Vector(0,0)

        self.max_voima= 1.0
        self.max_nopeus = 1.7
    

        

    def laske_uusi_paikka(self):
        
        self.edellinen_paikka = self.paikka
          
        if (self.voima.length() != 0):                  
                self.katkaise_vektori(self.voima, self.max_voima)

                kiihtyvyys = self.voima/self.massa
                
                nopeus = self.nopeus + kiihtyvyys                
                self.katkaise_vektori(nopeus,self.max_nopeus)
                
                #laske uusi paikkavektori
                uusi_paikka = nopeus+self.paikka
                #paivita katsomissuunta ja nopeus
                self.nopeus = nopeus     
                self.suunta = nopeus.normalize()
                          

        else:
                uusi_paikka = self.paikka
                

        #uusi paikka paivitetaan vasta kaytosmallissa
        return uusi_paikka
        
    

    def jarruta(self):
        """
        asettaa ohjausvoiman nollaan, eli pysayttaa yksilon valittomasti
        """
        
        self.voima -= self.voima
        
        



    def katkaise_vektori(self,vektori, maksimi):
        """
        lyhentaa parametria vektori jos se on pidempi kuin parametri maksimi.
        suunta pysyy ennallaan

        """
        if (vektori.length()>maksimi):
                vektori.scale_to_length(maksimi)

            






        

                



























                
