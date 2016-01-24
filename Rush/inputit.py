import pygame, pygame.font, pygame.event, pygame.draw, string,sys
from pygame.locals import *
from pygame.math import Vector2 as Vector

class Input(object):

    """
    KARTAN PARAMETRIT
    """
    RUUDUN_LEVEYS = 500 #500
    RUUDUN_KORKEUS = 400 #400
    SEINAN_LEVEYS = 3
    OVIAUKON_LEVEYS = 50
    KAYTAVAN_PITUUS = 40
    
    """
    RUUHKASIMULAATION PARAMETRIT
    """
    FPS = 70
     
    """
    KAYTOSMALLIN PARAMETRIT
    """
    KOHDE = Vector(RUUDUN_LEVEYS/2,-40)
    MIN_SEINAETAISYYS = 20
    MIN_KOHDEETAISYYS = 15
    MIN_VALIETAISYYS = 12
    MIN_JARRUTUSETAISYYS = 15
    KATSEALUE = 20.
  
        

    
    def pyyda_maara(self,kysymys,virheilmoitus):
            
        ruutu = pygame.display.set_mode((self.RUUDUN_LEVEYS,self.RUUDUN_KORKEUS))
        pygame.display.set_caption("Ruuhkasimulaatio")
        
        maara = int(self.kysy_kysymys(ruutu, kysymys,virheilmoitus))

        
        return maara

    def pyyda_merkki(self):
        while 1:
            event = pygame.event.poll()
            if event.type == KEYDOWN:
              return event.key
        
            if event.type == pygame.QUIT:
              sys.exit()                 
            else:
              pass
        
    
    def laatikko(self,ruutu,kysymys,virheilmoitus):
  
        fontobject = pygame.font.Font(None,18)
        pygame.draw.rect(ruutu, (0,0,0),((ruutu.get_width() / 2) - 100, \
                        (ruutu.get_height() / 2) - 10,200,20), 0)
        pygame.draw.rect(ruutu, (255,255,255),((ruutu.get_width() / 2) - 102, \
                        (ruutu.get_height() / 2) - 12,204,24), 1)

        if len(kysymys) != 0:
            ruutu.blit(fontobject.render(kysymys, 1, (255,255,255)),((ruutu.get_width() / 2) - 95, \
                      (ruutu.get_height() / 2) - 10))

        if (virheilmoitus):
            ruutu.blit(fontobject.render("Virheellinen arvo!", 1, (255,0,0)),((ruutu.get_width() / 2) - 95, \
                          (ruutu.get_height() / 2) - 50))
        pygame.display.flip()
    


    def kysy_kysymys(self,ruutu,kysymys,virheilmoitus):
  
        pygame.font.init()
        merkkijono = ""
  
        self.laatikko(ruutu, kysymys + ":  " + merkkijono,virheilmoitus)
        while 1:
            
            merkki = self.pyyda_merkki()    
            if merkki == K_BACKSPACE:
                merkkijono = merkkijono[0:-1]
            elif merkki == K_RETURN:
                break                
            elif merkki <= 127:          
                merkkijono+=chr(merkki)

            

            self.laatikko(ruutu, kysymys + ": " + merkkijono,virheilmoitus)
        return merkkijono

    
    







        

                



























                
