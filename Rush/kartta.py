import os,sys,math,pygame,pygame.mixer
from pygame.locals import *
from inputit import Input



class Kartta(object):

    
    def __init__(self):     
                        
        self.seinat = []
        ruudun_koko = (Input.RUUDUN_LEVEYS, Input.RUUDUN_KORKEUS)
        
        self.ruutu = pygame.display.set_mode(ruudun_koko)
                   
            
    def alusta_kartta(self,valiseinamaara):
        """
        Luo seinaobjektit ja lisaa ne self.seinat listaan
        """
                        

        y_coord = range(Input.KAYTAVAN_PITUUS,Input.RUUDUN_KORKEUS,Input.SEINAN_LEVEYS)
        x_coord = range(0,Input.RUUDUN_LEVEYS,Input.SEINAN_LEVEYS)
        
        
        #vasen seina
        for dy in y_coord:
            self.seinat.append(Seina((0,dy)))
        #oikea seina
        for dy in y_coord:
            self.seinat.append(Seina((Input.RUUDUN_LEVEYS-Input.SEINAN_LEVEYS,dy)))
        #alaseina
        for dx in x_coord:
            self.seinat.append(Seina((dx,Input.RUUDUN_KORKEUS-Input.SEINAN_LEVEYS)))

        x_coord = range(0,int(Input.RUUDUN_LEVEYS/2) - int(Input.OVIAUKON_LEVEYS/2),Input.SEINAN_LEVEYS)
        x_coord2 = range(int(Input.RUUDUN_LEVEYS/2)+ int(Input.OVIAUKON_LEVEYS/2),Input.RUUDUN_LEVEYS,Input.SEINAN_LEVEYS)

        #ylaseina
        for dx in x_coord:
            self.seinat.append(Seina((dx,Input.KAYTAVAN_PITUUS)))

        for dx in x_coord2:
            self.seinat.append(Seina((dx,Input.KAYTAVAN_PITUUS)))


        y_coord = range(0,int(Input.KAYTAVAN_PITUUS))
        vasen_x = int(Input.RUUDUN_LEVEYS/2)-int(Input.OVIAUKON_LEVEYS/2)
        oikea_x = int(Input.RUUDUN_LEVEYS/2)+int(Input.OVIAUKON_LEVEYS/2)
        for dy in y_coord:
            self.seinat.append(Seina((vasen_x,dy)))
            self.seinat.append(Seina((oikea_x,dy)))

        if (valiseinamaara != 0):

            self.luo_valiseinat(valiseinamaara)
        
        
            

    def luo_valiseinat(self,maara):
        """
        Lisaa parametrin maara verran valiseinia self.seinat listaan
        """

        if (maara >= 1):
            alku = int(Input.RUUDUN_LEVEYS/5.)
            loppu = int(Input.RUUDUN_LEVEYS*(2./5.))
            korkeudella = Input.RUUDUN_KORKEUS/2.
            for dx in range(alku,loppu,Input.SEINAN_LEVEYS):
                self.seinat.append(Seina((dx,korkeudella)))

        if (maara >= 2):
            alku = int(Input.RUUDUN_LEVEYS*(3./5.)) 
            loppu = int(Input.RUUDUN_LEVEYS*(4./5.)) 
            korkeudella = Input.RUUDUN_KORKEUS*(3./7.) 

            for dx in range(alku,loppu,Input.SEINAN_LEVEYS):
                self.seinat.append(Seina((dx,korkeudella)))
            
            
        if (maara >= 3):
            alku = int(Input.RUUDUN_LEVEYS*(4./7.)) 
            loppu = int(Input.RUUDUN_LEVEYS*(9./10.)) 
            korkeudella = Input.RUUDUN_KORKEUS*(2./3.) 

            for dx in range(alku,loppu,Input.SEINAN_LEVEYS):
                self.seinat.append(Seina((dx,korkeudella)))
        


    
    def piirra_yksilot(self,ruutu,yksilo_lista):
        """
        Piirtaa parametrin yksilo_lista:n kaikki yksilot ruudulle
        """

        for yksilo in yksilo_lista:

            paikka = (int(yksilo.paikka[0]),int(yksilo.paikka[1]))
            
            pygame.draw.circle(ruutu, (0,0,0), paikka,yksilo.size,yksilo.width) 



    def paivita_kartta(self,yksilo_lista):
            """
            Piirtaa yhden kierroksen yksilot ja seinat ruudulle
            """        
                                     
            pygame.display.set_caption("Ruuhkasimulaatio")
            piirra = True       #parametri joka antaa merkin ruuhkasimulaation lopettamisesta
                     
            for event in pygame.event.get():
            
                if event.type == pygame.QUIT:
                    piirra = False
                    return piirra            
                            
            self.ruutu.fill((255, 255, 255))

            for seina in self.seinat:
                pygame.draw.rect(self.ruutu, (0, 0, 255), seina.rect)                

            self.piirra_yksilot(self.ruutu,yksilo_lista)
            pygame.display.flip()

            return piirra
                




class Seina(object):
    """
    Lyhyt seinaluokka. Keskikohdan koordinaatteja kaytetaan Kaytosmallissa
    seinavoimien laskuun.
    """
    def __init__(self, pos):
        
        self.rect = pygame.Rect(pos[0], pos[1], Input.SEINAN_LEVEYS, Input.SEINAN_LEVEYS)
        self.keski_x = pos[0]
        self.keski_y = pos[1]






