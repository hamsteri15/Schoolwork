#include <SFML/Graphics.hpp>

#include <sstream>
#include "game.hpp"
#include "menu.hpp"

int main()
{	

	int numberOfPlayers;
	int winner = 0;
	
    
	do {

		numberOfPlayers = menu(winner);
		
		if(!numberOfPlayers) {
			break;
		}
		
		sf::RenderWindow window(sf::VideoMode(800, 600), "Artillery 3");
	    window.setPosition(sf::Vector2i(400, 50));
		Game game = Game(numberOfPlayers, window);
	    winner = game.newGame();

	    if(!winner){ 
	    	break;
	    }
    
    } while (1);

    return 0;
}
