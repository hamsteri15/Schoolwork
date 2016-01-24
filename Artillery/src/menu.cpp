#include <SFML/Graphics.hpp>


int menu(int winner)
{
    sf::RenderWindow window(sf::VideoMode(400, 400), "Menu");
    window.setPosition(sf::Vector2i(400, 50));
    
    sf::Font font;
    
    font.loadFromFile("../fonts/DejaVuSans-Bold.ttf");
    
    sf::Text text, gameOverText, textEsc;
    sf::FloatRect rect;
    text.setFont(font);
    gameOverText.setFont(font);
    textEsc.setFont(font);
    text.setCharacterSize(20);
    gameOverText.setCharacterSize(25);
    textEsc.setCharacterSize(15);

    if(winner) {
    	
    	std::string str  = "Game Over! \nPlayer " + std::to_string(winner) + " won!";
    	gameOverText.setString(str);
    	rect = gameOverText.getLocalBounds();
    	gameOverText.setOrigin(rect.width/2,rect.height/2);
    	gameOverText.setPosition(200,50);
    	gameOverText.setColor(sf::Color::Red);
    }

    text.setString("How many players? (2 - 5)");
    rect = text.getLocalBounds();
    text.setOrigin(rect.width/2,rect.height/2);
    text.setPosition(200, 150);
    text.setColor(sf::Color::Black);

    textEsc.setString("Press esc to quit");
    rect = textEsc.getLocalBounds();
    textEsc.setOrigin(rect.width/2, rect.height/2);
    textEsc.setPosition(200, 350);
    textEsc.setColor(sf::Color::Red);

    

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) {
                window.close();
                return 0;
            }
            if (event.type == sf::Event::KeyPressed) {
            	if (sf::Keyboard::Num2 <= event.key.code && event.key.code <= sf::Keyboard::Num5) {
            		window.close();
            		return (int)event.key.code - 26;
            	}
            	if (event.key.code == sf::Keyboard::Escape) {
            		window.close();
            		return 0;
            	}
            }
        }
	    window.clear(sf::Color(135,206,250));
	    window.draw(text);
	    window.draw(gameOverText);
	    window.draw(textEsc);
	    window.display();
    }

    return 0;
} 
