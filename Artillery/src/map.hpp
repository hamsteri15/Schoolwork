#ifndef map_hh
#define map_hh

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include <memory>
#include <SFML/Graphics.hpp>
#include "calc.hpp"



class Map
{
public:

    ////////////////////////////////////////////////////////////
    ///
    ///    \Class constructor.
    ///    \param sizee: the screen size vector
    ///
    ////////////////////////////////////////////////////////////
    Map(sf::Vector2u size);
    
    
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Default destructor.
    ///
    ////////////////////////////////////////////////////////////
    ~Map() = default;
    
    

    ////////////////////////////////////////////////////////////
    ///
    ///    \Randomly creates the initial land shape. Called
    ///     from constructor.
    ///    \Note that an image of the land must be stored to
    ///     get access to the color-picker functions in SFML.
    ///    
    ////////////////////////////////////////////////////////////
    sf::Image randomizeLand();
    
    
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Checks if a location is inside ground.
    ///    \param loc: location to be checked    
    ///     
    ////////////////////////////////////////////////////////////
    bool isLand(sf::Vector2u loc);



	////////////////////////////////////////////////////////////
	///
	///    \Checks if a location is outside the map.
	///    \param loc: location to be checked    
	///     
	////////////////////////////////////////////////////////////
	bool outOfBounds(sf::Vector2u loc);
    


    ////////////////////////////////////////////////////////////
    ///
    ///    \Makes a circular hole to the ground. Essentially this
    ///     this colors the map image with the background color.
    ///    \param loc: location of the hole    
    ///    \param r: radius of the hole
    ////////////////////////////////////////////////////////////
    void makeHole(sf::Vector2u loc, int r);
    


    //void updatePlayerLocs(std::vector<Player> players);

	
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Obvious getter functions. Nothing special here.
    ///
    ////////////////////////////////////////////////////////////
    sf::Vector2u getSize();
    //
    sf::Texture getTexture();

    void draw(sf::RenderWindow& window);
    
    
private:
    
    sf::Vector2u    m_size;
    sf::Texture     m_rockTexture;
    sf::Texture     m_airTexture;
    sf::Image       m_image;    
    sf::Texture     m_texture;
    
};


#endif
