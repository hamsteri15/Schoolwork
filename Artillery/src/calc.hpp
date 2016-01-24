#ifndef calc_hh
#define calc_hh

#include <vector>
#include <SFML/Graphics.hpp>
#include "parameters.hpp"
namespace Calculate
{
    ////////////////////////////////////////////////////////////
    ///
    ///    \Calculates the pixels inside a circle.
    ///    \param width: screen width
    ///    \param height: screen height
    ///    \param x: x-coordinate of the hole center
    ///    \param y: y-coordinate of the hole center
    ///    \param r: radius of the hole  
    ///
    ////////////////////////////////////////////////////////////
	std::vector<unsigned int> circlePixels(int width, int height, int x, int y, int r);
	
	
    ////////////////////////////////////////////////////////////
    ///
    ///    \Checks whether location is inside a circle 
    ///    \param loc: location to be checked
    ///    \param circleLoc: location of the circle
    ///    \param r: radius of the circle
    ///
    ////////////////////////////////////////////////////////////
    bool belongsToCircle(sf::Vector2u loc, sf::Vector2u circleLoc, int r);

	
	////////////////////////////////////////////////////////////
    ///
    ///    \Calculates the next location from the speed and current location
    ///     location of the object.
    ///    \param speed: current speed of the object
    ///    \param currentLocation: current location of the object
    ///
    ////////////////////////////////////////////////////////////
    sf::Vector2u nextLocation(sf::Vector2f speed, sf::Vector2u currentLocation);
    
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Calculates the initial speed from the mass and force
    ///     exerted on the object.
    ///    \param mass: mass of the object
    ///    \param energy: kinetic energy of the object
    ///    \param angle: angle between x-axis and the object mouse click
    ///
    ////////////////////////////////////////////////////////////
    sf::Vector2f initialSpeed(int mass, float energy, float angle);
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Calculates a new speed for the object based on
    ///     accelerations caused by wind and gravity.
    ///    \param speed: current speed of the object
    ///    \param wind: wind coefficient
    ///    \param mass: mass of the object
    ///
    ////////////////////////////////////////////////////////////
    sf::Vector2f newSpeed(sf::Vector2f speed, float wind, int mass);


    ////////////////////////////////////////////////////////////
    ///
    ///    \Calculates new wind
    ///    \param wind: reference to the game m_wind variable
    ///
    ////////////////////////////////////////////////////////////
    void wind(float& wind);


    ////////////////////////////////////////////////////////////
    ///
    ///    \Calculates the angle between the velocity vector and x-axis
    ///    \param bulletSpeed: velocity vector
    ///
    ////////////////////////////////////////////////////////////
    float angle(sf::Vector2f bulletSpeed);

}

#endif
