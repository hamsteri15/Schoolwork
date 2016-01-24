#include <math.h>
#include <vector>
#include "calc.hpp"
#include <iostream>


std::vector<unsigned int> Calculate::circlePixels(int width, int height, int x, int y, int r) {
	std::vector<unsigned int> pixelVector;

	int max_x, max_y, min_x, min_y;

    //Generate a small rectangle around point x,y. This is done to save calculation time.
	if ((x + r) > width) max_x = width;
	else max_x = x + r;

	if ((x - r) < 0) min_x = 0;
	else min_x = x - r;

	if ((y + r) > height) max_y = height;
	else max_y = y + r;

	if ((y - r) < 0) min_y = 0;
	else min_y = y - r;
	
    //Pick the pixels from a circle that lies inside the rectancle.
	int temp_x, temp_y;
	for (auto i = min_x; i < max_x; i++){
		for (auto j = min_y; j < max_y; j++){
			temp_x = i - x;
			temp_y = j - y;
			if ((pow(temp_x, 2) + pow(temp_y, 2)) <= pow(r,2)) {
				pixelVector.push_back(i);
				pixelVector.push_back(j);
			}

		}

	}


	return pixelVector;

}

bool Calculate::belongsToCircle(sf::Vector2u loc, sf::Vector2u circleLoc, int r)
{
    int temp_x = (int)circleLoc.x - (int)loc.x;
	int temp_y = (int)circleLoc.y - (int)loc.y;

    if ((pow(temp_x, 2) + pow(temp_y, 2)) < pow(r,2)) 
        return true;
    
    return false;

}



sf::Vector2u Calculate::nextLocation(sf::Vector2f speed, sf::Vector2u location)
{
	location.x = location.x + speed.x*param::dt;
    location.y = location.y + speed.y*param::dt;
    return location;
}


sf::Vector2f Calculate::initialSpeed(int mass, float energy, float angle)
{
    return (float)pow(2 * energy / mass, 0.5)*sf::Vector2f(cos(M_PI / 180 * angle),-sin(M_PI / 180 * angle));
}

sf::Vector2f Calculate::newSpeed(sf::Vector2f speed, float wind, int mass)
{
    speed.x += param::windEffect * wind * param::dt / mass;
	speed.y += param::dt*param::gravity; 
    return speed;
}

void Calculate::wind(float& wind)
{
	// Values between -10 - 10
	wind = -10 + static_cast <float> (rand() / (static_cast<float>(RAND_MAX)))*20;   
}

float Calculate::angle(sf::Vector2f bulletSpeed)
{
	 return (float)(180 / M_PI * atan2(-bulletSpeed.y, bulletSpeed.x));
}
















