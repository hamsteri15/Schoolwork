#include "map.hpp"
#include "calc.hpp"

 
Map::Map(sf::Vector2u size)
{

    m_size = size;
    
    //these textures are not used yet
    sf::Texture rock;
    rock.loadFromFile("../images/rock.jpg");
    m_rockTexture = rock;    
    
    sf::Texture air;
    air.loadFromFile("../images/air.jpg");    
    m_airTexture = air;
    
    //create the initial land shape as an image
    m_image = randomizeLand();

    //create a texture of the map_image for drawing
    sf::Texture temp;
    temp.loadFromImage(m_image);
    m_texture = temp;

}



sf::Image Map::randomizeLand()
{

    //create a vector of land heights at certain x-coordinate
    std::vector<unsigned int> landHeights;
	std::vector<int> coeff;
	unsigned int n = 2; //number of terms in series
	double pi = 3.14159;

	srand(time(NULL));
	coeff.push_back(rand() % 100 + m_size.y / 2);
	for (unsigned int j = 1; j < n; j++) {
		coeff.push_back(rand() % (300 - 50*j) - (150 - 25 * j)); //randomize the coefficients
	}

	double randomHeight;
	for (unsigned int i = 0; i < m_size.x; i++) {
		randomHeight = coeff[0];
		for (unsigned int j = 1; j < n; j++) {
			randomHeight += coeff[j] * sin(3*j*pi*i / m_size.x); // calculate height from fourier series 
		}

		landHeights.push_back(int(randomHeight));
	}
    
    //create a transparent image of screen size
    //this land image must be updated according to the
    //land shape updates.
    sf::Image mapImage;
    mapImage.create(m_size.x, m_size.y, sf::Color::Transparent);

    //color the image based on land heights
    for (unsigned int i=0; i < m_size.x; i++){
        for (unsigned int j=0; j < m_size.y; j++){
            
            if (j > landHeights[i]) 
                mapImage.setPixel(i, j, sf::Color(130 + std::rand() % 20, 80 + std::rand() % 20, 10 + std::rand() % 20));
                              
            else
                mapImage.setPixel(i, j, sf::Color(135, 206, 250));

        }

    }
    
    return mapImage;    

}

void Map::makeHole(sf::Vector2u loc, int r)
{
    //get a vector of pixels that represent the hole locations
    auto pixels = Calculate::circlePixels(m_size.x, m_size.y, loc.x, loc.y,r);
    
    //color the pixels with the background color
    for (unsigned int i = 0; i < pixels.size() - 1; i+=2){
        m_image.setPixel(pixels[i],pixels[i+1], sf::Color(135, 206, 250));
    }
    
    m_texture.loadFromImage(m_image);

}


bool Map::isLand(sf::Vector2u loc)
{
	if (!outOfBounds(loc)) return (m_image.getPixel(loc.x,loc.y) != sf::Color(135, 206, 250));
	return false;
}

sf::Vector2u Map::getSize()
{
    return m_size;
}

sf::Texture Map::getTexture()
{
    return m_texture;
}

bool Map::outOfBounds(sf::Vector2u loc)
{
    // Warning: unsigned will never be < 0
	return (/*loc.x < 0 ||*/ loc.x > m_size.x || /*loc.y < 0 ||*/ loc.y > m_size.y);
}


void Map::draw(sf::RenderWindow& window)
{
    auto tempVec = sf::Vector2f((float)m_size.x, (float)m_size.y);
    sf::RectangleShape rectangle(tempVec);
    rectangle.setTexture(&m_texture);

    window.draw(rectangle);
}



