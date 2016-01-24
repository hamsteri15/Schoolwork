#ifndef player_hh
#define player_hh

#include <SFML/Graphics.hpp>
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "map.hpp"
#include "weapon.hpp"
#include "calc.hpp"

class Weapon;
class Player
{
public:
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Class constructor.
    ///    \param number: player number     
    ///
    ////////////////////////////////////////////////////////////
    Player(int number);

    ////////////////////////////////////////////////////////////
    ///
    ///    \Class destructor. Frees the allocated memory for 
    ///      weapons.
    ///         
    ///
    ////////////////////////////////////////////////////////////
    ~Player() {
         for (auto it = m_weapons.begin(); it != m_weapons.end(); ++it){
            delete *it;
         }
        m_weapons.clear();
    }
            


    ////////////////////////////////////////////////////////////
    /// 
    ///    \Modifies the player location.
    ///    \param location: new location 
    ///
    ////////////////////////////////////////////////////////////
    void setLocation(sf::Vector2u location); 
    


    ////////////////////////////////////////////////////////////
    /// 
    ///    \Creates the weapons for the players. 
    ///    \Called from the constructor.
    ///
    ////////////////////////////////////////////////////////////
    std::vector<Weapon*> createWeapons();
    


	////////////////////////////////////////////////////////////
	/// 
	///    \Changes the equipped weapon. 
	///    \param weaponIndex: index of the weapon to be used.
	///
	////////////////////////////////////////////////////////////
	void changeWeapon(int weaponIndex);
   


	////////////////////////////////////////////////////////////
	/// 
	///    \Reduces the player's health by the amount of damage taken. 
	///    \param damage: amount of damage taken
	///
	////////////////////////////////////////////////////////////
	void reduceHealth(int damage);
    


    ////////////////////////////////////////////////////////////
	/// 
	///    \Prepares a shooting operation to be drawn from 
    ///     drawShot().
	///    \param aimed: vector, which indicates the initial shot power
	///    \param players: the players in the game (to check if shot hits one)
    ///    \param map: map of the game (to check if bullet hits ground)
    ///    \param wind: randomized wind coefficient
    ///
	////////////////////////////////////////////////////////////
    int shoot(sf::Vector2i aimed, std::vector<Player*> players, Map map, float wind);



    bool isShooting();
    ////////////////////////////////////////////////////////////
	/// 
	///    \Checks if the player is shooting.
    ///
	////////////////////////////////////////////////////////////



    bool isHit(sf::Vector2u loc);
    ////////////////////////////////////////////////////////////
	/// 
	///    \Checks if a location belongs to a rectangle drawn around
    ///     the player.
	///    \param loc: location to be checked
    ///
	////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////
	/// 
	///    \Draws the complete trajectory of the bullet. Special
    ///     effects are drawn from Game.
	///    \param window: a reference to the game window
    ///
	////////////////////////////////////////////////////////////
    void drawShot(sf::RenderWindow& window);



    ////////////////////////////////////////////////////////////
	/// 
	///    \Draws the player object and moves it to the ground
    ///     if necessary.
	///    \param window: a reference to the game window
    ///    \param map: a reference to the game map
    ///
	////////////////////////////////////////////////////////////
    void draw(sf::RenderWindow& window, Map& map);

    

    ////////////////////////////////////////////////////////////
    ///
    ///    \Obvious getter functions. Nothing special here.
    ///
    ////////////////////////////////////////////////////////////
    int getHealth() const;
    //
	sf::Color getColor() const;
    //
	sf::Texture getTexture() const;
    //
    sf::Vector2u getLastLocation() const;
    //    
    Weapon* getEquippedWeapon() const;
    //
    sf::Vector2u getLocation() const;
    //
    int getNumber() const;
private:
    
    sf::Vector2u            m_location;
    std::vector<Weapon*>     m_weapons;
    Weapon*                 m_equippedWeapon;
    int                     m_health;
    sf::Texture             m_texture;
    sf::Color               m_color;
    bool                    m_isShooting = false;
    int                     m_playerNumber;
    
    sf::Sprite                  m_bullet;
    std::vector<float>          m_bulletAngles;
    std::vector<sf::Vector2u>   m_bulletLocations;
    sf::Vector2u                m_lastLocation;
    int                         m_bulletCourseIndex = 0;
    int                         m_bulletCourseCount = 0;
    int                         m_targetTag;

    ////////////////////////////////////////////////////////////
	/// 
	///    \Prepares the shot trajectory. This is called from
    ///     shoot().
    ///    \param aimed: vector, which indicates the initial shot power
	///    \param players: the players in the game (to check if shot hits one)
    ///    \param map: map of the game (to check if bullet hits ground)
    ///    \param wind: randomized wind coefficient
	///
	////////////////////////////////////////////////////////////
    void setBulletLocs(sf::Vector2i aimed, std::vector<Player*> players, Map map, float wind);
    
};

#endif
