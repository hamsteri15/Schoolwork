#ifndef weapon_hh
#define weapon_hh

#include <SFML/Graphics.hpp>


class Weapon
{
public:

    
    
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Class constructor.
    ///    \param name: name of the weapon
    ///    \param radius: effect radius of the weapon
    ///    \param dam: damage exerted on the player when hit
    ///    \param weight: weight of the bullet (used in calculations)
    ///    \param ammoCount: number of bullets in the weapon
    ///
    ////////////////////////////////////////////////////////////
    Weapon(std::string name, int radius, int damage, int weight, int ammoCount) 
    : m_name(name)
    , m_effectRadius(radius)
    , m_damage(damage)
    , m_ammoWeight(weight)
    , m_ammoCount(ammoCount)
	, m_maxAmmo(ammoCount)
    , m_texture(setTexture())
    {}
    
    ////////////////////////////////////////////////////////////
    ///
    ///    \Default destructor.
    ///     
    ///
    ////////////////////////////////////////////////////////////
    ~Weapon() = default;

	////////////////////////////////////////////////////////////
	///
	///    \Increments ammo counter down by one.
	///
	////////////////////////////////////////////////////////////  
	void reduceAmmoCount();
    
	////////////////////////////////////////////////////////////
	///
	///    \Reloads the weapon.
	///
	////////////////////////////////////////////////////////////  
	void reload();

    ////////////////////////////////////////////////////////////
    ///
    ///    \Loads the weapon texture from file.
    ///
    ////////////////////////////////////////////////////////////        
    sf::Texture setTexture();

    ////////////////////////////////////////////////////////////
    ///     
    ///   \Obvious getter functions. Nothing special here.
    ///
    ////////////////////////////////////////////////////////////
	std::string getName() const;
	//
    int getEffectRadius() const;
    //
    int getDamage() const;
    //
    int getAmmoWeight() const;
    //
    int getAmmoCount() const;
    //
    sf::Texture& getTexture();

private:
    std::string m_name;
    int         m_effectRadius;
    int         m_damage;
    int         m_ammoWeight;
    int         m_ammoCount;
	int			m_maxAmmo;

    sf::Texture m_texture;
    
};

#endif
