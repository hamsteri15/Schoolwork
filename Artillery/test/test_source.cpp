#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>
#include <stdio.h> 
#include <stdlib.h>
#include "../src/calc.hpp"
#include "../src/weapon.hpp"
#include "../src/player.hpp"
#include "../src/map.hpp"
#include "../src/game.hpp"


TEST(Calculate, circlePixels) {
    std::cout << "------Testing calculate class-----" << std::endl <<std::endl;
	int r = 1;
    int x = 2; int y = 2;
    auto temp = Calculate::circlePixels(400,400,x,y,r);
    EXPECT_EQ(1, temp[0]);
    EXPECT_EQ(2, temp[1]);
    
}

TEST(Calculate, belongsToCircle) {
    
	sf::Vector2u p_loc(10,10);
    sf::Vector2u c_loc(10,10);
    EXPECT_EQ(true, Calculate::belongsToCircle(p_loc,c_loc,5));
    p_loc.x = 20;
    EXPECT_EQ(false, Calculate::belongsToCircle(p_loc,c_loc,5));    
}



TEST(Calculate, nextLocation) {

    sf::Vector2f speed(1,0);
    sf::Vector2u loc(1,1);
    auto temp = Calculate::nextLocation(speed,loc);
    EXPECT_EQ(loc.y, temp.y);
    EXPECT_EQ(int(loc.x+param::dt*speed.x), temp.x);
        
}


TEST(Calculate, initialSpeed) {

    int mass = 1;
    float energy = 1; float angle = 0;
    auto temp = Calculate::initialSpeed(mass,energy,angle);
    float correct =  pow(2*energy / mass, 0.5);
    EXPECT_EQ(correct, temp.x);
    EXPECT_EQ(0.f, temp.y);
    
}

TEST(Calculate, newSpeed) {

    int mass = 1;
    float wind = 3;
    sf::Vector2f speed(1,1);
    auto temp = Calculate::newSpeed(speed,wind,mass);
    speed.x += 3*param::windEffect*param::dt;     
    speed.y += param::dt*param::gravity;
    EXPECT_EQ(speed.x, temp.x);
    EXPECT_EQ(speed.y, temp.y);
    
}

TEST(Calculate, wind) {

    float wind = 1;
    Calculate::wind(wind);
    for (int i=0; i<10; i++){
        Calculate::wind(wind);
        EXPECT_EQ(true, wind <= 10);
        EXPECT_EQ(true, wind >= -10);
    }
    
}

TEST(Calculate, angle) {

    sf::Vector2f vec1(1,1);
    sf::Vector2f vec2(-1,-1);
    auto temp1 = Calculate::angle(vec1);
    auto temp2 = Calculate::angle(vec2);
    EXPECT_EQ((float)-45, temp1); //positive y-axis points down
    EXPECT_EQ((float)135, temp2);
    
}


TEST(Player, Constructor){

   std::cout << "------Testing player class-----" << std::endl <<std::endl;
   Player p = Player(1);
   sf::Vector2u correct = sf::Vector2u(0,0);
   EXPECT_EQ(correct, p.getLocation());
   EXPECT_EQ("Cannon", p.getEquippedWeapon()->getName());   

}

TEST(Player, setLocation){
   
   Player p = Player(1);
   p.setLocation(sf::Vector2u(10,10));
   EXPECT_EQ(sf::Vector2u(10,10), p.getLocation());
      

}

TEST(Player, createWeapons){
   
   Player p = Player(1);
   std::vector<Weapon*> wep = p.createWeapons();
   EXPECT_EQ("Cannon", wep[0]->getName());
   EXPECT_EQ("Railgun", wep[1]->getName());
   EXPECT_EQ("Teleport", wep[2]->getName());     

}

TEST(Player, changeWeapon){
   
   Player p = Player(1);
   p.changeWeapon(2);
   EXPECT_EQ("Teleport", p.getEquippedWeapon()->getName());      

}





TEST(Game, generalTests){

    //sf::RenderWindow window(sf::VideoMode(800, 600), "Artillery 3");
    

    //int players = 2;

    //Game game = Game(players,window);
    //game.newGame();
   
   

}
