
#include <iostream>
#include "structures.h"

//using namespace std;

using namespace solar_system;

int main(){
    /* 
    body Jupiter;
    Vector a(1, 1, 1);
    cout <<  a << endl;
    Vector Jupiter_vel(1, 4, 5);
    //Vector vel (1.0, 2.0, 3.0);
    cout << Jupiter.vel << endl;
    return 0;
    */
    //cout << solar_system::m1.name << "has mass" << solar_system::m1.mass << endl;
   body Jup;
   Jup.pos = {1, 4, 5};
   std::cout << "Jupiter position is"<< m9.pos.x <<std::endl;
   return 0;
}
