//-----------------------------------------------------------------------------//
// forward_euler3.C
// -----------------------------------------------------------------------------
#include  <iostream>
#include  <cmath>
#include  "include/integration.h"
using namespace std;
using namespace Two_Body_Algorithms;

int main(){
    Integrator body;
    body.calculate_single_body_acceleration(1);
    euler_forward();
    return 0;
}