//-----------------------------------------------------------------------------//
// forward_euler3.C
// -----------------------------------------------------------------------------
#include  <iostream>
#include  <cmath>
#include  "include/integration.h"
using namespace std;
using namespace Two_Body_Algorithms;

int main(){
    calculate_single_body_acceleration();
    euler_forward();
    return 0;
}