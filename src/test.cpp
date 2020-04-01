
#include <iostream>
#include "structures.h"
#include "integration.h"
#include "benchmark.h"
//using namespace std;

using namespace Integration_Algorithms;

int main(){

int value=0;
{

Timer timer;
for(int i=0; i<100000000; i++){
    value += i;
}
}
std::cout<< "En result:" << value << std::endl;
return 0;

}
