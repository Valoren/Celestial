/*
Simple chunk of data for solar system planets. Main purpose is for benchmark tests.
Units given in SI standard system(kg, m, m/sec)
*/

namespace solar_system{

constexpr double G = 6.673e-11;

body m0;
void setup_m0(){
//Sunextern body m0;
m0.name = "Sun";
m0.mass = 1.989e30;
m0.pos.x = 0.0;
m0.pos.y = 0.0;
m0.pos.z = 0.0;

m0.vel.x = 0.0;
m0.vel.y = 0.0;
m0.vel.z = 0.0;
m0.radius = 695700000;
}
 
//Mercury
body m1;
void setup_m1(){
m1.name = "Mercury";
m1.mass = 3.3011e23;
m1.pos.x = -46000000000;
m1.pos.y = 0.0;
m1.pos.z = 0.0;

m1.vel.x = 0.0;
m1.vel.y = -58980: 
m1.vel.z = 0.0;
m1.radius = 2439700;
}

//Venus
body m2;
void setup_m2(){
m2.name = "Venus";
m2.mass = 4.8675e24;
m2.pos.x = -107480000000;
m2.pos.y = 0.0;
m2.pos.z = 0.0;

m2.vel.x = 0.0;
m2.vel.y = -35260;
m2.vel.z = 0.0;
m2.radius = 6051800;
}

//Earth
body m3;
void setup_m3(){
m3.name = "Earth";
m3.mass = 5.972e24;
m3.pos.x = -147095000000;
m3.pos.y = 0.0;
m3.pos.z = 0.0;

m3.vel.x = 0.0; 
m3.vel.y = -30300;
m3.vel.z = 0.0;
m3.radius = 6371000;
}
 
//Mars
body m4;
void setup_m4(){
m4.name = "Mars";
m4.mass = 6.4171e23;
m4.pos.x = -206620000000;
m4.pos.y = 0.0;
m4.pos.z = 0.0;

m4.vel.x = 0.0;
m4.vel.y = -26500;
m4.vel.z = 0.0;
m4.radius = 3389500;
}

//Jupiter
body m5;
void setup_m5(){
m5.name = "Jupiter";
m5.mass = 1898.19e24
m5.pos.x = -740520000000;
m5.pos.y = 0.0;
m5.pos.z = 0.0;

m5.vel.x = 0.0;
m5.vel.y = -13720;
m5.vel.z = 0.0;
m5.radius = 71492000;
}

//Saturn
body m6;
void setup_m6(){
m6.name = "Saturn";
m6.mass = 568.34e24;
m6.pos.x = -1352550000000;
m6.pos.y = 0.0;
m6.pos.z = 0.0;

m6.vel.x = 0.0;
m6.vel.y = -10180;
m6.vel.z = 0.0;
m6.radius = 54364000;
}

//Uranus
body m7;
void setup_m7(){
m7.name = "Uranus";
m7.mass = 86.813e24;
m7.pos.x = -2741300000000;
m7.pos.y = 0.0;
m7.pos.z = 0.0;

m7.vel.x = 0.0;
m7.vel.y = -7110;
m7.vel.z = 0.0;
m7.radius = 24973000;
}

//Neptune
body m8;
void setup_m8(){
m8.name = "Neptune";
m8.mass = 102.413e24;
m8.pos.x = -4444450000000;
m8.pos.y = 0.0;
m8.pos.z = 0.0;

m8.vel.x = 0.0;
m8.vel.y = -5500;
m8.vel.z = 0.0;
m8.radius = 24341000;
}

//Pluto
body m9;
void setup_m9(){
m9.name = "Pluto";
m9.mass = 1.30900e22;
m9.pos.x = 5906376272000;
m9.pos.y = 0.0;
m9.pos.z = 0.0;

m9.vel.x = 0.0;
m9.vel.y = -4760;
m9.vel.z = 0.0;
m9.radius = 1188300; 
}
}