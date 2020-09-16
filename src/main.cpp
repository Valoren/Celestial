/*
 * main.cpp
 * 
 * Copyright 2019 Miquel Bernat Laporta i Granados 
 * <mlaportaigranados@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

/*
Matplotlib link instructions:
Compile with folowing flags:
-I/usr/include/python2.7 -lpython2.7
*/

#include <iostream>
#include <cmath>
#include <omp.h>
#include <functional> 
#include <map>          
#include <memory>       
#include <string>      
#include <sstream>     
#include <string_view>  
#include <variant>
#include "include/structures.h"
#include "include/integration.h"
#include "include/planet_data.h"
#include "include/menu.h"
#include "include/benchmark.h"
#include "parser.h"
#include "astro_constants.h"
//#include "include/astro_epochs.h"
//#include "matplotlibcpp.h"

//STANDARD INTEGRATOR TEMPLATE
template <typename Integrator>
void run_simulation(Integrator integrator, int iterations, int report_frequency)
{
    for (auto i = 0; i < iterations; i++)
    {
        if (i % report_frequency == 0)
            record_state(integrator.get_bodies());
        integrator.compute_gravity_step();
    }
    output_states(integrator.get_bodies());
}

//STANDARD PARSER TEMPLATE
template <class Opts>
struct CmdOpts : Opts
{
    using MyProp = std::variant<std::string Opts::*, int Opts::*, double Opts::*, bool Opts::*>;
    using MyArg = std::pair<std::string, MyProp>;

    ~CmdOpts() = default;

    Opts parse(int argc, const char* argv[])
    {
        std::vector<std::string_view> vargv(argv, argv+argc);
        for (int idx = 0; idx < argc; ++idx)
            for (auto& cbk : callbacks)
                cbk.second(idx, vargv);

        return static_cast<Opts>(*this);
    }

    static std::unique_ptr<CmdOpts> Create(std::initializer_list<MyArg> args)
    {
        auto cmdOpts = std::unique_ptr<CmdOpts>(new CmdOpts());
        for (auto arg : args) cmdOpts->register_callback(arg);
        return cmdOpts;
    }

private:
    using callback_t = std::function<void(int, const std::vector<std::string_view>&)>;
    std::map<std::string, callback_t> callbacks;

    CmdOpts() = default;
    CmdOpts(const CmdOpts&) = delete;
    CmdOpts(CmdOpts&&) = delete;
    CmdOpts& operator=(const CmdOpts&) = delete;
    CmdOpts& operator=(CmdOpts&&) = delete;

    auto register_callback(std::string name, MyProp prop)
    {
        callbacks[name] = [this, name, prop](int idx, const std::vector<std::string_view>& argv)
        {
            if (argv[idx] == name)
            {
                visit(
                    [this, idx, &argv](auto&& arg)
                    {
                        if (idx < argv.size() - 1)
                        {
                            std::stringstream value;
                            value << argv[idx+1];
                            value >> this->*arg;
                        }
                    },
                    prop);
            }
        };
    };

    auto register_callback(MyArg p) { return register_callback(p.first, p.second); }
};

int main(int argc, char *argv[]){
    
    std::vector<body> bodies;
    const char *algorithms[3] = {"Euler", "RK4", "Leapfrog"}; 


    //Using solar system data in planet_data.h for benchmarking
    bodies.push_back(solar_system::sun);
    bodies.push_back(solar_system::mercury);
    bodies.push_back(solar_system::venus);
    bodies.push_back(solar_system::earth);
    bodies.push_back(solar_system::mars);
    bodies.push_back(solar_system::saturn);
    bodies.push_back(solar_system::jupiter);
    bodies.push_back(solar_system::uranus);
    bodies.push_back(solar_system::neptune);
    bodies.push_back(solar_system::pluto);
    
    struct MyOpts
    {
        std::string AlgorithmOpt{}; //Selects integrator from aviable scheme.
        int intOpt{}; //Number of integratio steps
        double errorOpt{}; //Numerical integration error tolerance
        bool boolOpt{}; //True/False for Debug flag
        bool boolOpt_test{}; //True/False Solar system test flag
        std::string filenameOpt{}; //File name with system data
    };
    //{"-tol", &MyOpts::errorOpt}
    auto parser = CmdOpts<MyOpts>::Create({
        {"--integrator", &MyOpts::AlgorithmOpt },
        {"--steps", &MyOpts::intOpt},
        {"--DEBUG", &MyOpts::boolOpt},
        {"--test", &MyOpts::boolOpt},
        {"--error", &MyOpts::errorOpt},
        {"--file", &MyOpts::filenameOpt}});

    auto myopts = parser->parse(argc, argv);
    /*
    std::cout << "stringOpt = " << myopts.AlgorithmOpt << std::endl;
    std::cout << "intOpt = " << myopts.intOpt << std::endl;
    std::cout << "boolOpt = " << myopts.boolOpt << std::endl;
    */

   //Main program execution
   spawn_title();
   number_of_cores();
   //spawn_menu();
   //parse_file();
   //parse_data(argv[1]);
   {
        Timer timer;
        if(myopts.filenameOpt == ""){
            std::cout << "Initiating system with data file" << myopts.filenameOpt << std::endl;

        }
        if(myopts.AlgorithmOpt == "RK4"){
            Orbit_integration::RK4 orbit(bodies, 0.01);
            run_simulation(orbit, (int)myopts.intOpt, 1);
        }
        else if(myopts.AlgorithmOpt == "Euler"){
            Orbit_integration::Euler orbit(bodies, 0.01);
            run_simulation(orbit, (int)myopts.intOpt,1);
        }
        else{
            std::cout << "Non defined integrator" << std::endl;
        }
    }
    std::cout << "Execution terminated!!" << std::endl;

    return 0;
}
