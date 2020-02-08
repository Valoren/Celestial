#include <iostream>
#include <chrono>

class Timer
{

public:
      Timer();
      ~Timer();
     
      void Stop(){
          auto endTimepoint = std::chrono::high_resolution_clock::now();
          auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
          auto stop = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
          auto duration = stop - start;
          double ms = duration * 0.001;

          std::cout << duration << "us (" << ms << "ms)\n";
}
private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
};

Timer::~Timer(){

          Stop();
}

Timer::Timer(){

         m_StartTimepoint = std::chrono::high_resolution_clock::now();
}

int main(){
    int value = 0;
    {
        Timer timer1;
        for (int i=0; i<1000000; i++)
            value +=2;

    }
    std::cout << value << std::endl;

}