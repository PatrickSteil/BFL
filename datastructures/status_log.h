#ifndef STATUS_LOG_H
#define STATUS_LOG_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

class StatusLog {
 public:
  StatusLog(const std::string& msg) {
    std::cout << msg << " ... " << std::flush;
    startTime = Clock::now();
  }

  StatusLog(const char* msg) {
    std::cout << msg << " ... " << std::flush;
    startTime = Clock::now();
  }

  ~StatusLog() {
    auto endTime = Clock::now();
    auto duration =
        std::chrono::duration<double, std::milli>(endTime - startTime).count();
    std::cout << "done [" << duration << "ms]" << std::endl;
  }

 private:
  using Clock = std::chrono::steady_clock;
  Clock::time_point startTime;
};

#endif
