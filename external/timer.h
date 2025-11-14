/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#include <chrono>

class Timer {
 public:
  using Clock = std::chrono::steady_clock;
  using TimePoint = Clock::time_point;

  Timer() : running(false) {}

  inline void start() noexcept {
    startTime = Clock::now();
    running = true;
  }

  inline void restart() noexcept { start(); }

  inline void stop() noexcept {
    endTime = Clock::now();
    running = false;
  }

  inline void reset() noexcept { running = false; }

  inline double elapsedNanoseconds() const noexcept {
    return duration<std::chrono::nanoseconds>();
  }

  inline double elapsedMicroseconds() const noexcept {
    return duration<std::chrono::microseconds>();
  }

  inline double elapsedMilliseconds() const noexcept {
    return duration<std::chrono::milliseconds>();
  }

  inline double elapsedSeconds() const noexcept {
    return duration<std::chrono::duration<double>>();
  }

 private:
  template <typename Dur>
  inline double duration() const noexcept {
    if (running) {
      return std::chrono::duration_cast<Dur>(Clock::now() - startTime).count();
    } else {
      return std::chrono::duration_cast<Dur>(endTime - startTime).count();
    }
  }

 private:
  TimePoint startTime{};
  TimePoint endTime{};
  bool running;
};
