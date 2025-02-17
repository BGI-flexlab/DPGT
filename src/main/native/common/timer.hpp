/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
#ifndef DPGT_TIMER_H
#define DPGT_TIMER_H

#include <chrono>

/**
 * @brief A simple timer object for measuring time.
 *
 * This is a minimal class around a std::chrono::high_resolution_clock that
 * serves as a utility class for testing code.
 */
class Timer {
public:
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::nanoseconds ns;

    Timer() { Start(); }

    /**
    * @brief Starts a timer.
    */
    inline void Start() { start_time_ = clock::now(); }

    inline float NanoSeconds() {
        return static_cast<float>(
                std::chrono::duration_cast<ns>(clock::now() - start_time_).count());
    }

    /**
    * @brief Returns the elapsed time in milliseconds.
    */
    inline float MilliSeconds() { return NanoSeconds() / 1000000.f; }

    /**
    * @brief Returns the elapsed time in microseconds.
    */
    inline float MicroSeconds() { return NanoSeconds() / 1000.f; }

    /**
    * @brief Returns the elapsed time in seconds.
    */
    inline float Seconds() { return NanoSeconds() / 1000000000.f; }

protected:
    std::chrono::time_point<clock> start_time_;
};

#endif  // DPGT_TIMER_H