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
//
// Created by yangqi735 on 2021/3/5.
//

#ifndef INTERVAL_HPP
#define INTERVAL_HPP


template <typename TargetType, typename PosType>
class Interval {
public:
    Interval() = default;
    Interval(TargetType tid_, PosType start_, PosType end_):
        tid(tid_), start(start_), end(end_) {}
    Interval(const Interval &other) = default;
    Interval &operator=(const Interval &other) = default;
    Interval(Interval &&other) noexcept = default;
    Interval &operator=(Interval &&other) noexcept = default;
    virtual ~Interval() = default;

    TargetType tid = -1;
    PosType start = -1;
    PosType end = -1;
};


#endif //INTERVAL_HPP
