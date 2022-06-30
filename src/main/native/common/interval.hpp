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
