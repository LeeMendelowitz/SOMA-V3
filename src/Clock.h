#ifndef CLOCK_H
#define CLOCK_H

#include <ctime>

class Clock
{
    public:
    Clock() {
        start();
    }

    void start() {
        start_ = std::clock();
        lap_ = start_;
    }

    float lap() {
        clock_t now = std::clock();
        float delta = ((float) (now - lap_))/CLOCKS_PER_SEC;
        lap_ = now;
        return delta;
    }

    float elapsed() { 
        clock_t now = std::clock();
        float delta = ((float) (now - start_))/CLOCKS_PER_SEC;
        return delta;
    }

    private:
    clock_t start_;
    clock_t lap_;
};

#endif
