// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license
//
// performs random selections, without replacement
// uses a small exclude-list to hopefully increase cache-line usage

#ifndef PROJECT_NSELECTOR_H
#define PROJECT_NSELECTOR_H

#include "taus88.h"

namespace util
{

template<uint M>
struct NSelector
{
    const uint N;
    uint n[M];
    uint selected;

    NSelector(uint N_) : N(N_) {
        selected = 0;
    }

    void reset() { selected = 0; }

    uint select(Taus88 &fnRand)
    {
        if(selected == M)
            throw new std::runtime_error("too many selections");

        repick:
        uint a = fnRand() % N;
        for(int i=0; i<selected; i++)
        {
            if(a == n[i]) { goto repick; }
        }
        n[selected++] = a;
        return a;
    }
};

}

#endif //PROJECT_NSELECTOR_H
