// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license
//
// Homebrew OpenMP threadsafe Tausme88 PRNG
// http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps

#ifndef PSYCHICSNIFFLE_TAUS88_H
#define PSYCHICSNIFFLE_TAUS88_H

#include <omp.h>

namespace util {

/*
// instance only one of these in the main application thread:
Taus88State taus88State; // declare omp state
srand(int(time(NULL))); // seed the single-thread state
taus88State.seed(); // seed the omp state from the single-thread state
*/
struct Taus88State
{
    uint32_t *block = 0;

    Taus88State( const Taus88State& other ) = delete;
    Taus88State& operator=( Taus88State& other ) = delete;
    Taus88State& operator=( const Taus88State& other ) = delete;

#ifdef _OPENMP
    int ompMaxThreads() { return omp_get_max_threads(); }
    int ompThreadNum() { return omp_get_thread_num(); }
#else
    int ompMaxThreads() { return 1; }
    int ompThreadNum() { return 1; }
#endif

    void copyOut( uint32_t* stale )
    {
        int m = ompMaxThreads();
        int t = ompThreadNum();
        stale[0] = block[t*m+0];
        stale[1] = block[t*m+1];
        stale[2] = block[t*m+2];
        stale[3] = block[t*m+3];
    }

    void copyIn( uint32_t* dirty )
    {
        int m = ompMaxThreads();
        int t = ompThreadNum();
        block[t*m+0] = dirty[0];
        block[t*m+1] = dirty[1];
        block[t*m+2] = dirty[2];
        block[t*m+3] = dirty[3];
    }

    void seed()
    {
        int m = ompMaxThreads();
        for( int i=0; i<m*4; i++ )
            block[i] = ((uint32_t)rand() << 8) + (uint32_t)rand(); // 32 bits please
    }

    Taus88State()
    {
        int m = ompMaxThreads();
        block = new uint32_t[m*4];
    }
    virtual ~Taus88State()
    {
        delete[] block;
    }
};

/*
// instance one of these in each omp thread:
#pragma omp parallel // declare a parallel block
{
    Taus88 taus88(taus88State); // each thread will get it's own taus88 object, independently initialized
#pragma omp parallel for
    ...
} // when this scope closes, each taus88 object state is written back to the global omp state. no reseeding required.
 */
struct Taus88
{
    Taus88State& master;
    uint32_t state[4];

    // uncopyable and unassignable
    Taus88() = delete;
    Taus88( const Taus88& other ) = delete;
    Taus88& operator=( Taus88& other ) = delete;
    Taus88& operator=( const Taus88& other ) = delete;

    // initialize local state from master block
    Taus88( Taus88State &master_ ) : master(master_)
    {
        master.copyOut(state);
    }

    // restore local state to master block
    virtual ~Taus88()
    {
        master.copyIn(state);
    }

    // permute local state each time () is called
    uint32_t operator()()
    {
        state[3] = (((state[0] << 13) ^ state[0]) >> 19);
        state[0] = (((state[0] & 0xFFFFFFFE) << 12) ^ state[3]);
        state[3] = (((state[1] << 2) ^ state[1]) >> 25);
        state[1] = (((state[1] & 0xFFFFFFF8) << 4) ^ state[3]);
        state[3] = (((state[2] << 3) ^ state[2]) >> 11);
        state[2] = (((state[2] & 0xFFFFFFF0) << 17) ^ state[3]);
        return state[0] ^ state[1] ^ state[2];
    }
};

}

#endif //PSYCHICSNIFFLE_TAUS88_H
