// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license

#include <iostream>
#include <cstring>

#include "cpuinfo.h"
#include "schwefel/schwefel.h"

int main()
{
    srand(int(time(NULL)));

#if defined(NDEBUG)
    {
        int cores = EnumCores();
        printf("OpenMP using %d threads\n", cores);
        omp_set_num_threads(cores);
    }
#endif

    Schwefel<uint16_t, 20, 400>::Solve(10);

    return 0;
}
