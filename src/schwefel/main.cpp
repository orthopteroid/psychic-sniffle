#include <iostream>
#include <cstring>

#include "schwefel/schwefel.h"

int main()
{
    srand(int(time(NULL)));

    Schwefel<uint16_t, 20, 400>::Solve(10);

    return 0;
}
