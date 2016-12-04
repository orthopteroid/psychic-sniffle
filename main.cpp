//
// Created by orthopteroid on 02/12/16.
//

#include <iostream>
#include <cstring>

#include "sniffle.h"

#include "schwefel.h"

int main()
{
    srand(int(time(NULL)));

    Schwefel<uint16_t, 5, 400>::Solve(10);

    return 0;
}
