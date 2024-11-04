#include "model.h"
#include <cmath>
#include <compare>

int main()
{
    ld alpha = 0.1;
    ld l = 1;

    function<ld(ld,ld)> phi 
    {
        [](ld x, ld y)->ld
        {
            return sin(5*x) * sin(y);
        }
    };

    function<ld(ld,ld)> beta 
    {
        [](ld x, ld y)->ld
        {
            return 0.5;
        }
    };


}