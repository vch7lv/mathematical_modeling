#include "model.h"
#include <cmath>
#include <compare>
#include <iostream>

using std::cout;

int main() {

    double res = biv::makeCompoundSF(0, 1, [](long double x){ return x*x;});
    std::cout << res << "\n";
    return 0;
    ld alpha = 0.1;
    ld l = 1;

    function<ld(ld,ld)> dzdx
    {
        [](ld x, ld y)->ld
        {
            return 5 * sin(y) * cos(5*x);
        }
    };

    function<ld(ld,ld)> dzdy
    {
        [](ld x, ld y)->ld
        {
            return sin(5*x) * cos(y);
        }
    };

    function<ld(ld,ld)> beta 
    {
        [](ld x, ld y)->ld
        {
            return 0.5;
        }
    };

    Model model(dzdx, dzdy, beta, alpha, l);
    model.optimaze(50);

    cout << "J value: " << model.J_total_value << '\n';

    /*
        y:
        for (auto [x,y] : model.Points)
        {
            ...
        }
    */
}