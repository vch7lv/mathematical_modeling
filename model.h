#pragma once
#include "integration.h"


#include <functional>
#include <map>
#include <vector>

using std::function, std::map, std::vector;
using ld = long double;

ld integrate(function<ld(ld)> fun, ld x1, ld x2) {
    ld result = biv::makeCompoundSF(x1, x2, fun);
    
}


struct Point
{
    ld x;
    ld y;
};


struct Model
{
    Model(
        function<ld(ld,ld)> dzdx_, 
        function<ld(ld,ld)> dzdy_, 
        function<ld(ld,ld)> beta_,
        ld alpha_,
        ld l_
        );


    void optimaze(size_t iter);

    //y = ax + b
    ld J();
    ld J(ld x1, ld x2, ld a, ld b);

    //calculated y:
    ld J_total_value = -1;
    vector<Point> Points;

    //task parameters
    function<ld(ld,ld)> dzdx;
    function<ld(ld,ld)> dzdy;
    function<ld(ld,ld)> beta; 

    ld alpha;
    ld l;

    //method parameters
    static ld dx;
    static ld dy;
};
