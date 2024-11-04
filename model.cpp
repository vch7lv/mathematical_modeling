#include "model.h"
#include <cmath>



Model::Model(
    function<ld(ld,ld)> dzdx_, 
    function<ld(ld,ld)> dzdy_, 
    function<ld(ld,ld)> beta_,
    ld alpha_,
    ld l_
    ) :  dzdx(dzdx_), dzdy(dzdy_), beta(beta_), alpha(alpha), l(l_)
{
    int n = ceil(l/dx);

    Points.resize(n+1);

    for (int i = 0; i <= n; ++i)
    {
        Points[i] = {i * l/n, 0};
    }
}


void Model::optimaze(size_t iter) 
{
    int n = Points.size() - 1;

    for (int t = 0; t < iter; ++t)
    {
        J_total_value = 0;

        for (int i = 0; i < n; ++i)
        {
            ld best_dy1 = 0;
            ld best_dy2 = 0;

            ld J_best_local = 1/0;

            for (ld dy1 = -dy; dy1 <= dy; dy1 += dy)
            for (ld dy2 = -dy; dy2 <= dy; dy2 += dy)
            {
                if (i == 0 && dy1 != 0) continue;
                if (i == n-1 && dy2 != 0) continue;

                auto [x1,y1] = Points[i];
                auto [x2,y2] = Points[i+1];

                ld a = (y2 + dy2 - (y1 + dy1)) / (x2 - x1);
                ld b = y1 + dy1 - a * x1;

                ld J_cur = J(x1,x2,a,b);

                if (J_cur < J_best_local)
                {
                    best_dy1 = dy1;
                    best_dy2 = dy2;
                    J_best_local = J_cur;
                }
            }

            J_total_value += J_best_local;

            Points[i].y += best_dy1;
            Points[i+1].y += best_dy2;
        }
    }
}

ld Model::J()
{
    int n = Points.size() - 1;

    ld out = 0;

    for (int i = 0; i < n; ++i)
    {
        ld a = (Points[i+1].y - Points[i].y) / (Points[i+1].x - Points[i].x);
        ld b = Points[i].y - a * Points[i].x;

        out += J(Points[i].x, Points[i+1].x, a, b);
    }

    return out;
}

ld Model::J(ld x1, ld x2, ld a, ld b) 
{
    function<ld(ld)> f1
    {
        [a,b, dzdx = this->dzdx, dzdy = this->dzdy](ld x)->ld
        {
            ld y = a*x + b;
            ld dzdx_value =  dzdx(x,y) + dzdy(x,y) * a;

            return sqrt( 1 + pow(a, 2) + pow(dzdx_value, 2) );
        }
    };

    function<ld(ld)> f2
    {
        [a,b, dzdx = this->dzdx, dzdy = this->dzdy, beta = this->beta](ld x)->ld
        {
            ld y = a*x + b;
            ld dzdx_value =  dzdx(x,y) + dzdy(x,y) * a;

            return beta(x,y) * sqrt( 1 + pow(a, 2) + pow(dzdx_value, 2) );
        }
    };

    return alpha * 0.5 * pow(integrate(f1, x1, x2), 2) + integrate(f2, x1, x2);
}


ld Model::dx {0.02};
ld Model::dy {0.02};

