#include <cmath>

namespace hydro {

using namespace std;

// https://www.e-education.psu.edu/natureofgeoinfo/c7_p9.html
float CalcInterpolate( const float *chart, float Qp, float Qh)
{
    struct {
        float p, h, e;
    } I[3];
    float Dist[3] = {1e+10, 1e+10, 1e+10};
    const float *pf = chart;
    float minP = *(pf++), minH = *(pf++);
    float deltaP = *(pf++) - minP, deltaH = *(pf++) - minH;
    while (*pf > 0) // for all h
    {
        float h = *(pf++);
        while (*pf > 0) // for all p,e in h
        {
            float p = *(pf++);
            float e = *(pf++);

            // farthest-elimination
            int iFarthest = 0;
            if (Dist[1] > Dist[iFarthest]) iFarthest = 1;
            if (Dist[2] > Dist[iFarthest]) iFarthest = 2;

            // normalize the coordinate scales
            float Ndh = (Qh - h) / deltaH;
            float Ndp = (Qp - p) / deltaP;
            float d = sqrt(Ndh * Ndh + Ndp * Ndp);
            if (d < Dist[iFarthest]) {
                Dist[iFarthest] = d;
                I[iFarthest].p = p;
                I[iFarthest].h = h;
                I[iFarthest].e = e;
            }
        }
        pf++;
    }
    float N = 0, D = 0;
    for (int i = 0; i < 3; i++)
    {
        float p = I[i].p, h = I[i].h, e = I[i].e;

        float dh = (Qh - h);
        float dp = (Qp - p);
        float d = sqrt(dh * dh + dp * dp);
        N += d > .001 ? e / d : e;
        D += d > .001 ? 1 / d : 1;
    }
    return N / D;
}

// http://geomalgorithms.com/a03-_inclusion.html
bool CalcContains( const float* poly, float Qp, float Qh )
{
    int wn = 0;    // the winding number counter

    const float* pf = poly;
    while( true )
    {
        float x0 = *(pf+0), y0 = *(pf+1);
        float x1 = *(pf+2), y1 = *(pf+3);
        float marker = x1;
        if( marker < 0 )
        {
            x1 = *(poly+0); y1 = *(poly+1); // wrap
        }

        auto isLeft = [&] () -> float
        {
            return (x1 - x0) * (Qh - y0) - (Qp - x0) * (y1 - y0);
        };

        if( y0 <= Qh ) {
            if( y1 > Qh )
                if( isLeft() > 0.f ) // -ve side of upward vector
                    ++wn;
        } else {
            if( y1 <= Qh )
                if( isLeft() < 0.f ) // -ve side of downward vector
                    --wn;
        }
        if( marker < 0 ) break; // exit test after bookkeeping

        pf += 2;
    }
    return (wn != 0);
}

void CalcSpan( float& min, float& span, float* poly, float head )
{
    int wn = 0;    // the winding number counter
    float spanL = HUGE_VALF, spanR = HUGE_VALF;
    const float avgP = poly[0];
    float* pf = &(poly[1]);
    while( true )
    {
        float x0 = *(pf+0), y0 = *(pf+1);
        float x1 = *(pf+2), y1 = *(pf+3);
        float marker = x1;
        if( marker < 0 )
        {
            x1 = *(poly+0); y1 = *(poly+1); // wrap
        }

        // what side of the vector (x0,y0)-(x1,y1) is the origin of a positive ray (Qp,Qh)-(+inf,Qh)?
        // +ve means left side of vector, -ve means right side of vector (from vector's reference point)
        // assumes: 1st quadrant
        auto isLeft = [&] () -> float
        {
            return (x1 - x0) * (head - y0) - (avgP - x0) * (y1 - y0);
        };

        // find smallest span from Qp,Qh to Px,Py where:
        // P is determined by similar triangle law, and
        // span is distance from P to Q
        auto checkSpan = [&] (float& sp, float piy /* always positive */, float pdy /* always positive */) -> void
        {
            float p = pdy < .001f ? 0 : (piy / pdy);
            float D = fabs( avgP - (x0 + p * (x1 - x0)) );
            if( D < sp ) sp = D;
        };

        // handle L and R span checking for point in a concave poly
        // assumes: clockwise poly, 1st quadrant, interior point
        if( y0 <= head )
        {
            if( y1 > head ) // upward pointing, Q in range
            {
                if( isLeft() < 0.f ) // +ve side of upward vector
                    checkSpan( spanL, head - y0, y1 - y0 );
                else
                    ++wn;
            }
        }
        else
        {
            if( y1 <= head ) // downward pointing, Q in range
            {
                if( isLeft() < 0.f ) // -ve side of downward vector
                {
                    checkSpan( spanR, y0 - head, y0 - y1 );
                    --wn;
                }
            }
        }
        if( marker < 0 ) break; // exit test after bookkeeping

        pf += 2;
    }

    if( wn == 0 )
    {
        min = span = 0;
        return;
    }

    min = avgP - spanL;
    span = spanL + spanR;
}

}
