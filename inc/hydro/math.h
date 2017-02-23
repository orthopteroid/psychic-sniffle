#ifndef PSYCHICSNIFFLE_HYDRO_MATH_H
#define PSYCHICSNIFFLE_HYDRO_MATH_H

#include <sys/types.h>
#include <cmath>

namespace hydro {

using namespace std;

struct StatPosNeg
{
    float pos, neg;
    StatPosNeg() : pos(0), neg(0) {}
    inline void clear() { pos=neg=0.; }
    inline void inc( float v ) { if( v > 0 ) pos+=v; else neg+=-v; }
};

struct StatMinMax
{
    float min, max;
    StatMinMax() : min(HUGE_VALF), max(-HUGE_VALF) {}
    inline void clear() { min=HUGE_VALF; max=-HUGE_VALF; }
    inline void inc( float v ) { {if(v<min)min=v;} {if(v>max)max=v;} }
    inline float maximum() { return fabs(max) > fabs(min) ? max : min; }
};

struct StatStdDev
{
    float tot; uint num;
    StatStdDev() : tot(0), num(0) {}
    inline void clear() { tot=0.; num=0; }
    inline void inc( float v ) { tot += sqrt( v * v ); num++; }
    inline void incNZ( float v, const float tol = .01f ) { if( fabs( v ) > tol ) { inc( v ); } }
    inline float stdev() { return num == 0 ? 0 : tot / num; }
};

struct StatAvg
{
    float tot; uint num;
    StatAvg() : tot(0), num(0) {}
    inline void clear() { tot=0.; num=0; }
    inline void inc( float v ) { tot += v; num++; }
    inline void incGZ( float v, const float tol = .01f ) { if( v > tol ) { inc( v ); } }
    inline float avg() { return num == 0 ? 0 : tot / num; }
};

inline void Clamp(float &v, const float min, const float max)
{
    if( v< min ) v = min; else if( v > max ) v = max;
}

inline float CalcQ(float p, float e, float h, float units)
{
    return p / ((e / 100.f) * h * units);
}

inline float CalcE(float p, float q, float h, float units)
{
    return q < .1 ? 0.f : 100.f * p / (q * h * units);
}

// https://www.e-education.psu.edu/natureofgeoinfo/c7_p9.html
// uses inverse distce weighting
float CalcInterpolate( const float *chart, float Qp, float Qh );

// http://geomalgorithms.com/a03-_inclusion.html
bool CalcContains( const float* poly, float Qp, float Qh );

void CalcSpan( float& min, float& span, float* poly, float head );

}

#endif //PSYCHICSNIFFLE_HYDRO_MATH_H
