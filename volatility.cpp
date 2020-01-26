#include "volatility.hpp"

Volatility::Volatility(): v(0.3)
{
    
}

Volatility::Volatility(double vol): v(vol)
{
    
}

double Volatility::vol() const
{
    return v;
}
