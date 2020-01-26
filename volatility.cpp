#include "volatility.hpp"

Volatility::Volatility(): vol(0.3)
{
    
}

Volatility::Volatility(double vol): vol(vol)
{
    
}

double Volatility::vol() const
{
    return vol;
}
