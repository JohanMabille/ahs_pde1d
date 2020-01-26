#include "rate.hpp"

Rate::Rate(): r(0.0)
{
    
}

Rate::Rate(double r): r(r)
{
    
}

double Rate::rf() const
{
    return r;
}
