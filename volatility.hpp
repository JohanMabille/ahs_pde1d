#ifndef volatility_hpp
#define volatility_hpp


class Volatility
{
    
public:
    
    Volatility();
    Volatility(double vol);
    double vol() const;
    
private:
    
    const double vol;

};

#endif /* volatility_hpp */

