#ifndef rate_hpp
#define rate_hpp

class Rate
{
    
public:
    
    Rate();
    Rate(double r);
    double rf() const;
    
private:
    
    const double r;

};

#endif /* rate_hpp */
