
#ifndef solver_hpp
#define solver_hpp

#include <string>
#include "matrix.hpp"
#include "volatility.hpp"
#include "rate.hpp"

class Solver
{
// X*f(n) = X'*f(n+1) + u
// with X' = X(theta -1)
public:
    
    void matrixes(); // Matrixes for the linear equation
    Solver();
    Solver(double Spot, double T, double sigma, double r, double theta, size_t M, size_t N, double (*payoff)(double), std::vector<double> boundaries);

    std::vector<double> right_side(const std::vector<double> &f) const; // X'f(n+1) becomes X'f(n+1)+u to include Dirichlet boundariew
    void pricing();
    double print_price() const;
    void greeks(const std::vector <std::vector <double>> &f);
    std::vector <double> print_greeks() const;
    
private:
    
    const double Spot;
    const double T;
    Volatility sigma;
    Rate r;
    const double theta;
    double price;
    std::vector <double> _greeks;
    double (*payoff)(double); // pointer to function (need to define it outside the class)
    
    const size_t M; // size of space set
    const size_t N; // size of time set
    const double dx; // space step
    const double dt; // time step
    std::vector<double> space;
    std::vector<double> time;
    
    
    std::vector<double> boundaries; // Dirichlet boundaries
    mat_Tridiagonal X; // matrix X(θ)
    mat_Tridiagonal Xprime; // matrix X(θ-1)
    std::vector<double> u; // boundaries
    
};

std::vector<double> thomas_algo(const mat_Tridiagonal &M, const std::vector<double> &y); // Method of Gauss for a Tridiagonal matrix (Mx = y | M tridiagonal)

double payoff(double S); // call with strike 100 by default

#endif /* solver_hpp */

