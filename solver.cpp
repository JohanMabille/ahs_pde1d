#include <iostream>
#include "solver.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

void Solver::matrixes()
{
    // Implementation: second operands of x1 and x3 should have opposite signs
    // Design: what if we want to pass non constant volatility or rate?
    // Implementation: x * x is faster than std::pow(x, 2) and easier to read
    auto x1 = [&](double theta) { return -dt*theta*(pow(sigma.vol(),2)/(2*pow(dx,2)) + (pow(sigma.vol(),2))/(4*dx) - (r.rf())/(2*dx)); };
    auto x2 = [&](double theta) { return 1 + (pow(sigma.vol(),2)*theta*dt/pow(dx,2)) + r.rf()*theta*dt; };
    auto x3 = [&](double theta) { return -dt*theta*(pow(sigma.vol(),2)/(2*pow(dx,2)) + (pow(sigma.vol(),2))/(4*dx) - (r.rf())/(2*dx)); };
    

    std::vector <std::vector <double>> boundary({{x2(theta), x3(theta)}, {x1(theta), x2(theta)}});
    std::vector<double> interior({x1(theta), x2(theta), x3(theta)});
    X.fill_matrix(boundary, interior);
        
    std::vector <std::vector <double>> boundary_prime({{x2(theta-1), x3(theta-1)}, {x1(theta-1), x2(theta-1)}});
    std::vector<double> interior_prime({x1(theta-1), x2(theta-1), x3(theta-1)});
    Xprime.fill_matrix(boundary_prime, interior_prime);
        
    u[0] = boundaries[0]*(x1(theta-1)-x1(theta));
    u[M-2] = boundaries[1]*(x3(theta-1)-x3(theta));
    
}


Solver::Solver(double Spot, double T, double vol, double r, double theta, size_t M, size_t N, double (*payoff)(double), std::vector<double> boundaries): Spot(Spot), T(T), sigma(vol), r(r), theta(theta), price(0), _greeks(4, 0.), M(M), N(N), dx(10*sigma.vol()*sqrt(T)/M), dt(T/N), space(M+1, 10*sigma.vol()*sqrt(T)/M), time(N+1, dt), payoff(payoff), boundaries(boundaries), X(M-1), Xprime(M-1), u(M-1, 0)
{

    matrixes();
    
    space[0] = 0;
    time[0] = 0;

    std::partial_sum(space.begin(), space.end(), space.begin(), std::plus<double>());
    
    std::transform(space.begin(), space.end(), space.begin(), std::bind2nd(std::plus<double>(), log(Spot)-5*sigma.vol()*sqrt(T))); // {log(spot)-c, ..., log(spot), ..., log(spot)+c}
    
    std::partial_sum(time.begin(), time.end(), time.begin(), std::plus<double>()); // {0, dt, 2dt, ... _N*dt}
    

}

Solver::Solver(): Spot(100), T(2), sigma(0.3), r(0.05), theta(0.5), price(0.), _greeks(4, 0.), payoff(payoff), M(100), N(T*365), dx(10*sigma.vol()*sqrt(T)/M), dt((double)1/365), space(M+1, 10*sigma.vol()*sqrt(T)/M), time(N+1, dt), boundaries({0, Spot*exp(5*sigma.vol()*sqrt(T))-100}), X(M-1), Xprime(M-1), u(M-1, 0)
{
    matrixes();
    space[0] = 0;
    time[0] = 0;
    

    std::partial_sum(space.begin(), space.end(), space.begin(), std::plus<double>());
    std::transform(space.begin(), space.end(), space.begin(), std::bind2nd(std::plus<double>(), log(Spot)-5*sigma.vol()*sqrt(T)));
    
    std::partial_sum(time.begin(), time.end(), time.begin(), std::plus<double>());
    

    
}



std::vector<double> Solver::right_side(const std::vector<double> &f) const
{
    std::vector<double> res = Xprime*f;
    std::transform(res.begin(), res.end(), u.begin(), res.begin(), std::plus<double>());
    
    return res;
}

void Solver::pricing()
{

    std::vector<double> f(M-1, 0);
    std::cout << "space size = " << space.size() << std::endl;
    std::cout << "f size = " << f.size() << std::endl;
    std::transform(space.begin(), space.end(), f.begin(), [&](double arg1) { return (*payoff)(exp(arg1)); });
    
    // Implementation: you can handle boundary conditions more easily: given the system A X  = B,
    // fill A and B for i in (1, S-2). Then directly fill A(0, 0), A(0, 1), B(0), A(S-1, S-2), A(S-1, S-1) and
    // B(S-1) depending on the boundary conditions. For instance, for Dirichlet condition, you have:
    // A(0, 0) = A(S-1, S-1) = 1;
    // A(0, 1) = A(S-1, s-2) = 0;
    // B[0] = previous_sol[0] and B[S-1] = previous_sol[S -1]
    std::cout << "Segfault in the previous line becausespace and f have different sizes" << std::endl;
    for(std::size_t i=0; i<N-1; i++)
    {
        f = thomas_algo(X, right_side(f)); // compute f(n) knowing f(n+1)
    }
    
    std::vector <std::vector <double>> res(2); // store 2 price vector in time to compute the theta later
    res[1] = {f[(M/2)-2], f[(M/2)-1], f[(M/2)]}; // prices around the forward at day 1
    f = thomas_algo(X, right_side(f));
    res[0] = {f[(M/2)-2], f[(M/2)-1], f[(M/2)]}; // prices around the forward at day 0 (today)
    
    price = f[(M/2)-1];
    greeks(res);
}

double Solver::print_price() const
{
    return price;
}

void Solver::greeks(const std::vector <std::vector <double>> &f)
{
    double h1 = exp(space[(M/2)+2])-exp(space[(M/2)+1]);
    double h2 = exp(space[(M/2)+1])-exp(space[M/2]);
    
    _greeks[0] = (f[0][2]-f[0][0])/(h1+h2); // delta
    _greeks[1] = 2*(f[0][2]+f[0][0]-2*f[0][1]-_greeks[0]*(h1-h2))/(pow(h1,2)+pow(h2,2)); // gamma
    _greeks[2] = f[1][1]-f[0][1]; // theta
    _greeks[3] = _greeks[1]*pow(Spot,2)*sigma.vol()*T; // vega
}


std::vector <double> Solver::print_greeks() const
{
    return _greeks;
}


std::vector<double> thomas_algo(const mat_Tridiagonal &M, const std::vector<double> &y)
{
    size_t N = M.dim();
    std::vector<double> alpha(N, 0);
    std::vector<double> beta(N, 0);
    std::vector<double> x(N, 0);
    
    alpha[1] = -M(0,1)/M(0,0);
    beta[1] = y[0]/M(0,0);
    for (std::size_t i=1; i<N-1; i++)
    {
        alpha[i+1] = -M(i,i+1)/(M(i,i-1)*alpha[i]+M(i,i));
        beta[i+1] = (y[i]-beta[i]*M(i,i-1)) / (M(i,i-1)*alpha[i]+M(i,i));
    }
    
    x[N-1] = (y[N-1]-beta[N-1]*M(N-1,N-2)) / (M(N-1,N-2)*alpha[N-1]+M(N-1,N-1));
    for (std::size_t i=N-1; i>=1; i--)
    {
        x[i-1] = alpha[i]*x[i] + beta[i];
    }
    
    return x;
}

// Design: should accept strike as parameter
double payoff(double S)
{
    return std::max(S-100, 0.);
}



