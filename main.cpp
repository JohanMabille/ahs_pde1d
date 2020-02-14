#include <ctime>
#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "solver.hpp"
#include "closed_form.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

int main(int argc, const char * argv[]) {
    

    double bs_value = dauphine::bs_price(100, 100, 0.3, 2, true);
    std::cout << "Constructor" << std::endl;
    Solver call_test;
    std::cout << "Constructor done" << std::endl;
    call_test.pricing();
    double value = call_test.print_price();
    std::vector<double> greeks = call_test.print_greeks();
        
    std::cout << "Value as of today " << " " << value << std::endl;
    std::cout << "Error of pricing : " << " " << abs(bs_value-value) << std::endl << std::endl << std::endl;
    std::cout << "Greeks :" << std::endl;
    std::cout << "1) Delta : " << " " << greeks[0] << std::endl;
    std::cout << "2) Gamma : " << " " << greeks[1] << std::endl;
    std::cout << "3) Theta : " << " " << greeks[2] << std::endl;
    std::cout << "4) Vega : " << " "  << greeks[3] <<std::endl;


    
return 0;
}
