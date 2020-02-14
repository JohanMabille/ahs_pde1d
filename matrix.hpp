#ifndef matrix_hpp
#define matrix_hpp

#include <vector>

class mat_Tridiagonal
{

public:
    
    mat_Tridiagonal(size_t M);
    mat_Tridiagonal(mat_Tridiagonal const& temp); 
    ~mat_Tridiagonal(); 
    // Implementation: missing value semantic
    // the following will crash: mat_Tridiagonal m1(10), m2(45); m1 = m2
    std::size_t dim() const;
    void fill_matrix(std::vector <std::vector <double>> boundary, std::vector<double> interior);
    const double& operator()(std::size_t i, std::size_t j) const; 
    std::vector<double> col(std::size_t i) const;
    //mat_Tridiagonal& operator+=(const mat_Tridiagonal& rhs);


private:
    
    // Design: Why using double** instead of std::vector?
    const size_t M; // dimension
    double **mat; // data matrix
    
    // operator() is here to prevent this
    friend const std::vector<double> operator* (const mat_Tridiagonal&, const std::vector<double>&); // otherwise can't use mat its in private
    
};

std::ostream& operator<<(std::ostream& res, const mat_Tridiagonal &X); // print matrix
const std::vector<double> operator* (const mat_Tridiagonal &X, const std::vector<double> &y); // multiplication matrix/vector

#endif
