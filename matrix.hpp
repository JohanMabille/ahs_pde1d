#ifndef matrix_hpp
#define matrix_hpp

#include <vector>

class mat_Tridiagonal
{

public:
    
    mat_Tridiagonal(size_t N);
    mat_Tridiagonal(mat_Tridiagonal const& temp); 
    ~mat_Tridiagonal(); 
    std::size_t dim() const;
    const double& operator()(std::size_t i, std::size_t j) const; // A(i,j)
    std::vector<double> col(std::size_t i) const; 

private:
    
    const size_t N; // dimension
    double **mat; // data matrix
    
    
    
};

const std::vector<double> operator* (const mat_Tridiagonal &X, const std::vector<double> &y); // multiplication matrix/vector
std::ostream& operator<<(std::ostream& out, const mat_Tridiagonal& m); // print matrix


#endif
