#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>


std::ostream& operator<<(std::ostream& res, const mat_Tridiagonal &X)
{
    size_t N = X.dim();

    for(std::size_t i=0; i<N; i++)
    {
        for(std::size_t j=0; j<N; j++)
        {
            res << X(i, j) << ", ";
        }
        res << std::endl;
    }
    return res;
}

//mat_Tridiagonal& mat_Tridiagonal::operator+=(const mat_Tridiagonal& rhs)
    //{
        //for(std::size_t i = 0; i < rhs.dim(); ++i)
        //{
            //for(std::size_t j = 0; j < rhs.dim(); ++j)
            //{
                //mat[i * rhs.dim()+j] = mat[i * rhs.dim() + j] + rhs.mat[i * rhs.dim() + j];
            //}
        //}        

    //}

// Implementation: no need to return a const vector, since you return by value,
// let th user do what he wants with it
const std::vector<double> operator* (const mat_Tridiagonal &X, const std::vector<double> &y)
{
    std::vector<double> mult(X.dim(), 0);
    
    for (std::size_t i=0; i<X.dim(); i++)
    {
        for (std::size_t j=0; j<X.dim(); j++)
        {
            std::vector<double> store(X.dim());
            std::transform(X.mat[i], X.mat[i]+X.dim(), y.begin(), store.begin(), [](double val1, double val2) { return val1 * val2; });
            mult[i] = std::accumulate(store.begin(), store.end(), 0., std::plus<double>());
        }
    }
    return mult;
}



 mat_Tridiagonal::mat_Tridiagonal(size_t N): M(N)
{
    mat = new double *[M]; 
    for (std::size_t i=0; i<M; i++)
    {
        mat[i] = new double[M]; 
        std::fill(&mat[i][0], &mat[i][M], 0);
}
}



mat_Tridiagonal::mat_Tridiagonal(mat_Tridiagonal const &temp): M(temp.M)
{
        {
            mat = new double *[M];
            for (std::size_t i=0; i<M; i++)
            {
                mat[i] = new double[M];
                for (std::size_t j=0; j<M; j++)
                {
                    *(mat[i]+j) = temp(i,j);
                }
            }
        }
}





mat_Tridiagonal::~mat_Tridiagonal()
{
    for (std::size_t i=0; i<M; i++)
    {
        delete mat[i];
    }
    
    delete [] mat;
}



const double& mat_Tridiagonal::operator()(std::size_t i, std::size_t j) const
{
    return *(mat[i]+j);
}



void mat_Tridiagonal::fill_matrix(std::vector <std::vector <double>> boundary, std::vector<double> interior)
{
    for (std::size_t i=0; i<M; i++)
    {
        if(i==0 || i==M-1)
        {
            if(i==0)
            {
                *(mat[i]) = boundary[0][0];
                *(mat[i] + 1) = boundary[0][1];
            }
            if(i==M-1)
            {
                *(mat[i] + M-2) = boundary[1][0];
                *(mat[i] + M-1) = boundary[1][1];
            }
        }
        else
        {
            *(mat[i]+i-1) = interior[0];
            *(mat[i]+i) = interior[1];
            *(mat[i]+i+1) = interior[2];
        }
    }
}




std::size_t mat_Tridiagonal::dim() const
{
    return M;
}





std::vector<double> mat_Tridiagonal::col(std::size_t j) const
{
    std::vector<double> col(M);
    for(std::size_t i=0; i<M; i++)
    {
        col[i] = this->operator()(i,j);
    }
    return col;
}
    
    

