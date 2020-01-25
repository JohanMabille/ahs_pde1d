#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <numeric>
#include <string>
 
 mat_Tridiagonal::mat_Tridiagonal(size_t N): M(N)
{
    mat = new double *[M]; 
    for (std::size_t i=0; i<M; i++)
    {
        mat[i] = new double[M]; 
        std::fill(&mat[i][0], &mat[i][M], 0);
}
}



mat_Tridiagonal::mat_Tridiagonal(mat_Tridiagonal const&  temp): M(temp.M)
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
    std::vector<double> res(M);
    for(std::size_t i=0; i<M; i++)
    {
        res[i] = this->operator()(i,j);
    }
    return res;
}
    
    

