#include <iostream>
#include <complex>
#include <cmath>
#include "solver.h"

cSqMatrix::cSqMatrix(int dims, std::complex<double> *matrixData, std::complex<double> *vectorData)
{
    m_dims = dims;
    m_nElements = dims*dims;
    m_matrixData = new std::complex<double>[m_nElements];
    m_vectorData = new std::complex<double>[m_dims];
    for(int i = 0 ; i < m_nElements; i++)
    {
        m_matrixData[i] = matrixData[i];
    }
    for(int i = 0; i < m_dims; i++)
    {
        m_vectorData[i] = vectorData[i];
    }
    LU_decomposition(matrixData, vectorData, m_dims);
}

cSqMatrix::~cSqMatrix()
{
    if(m_matrixData != nullptr)
    {
        delete[] m_matrixData;
    }
    if(m_vectorData != nullptr)
    {
        delete[] m_vectorData;
    }
}

void cSqMatrix::LU_decomposition(std::complex<double> *pMatrix, std::complex<double> *pVector, int dim)
{
    std::complex<double> LUmatrix[dim][dim] = {0};
    std::complex<double> tempColumn[dim] = {0};
    std::complex<double> yVector[dim] = {0};
    std::complex<double> big = 0, temp;
    int bigIndx = 0;
    for(int j = 0; j < dim; j++)
    {
        for(int i = 0; i < dim; i++) // calculate column
        {
            if(i == 0)
            {
                tempColumn[i] = *(pMatrix + i*dim + j);
            } else if (i > 0 && i <=j)
            {
                temp = 0;
                for(int k = 0; k < i; k++)
                {
                    temp += LUmatrix[i][k] * tempColumn[k];
                }
                tempColumn[i] = *(pMatrix + i*dim + j) - temp;
            } else
            {
                temp = 0;
                for(int k = 0; k < j; k++)
                {
                    temp += LUmatrix[i][k] * tempColumn[k];
                }
                tempColumn[i] = *(pMatrix + i*dim + j) - temp;
            }
        }
        if(tempColumn[j] != static_cast<std::complex<double>>(0))
        {
            for(int i = 0; i < dim; i++)
            {
                if(i <= j)
                {
                    LUmatrix[i][j] = tempColumn[i];
                } else
                {
                    LUmatrix[i][j] = tempColumn[i]/tempColumn[j];
                }
            }
        } else
        {
            big = 0;
            for(int i = j; i < dim; i++)
            {
                if(abs(tempColumn[i]) > abs(big))
                {
                    big = tempColumn[i];
                    bigIndx = i;
                }
            }
            swap_two(tempColumn+j, tempColumn+bigIndx);
            swap_two(pVector+j, pVector+bigIndx);
            for(int i = 0; i < dim; i++)
            {
                swap_two(pMatrix + j*dim + i, pMatrix + bigIndx*dim + i);
                if(i <= j)
                {
                    LUmatrix[i][j] = tempColumn[i];
                }
                else
                {
                    LUmatrix[i][j] = tempColumn[i]/tempColumn[j];
                }
            }
        }
    }
    for(int i = 0; i < dim; i++) // calculate y vector
    {
        if(i == 0)
        {
            yVector[i] = *(pVector);
        } else
        {
            temp = 0;
            for(int j = 0; j < i; j++)
            {
                temp += LUmatrix[i][j]*yVector[j];
            }
            yVector[i] = (*(pVector+i) -temp);
        }
    }
    for(int i = dim-1; i >= 0; i--)
    {
        if(i == dim-1)
        {
            *(pVector+i) = yVector[i]/LUmatrix[i][i];
        } else
        {
            temp = 0;
            for(int j = i+1; j < dim; j++)
            {
                temp += LUmatrix[i][j]* (*(pVector+j));
            }
            *(pVector+i) = (yVector[i]-temp)/LUmatrix[i][i];
        }
    }
}

void cSqMatrix::swap_two(std::complex<double> *pt1, std::complex<double> *pt2)
{
    std::complex<double> temp = *pt1;
    *pt1 = *pt2;
    *pt2 = temp;
}


cEquations::cEquations(int dims, std::complex<double> *lambdas, std::complex<double> *energies, double g)
{
    eq_dims = dims;
    eq_g = g;
    eq_lambdaData = new std::complex<double>[dims];
    eq_energiesData = new std::complex<double>[dims];
    for(int i = 0; i < dims; i++)
    {
        eq_lambdaData[i] = lambdas[i];
        eq_energiesData[i] = energies[i];
    }
    newton_raphson(lambdas);
}

cEquations::~cEquations()
{
    if(eq_lambdaData != nullptr)
    {
        delete[] eq_lambdaData;
    }
    if(eq_energiesData != nullptr)
    {
        delete[] eq_energiesData;
    }
}

std::complex<double> cEquations::function_i(std::complex<double> *x_vector, int i_index)
{
    std::complex<double> sum = 0;
    for(int i =0; i < eq_dims; i++)
    {
        if(i != i_index)
        {
            sum+= (*(x_vector+i_index)-(*(x_vector+i)))/(eq_energiesData[i_index]-eq_energiesData[i]);
        }
    }
    return pow(*(x_vector+i_index),2)-(eq_g*sum)-(*(x_vector+i_index));
}

void cEquations::newton_raphson(std::complex<double> *init_guess)
{
    double PREC = 0.000001;
    std::complex<double> jacobian_matrix[eq_dims][eq_dims] = {0};
    std::complex<double> dx_vector[eq_dims] = {0};
    std::complex<double> delta_x[eq_dims] = {0};
    std::complex<double> temp_x_vector[eq_dims] = {0};
    std::complex<double> minus_f_vector[eq_dims] = {0};
    std::complex<double> difference_sqared = 0;
    int iteration = 1;
    do
    {
        //std::cout << iteration << '\n';
        for(int i = 0; i < eq_dims; i++)
        {
            // initializes temporary vector to work with
            temp_x_vector[i] = *(init_guess+i);
        }
        for(int i = 0; i < eq_dims; i++)
        {
            for(int j = 0 ; j < eq_dims; j++)
            {
                for(int k = 0; k < eq_dims; k++)
                {
                    // initializes x+dx vector
                    if(k == j)
                    {
                        dx_vector[k] = temp_x_vector[k] + PREC;
                    } else
                    {
                        dx_vector[k] = temp_x_vector[k];
                    }
                }
                // calculates jacobian matrix at a given initial guess point
                jacobian_matrix[i][j] = (function_i(dx_vector, i)-function_i(temp_x_vector,i))/PREC;

            }
        }
        //matrix_display(jacobian_matrix);
        for(int i = 0; i < eq_dims; i++)
        {
            // initializes b_vector
            minus_f_vector[i] = -function_i(init_guess, i);
        }
        cSqMatrix jaco_solver = cSqMatrix(eq_dims, &jacobian_matrix[0][0], &minus_f_vector[0]);
        for(int i = 0; i < eq_dims; i++)
        {
            delta_x[i] = minus_f_vector[i];
        }
        difference_sqared = 0; // resets "stopper"
        for(int i = 0; i < eq_dims; i++)
        {
            difference_sqared += pow(delta_x[i],2); // finds value that terminates loop
        }
        for(int i = 0; i < eq_dims; i++)
        {
            init_guess[i] += delta_x[i]; // changes initial values to a better guess
        }
        iteration ++;
    } while (sqrt(abs(difference_sqared)) > 0.00001); // loop terminates when desired precision is reached
}

cPolynomial::cPolynomial(int degree, std::complex<double> *lambdas, std::complex<double> *energies, double g)
{
    p_degree = degree;
    p_g = g;
    p_lambdaData = new std::complex<double>[degree];
    p_energiesData = new std::complex<double>[degree];
    p_coefficients = new std::complex<double>[degree+1];
    p_zeroes = new std::complex<double>[degree];
    for(int i = 0; i < p_degree; i++)
    {
        p_lambdaData[i] = lambdas[i];
        p_energiesData[i] = energies[i];
        p_zeroes[i] = 0;
    }
    coefficients();
    root_finder();
}

cPolynomial::~cPolynomial()
{
    if(p_lambdaData != nullptr)
    {
        delete[] p_lambdaData;
    }
    if(p_energiesData != nullptr)
    {
        delete[] p_energiesData;
    }
    if(p_coefficients != nullptr)
    {
        delete[] p_coefficients;
    }
    /*if(p_zeroes != nullptr)
    {
        delete[] p_zeroes;
    }*/
}

void cPolynomial::coefficients()
{
    std::complex<double> proper_lambda[p_degree] = {0};
    std::complex<double> a_matrix[p_degree][p_degree] = {0};
    std::complex<double> b_vector[p_degree] = {0};
    for(int i = 0; i < p_degree; i++)
    {
        proper_lambda[i] = p_lambdaData[i]/p_g; // changes from Lambda_j = gLambda(e_j) to Lambda(e_j)
    }
    for(int i = 0; i < p_degree; i++) // initializes set of equations for polynomial coefficients
    {
        for(int j = 0; j < p_degree; j++)
        {
            a_matrix[i][j] = static_cast<std::complex<double>>(j)*pow(p_energiesData[i],j-1) - proper_lambda[i]*pow(p_energiesData[i],j);
        }
    }
    for(int i = 0; i < p_degree; i++)
    {
        b_vector[i] = proper_lambda[i]*pow(p_energiesData[i],p_degree) - static_cast<std::complex<double>>(p_degree)*pow(p_energiesData[i],p_degree-1);
    }
    cSqMatrix coeff = cSqMatrix(p_degree, &a_matrix[0][0], &b_vector[0]);
    for(int i = p_degree; i >= 0; i--)
    {
        if(i == 0)
        {
            p_coefficients[i] = 1;
        } else
        {
            p_coefficients[i] = b_vector[p_degree-i];
        }
    }
}

std::complex<double> cPolynomial::value(std::complex<double> z_value)
{
    std::complex<double> result = 0;
    std::complex<double> multiplier = 1;
    for(int i = 0; i < p_degree; i++)
    {
        if(abs(p_zeroes[i]) != 0)
        {
            multiplier *= z_value - p_zeroes[i];
        }
    }
    // used to deflate polynomials
    for(int i = p_degree; i >= 0; i--)
    {
        result += (p_coefficients[i]*(pow(z_value,p_degree-i)))/multiplier;
    }
    return result;
}

std::complex<double> cPolynomial::first_der(std::complex<double> z_value)
{
    std::complex<double> result = 0;
    double prec = 0.00001;
    std::complex<double> delta_z = z_value + prec;
    result = (value(delta_z) - value(z_value))/prec;
    return result;
}

std::complex<double> cPolynomial::second_der(std::complex<double> z_value)
{
    std::complex<double> result = 0;
    double prec = 0.00001;
    std::complex<double> delta_z = z_value + prec;
    result = (first_der(delta_z) - first_der(z_value))/prec;
    return result;
}

void cPolynomial::root_finder()
{
    std::complex<double> init_z = 5;
    std::complex<double> a_change;
    std::complex<double> root;
    std::complex<double> G_val, H_val;
    // Numerical Recipes Laguerre Method
    for(int n = 0; n < p_degree; n++)
    {
        do
        {
            a_change = 0;
            G_val = first_der(init_z)/value(init_z);
            H_val = pow(G_val,2) - second_der(init_z)/value(init_z);
            root = sqrt(std::complex<double>((p_degree-1))*(std::complex<double>(p_degree)*H_val - pow(G_val,2)));
            if(abs(G_val + root) > abs(G_val - root))
            {
                a_change = (std::complex<double>(p_degree))/(G_val+root);
            }
            else
            {
                a_change = (std::complex<double>(p_degree))/(G_val-root);
            }
            init_z -= a_change;
        } while(abs(a_change) > 0.00001);
        p_zeroes[n] = init_z; // saves the result that will later be used to deflate polynomials
        init_z=5;
    }
}

std::complex<double> *cPolynomial::zeroes_index()
{
    return &p_zeroes[0];
}
