#ifndef SOLVER_H_
#define SOLVER_H_
#include <complex>

class cSqMatrix
{
private:
    int m_dims;
    int m_nElements;
    std::complex<double> *m_matrixData;
    std::complex<double> *m_vectorData;
public:
    cSqMatrix(int dims, std::complex<double> *matrixData, std::complex<double> *vecData);
    ~cSqMatrix();
    void LU_decomposition(std::complex<double> *pMatrix, std::complex<double> *pVector, int dim);
    void swap_two(std::complex<double> *pt1, std::complex<double> *pt2);
};

class cEquations
{
private:
    int eq_dims;
    double eq_g;
    std::complex<double> *eq_lambdaData;
    std::complex<double> *eq_energiesData;

public:
    cEquations(int dims, std::complex<double> *lambdas, std::complex<double> *energies, double g);
    ~cEquations();
    std::complex<double> function_i(std::complex<double> *x_vector, int i_index);
    void newton_raphson(std::complex<double> *init_guess);

};

class cPolynomial
{
private:
    double p_g;
    int p_degree;
    std::complex<double> *p_lambdaData;
    std::complex<double> *p_energiesData;
    std::complex<double> *p_coefficients;
    std::complex<double> *p_zeroes;
public:
    cPolynomial(int degree, std::complex<double> *lambdas, std::complex<double> *energies, double g);
    ~cPolynomial();
    void coefficients();
    std::complex<double> value(std::complex<double> z_value);
    std::complex<double> first_der(std::complex<double> z_value);
    std::complex<double> second_der(std::complex<double> z_value);
    void root_finder();
    std::complex<double> *zeroes_index();
};
#endif // SOLVER_H_
