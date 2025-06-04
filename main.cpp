#include <iostream>
#include <complex>
#include <iomanip>
#include <fstream>
#include "solver.h"

int n_energies, n_pairs, calc_type;

void init();
void energies_init(std::complex<double> *lambda, std::complex<double> *energy);

int main()
{
    std::cout << std::setprecision(6) << std::fixed;
    std::ofstream energiesRe, energiesIm, lambdas, energiesReCSV, energiesImCSV, lambdasCSV;
    energiesRe.open("energiesRe.txt");
    energiesIm.open("energiesIm.txt");
    lambdas.open("Lambdas.txt");
    energiesReCSV.open("energiesRe.csv");
    energiesImCSV.open("energiesIm.csv");
    lambdasCSV.open("lambdas.csv");
    double g = 0.005;
    init();
    std::complex<double> energies[n_energies] = {0};
    std::complex<double> lambda[n_energies] = {0};
    std::complex<double> pair_energies[n_pairs] = {0};
    energies_init(lambda, energies);
    energiesReCSV << "g";
    energiesImCSV << "g";
    lambdasCSV << "g";
    for(int i = 0; i < n_pairs; i++)
    {
        energiesReCSV << ",E" << i;
        energiesImCSV << ",E" << i;
    }
    energiesReCSV << '\n';
    energiesImCSV << '\n';
    for(int i = 0; i < n_energies; i++)
    {
        lambdasCSV << ",L" << i;
    }
    lambdasCSV << '\n';
    do
    {
        std::cout << "g = " << g << '\n';
        cEquations richardson = cEquations(n_energies, lambda, energies, g);
        cPolynomial polynomial = cPolynomial(n_pairs, lambda, energies, g);
        for(int i = 0; i < n_pairs; i++)
        {
            pair_energies[i] = *(polynomial.zeroes_index()+i);
        }
        energiesRe << g << '\t';
        energiesIm << g << '\t';
        lambdas << g << '\t';
        energiesReCSV << g;
        energiesImCSV << g;
        lambdasCSV << g;
        for(int i = 0; i < n_pairs; i++)
        {
            energiesRe << real(pair_energies[i]) << '\t';
            energiesReCSV << "," <<  real(pair_energies[i]);
            energiesIm << imag(pair_energies[i]) << '\t';
            energiesImCSV << "," << imag(pair_energies[i]);
        }
        for(int i = 0; i < n_energies; i++)
        {
            lambdas << lambda[i] << '\t';
            lambdasCSV << real(lambda[i]) << ",";
        }
        energiesRe << '\n';
        energiesIm << '\n';
        lambdas << '\n';
        energiesReCSV << '\n';
        energiesImCSV << '\n';
        lambdas << '\n';
        g+= 0.005;
    } while(g < 1.5);
    energiesIm.close();
    energiesRe.close();
    lambdas.close();
    energiesImCSV.close();
    energiesReCSV.close();
    lambdasCSV.close();
    std::cout << "FINISHED" << '\n';

    return 0;
}

void init()
{
    std::ifstream initFile;
    initFile.open("init.txt");
    initFile >> n_energies >> n_pairs >> calc_type;
    initFile.close();
}

void energies_init(std::complex<double> *lambda, std::complex<double> *energy)
{
    std::ifstream energiesFile;
    energiesFile.open("energies.txt");
    int temp;
    if(calc_type == 0)
    {
        int first_level, second_level;
        energiesFile >> first_level >> second_level;
        temp = 0;
        double first_im = -static_cast<double>(first_level)/2+0.5;
        for(double j = first_im; j <= -first_im; j++)
        {
            *(energy + temp) = std::complex<double>(-1, j);
            temp++;
        }
        double second_im = -static_cast<double>(second_level)/2+0.5;
        for(double j = second_im; j <= -second_im; j++)
        {
            *(energy+temp) = std::complex<double>(1,j);
            temp++;
        }
    } else if(calc_type == 1)
    {
        for(int i = 0; i < n_energies; i++)
        {
            *(energy+i) = std::complex<double>(2*i+1, 0);
        }
        for(int i = 0; i < n_energies; i++)
        {
            std::cout << *(energy+i) << '\t';
        }
        std::cout << '\n';
    } else if(calc_type == 2)
    {
        for(int i = 0; i < n_energies; i++)
        {
            *(energy+i) = std::complex<double>(0.2*(pow(i+1,2)), 0);
        }
        for(int i = 0; i < n_energies; i++)
        {
            std::cout << *(energy+i) << '\t';
        }
        std::cout << '\n';
    } else if(calc_type == 3)
    {
        for(int i = 0; i < n_energies; i++)
        {
            *(energy+i) = std::complex<double>(-5/(pow(i+1,2)), 0);
        }
        for(int i = 0; i < n_energies; i++)
        {
            std::cout << *(energy+i) << '\t';
        }
        std::cout << '\n';
    } else if(calc_type == 4)
    {
        int first_level, second_level, third_level, fourth_level, fifth_level;
        energiesFile >> first_level >> second_level >> third_level >> fourth_level >> fifth_level;
        temp = 0;
        double first_im = -static_cast<double>(first_level)/2+0.5;
        for(double j = first_im; j <= -first_im; j++)
        {
            *(energy + temp) = std::complex<double>(2, j);
            temp++;
        }
        double second_im = -static_cast<double>(second_level)/2+0.5;
        for(double j = second_im; j <= -second_im; j++)
        {
            *(energy+temp) = std::complex<double>(5,j);
            temp++;
        }
        double third_im = -static_cast<double>(third_level)/2+0.5;
        for(double j = third_im; j <= -third_im; j++)
        {
            *(energy + temp) = std::complex<double>(8, j);
            temp++;
        }
        double fourth_im = -static_cast<double>(fourth_level)/2+0.5;
        for(double j = fourth_im; j <= -fourth_im; j++)
        {
            *(energy + temp) = std::complex<double>(10, j);
            temp++;
        }
        double fifth_im = -static_cast<double>(fifth_level)/2+0.5;
        for(double j = fifth_im; j <= -fifth_im; j++)
        {
            *(energy + temp) = std::complex<double>(13, j);
            temp++;
        }
        for(int i = 0; i < n_energies; i++)
        {
            std::cout << *(energy+i) << '\t';
        }
        std::cout << '\n';
    }else if(calc_type == 5)
    {
        int first_level, second_level, third_level, fourth_level, fifth_level;
        energiesFile >> first_level >> second_level >> third_level >> fourth_level >> fifth_level;
        temp = 0;
        double first_im = -static_cast<double>(first_level)/2+0.5;
        for(double j = first_im; j <= -first_im; j++)
        {
            *(energy + temp) = std::complex<double>(2, j);
            temp++;
        }
        double second_im = -static_cast<double>(second_level)/2+0.5;
        for(double j = second_im; j <= -second_im; j++)
        {
            *(energy+temp) = std::complex<double>(4,j);
            temp++;
        }
        double third_im = -static_cast<double>(third_level)/2+0.5;
        for(double j = third_im; j <= -third_im; j++)
        {
            *(energy + temp) = std::complex<double>(6, j);
            temp++;
        }
        double fourth_im = -static_cast<double>(fourth_level)/2+0.5;
        for(double j = fourth_im; j <= -fourth_im; j++)
        {
            *(energy + temp) = std::complex<double>(8, j);
            temp++;
        }
        double fifth_im = -static_cast<double>(fifth_level)/2+0.5;
        for(double j = fifth_im; j <= -fifth_im; j++)
        {
            *(energy + temp) = std::complex<double>(10, j);
            temp++;
        }
        for(int i = 0; i < n_energies; i++)
        {
            std::cout << *(energy+i) << '\t';
        }
        std::cout << '\n';
    }
    for(int i = 0; i < n_energies; i++)
    {
        energiesFile >> *(lambda + i);
    }
    energiesFile.close();
}

