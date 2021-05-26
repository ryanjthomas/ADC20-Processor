#ifndef FFT_H
#define FFT_H

#include <complex>
#include <iostream>
#include <cmath>
#include <valarray>
 
const double PI = 3.141592653589793238460;
 
typedef std::complex<long double> Complex;
typedef std::valarray<Complex> CArray;


void fft(CArray &x);
void ifft(CArray &x);

#endif
