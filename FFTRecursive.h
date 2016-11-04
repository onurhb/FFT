//
// Created by Onur on 02.10.2016.
//

#ifndef FFT_FFTRECURSIVE_H
#define FFT_FFTRECURSIVE_H


#include <complex>

#define PI 3.141592653589793238462
class FFTRecursive {
public:
    std::complex<double>* doFFT(std::complex<double>* samples, unsigned int size);
};


#endif //FFT_FFTRECURSIVE_H
