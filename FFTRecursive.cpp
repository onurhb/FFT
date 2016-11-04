//
// Created by Onur on 02.10.2016.
//

#include "FFTRecursive.h"


/*
* Coution! You need to understand how FFT works!
* For that you have to understand DFT, and why FFt is faster.
* If you understand DFT you'll notice FFT is just DFT, but uses symetry and shortcuts to be faster.
* FFT O(N LogN) : DFT O(N^2)
* */
std::complex<double>* FFTRecursive::doFFT(std::complex<double>* samples, unsigned int size){

    unsigned int N = size; // - Size of samples

    if(N==1) return samples;   // - Stop recursion if sample size is one

    int M = N / 2;

    std::complex<double> * Xeven = new std::complex<double>[M]; // - Half of evens
    std::complex<double> * Xodd = new std::complex<double>[M];  // - Half of odds

    for (int i = 0; i < M; ++i) {
        Xeven[i] = samples[2 * i];       // - Push even samples
        Xodd[i] = samples[2 * i + 1];    // - Push odd samples
    }

    std::complex<double> * Feven; // - Half of transformed evens
    std::complex<double> * Fodd;  // - Half of transformed odds

    Feven = doFFT(Xeven, M);                     // - We should end up with all samples
    Fodd = doFFT(Xodd, M);                       // - We should end up with all samples

    std::complex<double> * Fbinds = new std::complex<double>[N]; // - We have N freq bins

    // Calculate frequency bins using symetry
    for (int k = 0; k != N / 2; k ++) {
        std::complex<double> w = std::polar(1.0, -2.0 * PI * k/N) * Fodd[k];
        Fbinds[k] = Feven[k] + w;
        Fbinds[k + N/ 2] = Feven[k] - w;    // - Using symetry
    }
    return  Fbinds;

}