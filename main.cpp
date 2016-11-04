#include <iostream>
#include <complex>
#include <vector>
#include "FFT.h"
#include <ctime>
#include <cstdio>


// -printf("%f\n", std::abs(freqbins[j] * 2.0) * 1.0 / N);
int main() {

    // - 1Hz, amp = 1
    int N = 8;
    int sampleRate = 8; // - Samples each sec, needed when we are going to print out spectrum
    double rawSamples[N] = {0.0, 0.707, 1.0, 0.707, 0.0, -0.707, -1.0, -0.707};

    /*
    // - 1Hz, amp = 1
    unsigned int N = 4;
    double rawSamples[N] = {0.0, 1.0, 0.0, -1.0};
    */

    std::complex<double>* samples = new std::complex<double>[N];

    for (int i = 0; i < N; ++i) {
        samples[i] = rawSamples[i];
    }

    Algorithm::FFT fft(N, Algorithm::NO_WINDOW);

    clock_t begin = clock();
    std::complex<double>* frequencyBins = fft.forwardFFT(samples);

    // - Printing out frequency bins!
    for (int i = 0; i < N; ++i) {
        std::cout <<  abs(frequencyBins[i]) * 2 / N << std::endl;
    }

    printf("------\n");

    // - Getting the samples back from the frequency bins!
    samples = fft.inverseFFT(frequencyBins);
    for (int i = 0; i < N; ++i) {
        std::cout <<  samples[i].real() << std::endl;
    }


    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Time: %f", elapsed_secs);
    return 0;
}