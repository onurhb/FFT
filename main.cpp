#include <complex>
#include "FFT.h"
#include <ctime>
#include <iostream>
#include <kiss_fft.h>



int main() {

    clock_t begin, end;
    double elapsed_secs;
    using namespace Algorithm;

    int N = 8;
    std::vector<std::complex<float>> spectrum;

    int sampleSize = 2048;
    int rounds = 10000000;
    std::cout << "running " << sampleSize << " samples" << std::endl;

    for (int j = 0; j < sampleSize; ++j) {
        spectrum.push_back(10.0f);
    }


    Algorithm::FFT fft(N);
    std::vector<std::complex<float>> frequencyBins;
    begin = clock();
    for (int i = 0; i < rounds; ++i) {
        frequencyBins = fft.forwardFFT(spectrum);
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Time FFT: %f\n", elapsed_secs);



    kiss_fft_cfg fwd = kiss_fft_alloc(N,0,NULL,NULL);
    begin = clock();
    for (int k = 0; k < rounds; ++k) {
        kiss_fft(fwd,(kiss_fft_cpx*) &spectrum[0],(kiss_fft_cpx*) &frequencyBins[0]);
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Time KissFFT: %f\n", elapsed_secs);


    kiss_fft_free(fwd);
    return 0;
}