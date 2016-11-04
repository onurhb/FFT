
// ---------------- INCLUDES
#include <cstdio>
#include <iostream>
#include "FFT.h"
// ----------------

namespace Algorithm {

    FFT::FFT(int N, WINDOW windowFunction) {
        initialize(N, windowFunction);
    }

    FFT::~FFT() {
        delete container.left;
        delete container.right;
    }

    /*
     * This should be called if you decide to change N or windowing function
     */
    void FFT::initialize(int N, WINDOW windowFunction) {
        this->N = N;
        this->windowFunction = windowFunction;
        if (!validSize()) {
            printf("%s", "size of Samples needs to be power of 2.\n");
            return;
        };

        // - Initialize both left and right container
        container.left = new std::complex<double>[N];
        container.right = new std::complex<double>[N];
    }

    /*
     * This method will convert samples to time domain frequency container
     * Reference            : https://en.wikipedia.org/wiki/Butterfly_diagram
     * @param samples       : raw samples
     * @return              : time domain frequency container
     */
    std::complex<double> *FFT::forwardFFT(std::complex<double> *samples) {

        std::complex<double> *frquencyBins = container.left;
        // - Reversing and using windowing functions, pre-processed before butterfly algorithm is executed
        if (this->windowFunction == HANNING_WINDOW)
            bizarreOrder(frquencyBins, hanningWindow(samples));    // - Reversing hanned samples
        else if (this->windowFunction == HARRIS_WINDOW)
            bizarreOrder(frquencyBins, harrisWindow(samples));
        else
            bizarreOrder(frquencyBins, samples);

        // - Initialy, the butterfly size if 1 and next butterfly is at offset 2
        int b_size = 1;
        int offset = 2;

        // - Butterfly algorithm, N = 8 : 3 stages
        for (int stage = 0; stage < getStages(); ++stage) {
            for (int k = 0; k < N; k++) {
                // - b_size = 0 & k gives even numbers, b_size = 2 includes 2 and skips 2
                if (!(b_size & k)) {
                    // - Calculate twiddle factor, using symetry and relations
                    std::complex<double> W = getTwiddle(offset, k % offset);

                    // - Temporary store upper
                    std::complex<double> temp = frquencyBins[k + b_size] * W;

                    // - Bottom is = upper - bottom * W
                    frquencyBins[k + b_size] = frquencyBins[k] - temp;
                    // - Upper is = upper + bottom * W
                    frquencyBins[k] += temp;

                }
            }

            b_size *= 2;
            offset *= 2;
        }

        return frquencyBins;
    }


    /*
     * This method will convert samples to time domain frequency container for two channels parallell
     * Reference            : https://en.wikipedia.org/wiki/Butterfly_diagram
     * @param left          : raw samples of first sample array
     * @param right          : raw samples of secound sample array
     * @return              : time domain frequency container for left and right
     */
    FFT::CONTAINER FFT::forwardFFT_2D(std::complex<double> *left, std::complex<double> *right) {

        // - Pointer to left and right container that we want to fill
        std::complex<double> *leftSpectrum = container.left;
        std::complex<double> *rightSpectrum = container.right;

        // - Reversing and using windowing functions, pre-processed before butterfly algorithm is executed
        if (this->windowFunction == HANNING_WINDOW) {
            bizarreOrder(leftSpectrum, hanningWindow(left));
            bizarreOrder(rightSpectrum, hanningWindow(right));
        } else if (this->windowFunction == HARRIS_WINDOW) {
            bizarreOrder(leftSpectrum, harrisWindow(left));
            bizarreOrder(rightSpectrum, harrisWindow(right));
        } else {
            bizarreOrder(leftSpectrum, left);
            bizarreOrder(rightSpectrum, right);
        }

        // - Initialy, the butterfly size is 1 and next butterfly is at offset 2
        int b_size = 1;
        int offset = 2;

        // - Butterfly algorithm, N = 8 : 3 stages
        for (int stage = 0; stage < getStages(); ++stage) {
            for (int k = 0; k < N; k++) {
                // - b_size = 0 & k gives even numbers, b_size = 2 includes 2 and skips 2
                if (!(b_size & k)) {
                    // - Calculate twiddle factor, using symetry and relations
                    std::complex<double> W = getTwiddle(offset, k % offset);

                    // - Temporary store upper
                    std::complex<double> tempLeft = leftSpectrum[k + b_size] * W;
                    std::complex<double> tempRight = rightSpectrum[k + b_size] * W;

                    // - Bottom is = upper - bottom * W
                    leftSpectrum[k + b_size] = leftSpectrum[k] - tempLeft;
                    rightSpectrum[k + b_size] = rightSpectrum[k] - tempRight;
                    // - Upper is = upper + bottom * W
                    leftSpectrum[k] += tempLeft;
                    rightSpectrum[k] += tempRight;

                }
            }

            b_size *= 2;
            offset *= 2;
        }
        return container;
    }

    /*
     * This method will transform time domain frequency bins to ectual samples
     * This technique uses the FFT algorithm to calculate the inverse FFT
     * Reference            : http://www.originlab.com/doc/Origin-Help/IFFT
     * @param samples       : time domain frequency to inverse
     * @return              : raw samples
     */
    std::complex<double> *FFT::inverseFFT(std::complex<double> *input) {
        // - First we find conjugate of each frequency, then we store in left
        for (int i = 0; i < N; ++i) {
            //container.left[i] = std::conj(input[i]); // <--- does not work? (because container.left = input)
            container.right[i] = std::conj(input[i]);  // <--- works
        }

        // - Using FFT on the conjugated values
        container.right = forwardFFT(container.right);

        // Applying conjugate again to reverse back
        for (int i = 0; i < N; ++i) {
            container.right[i] = std::conj(container.right[i]);
            container.right[i] /= N;
        }

        return container.right;
    }

    /*
     * This method will transform time domain frequency bins to ectual samples
     * This technique uses the FFT algorithm to calculate the inverse FFT
     * Reference            : http://www.originlab.com/doc/Origin-Help/IFFT
     * @param samples       : time domain frequency to inverse
     * @return              : raw samples
     */
    FFT::CONTAINER FFT::inverseFFT_2D(std::complex<double> *left, std::complex<double> *right) {
        // - First we find conjugate of each frequency, then we store in left
        for (int i = 0; i < N; ++i) {
            container.right[i] = std::conj(left[i]);
            container.left[i] = std::conj(right[i]);
        }

        // - Using FFT on the conjugated values
        forwardFFT_2D(container.right, container.left);

        // Applying conjugate again to reverse back
        for (int i = 0; i < N; ++i) {
            container.right[i] = std::conj(container.right[i]);
            container.right[i] /= N;
            container.left[i] = std::conj(container.left[i]);
            container.left[i] /= N;
        }

        return container;
    }

    /*
     * This method will execute "windowing function" on the samples
     * Reference        : http://www.ni.com/white-paper/4844/en/
     * @param samples   : samples to manipulate
     * @return          : manipulated samples
     */
    std::complex<double> *FFT::harrisWindow(std::complex<double> *samples) const {
        double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
        std::complex<double> temp;
        for (int i = 0; i < N; i++) {
            temp = 2 * PI * samples[i] / double(N);
            samples[i] *= (a0 - (a1 * cos(temp)) + (a2 * cos(2.0 * temp)) - (a3 * cos(3.0 * temp)));
        }

        return samples;
    }

    /*
     * This method will execute "windowing function" on the samples
     * Reference            : http://www.ni.com/white-paper/4844/en/
     * @param samples       : samples to manipulate
     * @return              : manipulated samples
     */
    std::complex<double> *FFT::hanningWindow(std::complex<double> *samples) const {
        for (int i = 0; i < N; i++) {
            double multiplier = 0.5 * (1 - cos(2 * PI * i / N));
            samples[i] *= multiplier;
        }
        return samples;
    }

    /*
     * Re-orders an array of samples by inversing indexes
     * Reference            : http://www.dspguide.com/ch12/2.htm
     * @param container     : reversed array will be placed into this array
     * @param samples       : array of samples to process
     * @return              : pointer to target
     */
    void FFT::bizarreOrder(std::complex<double> *container, std::complex<double> *samples) const {
        // - Target should be size of atleast N!
        for (int i = 0; i < N; i++) {
            // - Index to swap
            int indexSwap = reverseBit(i);
            container[i] = samples[indexSwap];
        }
    }

    /*
     * Calculates the twiddlefactor, e^(-i2PIn)/N
     * Reference            : https://en.wikipedia.org/wiki/Twiddle_factor
     * @param N             : number of samples, same N as in e^(-i2PIn)/N
     * @param n             : sample index, same n as in e^(-i2PIn)/N
     */
    std::complex<double> FFT::getTwiddle(int N, int n) const {
        return std::polar(1.0, -2.0 * PI * n / N);
    }


    /*
     * Returns number of required stages in order to perform FFT, same as log(N)/log(2)
     * FFT should be initilized with N, if N is 8, number of stages will be 3.
     * @return          : number of stages
     */
    int FFT::getStages() const {
        int p = 1, i = 1;;
        // - Same as taking power of two
        while ((p <<= 1) < N) i++;
        return i;
    }

    /*
     * Checks if N is power of two
     * @return          : true if power of two
     */
    bool FFT::validSize() const {
        // - We can decide if a number is power of two by AND'ing with a odd number
        return !(N & (N - 1));      // - Should give 1 if power of two
    }

    /*
     * Reverse a number by reversing bit
     * Handles out of index range by using number of stages
     */
    int FFT::reverseBit(int bit) const {
        int stages = getStages();
        int r = 0;
        // - Making sure that x-amount of bits are included when we perform bit-reversing
        for (int i = 0; i < stages; ++i) {
            int b = bit & 1;
            b <<= (stages - 1) - i;
            r |= b;
            bit >>= 1;
        }

        return r;
    }

}
