
// ---------------- INCLUDES
#include "FFT.h"
#include <iostream>
// ----------------

namespace Algorithm {

    FFT::FFT(int N) {
        this->N = N;
        if (!validSize()) {
            printf("%s", "size of samples needs to be power of 2.\n");
            return;
        };
        // - Initialize vectors
        this->container.left = cpx_vector(N, 0.0f);
        this->container.right = cpx_vector(N, 0.0f);
    }


    /*
     * This method will convert samples to time domain frequency container
     * Reference            : https://en.wikipedia.org/wiki/Butterfly_diagram
     * @param samples       : raw samples
     * @return              : time domain frequency container
     */
    cpx_vector FFT::forwardFFT(cpx_vector &samples) {

#ifdef HANNING_WINDOW
        container.left = bizarreOrder(hanningWindow(samples));    // - Reversing hanned samples
#endif
#ifndef HANNING_WINDOW
        container.left = bizarreOrder(samples);
#endif

        // - Initialy, the butterfly size if 1 and next butterfly is at offset 2
        int b_size = 1;
        int offset = 2;

        // - Butterfly algorithm, N = 8 : 3 stages
        for (int stage = 0; stage < getStages(); ++stage) {
            for (int k = 0; k < N; k++) {
                // - b_size = 0 & k gives even numbers, b_size = 2 includes 2 and skips 2
                if (!(b_size & k)) {
                    // - Calculate twiddle factor, using symetry and relations
                    std::complex<float> W = getTwiddle(offset, k % offset);

                    // - Temporary store upper
                    std::complex<float> temp = container.left[k + b_size] * W;

                    // - Bottom is = upper - bottom * W
                    container.left[k + b_size] = container.left[k] - temp;
                    // - Upper is = upper + bottom * W
                    container.left[k] += temp;

                }
            }

            b_size *= 2;
            offset *= 2;
        }

        return container.left;
    }


    /*
     * This method will convert samples to time domain frequency container for two channels parallell
     * Reference            : https://en.wikipedia.org/wiki/Butterfly_diagram
     * @param left          : raw samples of first sample array
     * @param right          : raw samples of secound sample array
     * @return              : time domain frequency container for left and right
     */
    FFT::CONTAINER &FFT::forwardFFT_2D(cpx_vector &left, cpx_vector &right) {

        // - Reversing and using windowing functions, pre-processed before butterfly algorithm is executed
#ifdef HANNING_WINDOW
        container.left = bizarreOrder(hanningWindow(left));
            container.right = bizarreOrder(hanningWindow(right));
#endif
#ifndef HANNING_WINDOW
        container.left = bizarreOrder(right);
        container.left = bizarreOrder(left);
#endif

        // - Initialy, the butterfly size is 1 and next butterfly is at offset 2
        int b_size = 1;
        int offset = 2;

        // - Butterfly algorithm, N = 8 : 3 stages
        for (int stage = 0; stage < getStages(); ++stage) {
            for (int k = 0; k < N; k++) {
                // - b_size = 0 & k gives even numbers, b_size = 2 includes 2 and skips 2
                if (!(b_size & k)) {
                    // - Calculate twiddle factor, using symetry and relations
                    std::complex<float> W = getTwiddle(offset, k % offset);

                    // - Temporary store upper
                    std::complex<float> tempLeft = container.left[k + b_size] * W;
                    std::complex<float> tempRight = container.right[k + b_size] * W;

                    // - Bottom is = upper - bottom * W
                    container.left[k + b_size] = container.left[k] - tempLeft;
                    container.right[k + b_size] = container.right[k] - tempRight;
                    // - Upper is = upper + bottom * W
                    container.left[k] += tempLeft;
                    container.right[k] += tempRight;

                }
            }

            b_size *= 2;
            offset *= 2;
        }
        return container;
    }

    /*
     * This method will transform time domain frequency bins to ectual samples
     * Note that the parameter is not a reference, this is so we don't override the original forwarded fft
     * This technique uses the FFT algorithm to calculate the inverse FFT
     * Reference            : http://www.originlab.com/doc/Origin-Help/IFFT
     * @param samples       : time domain frequency to inverse
     * @return              : raw samples
     */
    cpx_vector FFT::inverseFFT(cpx_vector input) {
        // - First we find conjugate of each frequency
        for (int i = 0; i < N; ++i) {
            input[i] = std::conj(input[i]);
        }

        // - Using FFT on the conjugated values
        input = forwardFFT(input);

        // Applying conjugate again to reverse back
        for (int i = 0; i < N; ++i) {
            input[i] = std::conj(input[i]);
            input[i] /= N;
        }

        return input;
    }

    /*
     * This method will transform time domain frequency bins to ectual samples
     * This technique uses the FFT algorithm to calculate the inverse FFT
     * Reference            : http://www.originlab.com/doc/Origin-Help/IFFT
     * @param samples       : time domain frequency to inverse
     * @return              : raw samples
     */
    FFT::CONTAINER &FFT::inverseFFT_2D(cpx_vector &left, cpx_vector &right) {
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
     * Reference            : http://www.ni.com/white-paper/4844/en/
     * @param samples       : samples to manipulate
     * @return              : manipulated samples
     */
    cpx_vector &FFT::hanningWindow(cpx_vector &samples) const {
        for (float i = 0.0f; i < N; i++) {
            samples[i] *=  0.5f * (1.0f - cos(2.0f * PI * i / N));
        }
        return samples;
    }

    cpx_vector &FFT::hammingWindow(cpx_vector &samples) const{
        for (int i = 0; i < N; ++i) {
            samples[i] *= 0.54f - 0.46f * (float) cos((2 * PI * i) / (N - 1));
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
    cpx_vector& FFT::bizarreOrder(cpx_vector &fbins)  {
        // - Target should be size of atleast N!
        for (int i = 0; i < N; i++) {
            // - Index to swap
            int indexSwap = reverseBit(i);
            container.left[i] = fbins[indexSwap];
        }

        return container.left;
    }

    /*
     * Calculates the twiddlefactor, e^(-i2PIn)/N
     * Reference            : https://en.wikipedia.org/wiki/Twiddle_factor
     * @param N             : number of samples, same N as in e^(-i2PIn)/N
     * @param n             : sample index, same n as in e^(-i2PIn)/N
     */
    std::complex<float> FFT::getTwiddle(int N, int n) const {
        return std::polar<float>(1.0f, -2.0f * PI * n / N);
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
