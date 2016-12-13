/**
 * USAGE:
 * - Initialize FFT with size of N that you're planning to use. If you have 1024 samples that you want to analyze
 * N will be 1024. Specify windowing function that you want FFT to use, currently supporting windows specified in WINDOW enum
 * - Example: FFT fft(N, NO_WINDOW);
 * - FFT is now ready to rock. All you need to do is executing forward, inverse or two dimensional FFT.
 * - Example:
 * Let's say left has raw samples converted to complex numbers
 *      std::complex<double>* left  = new std::complex<double>[N];
 * Executiong forward FFT on left
 *      std::complex<double>* spectrum = fft.forwardFFT(left);
 * Variable spectrum consists of frequency spectrum
 * Verify with
 *  for (int j = 0; j < N; ++j) {
 *       std::cout << spectrum[j] << std::endl;
 *   }
 *
 *   - If you decide to change size of FFT or you want different windowing, simply initialize FFT with initialize()
 */

// ------------ Includes
#include <complex>
#include <vector>
// ------------ Constants
#define PI 3.14159265359
// -----------

namespace Algorithm {

    // - Write less, save time
    typedef std::vector<std::complex<float>> cpx_vector;

    class FFT {
        // - N specifies size of FFT, and windowsFunction is currently used windowing
        int N;

        // - If 1D FFT/IFFT is used, only left will be used
        // - If 2D FFT/IFFT is used, both will be used
        struct CONTAINER {
            cpx_vector left;     // - Actual frequency domain for first channel
            cpx_vector right;    // - Actual frequency domain for secound channel
        } container;


    public:
        /*
         * We should tell FFT size and what kind of windowing we want
         */
        FFT(int N);

        /*
         * Forward FFT will transform samples to frequency domain
         */
        cpx_vector forwardFFT(cpx_vector& samples);

        /*
         * Forward 2 dimensional FFT, use this if you want to execute FFT on two sample arrays parallell
         * Recommended if you have two channels that you want to analyze
         */
        CONTAINER& forwardFFT_2D(cpx_vector& left, cpx_vector& right);

        /*
         * Inverse FFT will transform frequency domain back to samples
         */
        cpx_vector inverseFFT(cpx_vector input);

        /*
         * Inverse 2 dimensional FFT, use this if you want to execute IFFT on two sample arrays parallell
         * Recommended if you have two channels that you want to convert back to raw samples
         */
        CONTAINER& inverseFFT_2D(cpx_vector& left, cpx_vector& right);

        /*
         * Windowing functions
         */
        cpx_vector& hanningWindow(cpx_vector& samples) const;

        cpx_vector &hammingWindow(cpx_vector &samples) const;

    private:

        /*
        * @param N : N in twiddle factor, @param n : n in twiddle factor , @return twiddle factor
        */
        std::complex<float> getTwiddle(int N, int n) const;

        /*
        * @param complexSamples : samples to reverse order by reversed bit index
        * @param N : size of complexSamples
        * @return array of std::complex<double> ordered by reverse bit index
        */
        cpx_vector& bizarreOrder(cpx_vector& fbins) ;

        /*
         * Reverses a number (it's bits)
         */
        int reverseBit(int b) const;

        /*
        * @param N : size of samples, @return total steps to execute in order to calculate FFT
        */
        int getStages() const;

        /*
        * FFT works best with power of two size, zero padding might be used to obtain need if it is not the case
        * @return true/false if size of N is power of 2
        */
        bool validSize() const;
    };
}