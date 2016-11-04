/*
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

// ------------ Constants
#define PI 3.14159265359
// -----------

namespace Algorithm {

    enum WINDOW {
        NO_WINDOW,
        HANNING_WINDOW,
        HARRIS_WINDOW
    };

    class FFT {
        // - If 1D FFT/IFFT is used, only left will be set
        // - If 2D FFT/IFFT is used, both will be set
        struct CONTAINER {
            std::complex<double> *left;     // - Actual frequency domain for first channel
            std::complex<double> *right;    // - Actual frequency domain for secound channel
        } container;

        // - N specifies size of FFT, and windowsFunction is currently used windowing
        int N;
        WINDOW windowFunction;

    public:
        /*
         * We should tell FFT size and what kind of windowing we want
         */
        FFT(int N, WINDOW windowFunction);

        ~FFT();

        /*
         * FFT should be initialized with number of samples and a window function (use macros)
         */
        void initialize(int N, WINDOW windowFunction);

        /*
         * Forward FFT will transform samples to frequency domain
         */
        std::complex<double> *forwardFFT(std::complex<double> *samples);

        /*
         * Forward 2 dimensional FFT, use this if you want to execute FFT on two sample arrays parallell
         * Recommended if you have two channels that you want to analyze
         */
        CONTAINER forwardFFT_2D(std::complex<double> *left, std::complex<double> *right);

        /*
         * Inverse FFT will transform frequency domain back to samples
         */
        std::complex<double> *inverseFFT(std::complex<double> *input);

        /*
         * Inverse 2 dimensional FFT, use this if you want to execute IFFT on two sample arrays parallell
         * Recommended if you have two channels that you want to convert back to raw samples
         */
        CONTAINER inverseFFT_2D(std::complex<double> *left, std::complex<double> *right);

    private:

        /*
        * @param N : N in twiddle factor, @param n : n in twiddle factor , @return twiddle factor
        */
        std::complex<double> getTwiddle(int N, int n) const;

        /*
        * @param complexSamples : Samples to reverse order by reversed bit index
        * @param N : size of complexSamples
        * @return array of std::complex<double> ordered by reverse bit index
        */
        void bizarreOrder(std::complex<double> *container, std::complex<double> *samples) const;

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

        /*
         * Windowing functions
         */
        std::complex<double> *harrisWindow(std::complex<double> *samples) const;

        std::complex<double> *hanningWindow(std::complex<double> *samples) const;
    };
}