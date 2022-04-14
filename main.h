#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

float* SpectrumAnalyst(int alg, int windowFunc, int windowSize, double rate, const float *data, int dataLen);
float* NewWindowFunc(int whichFunction, int NumSamplesIn, int extraSample, float *in);
float* WindowFunc(int whichFunction, int NumSamples, float *in);
float* PowerSpectrum(int NumSamples, const float *In, float *Out);

enum eWindowFunctions
{
   eWinFuncRectangular,
   eWinFuncBartlett,
   eWinFuncHamming,
   eWinFuncHann,
   eWinFuncBlackman,
   eWinFuncBlackmanHarris,
   eWinFuncWelch,
   eWinFuncGaussian25,
   eWinFuncGaussian35,
   eWinFuncGaussian45,
   eWinFuncCount
};

enum Algorithm
{
  Spectrum,
  Autocorrelation,
  CubeRootAutocorrelation,
  EnhancedAutocorrelation,
  Cepstrum,

  NumAlgorithms
};

typedef struct {
    int Points;
    int *BitReversed;
    float *SinTable;
} FFTParam;

typedef struct {
    float *RealPart;
    float *ImagPart;
} FFTOutput;

typedef struct {
    int length;
    float *A;
}ArrayInfo;

float* RealFFTf(float *buffer, const FFTParam *h);
FFTParam InitializeFFT(int fftlen);
FFTOutput RealFFT(int NumSamples, const float *RealIn, float *RealOut, float *ImagOut);
void InverseRealFFTf(float *buffer, const FFTParam *h);
float * InverseRealFFT(int NumSamples, const float *RealIn, const float *ImagIn, float *RealOut, float * bffInit);
float* ReorderToTime(const FFTParam *hFFT, const float *buffer, float *TimeOut);

#endif // MAIN_H_INCLUDED
