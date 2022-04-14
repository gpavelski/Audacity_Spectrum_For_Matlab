#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"

#define WinSize 512

int main()
{
    char * strncpy();
    char ch;
    int i = 0;
    int j = 0;
    int linesCount=0;
    float *Spec;
    enum Algorithm alg = Spectrum;
    enum eWindowFunctions winType = eWinFuncRectangular;
    int mRate = 44100;

    FILE *fp;
    fp =fopen ("sample.txt", "r");
    char a[21];
    char b[22];

    // Count the number of lines in the file
    while((ch=fgetc(fp))!=EOF) {
        if(ch=='\n')
            linesCount++;
    }
    fclose(fp);
    fp =fopen ("sample.txt", "r");
    float data[linesCount+1]; // Initialize the array

    for (i = 0; i <= linesCount; i++){
        j = 0;
        while ((ch=getc(fp))!='\n'){
            a[j] = ch;
            if (i == linesCount && j > 20) // End of File Reached
                break;
            j++;
        }

        strncpy(b, a, j);
        b[j] = '\0';
        data[i] = atof(b);
    }
    fclose(fp);



    // Here starts the program

    int dataLen;
    dataLen = sizeof data / sizeof data[0];

    Spec = SpectrumAnalyst(alg, winType, WinSize, mRate, data, dataLen);

    // Outputs data to txt file

    FILE *f = fopen("output.txt", "w");
    if (alg < 4)
    {
        fprintf(f,"Frequency (Hz)\tLevel (dB)\n");
        for (i = 1; i < WinSize/2; i++)
        {
            fprintf(f, "%.6f \t %.18f\n", (float)i*mRate/WinSize, Spec[i]);
        }
    }
    else
    {
        fprintf(f,"Lag (seconds)\tFrequency (Hz)\tLevel\n");
        for (i = 1; i < WinSize/2; i++)
        {
            fprintf(f, "%.6f \t %.6f \t %.18f\n", (float)i/mRate, (float)mRate/i, Spec[i]);
        }
    }
    fclose(f);


    return 0;
}

float* SpectrumAnalyst(int alg, int windowFunc, int windowSize, double rate, const float *data, int dataLen)
{
    double mRate = rate;
    int mWindowSize = windowSize;
    int half = mWindowSize / 2;
    static float mProcessed[WinSize];
    float in[mWindowSize];
    float out[mWindowSize];
    float out2[mWindowSize];
    float win[mWindowSize];
    float ipFFTinit[mWindowSize];
    int i;
    float *pwin;
    float *pout;
    FFTOutput r;

    for (i = 0; i < mWindowSize; i++)
    {
        mProcessed[i] = 0.0f;
        win[i] = 1.0f;
        ipFFTinit[i] = 0.0f;
    }

    pwin = WindowFunc(windowFunc, mWindowSize, win);

    for (i = 0; i < mWindowSize; i++)
    {
        win[i] = pwin[i];
    }

    // Scale window such that an amplitude of 1.0 in the time domain
   // shows an amplitude of 0dB in the frequency domain
    double wss = 0;
    for (i = 0; i<mWindowSize; i++){
        wss += win[i];
    }
    if(wss > 0)
        wss = 4.0 / (wss*wss);
    else
        wss = 1.0;

    int start = 0;
    int windows = 0;

    while (start + mWindowSize <= dataLen)
    {
        for (i = 0; i < mWindowSize; i++) {
            in[i] = win[i] * data[start + i];
        }

        if (alg == Spectrum)
        {
            pout = PowerSpectrum(mWindowSize, in, out);

            for (i = 0; i < mWindowSize; i++)
            {
                out[i] = pout[i];

            }
            for (i = 0; i < half; i++)
            {
                mProcessed[i] += out[i];
            }
        }

        else if (alg == Autocorrelation || alg == CubeRootAutocorrelation || alg == EnhancedAutocorrelation)
        {
            r = RealFFT(mWindowSize, in, out, out2);
            for (i = 0; i < mWindowSize; i++)
            {
                out[i] = r.RealPart[i];
                out2[i] = r.ImagPart[i];
                in[i] = (out[i] * out[i]) + (out2[i] * out2[i]);
            }
            if (alg == Autocorrelation)
            {
               for (i = 0; i < mWindowSize; i++)
               {
                   in[i] = sqrt(in[i]);
               }
            }

            if (alg == CubeRootAutocorrelation ||
                alg == EnhancedAutocorrelation)
            {

                // Tolonen and Karjalainen recommend taking the cube root
               // of the power, instead of the square root
               for (i = 0; i < mWindowSize; i++)
               {
                   in[i] = pow(in[i], 1.0f / 3.0f);
               }
            }
           // Take FFT
            r = RealFFT(mWindowSize, in, out, out2);
            for (i = 0; i < mWindowSize; i++)
            {
                out[i] = r.RealPart[i];
                out2[i] = r.ImagPart[i];

            }
            // Take real part of result
            for (i = 0; i < half; i++)
            {
                mProcessed[i] += out[i];
            }
        }
        else if (alg == Cepstrum)
        {
            r = RealFFT(mWindowSize, in, out, out2);
            for (i = 0; i < mWindowSize; i++)
            {
                out[i] = r.RealPart[i];
                out2[i] = r.ImagPart[i];
            }

            float power;
            float minpower = 1e-20*mWindowSize*mWindowSize;
            for (i = 0; i < mWindowSize; i++)
            {
                power = (out[i] * out[i]) + (out2[i] * out2[i]);
                if(power < minpower)
                    in[i] = log(minpower);
                else
                    in[i] = log(power);
            }
            pout = InverseRealFFT(mWindowSize, in, NULL, out, ipFFTinit);

            for (i = 0; i < half; i++)
            {
                mProcessed[i] += pout[i];
            }
        }
        start += half;
        windows++;
    }

    float mYMin = 1000000, mYMax = -1000000;
    double scale;
    if (alg == Spectrum)
    {
        // Convert to decibels
        mYMin = 1000000.;
        mYMax = -1000000.;
        scale = wss / (double)windows;

        for (i = 0; i < half; i++)
        {
            mProcessed[i] = 10 * log10(mProcessed[i] * scale);
            if(mProcessed[i] > mYMax)
            {
                mYMax = mProcessed[i];
            }
            else if(mProcessed[i] < mYMin)
            {
                mYMin = mProcessed[i];
            }
        }
    }
    else if (alg == Autocorrelation || alg == CubeRootAutocorrelation)
    {
        for (i = 0; i < half; i++)
        {
            mProcessed[i] = mProcessed[i] / windows;
        }
        // Find min/max
        mYMin = mProcessed[0];
        mYMax = mProcessed[0];
        for (i = 1; i < half; i++) {
            if (mProcessed[i] > mYMax)
            {
                mYMax = mProcessed[i];
            }
            else if (mProcessed[i] < mYMin)
            {
                mYMin = mProcessed[i];
            }
        }
    }
    else if (alg == EnhancedAutocorrelation)
    {
        for (i = 0; i < half; i++)
        {
            mProcessed[i] = mProcessed[i] / windows;
        }

      // Peak Pruning as described by Tolonen and Karjalainen, 2000

      // Clip at zero, copy to temp array
      for (i = 0; i < half; i++) {
         if (mProcessed[i] < 0.0)
         {
            mProcessed[i] = (float)(0.0);
         }
         out[i] = mProcessed[i];
      }

      // Subtract a time-doubled signal (linearly interp.) from the original
      // (clipped) signal
      for (i = 0; i < half; i++)
      {
         if ((i % 2) == 0)
         {
             mProcessed[i] -= out[i / 2];
         }
         else
         {
            mProcessed[i] -= ((out[i / 2] + out[i / 2 + 1]) / 2);
         }
      }
      // Clip at zero again
      for (i = 0; i < half; i++)
      {
         if (mProcessed[i] < 0.0)
         {
             mProcessed[i] = (float)(0.0);
         }
      }

      // Find NEW min/max
      mYMin = mProcessed[0];
      mYMax = mProcessed[0];
      for (i = 1; i < half; i++)
         if (mProcessed[i] > mYMax)
            mYMax = mProcessed[i];
         else if (mProcessed[i] < mYMin)
            mYMin = mProcessed[i];
    }
    else if (alg == Cepstrum)
    {
        for (i = 0; i < half; i++)
        {
            mProcessed[i] = mProcessed[i] / windows;

        }

      // Find min/max, ignoring first and last few values
      {
         int ignore = 4;
         mYMin = mProcessed[ignore];
         mYMax = mProcessed[ignore];
         for (i = ignore + 1; i + ignore < half; i++)
            if (mProcessed[i] > mYMax)
               mYMax = mProcessed[i];
            else if (mProcessed[i] < mYMin)
               mYMin = mProcessed[i];
      }
    }

//    if (pYMin)
//        *pYMin = mYMin;
//    if (pYMax)
//        *pYMax = mYMax;
   return mProcessed;
}

float* WindowFunc(int whichFunction, int NumSamples, float *in)
{
   int extraSample = 0;
   switch (whichFunction)
   {
        case eWinFuncHamming:
        case eWinFuncHann:
        case eWinFuncBlackman:
        case eWinFuncBlackmanHarris:
            extraSample = 1;
            break;
        default:
            break;
        case eWinFuncBartlett:
          // PRL:  Do nothing here either
          // But I want to comment that the old function did this case
          // wrong in the second half of the array, in case NumSamples was odd
          // but I think that never happened, so I am not bothering to preserve that
            break;
   }
   in = NewWindowFunc(whichFunction, NumSamples, extraSample, in);

   return in;
}

float* NewWindowFunc(int whichFunction, int NumSamplesIn, int extraSample, float *in)
{
   int NumSamples = (int)NumSamplesIn;

   switch (whichFunction) {
   default:
      break;
   case eWinFuncRectangular:
      // Multiply all by 1.0f -- do nothing
      break;

   case eWinFuncBartlett:
   {
      // Bartlett (triangular) window
      const int nPairs = (NumSamples - 1) / 2; // whether even or odd NumSamples, this is correct
      const float denom = NumSamples / 2.0f;
      in[0] = 0.0f;
      for (int ii = 1;
           ii <= nPairs; // Yes, <=
           ++ii) {
         const float value = ii / denom;

         in[ii] *= value;
         in[NumSamples - ii] *= value;
      }
      // When NumSamples is even, in[half] should be multiplied by 1.0, so unchanged
      // When odd, the value of 1.0 is not reached
   }
      break;
   case eWinFuncHamming:
   {
      // Hamming
      const double multiplier = 2 * M_PI / NumSamples;
      static const double coeff0 = 0.54, coeff1 = -0.46;
      for (int ii = 0; ii < NumSamples; ++ii)
      {
         in[ii] *= coeff0 + coeff1 * cos(ii * multiplier);
      }
   }
      break;
   case eWinFuncHann:
   {
      // Hann
      const double multiplier = 2 * M_PI / NumSamples;
      static const double coeff0 = 0.5, coeff1 = -0.5;
      for (int ii = 0; ii < NumSamples; ++ii)
         in[ii] *= coeff0 + coeff1 * cos(ii * multiplier);
   }
      break;
   case eWinFuncBlackman:
   {
      // Blackman
      const double multiplier = 2 * M_PI / NumSamples;
      const double multiplier2 = 2 * multiplier;
      static const double coeff0 = 0.42, coeff1 = -0.5, coeff2 = 0.08;
      for (int ii = 0; ii < NumSamples; ++ii)
         in[ii] *= coeff0 + coeff1 * cos(ii * multiplier) + coeff2 * cos(ii * multiplier2);
   }
      break;
   case eWinFuncBlackmanHarris:
   {
      // Blackman-Harris
      const double multiplier = 2 * M_PI / NumSamples;
      const double multiplier2 = 2 * multiplier;
      const double multiplier3 = 3 * multiplier;
      static const double coeff0 = 0.35875, coeff1 = -0.48829, coeff2 = 0.14128, coeff3 = -0.01168;
      for (int ii = 0; ii < NumSamples; ++ii)
         in[ii] *= coeff0 + coeff1 * cos(ii * multiplier) + coeff2 * cos(ii * multiplier2) + coeff3 * cos(ii * multiplier3);
   }
      break;
   case eWinFuncWelch:
   {
      // Welch
      const float N = NumSamples;
      for (int ii = 0; ii < NumSamples; ++ii) {
         const float iOverN = ii / N;
         in[ii] *= 4 * iOverN * (1 - iOverN);
      }
   }
      break;
   case eWinFuncGaussian25:
   {
      // Gaussian (a=2.5)
      // Precalculate some values, and simplify the fmla to try and reduce overhead
      static const double A = -2 * 2.5*2.5;
      const float N = NumSamples;
      for (int ii = 0; ii < NumSamples; ++ii) {
         const float iOverN = ii / N;
         // full
         // in[ii] *= exp(-0.5*(A*((ii-NumSamples/2)/NumSamples/2))*(A*((ii-NumSamples/2)/NumSamples/2)));
         // reduced
         in[ii] *= exp(A * (0.25 + (iOverN * iOverN) - iOverN));
      }
   }
      break;
   case eWinFuncGaussian35:
   {
      // Gaussian (a=3.5)
      static const double A = -2 * 3.5*3.5;
      const float N = NumSamples;
      for (int ii = 0; ii < NumSamples; ++ii) {
         const float iOverN = ii / N;
         in[ii] *= exp(A * (0.25 + (iOverN * iOverN) - iOverN));
      }
   }
      break;
   case eWinFuncGaussian45:
   {
      // Gaussian (a=4.5)
      static const double A = -2 * 4.5*4.5;
      const float N = NumSamples;
      for (int ii = 0; ii < NumSamples; ++ii) {
         const float iOverN = ii / N;
         in[ii] *= exp(A * (0.25 + (iOverN * iOverN) - iOverN));
      }
   }
      break;
   }

   if (extraSample && whichFunction != eWinFuncRectangular) {
      double value = 0.0;
      switch (whichFunction) {
      case eWinFuncHamming:
         value = 0.08;
         break;
      case eWinFuncGaussian25:
         value = exp(-2 * 2.5 * 2.5 * 0.25);
         break;
      case eWinFuncGaussian35:
         value = exp(-2 * 3.5 * 3.5 * 0.25);
         break;
      case eWinFuncGaussian45:
         value = exp(-2 * 4.5 * 4.5 * 0.25);
         break;
      default:
         break;
      }
      in[NumSamples] = 1;
      in[NumSamples] *= value;
   }
   return in;
}

float* PowerSpectrum(int NumSamples, const float *In, float *Out)
{
    FFTParam hFFT = InitializeFFT(NumSamples);
    FFTParam *p;
    p = &hFFT;
    int i = 0;
    float *buffer;
    float *pFFT;


    pFFT = In;

    buffer = RealFFTf(pFFT, p);

    for (i = 1; i < NumSamples/2; i++)
    {
        Out[i]= (buffer[ hFFT.BitReversed[i] ]*buffer[ hFFT.BitReversed[i] ]) + (buffer[ hFFT.BitReversed[i] +1 ]*buffer[ hFFT.BitReversed[i]+1 ]);
    }


   // Handle the (real-only) DC and Fs/2 bins
   Out[0] = buffer[0]*buffer[0];
   Out[NumSamples / 2] = buffer[1]*buffer[1];

   return Out;
}

FFTParam InitializeFFT(int fftlen)
{
    int temp;
    FFTParam h;
    h.Points = fftlen / 2;
    float temparray[2*h.Points];
    int temparray2[h.Points];
    int i;

    for (int i = 0; i < 2*h.Points; i++) {
        temparray[i] = 0;
    }
    for (i = 0; i < h.Points; i++) {
        temparray2[i] = 0;
    }

    for(i = 0; i < h.Points; i++)
    {
        temp = 0;
        for(int mask = h.Points / 2; mask > 0; mask >>= 1){
            temp = (temp >> 1) + (i & mask ? h.Points : 0);
        }
        temparray2[i] = temp;
    }
    for(i = 0; i < h.Points; i++)
    {
        temparray[temparray2[i]  ]=(float)-sin(2*M_PI*i/(2*h.Points));
        temparray[temparray2[i]+1]=(float)-cos(2*M_PI*i/(2*h.Points));
    }

    h.SinTable = temparray;
    h.BitReversed = temparray2;

    return h;
}

float* RealFFTf(float *buffer, const FFTParam *h)
{
    float *A,*B;
    const float *sptr;
    const float *endptr1,*endptr2;
    const int *br1,*br2;
    float HRplus,HRminus,HIplus,HIminus;
    float v1,v2,sin,cos;

    int ButterfliesPerGroup = h->Points/2;

   /*
   *  Butterfly:
   *     Ain-----Aout
   *         \ /
   *         / \
   *     Bin-----Bout
   */

   endptr1 = buffer + h->Points * 2;

   while(ButterfliesPerGroup > 0)
   {
      A = buffer;
      B = buffer + ButterfliesPerGroup * 2;
      sptr = h->SinTable;

      while(A < endptr1)
      {
         sin = *sptr;
         cos = *(sptr+1);
         endptr2 = B;
         while(A < endptr2)
         {
            v1 = *B * cos + *(B + 1) * sin;
            v2 = *B * sin - *(B + 1) * cos;
            *B = (*A + v1);
            *(A++) = *(B++) - 2 * v1;
            *B = (*A - v2);
            *(A++) = *(B++) + 2 * v2;
         }
         A = B;
         B += ButterfliesPerGroup * 2;
         sptr += 2;
      }
      ButterfliesPerGroup >>= 1;
   }

   /* Massage output to get the output for a real input sequence. */
   br1 = h->BitReversed + 1;
   br2 = h->BitReversed + h->Points - 1;


   while(br1<br2)
   {
      sin=h->SinTable[*br1];
      cos=h->SinTable[*br1+1];
      A=buffer+*br1;
      B=buffer+*br2;
      HRplus = (HRminus = *A     - *B    ) + (*B     * 2);
      HIplus = (HIminus = *(A+1) - *(B+1)) + (*(B+1) * 2);
      v1 = (sin*HRminus - cos*HIplus);
      v2 = (cos*HRminus + sin*HIplus);
      *A = (HRplus  + v1) * (float)0.5;
      *B = *A - v1;
      *(A+1) = (HIminus + v2) * (float)0.5;
      *(B+1) = *(A+1) - HIminus;

      br1++;
      br2--;
   }
   /* Handle the center bin (just need a conjugate) */
   A=buffer+*br1+1;
   *A=-*A;
   /* Handle DC bin separately - and ignore the Fs/2 bin
   buffer[0]+=buffer[1];
   buffer[1]=(fft_type)0;*/
   /* Handle DC and Fs/2 bins separately */
   /* Put the Fs/2 value into the imaginary part of the DC bin */
   v1=buffer[0]-buffer[1];
   buffer[0]+=buffer[1];
   buffer[1]=v1;

   return buffer;
}

FFTOutput RealFFT(int NumSamples, const float *RealIn, float *RealOut, float *ImagOut)
{
   FFTParam hFFT = InitializeFFT(NumSamples);
   FFTParam *p;
   p = &hFFT;
   float *pFFT;
   int i = 0;
   float *buffer;
   FFTOutput Result;

   // Copy the data into the processing buffer
   pFFT = RealIn;

   // Perform the FFT
   buffer = RealFFTf(pFFT, p);

   // Copy the data into the real and imaginary outputs
   for (i = 1; i<(NumSamples / 2); i++)
    {
        RealOut[i]=buffer[hFFT.BitReversed[i]  ];
        ImagOut[i]=buffer[hFFT.BitReversed[i]+1];
    }

   // Handle the (real-only) DC and Fs/2 bins
   RealOut[0] = buffer[0];
   RealOut[NumSamples / 2] = buffer[1];
   ImagOut[0] = ImagOut[NumSamples / 2] = 0;
   // Fill in the upper half using symmetry properties
   for(i = NumSamples / 2 + 1; i < NumSamples; i++) {
      RealOut[i] =  RealOut[NumSamples-i];
      ImagOut[i] = -ImagOut[NumSamples-i];
   }
   Result.RealPart = RealOut;
   Result.ImagPart = ImagOut;

    return Result;
}

float * InverseRealFFT(int NumSamples, const float *RealIn, const float *ImagIn, float *RealOut, float * bffInit)
{
    FFTParam hFFT = InitializeFFT(NumSamples);
    int i = 0;
    FFTParam *ptrFFT;
    ptrFFT = &hFFT;
    float *ipFFT;

    ipFFT = bffInit;

   // Copy the data into the processing buffer
   for (i = 0; i < (NumSamples / 2); i++)
   {
       ipFFT[2*i  ] = RealIn[i];
   }
   if(ImagIn == NULL) {
      for (i = 0; i < (NumSamples / 2); i++)
      {
          ipFFT[2*i+1] = 0;
      }
   }
   else {
      for (i = 0; i < (NumSamples / 2); i++)
         ipFFT[2*i+1] = ImagIn[i];
   }
   // Put the fs/2 component in the imaginary part of the DC bin
   ipFFT[1] = RealIn[NumSamples / 2];

   // Perform the FFT
   InverseRealFFTf(ipFFT, ptrFFT);

   // Copy the data to the (purely real) output buffer
   RealOut = ReorderToTime(ptrFFT, ipFFT, RealOut);

   return RealOut;
}

void InverseRealFFTf(float *buffer, const FFTParam *h)
{
   float *A,*B;
   const float *sptr;
   const float *endptr1,*endptr2;
   const int *br1;
   float HRplus,HRminus,HIplus,HIminus;
   float v1,v2,sin,cos;

   int ButterfliesPerGroup = h->Points / 2;

   /* Massage input to get the input for a real output sequence. */
   A = buffer + 2;
   B = buffer + h->Points * 2 - 2;
   br1 = h->BitReversed + 1;

   while(A<B)
   {
      sin=h->SinTable[*br1];
      cos=h->SinTable[*br1+1];
      HRplus = (HRminus = *A     - *B    ) + (*B     * 2);
      HIplus = (HIminus = *(A+1) - *(B+1)) + (*(B+1) * 2);
      v1 = (sin*HRminus + cos*HIplus);
      v2 = (cos*HRminus - sin*HIplus);
      *A = (HRplus  + v1) * (float)0.5;
      *B = *A - v1;
      *(A+1) = (HIminus - v2) * (float)0.5;
      *(B+1) = *(A+1) - HIminus;

      A+=2;
      B-=2;
      br1++;
   }

   /* Handle center bin (just need conjugate) */
   *(A+1)=-*(A+1);

   /* Handle DC bin separately - this ignores any Fs/2 component
   buffer[1]=buffer[0]=buffer[0]/2;*/
   /* Handle DC and Fs/2 bins specially */
   /* The DC bin is passed in as the real part of the DC complex value */
   /* The Fs/2 bin is passed in as the imaginary part of the DC complex value */
   /* (v1+v2) = buffer[0] == the DC component */
   /* (v1-v2) = buffer[1] == the Fs/2 component */
   v1=0.5f*(buffer[0]+buffer[1]);
   v2=0.5f*(buffer[0]-buffer[1]);

   buffer[0]=v1;
   buffer[1]=v2;

   /*
   *  Butterfly:
   *     Ain-----Aout
   *         \ /
   *         / \
   *     Bin-----Bout
   */

   endptr1 = buffer + h->Points * 2;

   while(ButterfliesPerGroup > 0)
   {
      A = buffer;
      B = buffer + ButterfliesPerGroup * 2;
      sptr = h->SinTable;

      while(A < endptr1)
      {
         sin = *(sptr++);
         cos = *(sptr++);
         endptr2 = B;

         while(A < endptr2)
         {
            v1 = *B * cos - *(B + 1) * sin;
            v2 = *B * sin + *(B + 1) * cos;
            *B = (*A + v1) * (float)0.5;
            *(A++) = *(B++) - v1;
            *B = (*A + v2) * (float)0.5;
            *(A++) = *(B++) - v2;
         }
         A = B;
         B += ButterfliesPerGroup * 2;
      }
      ButterfliesPerGroup >>= 1;
   }
}

float* ReorderToTime(const FFTParam *hFFT, const float *buffer, float *TimeOut)
{
   // Copy the data into the real outputs
   for(int i = 0; i < hFFT->Points; i++) {
      TimeOut[i*2  ]=buffer[hFFT->BitReversed[i]  ];
      TimeOut[i*2+1]=buffer[hFFT->BitReversed[i]+1];
   }

   return TimeOut;
}
