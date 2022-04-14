# Audacity_Spectrum_For_Matlab
Here I translated some of the Audacity Spectrum Analysing tools into Matlab code. The results are pretty similar.

The AudacitySpectrum file imports a waveform as input and outputs the corresponding Spectrum according to the chosen settings.
The logic is based on Audacity SpectrumAnalyst function.

Tested myself a few cases, it is useful, but of course the original C++ version is a way faster.

I also provided my C implementation of the algorithms. In this case, it receives as input a .txt file corresponding to the audio samples. 
Output is another .txt file with the Spectrum information.
