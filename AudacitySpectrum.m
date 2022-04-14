%% Program Init %% 

%clc; close all; clear all;

%% Read Waveform

filename = 'C:\Users\Charles\Documents\MATLAB\Wav Files\Test.wav';

%% Options Input

alg = 4;
windowFunc = 8;
windowSize = 2^9;
LogPlot = 0;

%% Spectrum Computation

[x,fs] = audioread(filename);

%% Extract Data
if size(x,2) == 2
    data = x(:,1) + x(:,2); % Comment this line if the file is mono
else
    data = x(:,1);
end

Algorithm =   ["Spectrum", "Autocorrelation", "CubeRootAutocorrelation", "EnhancedAutocorrelation", "Cepstrum"];

eWindowFunctions = ["eWinFuncRectangular","eWinFuncBartlett","eWinFuncHamming", "eWinFuncHann", ...
   "eWinFuncBlackman", "eWinFuncBlackmanHarris", "eWinFuncWelch", "eWinFuncGaussian25", "eWinFuncGaussian35",...
   "eWinFuncGaussian45"];
tic;
[Points, SinTable, BitReversed] = InitializeFFT(windowSize);
Spec = SpectrumAnalyst(Algorithm(alg), eWindowFunctions(windowFunc), windowSize, fs, data', size(data,1), Points, SinTable, BitReversed);
toc;
DisplaySpectrum(alg, LogPlot, fs, windowSize,Spec);


%% Function Definition

function DisplaySpectrum(alg, LogPlot, fs, windowSize,Spec)

    figure(1);
    if alg == 1
       f = linspace(fs/windowSize,fs/2, windowSize/2);
       if LogPlot == 1
           semilogx(f(1,1:end-1),Spec(1,2:end));
       else
           plot(f(1,1:end-1),Spec(1,2:end));
           grid on;
       end
    else
        t = linspace(1/fs, windowSize/(2*fs), windowSize/2);
        f = 1./t;
        
        fmin = f(1,end);
        fmax = f(1,1);
        ifmin = find(abs(fmin-f') == min(abs(fmin-f'))); % Compute the index related to the minimum frequency
        ifmax = find(abs(fmax-f') == min(abs(fmax-f'))); % Compute the index related to the maximum frequency

        %title('Short-Time Fourier Transform Representation of the Signal'); %Define the title for the figure
        xlabel('Note'); % Define the label for the x-axis of the figure
        ylabel('Amplitude');   %Define the label for the y-axis of the figure

        hold on; % Hold the plot
        notes = ["sil";...
                                               "C0";"C#0";"D0";"Eb0";"E0";"F0"; "F#0";...
                 "G0"; "Ab0"; "A0"; "Bb0";"B0";"C1";"C#1";"D1";"Eb1";"E1";"F1"; "#F#1";...
                 "G1"; "Ab1"; "A1"; "Bb1";"B1";"C2";"C#2";"D2";"Eb2";"E2";"F2"; "F#2";...
                 "G2"; "Ab2"; "A2"; "Bb2";"B2";"C3";"C#3";"D3";"Eb3";"E3";"F3"; "F#3";...
                 "G3"; "Ab3"; "A3"; "Bb3";"B3";"C4";"C#4";"D4";"Eb4";"E4";"F4"; "F#4";...
                 "G4"; "Ab4"; "A4"; "Bb4";"B4";"C5";"C#5";"D5";"Eb5";"E5";"F5"; "F#5";...
                 "G5"; "Ab5"; "A5"; "Bb5";"B5";"C6";"C#6";"D6";"Eb6";"E6";"F6"; "F#6"; ...
                 "G6"; "Ab6"; "A6"; "Bb6";"B6";"C7";"C#7";"D7";"Eb7";"E7";"F7"; "F#7";...
                 "G7"; "Ab7"; "A7"; "Bb7";"B7";"C8";"C#8";"D8";"Eb8";"E8";"F8"; "F#8";...
                 "G8"; "Ab8"; "A8"; "Bb8";"B8"]; % Define the names of the musical notes
             %1       1        1           1        1        2        2          2 
        fnotes = [0;...
                                                            16.35  ; 17.32   ; 18.35   ; 19.45  ; 20.60  ; 21.83   ; 23.12   ; ...
            24.50  ; 25.96   ; 27.50  ; 29.14   ; 30.87   ; 32.7   ; 34.65   ; 36.71   ; 38.89  ; 41.20  ; 43.65   ; 46.25   ; ...
            49     ; 51.91   ; 55     ; 58.27   ; 61.74   ; 65.41  ; 69.3    ; 73.42   ; 77.78  ; 82.41  ; 87.31   ; 92.5    ; ... 
            98     ; 103.83  ; 110    ; 116.54  ; 123.47  ; 130.81 ; 138.59  ; 146.83  ; 155.56 ; 164.81 ; 174.61  ; 185     ; ...
            196.00 ; 207.65  ; 220.00 ; 233.08  ; 246.94  ; 261.63 ; 277.18  ; 293.66  ; 311.13 ; 329.63 ; 349.23  ; 369.99  ; ...
            392.00 ; 415.30  ; 440.00 ; 466.16  ; 493.88  ; 523.25 ; 554.37  ;  587.33 ; 622.25 ; 659.25 ; 698.46  ; 739.99  ; ...
            783.99 ; 830.61  ; 880.00 ; 932.33  ; 987.77  ; 1046.50; 1108.73 ; 1174.66 ; 1244.51; 1318.51; 1396.91 ; 1479.98 ; ...
            1567.98 ; 1661.22; 1760.00; 1864.66 ; 1975.53; 2093.00 ; 2217.46 ; 2349.32; 2489.02 ; 2637.02 ; 2793.83; 2959.96 ;...
            3135.96 ; 3322.44; 3520.00; 3729.31 ; 3951.07; 4186.01 ; 4434.92 ; 4698.63; 4978.03 ; 5274.04 ; 5587.65; 5919.91 ;...
            6271.93 ; 6644.88; 7040.00; 7458.62; 7902.13];
        %     G        Ab        A        Bb        B        C          C#       D         Eb       E         F         F#

        k1 = f(ifmin);
        ind1 = find(abs(k1-fnotes') == min(abs(k1-fnotes')));
        k2 = f(ifmax);
        ind2 = find(abs(k2-fnotes') == min(abs(k2-fnotes')));

        % Plot the Spectrum
        s = area(t(1,1:end-1),Spec(1,2:end),'FaceColor', [0.784 0.196 0.588], 'EdgeColor', [0.9 0.39 0.9]);
        dtt = s.DataTipTemplate; % Get the DataTip information of this plot
        dtt.DataTipRows(1).Label = 'Freq'; % Define the first label as Time
        dtt.DataTipRows(1).Value = 1./s.XData;
        dtt.DataTipRows(2).Label = 'Amplitude';
        
        % Plot the notes lines
        for i = ind1:ind2
            s = plot(1/fnotes(i).*ones(1, 10), linspace(min(Spec(:)),max(Spec(:)),10) , 'Color', [0.8, 0.8, 0.8]);
            dtt = s.DataTipTemplate; % Get the DataTip information of this plot
            dtt.DataTipRows(1).Label = 'Freq'; % Define the first label as Time
            dtt.DataTipRows(1).Value = fnotes(i).*ones(1, 10);
            dtt.DataTipRows(2).Label = 'Note'; % Create a third DataTip called 'Note'
            dtt.DataTipRows(2).Value = repmat(notes(i),1,10); % Attribute the correct note to that DataTip
        end
        set(gca, 'YGrid', 'on', 'XGrid', 'off'); %Y-Axis Grid on
        hold off; % Release the plot
        xticks(1./fnotes(ind2:-1:ind1));
        xticklabels(notes(ind2:-1:ind1));
        axis([t(1), t(end), min(Spec(:)), max(Spec(:))]);

    end
end

function Spec = SpectrumAnalyst(alg, windowFunc, windowSize, ~, data,dataLen, Points, SinTable, BitReversed)
    
    mWindowSize = windowSize;
    half = mWindowSize / 2;
    mProcessed = zeros([1,mWindowSize/2]);
    in =  zeros([1,mWindowSize]);
    out =  zeros([1,mWindowSize]);
    win= ones([1,mWindowSize]);
    
    win = WindowFunc(windowFunc, mWindowSize, win);

    wss = sum(win(1,:));
    
    if(wss > 0)
        wss = 4 / (wss*wss);
    else
        wss = 1;
    end
    
    start = 0;
    windows = 0;
    
    while (start + mWindowSize <= dataLen)
       
       in(1,1:mWindowSize) = win(1,1:mWindowSize) .* data(1,start + 1:start+mWindowSize);

       if (alg == "Spectrum")
          out = PowerSpectrum(mWindowSize, in, Points, SinTable, BitReversed);

          mProcessed(1,1:half) = mProcessed(1,1:half) + out(1,1:half);

          
       elseif (alg == "Autocorrelation" || alg == "CubeRootAutocorrelation" || alg == "EnhancedAutocorrelation")
           
           [RealPart, ImagPart] = RealFFT(mWindowSize, in, Points, SinTable, BitReversed); 
           
           in = RealPart.*RealPart + ImagPart.*ImagPart;
               
           if (alg == "Autocorrelation")
               in = sqrt(in);
           end
           if (alg == "CubeRootAutocorrelation" || alg == "EnhancedAutocorrelation")
               in = in.^(1/3);
           end

           [RealPart, ~] = RealFFT(mWindowSize, in, Points, SinTable, BitReversed);
           
           mProcessed(1,1:mWindowSize/2) = mProcessed(1,1:mWindowSize/2) + RealPart(1,1:mWindowSize/2);

       elseif (alg == "Cepstrum")
           
           [RealPart, ImagPart] = RealFFT(mWindowSize, in, Points, SinTable, BitReversed);
           minpower = 1e-20*mWindowSize*mWindowSize;

            power = RealPart.*RealPart + ImagPart.*ImagPart;
            in(1,power < minpower) = log(minpower);
            in(1, power >= minpower) = log(power);
     
           out = InverseRealFFT(mWindowSize, in, 0, Points, SinTable, BitReversed);
           
           mProcessed(1,1:half) = mProcessed(1,1:half) + out(1,1:half);
       end
       start = start + half;
       windows = windows+1;
    end
    
     mYMin = 1000000;
     mYMax = -1000000;
     if alg == "Spectrum"
        scale = wss / windows;

        for i = 1:half
            mProcessed(1,i) = 10 * log10(mProcessed(1,i) * scale);
            if mProcessed(1,i) > mYMax
                mYMax = mProcessed(1,i);
            elseif mProcessed(1,i) < mYMin
                mYMin = mProcessed(1,i);
            end
        end
        
     elseif (alg == "Autocorrelation" || alg == "CubeRootAutocorrelation")

         mProcessed(1,1:half) = mProcessed(1,1:half) ./ windows;
         
         mYMin = mProcessed(1,1);
         mYMax = mProcessed(1,1);
         for i = 2:half
            if mProcessed(1,i) > mYMax
               mYMax = mProcessed(1,i); 
            end
            if mProcessed(1,i) < mYMin
                mYMin = mProcessed(1,i);
            end
         end
         
     elseif alg == "EnhancedAutocorrelation"
       
         mProcessed(1,1:half) = mProcessed(1,1:half) ./ windows;
         
         for i = 1:half
             if mProcessed(1,i) < 0
                 mProcessed(1,i) = 0;
             end
             out(1,i) = mProcessed(1,i);
         end
         
         for i = 2:half
            if rem(i,2) ~= 0
                mProcessed(1,i) = out(1,i) - out(1, ceil(i / 2)); 
            else
                mProcessed(1,i) = out(1,i) - ((out(1,floor(i / 2)) + out(1, floor(i / 2) + 1)) / 2);
            end
         end
         mProcessed(1,1) = 0;
         for i = 1:half
            if mProcessed(1,i) < 0
               mProcessed(1,i) = 0; 
            end
         end
         
         mYMin = mProcessed(1,1);
         mYMax = mProcessed(1,1);
         for i = 2:half
             if (mProcessed(1,i) > mYMax)
                 mYMax = mProcessed(1,i);
             elseif (mProcessed(1,i) < mYMin)
                 mYMin = mProcessed(1,i);
             end
         end
     elseif alg == "Cepstrum"
        
         mProcessed(1,1:half) = mProcessed(1,1:half) ./ windows;

         ignore = 4;
         mYMin = mProcessed(1,ignore);
         mYMax = mProcessed(1,ignore);
         i = ignore+1;
         while i + ignore < half
             if mProcessed(1,i) > mYMax
                 mYMax = mProcessed(1,i);
             elseif mProcessed(1,i) < mYMin
                 mYMin = mProcessed(1,i);
             end
            i = i+1; 
         end
         
     end
    Spec = mProcessed;
end

function Out = PowerSpectrum(NumSamples, In, Points, SinTable, BitReversed)
    
    pFFT = In;
    buffer = RealFFTf(pFFT,Points,SinTable,BitReversed);

    Out = zeros([1, NumSamples/2]);
    for i = 2:NumSamples/2
        Out(1,i)= buffer(1,BitReversed(1,i)+1)*buffer(1,BitReversed(1,i)+1) + buffer(1,BitReversed(1,i)+2)*buffer(1,BitReversed(1,i)+2);
    end

   Out(1,1) = buffer(1,1)*buffer(1,1);
   Out(1,NumSamples / 2 +1) = buffer(1,2)*buffer(1,2); 
end

function [Points, SinTable, BitReversed] = InitializeFFT(fftlen)
    
    Points = fftlen / 2;
    SinTable = zeros(1,2*Points);
    BitReversed = zeros(1,Points);
    
    for i = 0:Points-1
       temp = 0;
       mask = Points/2;
       while mask > 0
          if bitand(i,mask) > 0
              temp = bitshift(temp,-1) + Points;
          else
              temp = bitshift(temp,-1);
          end
          mask = bitshift(mask,-1);
       end
       BitReversed(1,i+1) = temp;
    end
    
   for i = 1:Points
      SinTable(1,BitReversed(1,i)+1  )=-sin(2*pi*(i-1)/(2*Points));
      SinTable(1,BitReversed(1,i)+2)=-cos(2*pi*(i-1)/(2*Points));
   end
end

function buffer = RealFFTf(buffer, Points, SinTable, BitReversed)
    
    ButterfliesPerGroup = Points/2;
    
    while(ButterfliesPerGroup > 0)
       A = buffer;
       B = buffer(1,2*floor(ButterfliesPerGroup)+1:length(buffer));
       sptr = SinTable;
       
       for i = 1:Points/(2*ButterfliesPerGroup)
           sin = sptr(1,1);
           cos = sptr(1,2);
           for j = 1:ButterfliesPerGroup
               v1 = B(1,1) * cos + B(1,2) * sin;
               v2 = B(1,1) * sin - B(1,2) * cos;
               B(1,1) = A(1,1) + v1;
               buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j-1) = B(1,1);
               buffer(1,(i-1)*4*ButterfliesPerGroup + 2*j-1) = buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j-1)-2*v1;
               A = circshift(A,[1 -1]);
               B = circshift(B,[1 -1]);
               A(end-2*j+2) = 0;
               B(end-2*j+2) = 0;
               B(1,1) = A(1,1) - v2;
               buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j) = B(1,1);
               buffer(1,(i-1)*4*ButterfliesPerGroup + 2*j) = buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j)+2*v2;
               A = circshift(A,[1 -1]);
               B = circshift(B,[1 -1]);
               A(end-2*j+1) = 0;
               B(end-2*j+1) = 0;
           end
         A = B;
         B = buffer(1,(6+4*(i-1))*floor(ButterfliesPerGroup) + 1:length(buffer));
         sptr = circshift(sptr,[1 -2]);
       end
       ButterfliesPerGroup = floor(ButterfliesPerGroup/2);
    end
    
    for i = 1:floor(Points/2)-1
        br1 = BitReversed(1,1+i);
        br2 = BitReversed(1,end+1-i);
        sin = SinTable(br1+1);
        cos = SinTable(br1+2);
        A = buffer(1,br1+1:end);
        B = buffer(1,br2+1:end);
        HRminus = A(1,1) - B(1,1);
        HRplus = HRminus + 2*B(1,1);
        HIminus = A(1,2) - B(1,2);
        HIplus = HIminus + 2*B(1,2);
        v1 = (sin*HRminus - cos*HIplus);
        v2 = (cos*HRminus + sin*HIplus);
        A(1,1) = 0.5 * (HRplus + v1);
        B(1,1) = A(1,1) - v1;
        A(1,2) = 0.5 * (HIminus + v2);
        B(1,2) = A(1,2) - HIminus;
        buffer(1,br1+1:br1+2) = A(1,1:2);
        buffer(1,br2+1:br2+2) = B(1,1:2); 
        
    end
    buffer(1,BitReversed(1,2+i)+2) = - buffer(1,BitReversed(1,2+i)+2);

   v1=buffer(1,1)-buffer(1,2);
   buffer(1,1)= buffer(1,1)+ buffer(1,2);
   buffer(1,2)=v1;
end

function in = WindowFunc(whichFunction, NumSamples, in)

    extraSample = 0;
    switch (whichFunction)
        case "eWinFuncHamming"
            extraSample = 1;
        case "eWinFuncHann"
            extraSample = 1;
        case "eWinFuncBlackman"
            extraSample = 1;
        case "eWinFuncBlackmanHarris"
            extraSample = 1;
        otherwise
    end
    
    in = NewWindowFunc(whichFunction, NumSamples, extraSample, in);
end

function in = NewWindowFunc(whichFunction, NumSamplesIn, extraSample, in)
    
    NumSamples = NumSamplesIn;
    
    switch (whichFunction)
        case "eWinFuncRectangular"
        case "eWinFuncBartlett"
            nPairs = floor((NumSamples - 1) / 2);
            denom = NumSamples / 2;
            
            in(1,1) = 0;
            for ii = 1:nPairs
                value = ii / denom;
                in(1,ii+1) = in(1,ii+1) * value;
                in(1,NumSamples + 1 - ii) = in(1,NumSamples + 1 - ii) * value;
            end
        case "eWinFuncHamming"
            multiplier = 2 * pi / NumSamples;
            coeff0 = 0.54;
            coeff1 = -0.46;
            for ii = 1:NumSamples
                in(1,ii) = in(1,ii) * (coeff0 + coeff1 * cos((ii - 1) * multiplier));
            end
        case "eWinFuncHann"
            multiplier = 2 * pi / NumSamples;
            coeff0 = 0.5;
            coeff1 = -0.5;
            for ii = 1:NumSamples
                in(1,ii) = in(1,ii)* (coeff0 + coeff1 * cos((ii - 1) * multiplier));
            end
        case "eWinFuncBlackman"
            multiplier = 2 * pi / NumSamples;
            multiplier2 = 2 * multiplier;
            coeff0 = 0.42; 
            coeff1 = -0.5; 
            coeff2 = 0.08;
            for ii = 1:NumSamples
                in(1,ii) = in(1,ii)* ( coeff0 + coeff1 * cos((ii - 1) * multiplier) + coeff2 * cos((ii - 1) * multiplier2));
            end
        case "eWinFuncBlackmanHarris"
            multiplier = 2 * pi / NumSamples;
            multiplier2 = 2 * multiplier;
            multiplier3 = 3 * multiplier;
            coeff0 = 0.35875; 
            coeff1 = -0.48829;
            coeff2 = 0.14128; 
            coeff3 = -0.01168;
            for ii = 1:NumSamples
                in(1,ii) = in(1,ii) * ( coeff0 + coeff1 * cos((ii - 1) * multiplier) + coeff2 * cos((ii - 1) * multiplier2) + coeff3 * cos((ii - 1) * multiplier3));
            end
        case "eWinFuncWelch"
            N = NumSamples;
            for ii = 1:NumSamples
                iOverN = (ii-1) / N;
                in(1,ii) = in(1,ii) * (4 * iOverN * (1 - iOverN));
            end
        case "eWinFuncGaussian25"
              A = -2 * 2.5 * 2.5;
              N = NumSamples;
              for ii = 1:NumSamples
                 iOverN = (ii-1) / N;
                 in(1,ii) =  in(1,ii) * exp(A * (0.25 + (iOverN * iOverN) - iOverN));
              end
        case "eWinFuncGaussian35"
            A = -2 * 3.5*3.5;
            N = NumSamples;
            for ii = 1:NumSamples
                iOverN = (ii-1) / N;
                in(1,ii) = in(1,ii) * exp(A * (0.25 + (iOverN * iOverN) - iOverN));
            end
        case "eWinFuncGaussian45"
            A = -2 * 4.5*4.5;
            N = NumSamples;
            for ii = 1:NumSamples
                iOverN = (ii-1) / N;
                in(1,ii) = in(1,ii) * exp(A * (0.25 + (iOverN * iOverN) - iOverN));
            end
        otherwise
    end
    
    if extraSample && whichFunction ~= "eWinFuncRectangular"
        value = 0;
        switch whichFunction
            case "eWinFuncHamming"
                value = 0.08;
            case "eWinFuncGaussian25"
                value = exp(-2 * 2.5 * 2.5 * 0.25);
            case "eWinFuncGaussian35"
                value = exp(-2 * 3.5 * 3.5 * 0.25);
            case "eWinFuncGaussian45"
                value = exp(-2 * 4.5 * 4.5 * 0.25);
            otherwise
        end
      %in(1,NumSamples) = in(1,NumSamples) * value;
    end
end

function [RealOut, ImagOut] = RealFFT(NumSamples, RealIn, Points, SinTable, BitReversed)

    buffer = RealFFTf(RealIn, Points, SinTable, BitReversed);

    RealOut = zeros([1, NumSamples]);
    ImagOut = zeros([1, NumSamples]);
    
    RealOut(1,2:NumSamples/2) = buffer(1,BitReversed(1,2:NumSamples/2)+1);
    ImagOut(1,2:NumSamples/2) = buffer(1,BitReversed(1,2:NumSamples/2)+2);

    RealOut(1,1) = buffer(1,1);
    RealOut(NumSamples / 2 + 1) = buffer(1,2);
    ImagOut(1,1) = 0;
    ImagOut(1, NumSamples / 2 + 1) = 0;
    
    for i = NumSamples / 2+2:NumSamples
        RealOut(1,i) =  RealOut(1,NumSamples+2-i);
        ImagOut(1,i) = -ImagOut(1,NumSamples+2-i);
    end
end

function TimeOut = InverseRealFFT(NumSamples, RealIn, ImagIn, Points, SinTable, BitReversed)

    in = zeros([1, NumSamples]);
    for i = 1:NumSamples/2
       in(1,2*i-1) = RealIn(1,i);
    end
    if ImagIn ~= 0
        for i = 1: NumSamples/2-1
            in(1,2*i) = ImagIn(1,i);
        end
    end
    in(1,2) = RealIn(1, NumSamples/2+1);

    RealOut = InverseRealFFTf(in, Points, SinTable, BitReversed);
    
    TimeOut = ReorderToTime(Points, SinTable, BitReversed, RealOut);   
end

function buffer = InverseRealFFTf(buffer, Points, SinTable, BitReversed)
    
    ButterfliesPerGroup = Points / 2;
    
    A = buffer(1,3:end);
    B = buffer(1,2*Points-1:end);
    br1 = BitReversed(1,2);

    for i = 1:floor(Points/2)-1
        sin = SinTable(1, br1+1);
        cos = SinTable(1, br1+2);
        HRminus = A(1,1) - B(1,1);
        HRplus = HRminus + 2*B(1,1);
        HIminus = A(1,2) - B(1,2);
        HIplus = HIminus + 2*B(1,2);
        v1 = (sin*HRminus + cos*HIplus);
        v2 = (cos*HRminus - sin*HIplus);
        A(1,1) = 0.5*(HRplus + v1);
        B(1,1) = A(1,1) - v1;
        A(1,2) = 0.5*(HIminus-v2);
        B(1,2) = A(1,2) - HIminus;
        buffer(1,3+2*(i-1):4+2*(i-1)) = A(1,1:2);
        buffer(1,2*Points-1-2*(i-1):2*Points-2*(i-1)) = B(1,1:2); 

        A = buffer(1,3+2*i:end);
        B = buffer(2*Points-1-2*i:end);
        br1 = BitReversed(1,2+i);
    end
    
    buffer(1,Points+2) = -buffer(1,Points+2);
    
    v1 = 0.5*(buffer(1,1) + buffer(1,2));
    v2 = 0.5*(buffer(1,1) - buffer(1,2));
    buffer(1,1) = v1;
    buffer(1,2) = v2;
    
    while(ButterfliesPerGroup > 0)
       A = buffer;
       B = buffer(1,2*floor(ButterfliesPerGroup)+1:length(buffer));
       sptr = SinTable;
       
       for i = 1:Points/(2*ButterfliesPerGroup)
           sin = sptr(1,2*(i-1)+1);
           cos = sptr(1,2*i);
           for j = 1:ButterfliesPerGroup
               v1 = B(1,1) * cos - B(1,2) * sin;
               v2 = B(1,1) * sin + B(1,2) * cos;
               B(1,1) = 0.5*(A(1,1) + v1);
               buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j-1) = B(1,1);
               buffer(1,(i-1)*4*ButterfliesPerGroup + 2*j-1) = buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j-1)-v1;
               A = circshift(A,[1 -1]);
               B = circshift(B,[1 -1]);
               A(end-2*j+2) = 0;
               B(end-2*j+2) = 0;
               B(1,1) = 0.5*(A(1,1) + v2);
               buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j) = B(1,1);
               buffer(1,(i-1)*4*ButterfliesPerGroup + 2*j) = buffer(1,2*(2*i-1)*ButterfliesPerGroup + 2*j)-v2;
               A = circshift(A,[1 -1]);
               B = circshift(B,[1 -1]);
               A(end-2*j+1) = 0;
               B(end-2*j+1) = 0;
           end
         A = B;
         B = buffer(1,(6+4*(i-1))*floor(ButterfliesPerGroup) + 1:length(buffer));
       end
       ButterfliesPerGroup = floor(ButterfliesPerGroup/2);
    end

    
end

function TimeOut = ReorderToTime(Points, ~, BitReversed, buffer)
    TimeOut = zeros([1,Points*2]);
   for i = 1:Points
      TimeOut(1,2*i-1)=buffer(1,BitReversed(1,i)+1 );
      TimeOut(1,2*i)=buffer(BitReversed(1,i)+2);
   end
end