
%%%% Read In Data %%%
clearvars
NEOreadin
%%
clearvars -except complexData method
%% Spectral Data

specData = real(complexData);

% Extract FWHM from abs. spectrum
ftSpec = fftshift(fft(specData));
ftSpec = abs(ftSpec) ./ max(abs(ftSpec));

lenSpec = length(ftSpec);

specWidth = method.PVM_SpecSWH; % Hz
%w0 = getValues.get_w0(21.1);
%w0 = method.PVM_FrqRef(1) * 1e6; 

df = specWidth/(lenSpec-1);

ftXrange = -specWidth/2:df:specWidth/2; % Hz

[maxVal, peakIndex] = max(ftSpec);

ppmRange = ftXrange ./ (method.PVM_FrqRef(1)); %ppm


% Does usially not take values exactly at 0.5 ... need to interpolate/fit
% Left of the peak
leftIndex = find(ftSpec(1:peakIndex) <= 0.5, 1, 'last');

% Right of the peak
rightIndex = find(ftSpec(peakIndex:end) <= 0.5, 1, 'first') + peakIndex - 1;

% Step 4: Calculate FWH
FWHM = ftXrange(rightIndex) - ftXrange(leftIndex);

figure();
plot(ppmRange, ftSpec, 'LineWidth',2);
%xline(w0 , '--', 'DisplayName', 'Larmor Frequency');
xlabel("Frequency [MHz] / ppm");
ylabel("Signal norm.");
title(sprintf('FWHM: %0.3f Hz', FWHM));
set(gca,'FontSize',14);
set(gca, 'fontname', 'Palatino');
set(gca,'linewidth',1.5);
legend();
grid on;

%% Prepare Data for Fit

fitVectorSPcorr = SPGibbsCorr(specData);

if fitVectorSPcorr(1) < 0
    fitVectorSPcorr = fitVectorSPcorr .* (-1);
end

fitVectorSP = fitVectorSPcorr ./ max(fitVectorSPcorr);
%fitVectorSP  = specData(78:end) ./ max(specData(78:end));

fitVectorX = (0:1:size(fitVectorSP,1)-1) .* 0.1;

if size(fitVectorX,1) ~= size(fitVectorSP,1)
    fitVectorSP = fitVectorSP';
end



%% Fit

sqStopIdx = 1800;

funSQ = @(x,xdata)  x(1)*exp(-xdata/x(2)) + x(3)*exp(-xdata/x(4))+ x(5);

x0 = [0.4 33 0.6 5 0 2 0];
x0low =  [0  0 0  0 -0.1  0 -0.1];
x0high = [2 70 2 50  0.1 20  0.1];


[X,~,residual,~,~,lambda,jacobian] = lsqcurvefit(funSQ, x0, fitVectorX, fitVectorSP, x0low, x0high);
conf = nlparci(X,residual,'jacobian',jacobian);


if X(2) < X(4) % change values if T2s smaller than T2f
    print("Changing T2 fit values");
    T2S = X(4);
    stdT2S = max([X(4)-conf(4,1) conf(4,2)-X(4)]);
    T2F = X(2);
    stdT2F = max([X(2)-conf(2,1) conf(2,2)-X(2)]);
    ASQF = X(1);
    ASQS = X(3);
    ATQ = X(6);
    DC = X(7);
    DCSQ = X(5);
else
    T2S = X(2);
    stdT2S = max([X(2)-conf(2,1) conf(2,2)-X(2)]);
    T2F = X(4);
    stdT2F = max([X(4)-conf(4,1) conf(4,2)-X(4)]);
    ATQ = X(6);
    DC = X(7);
    ASQS = X(1);
    ASQF = X(3);
    DCSQ = X(5);
    

end


%% Plot

figure();hold on;
fplot(@(x) ASQS.* exp(-x./T2S) + ASQF.*exp(-x./T2F)+ DCSQ, [fitVectorX(1) 100], 'linewidth',2);hold on;
fplot(@(x)(exp(-x./T2S) - exp(-x./T2F))+ DCSQ, [fitVectorX(1) 100],'LineWidth',2);
plot(fitVectorX(1:1000),fitVectorSP(1:1000));
%xline(evoTime, 'linewidth',2, Color=[0.9290 0.6940 0.1250]);
xlabel("Acquisition Time / [ms]")
ylabel("SQ Signal normalized");
title(sprintf("Single Pulse: T_{2s}: %0.2f ms, T_{2f}: %0.2f ms",T2S, T2F));
set(gca,'FontSize',14);
set(gca, 'fontname', 'Palatino');
set(gca,'linewidth',2);
grid on;
legend("Fit","TQ extraction", "Data");
