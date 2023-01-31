function [PowerSpec,Peak,Med,Mean,picMask,PeakPhase,mask,nframe,PeakPos] = fft_analysis_simple(data,CBF,varargin)
% adapted from the code written by Christa Ringers, Jan Niklas Hansen,...
% Benjamin Friedrich and Nathalie Jurisch-Yaksi 
% Reference is: Ringers et al, BioRXiv, 2021, https://doi.org/10.1101/2021.11.23.469646 


% sourceP = Path to the folder where the data is located
% name    = Name of the file containing the data
% Fs      = Frequency of acquisition
% power   = Do you want to retrieve the entire Powerspectrum, true or false.
% targetP = Path to the folder were the data should be saved

% Some information:
% The mask is automatically generated
% Peak frequency is calculated as the max value of the primary frequency map



%% Number of frames as small prime factors (not only 2^L)
nframe0=size(data,3); % number of frames should only have small prime factors (2^L best, but not necessary)
nframe=nframe0;
while max(factor(nframe))>5
    nframe=nframe-1;
end

% pixel dimensions of each frame
xmax=size(data,1);
ymax=size(data,2);

% A raw frame for further reference
raw10 = data(:,:,10);

%% Fourier transform ------------------------------------------------------
fprintf('Perform the FFT\n')

% CBF.w_min=15; % [Hz] lower frequency cut-off
iw_min=1+round(CBF.w_min/CBF.Fs*nframe); % corresponding frequency index
nyquist=floor(nframe/2) + 1; % Nyquist frequency OBS: I think nyquist needs a +1

% Allocate memory
PeakPos   = nan(xmax,ymax); % frequency index of peak position
PeakPower = nan(xmax,ymax); % power at peak position
PeakPhase = nan(xmax,ymax); % phase at peak position


PowerSpec = nan(xmax,ymax,nyquist); % entire powerspectrum
tic
for y=1:ymax
    
    % Fourier transform (vectorized for entire pixel column; entire movie would be too large)
    data_fft=squeeze(fft(single(data(:,y,1:nframe)),[],3) );
    
    % Power spectrum
    PowerSpec(:,y,1:nyquist)=abs(data_fft(:,1:nyquist));
    
    % Locate peak in power spectrum: first estimate
    [maxval,maxpos]=max(PowerSpec(:,y,iw_min:end),[],3);
    maxpos=maxpos+iw_min-1;
    
    % Store results
    PeakPos(:,y)   = maxpos;
    PeakPower(:,y) = maxval;
    for x=1:xmax
        PeakPhase(x,y)=angle(data_fft(x,maxpos(x)));
    end
end
toc


%% Retrieving some data

% Translate primary frequency indices (PeakPos) in the power spectra to
% actual frequencies in Hertz.
%pic = (PeakPos-1).*100./nframe;
pic = ((PeakPos-1).*CBF.Fs)./nframe;

% Create mask
if isempty(varargin)
    picSD=jnh_FreqFilterV1(pic,3,CBF.Fs);
    mask = create_mask(picSD);
    mask(mask == 0) = NaN;
else
    mask = varargin{1};
end

% Amplitude
ampl = 2*(PeakPower/ nframe); % Amplitude is just power rescaled/divided by the length of the original time-domain signal

% calculate values from masked CBF images
picMask = pic.*mask;
phase = PeakPhase.*mask;
amplMask = ampl.*mask;
Med = nanmedian(picMask(:));
Mean = nanmean(picMask(:));

% Calculate the average powerspectrum and identify its Peak frequency
MeanPowerSpec = nanmean(reshape(PowerSpec.*mask,[xmax*ymax,nyquist]));
[ ~, Peak] = max(MeanPowerSpec(iw_min:nyquist));
Peak = ((Peak+iw_min-1)/ nframe)*CBF.Fs;

% Mask boundary
mask_boundary=~edge(mask,'sobel',[],'nothinning');

%% Save data -------------------------------------------------------------
% save(fullfile(CBF.targetP,[CBF.name,'_fft.mat']),'PeakPos','PeakPhase','PeakPower', ...
%     'mask','nframe','raw10','PowerSpec', 'MeanPowerSpec', 'PSD','-v7.3');

%% Create figures  --------------------------------------------------------

figure; clf
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8])

% Raw data frame
ax1=subplot('Position', [0.05 0.52 0.45 0.4]);
imagesc(raw10);
colormap(ax1,'gray'), c=colorbar;
c.Label.String = 'pixel intensity';
title(ax1,'raw image')
box off; axis off; axis image
axis image
hold on
line([10 60/CBF.spatres],[CBF.x-10 CBF.x-10],'Color', 'k', 'LineWidth' , 2 );
text (30, CBF.x-50, '50µm')


% Primary frequency peak [raw]
ax2=subplot('Position', [0.05 0.05 0.45 0.4]);
imagesc(pic.*mask_boundary);
colormap(ax2,'jet'), caxis(ax2,CBF.caxis);
c=colorbar; c.Label.String = 'CBF (Hz)';
title(ax2,'CBF raw')
box off; axis off; axis image
hold on
line([10 60/CBF.spatres],[CBF.x-10 CBF.x-10],'Color', 'k', 'LineWidth' , 2 );


% Primary frequency peak [mask]
ax3=subplot('Position', [0.52 0.05 0.45 0.4]);
AlphDat =double(~isnan(picMask));
imagesc(picMask, 'AlphaData', AlphDat), colormap(ax3,jet);
caxis(ax3,CBF.caxis), c=colorbar; c.Label.String = 'CBF (Hz)';
title(ax3,'CBF masked')
box off; axis off; axis image
hold on
line([10 60/CBF.spatres],[CBF.x-10 CBF.x-10],'Color', 'k', 'LineWidth' , 2 );

% Freq spectrum mean
ax4=subplot('Position', [0.57 0.55 0.14 0.32]);
line ([Peak, Peak],[min(log10(MeanPowerSpec(2:end))*10),max(log10(MeanPowerSpec(2:end))*10-5)], 'LineWidth', 2,'Color', 'r')
hold on
plot((2:length(MeanPowerSpec))*(CBF.Fs/nframe),log10(MeanPowerSpec(2:end))*10,'k', 'linewidth', 2) 
ylabel('Power/frequency (dB/Hz)')
xlabel ('Frequency (Hz)')
xlim([0 CBF.Fs/2])
box off; grid off
text(5 ,max(log10(MeanPowerSpec(2:end))*10), sprintf('Peak Freq = %0.2f', Peak));
title ('Mean Power Spectrum')

% Histogram picMask
ax5=subplot('Position', [0.82 0.55 0.14 0.32]);
[N,edges] = histcounts(picMask(:),100, 'Normalization' ,'probability');
plot (edges(1:end-1)+ (edges(1)-edges(2))/2, N,'k', 'linewidth', 2)
hold on
line ([Med, Med],[0,max(N)], 'LineWidth', 2,'Color', 'r')
text(Med-10, max(N)+0.05*max(N), sprintf('CBF median = %0.2f', Med))
xlabel ('Frequency (Hz)')
ylabel('probability')
xlim(CBF.caxis); 
ylim([0 max(N)+0.1*max(N)]); 
box off; grid off
title(ax5, 'Histogram of CBF');

% Put a title without the underscore effect
suptitle(insertBefore(CBF.name,'_','\'));


%% Export nicely
saveas(gcf,[CBF.SourceP, CBF.name,'.png'])
close(gcf)
