%% Master Analysis file
% All analysis performed on light transmission recordings of ciliary
% beating in the zebrafish adult telencephalon

%% Step 0: If you need to align your recording go to Step 1, otherwise proceed to Step 2 

%% Step 1: Align recordings if there is a bit of x,y drift during the recording
% the alignement code will not allow to correct large drifts

% Define the path to the raw recording
[CBF.name,CBF.SourceP] = uigetfile;
cd(CBF.SourceP); 

% load the raw recording
load([CBF.SourceP, CBF.name]);

% Input the frequency of acquisition
CBF.Fs = 102.995; % Frequency of acquisition; you find it in place 8 of the metadata, no need for all the digits

CBF.aligned = 'true';

% Input the amount of data to be loaded
CBF.duration=10; %in seconds
data = double(data(:,:,1:CBF.duration*CBF.Fs)); %only take CBF.duration sec.

% Align the data 
[aligned] = align_stacks_simple(data); % Sometimes aligning is suboptimal

% Transform in 8bit to reduce datasize and save it as a variable called
% data
clear data
data=uint8(aligned); 

% Save the results
save([CBF.SourceP, CBF.name, '_aligned.mat'], 'data', '-v7.3');
save([CBF.SourceP, CBF.name, '_aligned.mat.CBF_parameters.mat'], 'CBF');


%% Step 2: Import an aligned or raw recording and set some variables

% Define the path to the raw or aligned recording
[CBF.name,CBF.SourceP] = uigetfile;

% Load the file
load([CBF.SourceP, CBF.name]);

% Load the aligned file or a previously analysed file, and transform the data in double
% These are the only cases where a CBF_parameters file exists

if exist([CBF.SourceP, CBF.name, '.CBF_parameters.mat'])
   load([CBF.SourceP, CBF.name, '.CBF_parameters.mat']);
   data = double(data);
   
% If you load a new file that was not aligned, you should enter some extra information such as the frequency of aquisition and
% the duration of the file to be analysed in sec
else 
    prompt = 'What is the CBF (enter the value and press enter)? ';
    CBF.Fs = input(prompt);
    prompt = 'What duration of the recording would you like to analyse in sec (enter the value and press enter)? ';
    CBF.duration = input(prompt);
    CBF.aligned = 'false';
    data=double(data);
end

%% Step 3: Defining variables
 
% Define variables and save them in one common structure
CBF.x = size(data,1); 
CBF.y = size(data,2); 
CBF.w_min = 10; % Lower frequency cutoff 
CBF.caxis = [10 40]; % Upper and lower bound for any frequency plot. 
CBF.spatres = 0.314;  % [um/pixel] spatial resolution, depends on your microscope

% Save the results
save([CBF.SourceP, CBF.name, '.CBF_parameters.mat'], 'CBF');

%% Step 4: Perform the fast Fourier transform 

% Fourier Transform (~20s to run the fft & ~200s with plotting)
[PowerSpec,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe,PeakPos] = fft_analysis_simple(data,CBF);

% Save the results
save([CBF.SourceP, CBF.name, '.CBF_parameters.mat'], 'CBF');

% in the results, you will find
% CBF.picSD is the masked CFP map
% CBF.mask is the mask
% you can use this matrix to plot the masked CBF map
% Note that the Power spectrum is not saved in the CBF structure because it
% is too big. If you wish to look at the power spectrum of all individual
% pixel, this is saved in the PowerSpec matrix. Each data point in the 3rd
% dimension corresponds to one frequency. The range of frequency is given by (1:size(PowerSpec,3))*(CBF.Fs/CBF.nframe)
