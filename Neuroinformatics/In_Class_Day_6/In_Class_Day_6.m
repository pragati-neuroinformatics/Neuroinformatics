% 1a.) Inverted U- Kernel 
Inverted_U_kernel= [0.2, 0.4, 0.6, 0.4, 0.2];
subplot(4,1,1)
plot(Inverted_U_kernel ,'.-')
set(gca,'xlim',[0 6],'ylim',[0 1])
xlabel('Time Points');
ylabel('Kernel Value');
title('Inverted U Kernel');

%To Convolve the Inverted_U_ kernels with 50 time-points of EEG data
% from one electrode - use the recording of the Fcz electrode for the first trial from sampleEEGdata.

% Load EEG data
load sampleEEGdata.mat
data = 'sampleEEGdata.mat';

% Specify the label of the channel to plot
which_channel_to_plot = 'fcz';

% Find the index (channel number) of that label
channel_index = find(strcmpi(which_channel_to_plot, {EEG.chanlocs.labels}));
eeg_data = squeeze(EEG.data(channel_index, 1:50, 1));
mean_data = mean(eeg_data);
sd_data = std(eeg_data);
eeg_data_zscored = (eeg_data - mean_data)/sd_data;

% Plot EEG data
subplot(4,1,2)
time_range = [1:50]; % in ms
plot(time_range, eeg_data_zscored);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['EEG Data - Electrode ', which_channel_to_plot, ', Trial 1']);
%Plotting using convolved function (inbuilt)
matlab_conv_result_1 = conv(eeg_data_zscored,Inverted_U_kernel,'same');
% plotting  the result of Inverted_U_convolution (inbuilt)
subplot(4,1,3)
time_range = [1:50]; % in ms
plot(time_range, matlab_conv_result_1);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['convolved-EEG Data - Electrode (done using conv function_using Inverted_U_kernel) ', which_channel_to_plot, ', Trial 1']);

% using manual custom code convolution (for_inverted_U_kernel)
dat4conv = [zeros(1,length(Inverted_U_kernel)-1) eeg_data_zscored zeros(1,length(Inverted_U_kernel)-1) ]; % padding zeroes on both sides which depends on the length of kernel and data points. 

% running convolution
for ti=1:numel(dat4conv) - length(Inverted_U_kernel) + 1
    convolution_result(ti) = sum(dat4conv(ti:ti+length(Inverted_U_kernel)-1).*Inverted_U_kernel); 
end

convolution_result_trimmed = convolution_result(3:52);
% Plotting convolve data
subplot(4,1,4)
time_range = [1:50]; % in ms
plot(time_range, convolution_result_trimmed, 'g');
xlabel('Time (ms)');
ylabel('Amplitude');
title(['convolved-EEG Data - Electrode (done manually_for inverted_U_kernel) ', which_channel_to_plot, ', Trial 1']);


%% Using_ Decay function_kernel 
kernel_decay= [1, 0.6, 0.4, 0.30,.2];
figure;
subplot (4,1,1)
plot(kernel_decay ,'.-')
set(gca,'xlim',[0 6],'ylim',[0 1])
xlabel('Time Points');
ylabel('Kernel Value');
title('Decay Function Kernel');

% 2.) To Convolve the decay_function kernels with 50 time-points of EEG data
% from one electrode - use the recording of the Fcz electrode for the first trial from sampleEEGdata.

% Load EEG data
load sampleEEGdata.mat
data = 'sampleEEGdata.mat';

% Specify the label of the channel to plot
which_channel_to_plot = 'fcz';

% Find the index (channel number) of that label
channel_index = find(strcmpi(which_channel_to_plot, {EEG.chanlocs.labels}));
eeg_data = squeeze(EEG.data(channel_index, 1:50, 1));
mean_data = mean(eeg_data);
sd_data = std(eeg_data);
eeg_data_zscored = (eeg_data - mean_data)/sd_data;

% Plot EEG data
subplot(4,1,2)
time_range = [1:50]; % in ms
plot(time_range, eeg_data_zscored);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['EEG Data - Electrode ', which_channel_to_plot, ', Trial 1']);

%Plotting using convolved function (inbuilt)
matlab_conv_result_2 = conv(eeg_data_zscored,kernel_decay,'same');


% plotting  the result of Decay_funciton_kernel (inbuilt)
subplot(4,1,3)
time_range = [1:50]; % in ms
plot(time_range, matlab_conv_result_2);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['convolved-EEG Data - Electrode (done using conv function_using decay function) ', which_channel_to_plot, ', Trial 1']);

% using manual custom code convolution (using kernel_decay)
dat4conv = [zeros(1,length(kernel_decay)-1) eeg_data_zscored zeros(1,length(kernel_decay)-1) ]; % padding zeroes on both sides which depends on the length of kernel and data points. 

% running convolution
for ti=1:numel(dat4conv) - length(kernel_decay) + 1
    convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel_decay)-1).*kernel_decay); 
end

convolution_result_trimmed = convolution_result(3:52);
% Plotting convolve data
subplot(4,1,4)
time_range = [1:50]; % in ms
plot(time_range, convolution_result_trimmed, 'g');
xlabel('Time (ms)');
ylabel('Amplitude');
title(['convolved-EEG Data - Electrode (done manually_for kernel_decay) ', which_channel_to_plot, ', Trial 1']);

%% 3.) Based on visual inspection, what is the effect of convolving the EEG data with these two kernels? 


%% 4.) To perform convolution using DTFT, FFT and iFFT


% Find the index (channel number) of that label
eeg_data_all_data_points = squeeze(EEG.data(channel_index, 1:640, 1));
mean_data = mean(eeg_data_all_data_points);
sd_data = std(eeg_data_all_data_points);
eeg_data_zscored_all_data_points = (eeg_data_all_data_points - mean_data)/sd_data;

% Plot EEG data
figure;
subplot(6,1,1)
time_range = [1:640]; % in ms
plot(time_range, eeg_data_zscored_all_data_points);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['EEG Data - Electrode- all_time_points', which_channel_to_plot, ', Trial 1']);

N = length(eeg_data_all_data_points);
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform of data
for fi=1:N 
    % create sine wave
    sine_wave = exp(-1i*2*pi*(fi-1).*time);
    % compute dot product between sine wave and data 
    fourier_data(fi) = sum(sine_wave.*eeg_data_all_data_points);
end

N = length(Inverted_U_kernel);
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1

% Fourier transform of inverted kernel
for fj=1:N
    % create sine wave
    sine_wave = exp(-1i*2*pi*(fj-1).*time);
    % compute dot product between sine wave and data 
    fourier_inv_kernel(fj) = sum(sine_wave.*Inverted_U_kernel);
end

product_of_fouriers = 

%% Performing FFT and IFFT 

load sampleEEGdata.mat
data = 'sampleEEGdata.mat';
% extracting_eeg_data_points
which_channel_to_plot = 'fcz';
channel_index = find(strcmpi(which_channel_to_plot, {EEG.chanlocs.labels}));
eeg_data_all_data_points = squeeze(EEG.data(channel_index, 1:640, 1));

Inverted_U_kernel= [0.2 0.4 0.6 0.4 0.2];

kernel_decay= [1, 0.6, 0.4, 0.30, .2];

fft_eeg_data = fft (eeg_data_all_data_points, length (eeg_data_all_data_points) + length (Inverted_U_kernel)-1);
fft_kernel = fft (Inverted_U_kernel,length (eeg_data_all_data_points) + length (Inverted_U_kernel)-1); 
fft_product= fft_eeg_data.* fft_kernel;
conv_inverse_output = ifft(fft_product);
convo_inverse_output = conv_inverse_output(3:642); 

tic
for k=1:1000
    fft_eeg_data = fft (eeg_data_all_data_points, length (eeg_data_all_data_points) + length (Inverted_U_kernel)-1);
    fft_kernel = fft (Inverted_U_kernel,length (eeg_data_all_data_points) + length (Inverted_U_kernel)-1); 
    fft_product= fft_eeg_data.* fft_kernel;
    conv_inverse_output = ifft(fft_product);
    convo_inverse_output = conv_inverse_output(3:642); 

end
toc

auto_conv = conv(eeg_data_all_data_points, Inverted_U_kernel, 'same');

figure;
plot(convo_inverse_output)
hold on 
plot(auto_conv)
legend
