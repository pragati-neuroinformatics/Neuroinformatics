% 1.) To create a family of Morlet wavelets ranging in frequency from 2 Hz to 30 Hz.
load sampleEEGdata.mat

fs = EEG.srate; 
t = -1:1/fs:1; 
frequencies = [2, 9, 16, 23, 30]; %Frequencies from 2:30 

% Plotting the morlet wavelet 
figure;
for i = 1:length(frequencies)
    f0 = frequencies(i);  % Current frequency
    
    % Morlet wavelet equation
    n_cycles = 6;
    s = n_cycles/(2*pi*f0);
    morlet_wavelet(i,:) = (1/(sqrt(s*sqrt(pi)))) * (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / (2*s^2)));
    % morlet_wavelet(i,:) = (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / 2)); % old construction w/o std devn

    % Plot the Morlet wavelet
    subplot(length(frequencies), 1, i);
    plot(t, real(morlet_wavelet(i,:)), 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    title(['Morlet Wavelet - Frequency ', num2str(f0)]);
    grid on;
end

%% Convolve each Morlet wavelet with EEG data from all trials from the electrode 

data = 'sampleEEGdata.mat';

channel = 'fcz';

channel_index = find(strcmpi(channel, {EEG.chanlocs.labels}));

eeg_data_all_trials = squeeze(EEG.data(channel_index, 1:640, :))';


% Convolve each frequency's data
for i = 1:length(frequencies) % Loop through freqs
    for l = 1:size(eeg_data_all_trials, 1) % Loop through trials
        convol_result(l, :) = real(conv(eeg_data_all_trials(l, :), morlet_wavelet(i, :), 'same'));
    end
    convol_results{i} = convol_result;
end

sample_convolve_data = convol_results{4}(66,:);  % 23-Hz data chosen 

% Y = fft(sample_convolve_data);
% Z = abs (Y);
% figure;
% plot(Z);% the data has been band pass filtered at ~15-26 Hz with peak pass at 20Hz

n = length(sample_convolve_data);  % Number of samples
freq_fourier = linspace(0, EEG.srate/2, n/2 + 1);  % Frequency axis
fft_out = fft(sample_convolve_data);
magnitude_spectrum = abs(fft_out(1:n/2 + 1));
figure
plot(freq_fourier, magnitude_spectrum);
%% Averaging the result(s) of convolution over all trials and  and plotting  an ERP corresponding to each wavelet frequency.

% 
% for j = 1:5
%     for l = 1:size(eeg_data_all_trials,1)
%          same_freq_store{j}(l,:) = convol_result{l}(1,:);
%     end
% 
%     mean_result = mean(same_freq_store{j},1);
%     figure;
%     plot(real(mean_result));
%     disp(mean_result);
% 
% end


% Compute the broadband ERP (without convolution)
broadband_ERP = mean(eeg_data_all_trials, 1);

% Plotting
% figure;


mean_results = cell(1, length(frequencies));
% Plot each frequency's ERP in its own subplot
for i = 1:length(frequencies)
    subplot(length(frequencies)+1, 1, i);
    mean_result = mean(convol_results{i}, 1);
    mean_results{i} = mean_result;
    plot(EEG.times(193:513),mean_result(193:513));
    title(['ERP for frequency: ' num2str(frequencies(i)) ' Hz']);
    xlabel('Time');
    ylabel('Amplitude');
end

% Plotting the broadband ERP
subplot(length(frequencies)+1, 1, length(frequencies)+1);
plot(EEG.times(193:513),broadband_ERP(193:513));
title('Broadband ERP (without convolution)');
xlabel('Time');
ylabel('Amplitude');

% Reasons when we compare the wavelet convolved ERP's with Broadband comparison
%Wavelet_convolved ERPs exhibit activity at specific frequency and we can
%observe that there is some brain activity going at around 300ms for
%9,16,23 and 30 Hz, whereas in the broadband ERP without convolution it 
%provides a summary of the overall brain acitivity without any specific
% frequencies role in the brain acivity but overall we can say there is
% some brain activity around 300ms.
