% 1.) To create a family of Morlet wavelets ranging in frequency from 2 Hz to 30 Hz.
fs = 1000; 
t = -2:1/fs:2; 
frequencies = [2, 9, 16, 23, 30]; %Frequencies from 2:30 

% Plotting the morlet wavelet 
figure;
for i = 1:length(frequencies)
    f0 = frequencies(i);  % Current frequency
    
    % Morlet wavelet equation
   % morlet_wavelet(i,:) = (pi^(-1/4)) .* exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / 2);
    morlet_wavelet(i,:) = (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / 2));
    % Plot the Morlet wavelet
    subplot(length(frequencies), 1, i);
    plot(t, real(morlet_wavelet(i,:)), 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    title(['Morlet Wavelet - Frequency ', num2str(f0)]);
    grid on;
end

%% Convolve each Morlet wavelet with EEG data from all trials from the electrode 

load sampleEEGdata.mat
data = 'sampleEEGdata.mat';

channel = 'fcz';

channel_index = find(strcmpi(channel, {EEG.chanlocs.labels}));

eeg_data_all_trials = squeeze(EEG.data(channel_index, 1:640, :))';


for i = 1:length(frequencies)

   for j = 1:size(eeg_data_all_trials,1)
        eeg_data = eeg_data_all_trials(j,:);

        convol_result{j}(i,:) = real(conv(eeg_data, morlet_wavelet(i,:) , 'same'));
            
    end
    
end

sample_convolve_data = convol_result{66}(5,:);  % 16-Hz data chosen 
Y = fft(sample_convolve_data);
Z = abs (Y);
figure;
plot(Z);% the data has been band pass filtered at ~15-26 Hz with peak pass at 20Hz
