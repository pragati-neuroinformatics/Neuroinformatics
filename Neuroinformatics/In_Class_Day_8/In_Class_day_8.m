% 1.) To create a family of Morlet wavelets ranging in frequency from 2 Hz to 30 Hz.
clear;
choice_trial = 88;
load sampleEEGdata.mat

fs = EEG.srate; 
t = -1:1/fs:1; 
f_min = 2;
f_max = 30;
steps = 5;
frequencies = linspace(f_min,f_max,steps); %linspace
% frequencies = logspace(log10(f_min),log10(f_max),steps); %logspaced
%%
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
    % subplot(length(frequencies), 1, i);
    plot(t, real(morlet_wavelet(i,:)), 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    % title(['Morlet Wavelet - Frequency ', num2str(f0)]);
    grid on;
    hold on;
    legend;
end
hold off

%% Convolve each Morlet wavelet with EEG data from all trials from the electrode 

data = 'sampleEEGdata.mat';

eeg_data_all_electrodes = squeeze(EEG.data(:, :, choice_trial));

for i = 1:length(frequencies)
   for j = 1:size(eeg_data_all_electrodes,1)
        eeg_data = eeg_data_all_electrodes(j,:);
        convol_result{j}(i,:) = conv(eeg_data, morlet_wavelet(i,:) , 'same');           
   end    
end

%% topographical plots of power and phase at (a) 180 ms and (b) 360 ms at all the 5 frequencies.

power_phase_matrix = zeros(EEG.pnts, length(frequencies), size(eeg_data_all_electrodes,1), 2);
for i = 1:length(frequencies)
    for electrode = 1:size(eeg_data_all_electrodes,1)
        convolution_result = squeeze(convol_result{electrode}(i,:));
        % Extract power and phase
        power_phase_matrix(:, i, electrode, 1) = abs(convolution_result) .* abs(convolution_result) ; % Power
        power_phase_matrix(:, i, electrode, 2) = angle(convolution_result); % Phase
    end
end
time_points = [180,360];

for k = 1:length(time_points)   

    t = time_points(k); 
    times_delta = abs(EEG.times - t);
    [~, index] = min(times_delta);
    
    for i = 1:length(frequencies)
        figure;
        
        subplot(2, 1, 1);
        topoplot(squeeze(power_phase_matrix(index, i, :, 1)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'map');
        title(['Power at ' num2str(frequencies(i)) 'Hz, t = ' num2str(t) 'ms']);
        
        subplot(2, 1, 2);
        topoplot(squeeze(power_phase_matrix(index, i, :, 2)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'map');
        title(['Phase at ' num2str(frequencies(i)) 'Hz, t = ' num2str(t) 'ms']);
    end
end


%%
%If we change the number of cycles in our  construction of wavelets 

% load sampleEEGdata.mat
% % Morlet wavelet equation
% for i = 1:length(frequencies)
%     f0 = frequencies(i);  % Current frequency
%     
%     % Morlet wavelet equation
%     n_cycles_change = [6, 30];
%     for g = 1: length (n_cycles_change)
%       w= n_cycles_change(g);
%       s = w/(2*pi*f0);
%       morlet_wavelet{(i,:) = (1/(sqrt(s*sqrt(pi)))) * (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / (2*s^2)));
%     end
% end

 
% for i = 1:length(frequencies)
%     for electrode = 1:size(eeg_data_all_electrodes,1)
%         convolution_result = squeeze(convol_result{electrode}(i,:));
%         % Extract power and phase
%         power_phase_matrix(:, i, electrode, 1) = abs(convolution_result) .* abs(convolution_result) ; % Power
% 
%  power_phase_matrix(:, i, electrode, 1)
% % 
% 
% 
% 
% 
% t2 = 360;
% times_delta = abs(EEG.times - t2);
% [~, index] = min(times_delta);
% 
% for i = 1:length(frequencies)
%     figure;
%     
%     subplot(2, 1, 1);
%     topoplot(squeeze(power_phase_matrix(index, i, :, 1)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'map');
%     title(['Power at ' num2str(frequencies(i)) 'Hz, t=360ms']);
%     
%     % Phase at 180ms
%     subplot(2, 1, 2);
%     topoplot(squeeze(power_phase_matrix(index, i, :, 2)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'map');
%     title(['Phase at ' num2str(frequencies(i)) 'Hz, t=360ms']);
% end
% %% If we change the number of cycles in our  construction of wavelets 
% 
% load sampleEEGdata.mat
% 
% % fs = EEG.srate; 
% % t = -1:1/fs:1; 
% % frequencies = [2, 9, 16, 23, 30]; %Frequencies from 2:30 
% % 
% % % Plotting the morlet wavelet 
% % figure;
% % for i = 1:length(frequencies)
% %     f0 = frequencies(i);  % Current frequency
% %     
% %     % Morlet wavelet equation
%     n_cycles = 30;% using of 3 Hz wavelet 
%     s = n_cycles/(2*pi*f0);
%     morlet_wavelet(i,:) = (1/(sqrt(s*sqrt(pi)))) * (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / (2*s^2)));
%     % morlet_wavelet(i,:) = (exp(1i * 2 * pi * f0 .* t) .* exp(-t.^2 / 2)); % old construction w/o std devn
% 
%     % Plot the Morlet wavelet
%     subplot(length(frequencies), 1, i);
%     plot(t, real(morlet_wavelet(i,:)), 'LineWidth', 1.5);
%     xlabel('Time');
%     ylabel('Amplitude');
%     title(['Morlet Wavelet - Frequency ', num2str(f0)]);
%     grid on;
% end
% %%
% 
% data = 'sampleEEGdata.mat';
% 
% channel = 'fcz';
% 
% channel_index = find(strcmpi(channel, {EEG.chanlocs.labels}));
% 
% eeg_data_all_trials = squeeze(EEG.data(channel_index, 1:640, :))';
% 
% 
% for i = 1:length(frequencies)
% 
%    for j = 1:size(eeg_data_all_trials,1)
%         eeg_data = eeg_data_all_trials(j,:);
% 
%         convol_result{j}(i,:) = real(conv(eeg_data, morlet_wavelet(i,:) , 'same'));
%             
%     end
%     
% end
% 
% sample_convolve_data = convol_result{66}(4,:);  % 23-Hz data chosen 
% 
% % Y = fft(sample_convolve_data);
% % Z = abs (Y);
% % figure;
% % plot(Z);% the data has been band pass filtered at ~15-26 Hz with peak pass at 20Hz
% 
% n = length(sample_convolve_data);  % Number of samples
% freq_fourier = linspace(0, EEG.srate/2, n/2 + 1);  % Frequency axis
% fft_out = fft(sample_convolve_data);
% magnitude_spectrum = abs(fft_out(1:n/2 + 1));
% figure
% plot(freq_fourier, magnitude_spectrum);
% %%
% 
% data = 'sampleEEGdata.mat';
% 
% channel = 'fcz';
% 
% channel_index = 64;
% 
% eeg_data_all_trials = squeeze(EEG.data(channel_index, 1:640, :))';
% 
% power_phase_matrix = zeros(EEG.pnts, length(frequencies), channel_index, 2);
% for i = 1:length(frequencies)
%     for electrode = 1:channel_index
%         convolution_result = squeeze(convol_result{electrode}(i,:));
%         
%         % Extract power and phase
%         power_phase_matrix(:, i, electrode, 1) = abs(convolution_result) .* abs(convolution_result) ; % Power
%         power_phase_matrix(:, i, electrode, 2) = angle(convolution_result); % Phase
%     end
% end
% 
% 
% % t = 180ms
% tm_pnt = nearest(180 * EEG.srate / 1000);
% 
% 
% for i = 1:length(frequencies)
%     figure;
%     
%     subplot(2, 1, 1);
%     topoplot(squeeze(power_phase_matrix(tm_pnt, i, :, 1)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'fill', 'numcontour', 8);
%     title(['Power at ' num2str(frequencies(i)) 'Hz, t=180ms']);
%     
%     subplot(2, 1, 2);
%     topoplot(squeeze(power_phase_matrix(tm_pnt, i, :, 2)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'fill', 'numcontour', 8);
%     title(['Phase at ' num2str(frequencies(i)) 'Hz, t=180ms']);
% end
% 
% 
% % t = 360ms
% tm_pnt = nearest(360 * EEG.srate / 1000);
% 
% 
% for i = 1:length(frequencies)
%     figure;
%     
%     subplot(2, 1, 1);
%     topoplot(squeeze(power_phase_matrix(tm_pnt, i, :, 1)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'fill', 'numcontour', 8);
%     title(['Power at ' num2str(frequencies(i)) 'Hz, t=360ms']);
%     
%     % Phase at 360ms 
%     subplot(2, 1, 2);
%     topoplot(squeeze(power_phase_matrix(tm_pnt, i, :, 2)), EEG.chanlocs, 'electrodes', 'off', 'colormap', 'jet', 'style', 'fill', 'numcontour', 8);
%     title(['Phase at ' num2str(frequencies(i)) 'Hz, t=360ms']);
% end
