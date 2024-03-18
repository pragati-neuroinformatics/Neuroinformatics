% Performing wavelet convolution 
% Load EEG data 
load sampleEEGdata

% Select electrode
chan2use = 'fcz';

% Frequency of interest
freq_of_interest = 10; % For example, 10 Hz

% Wavelet parameters
time = -1:1/EEG.srate:1;
s = 5 / (2*pi*freq_of_interest); % "5" for the number of wavelet cycles

% Wavelet and data sizes
n_wavelet = length(time);
n_data = EEG.pnts * EEG.trials;
n_convolution = n_wavelet + n_data - 1;
n_conv_pow2 = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet - 1) / 2;

% To extract FFT of the data for the selected channel
chanidx = strcmpi(chan2use, {EEG.chanlocs.labels});
eegfft = fft(reshape(EEG.data(chanidx,:,:), 1, EEG.pnts*EEG.trials), n_conv_pow2);

% Create the wavelet
wavelet = fft(sqrt(1/(s*sqrt(pi))) * exp(2*1i*pi*freq_of_interest.*time) .* exp(-time.^2./(2*s^2)), n_conv_pow2);

% Perform convolution
eegconv = ifft(wavelet.*eegfft);
eegconv = eegconv(1:n_convolution);
eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);

% Reshape and compute power
eegpower_time = abs(reshape(eegconv, EEG.pnts, EEG.trials)).^2;
mean_power = mean(eegpower_time, 2);

mu = mean(mean_power);
sigma = std(mean_power);

% Apply z-score normalization
zScoreNormalizedSignal = (mean_power - mu) / sigma;
figure
subplot(311)
plot(EEG.times, eegpower_time);
xlabel('Time (ms)');
ylabel('Power (dB)');
title(sprintf('Power over time at %s electrode, %d Hz', chan2use, freq_of_interest));
xlim([EEG.times(1), EEG.times(end)]);

subplot(312)
plot(EEG.times, mean_power);
xlabel('Time (ms)');
ylabel('Power (uV^2)'); % Adjust the unit based on your data
title(sprintf('Trial-Averaged Power at %d Hz for Electrode %s', freq_of_interest, chan2use));
xlim([EEG.times(1), EEG.times(end)]); % Adjust the x-axis limits to the EEG time range


subplot(313)
plot(EEG.times, zScoreNormalizedSignal)
xlabel('Time (ms)')
ylabel('Z-score normalized amplitude')
title(sprintf('Z-score Normalized Trial-Averaged Signal for Wavelet Convolution at Electrode FCz'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing Hilbert Transform 
% Load sample EEG data
load sampleEEGdata

% Define the electrode of interest
chan2use = 'FCz'; % For example, FCz
chanidx = strcmpi({EEG.chanlocs.labels}, chan2use);

% Filter specifications
nyquist = EEG.srate / 2;
lower_filter_bound = 4; % Lower frequency bound in Hz
upper_filter_bound = 10; % Upper frequency bound in Hz
transition_width = 0.2; % Transition width
filter_order = round(3 * (EEG.srate / lower_filter_bound));

% Design the filter
ffrequencies = [0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist] / nyquist;
idealresponse = [0 0 1 1 0 0];
filterweights = firls(filter_order, ffrequencies, idealresponse);

% Pre-allocate the matrix for filtered data
filtered_data = zeros(size(EEG.data, 1), EEG.pnts, EEG.trials);

% Apply the filter to the data
for trial = 1:EEG.trials
    for chani = 1:EEG.nbchan
        filtered_data(chani, :, trial) = filtfilt(filterweights, 1, double(EEG.data(chani, :, trial)));
    end
end

% Compute the Hilbert transform for the filtered data at the specified electrode
hilbert_transformed_data = hilbert(squeeze(filtered_data(chanidx, :, :)).').';

% Define the frequency of interest
freq_of_interest = 10; % For example, 10 Hz

% Find the index of the frequency of interest
freq_idx = round(freq_of_interest / (EEG.srate / EEG.pnts));

% Extract the power at the frequency of interest for each time point across all trials
power_at_freq = abs(hilbert_transformed_data).^2;

mu = mean(power_at_freq);
sigma = std(power_at_freq);

% Apply z-score normalization
zScoreNormalizedSignal2 = (power_at_freq - mu) / sigma;


figure;
subplot(311)
plot(EEG.times, power_at_freq(:,:,1)); 
xlabel('Time (ms)');
ylabel('Power (\muV^2)');
title(['Power over time at ' chan2use ' electrode Frequency of interest : ' num2str(freq_of_interest) ' Hz ']);

subplot(312)
plot(EEG.times, mean(power_at_freq(:,:,1),2)); 
xlabel('Time (ms)');
ylabel('Power (\muV^2)');
title(['Trial Averaged Power over time at ' chan2use ' electrode Frequency of interest : ' num2str(freq_of_interest) ' Hz ']);


subplot(313)
plot(EEG.times, zScoreNormalizedSignal2)
xlabel('Time (ms)')
ylabel('Z-score normalized amplitude')
title(sprintf('Z-score Normalized Trial-Averaged Signal at Electrode FCz'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing STFT

% Load EEG data
load sampleEEGdata.mat

% Select electrode of interest
electrode = 'FCz';
chanIdx = find(strcmpi({EEG.chanlocs.labels}, electrode));

% Frequency of interest
targetFrequency = 10; % in Hz

% Parameters for STFT
windowSize = round(EEG.srate / 5); % Window size for the spectrogram (e.g., 0.5 second)
noverlap = round(windowSize * 0.9); % Overlap between windows (e.g., 90% overlap)
nfft = 2^nextpow2(windowSize); % Number of FFT points

% Pre-allocate for efficiency
powerOverTime = zeros(EEG.pnts, EEG.trials);

figure;
subplot(311);
% Loop through each trial
for trial = 1:EEG.trials
    % Extract data for current trial at the selected electrode
    data = EEG.data(chanIdx,:,trial);
    
    % Use the spectrogram function to compute the STFT and power spectrum
    [s,f,t,p] = spectrogram(data, windowSize, noverlap, nfft, EEG.srate);
    
    % Find the closest frequency index to the target frequency
    [~, freqIndex] = min(abs(f - targetFrequency));
    
    % Extract power at the target frequency and interpolate to match EEG time points
    powerAtFreq = p(freqIndex, :);
    powerInterpolated = interp1(t, powerAtFreq, linspace(t(1), t(end), EEG.pnts), 'pchip');
    
    % Store the interpolated power for the current trial
    powerOverTime(:, trial) = powerInterpolated;

    plot(EEG.times,powerOverTime);
    xlabel('Time (ms)');
    ylabel('Power (\muV^2)');
    title([' Power over Time for all trials at ' electrode ' Electrode, ' num2str(targetFrequency) ' Hz Frequency']);
    hold on;
end

% Calculate the average power over all trials
meanPowerOverTime = mean(powerOverTime, 2);

% Plotting
% figure;
subplot(312);
plot(EEG.times, meanPowerOverTime);
xlabel('Time (ms)');
ylabel('Power (\muV^2)');
title(['Average Power over Time at ' electrode ' Electrode, ' num2str(targetFrequency) ' Hz Frequency']);

mu = mean(meanPowerOverTime);
sigma = std(meanPowerOverTime);

% Apply z-score normalization
zScoreNormalizedSignal3 = (meanPowerOverTime - mu) / sigma;

subplot(313)
plot(EEG.times, zScoreNormalizedSignal3)
xlabel('Time (ms)')
ylabel('Z-score normalized amplitude')
title(sprintf('Z-score Normalized Trial-Averaged Signal at Electrode %s', EEG.chanlocs(chanIdx).labels))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing Multi-taper method 

% Load sample EEG data
load sampleEEGdata

% Define the electrode of interest
channel2plot = 'FCz'; % For example, FCz

% Define the frequency of interest
frequency2plot = 10; % in Hz

% Define time window
timewin = 400; % in ms

% Define time points
times2save = -1000:50:1500;

% Convert time points to indices
times2saveidx = dsearchn(EEG.times', times2save');

% Convert time window to indices
timewinidx = round(timewin / (1000 / EEG.srate));

% Define multitaper parameters
nw_product = 3; % Determines the frequency smoothing, given a specified time window
tapers = dpss(timewinidx, nw_product);
f = linspace(0, EEG.srate / 2, floor(timewinidx / 2) + 1);

% Find channel index
chanidx = strcmpi(channel2plot, {EEG.chanlocs.labels});

% Initialize output matrix
multitaper_tf = zeros(floor(timewinidx / 2) + 1, length(times2save));

% Loop through time bins
for ti = 1:length(times2saveidx)
    
    % Initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx / 2) + 1, 1);
    
    % Loop through tapers
    for tapi = 1:size(tapers, 2)
        
        % Window and taper data, and get power spectrum
        % Extract EEG data within the time window for all trials
        time_window = times2saveidx(ti) - floor(timewinidx / 2) + 1 : times2saveidx(ti) + ceil(timewinidx / 2);
        time_window = max(1, time_window); % Ensure time window indices are within range
        time_window = min(time_window, EEG.pnts); % Ensure time window indices are within range
        
        % Extract EEG data within the time window for the specified channel and all trials
        data = EEG.data(chanidx, time_window, :); % No need for squeeze
        
        % Replicate tapers along the third dimension to match the number of trials
        tapers_replicated = repmat(tapers(:, tapi), [1 1 size(data, 3)]);
        
        % Perform element-wise multiplication
        data_tapered = bsxfun(@times, data, tapers_replicated);
        
        % Compute power spectrum
        pow = fft(data_tapered, timewinidx) / timewinidx;
        pow = pow(1:floor(timewinidx / 2) + 1, :, :);
        
        % Average power across trials and tapers
        taperpow = taperpow + mean(mean(abs(pow).^2, 3), 2);
    end
    
    % Store multitaper power for this time point
    multitaper_tf(:, ti) = taperpow / size(tapers, 2);
end

% Plot full TF map
figure
subplot(211)
plot(times2save, multitaper_tf)
title(['Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz'])
xlabel('Time (ms)'), ylabel('Power')
axis square
xlim([times2save(1) times2save(end)])


subplot(212)
plot(times2save, mean(multitaper_tf,1))
title(['Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz'])
xlabel('Time (ms)'), ylabel('Power')
axis square
xlim([times2save(1) times2save(end)])


figure;
plot(EEG.times, zScoreNormalizedSignal);
hold on;
plot(EEG.times, zScoreNormalizedSignal2);
hold on;
plot(EEG.times, zScoreNormalizedSignal3)
xlabel('Time (ms)')
ylabel('Z-score normalized amplitude')
title(sprintf('Z-score Normalized Trial-Averaged Signal at Electrode %s', EEG.chanlocs(chanIdx).labels))

legend('Wavelet Convolution' , 'Hilbert Filter','STFT');

%% How visually similar are the results from these three methods? 
% If the results from the three methods are different, how are they different,
% and what parameters do you think you could change in the three methods to make the results look more or less similar?
% The three methods for computing power used here are  wavelet convolution,
% filter-hilbert, Short term Fourier Transform (STFT) and Multitaper
% method. 
% Wavelet convolution often provides better control over the trade-off between time and 
% frequency resolution due to its multiplicative scaling of the wavelet with frequency.
% STFT has a fixed window size, leading to a constant time-frequency resolution across all frequencies. 
% Hilbert transform is primarily used to extract the amplitude envelope and instantaneous phase of signals within a specific frequency band, 
% filtered prior to the transform, and doesn't inherently provide a time-frequency representation.
% The Hilbert transform and wavelet convolution might handle edge artifacts differently compared to STFT, 
% depending on how the data is preprocessed and how the convolution is implemented.
% Here, in the Figure 1 and 3 we can observe after (z-scoring) wavelet convolution
% and STFT shows the similar trend. But, in case of Figure 2 we are getting
% a different plot which might be due to the extraction of amplitude
% envelope and its instantaneous phase of signals.  Also, the z scored
% normalized Trial averaged signal shows almost similar plot among the 3
% metric (Wavelet convolution, hilbert transform and the STFT) after
% changing the parameters such as number of cycles (wavelet convolution),
% filter order (Hilbert transform), and modifying the window length (STFT).
% Parameters to Adjust:
% Wavelet Convolution: Adjusting the number of cycles in the wavelet can make the results more similar or different by changing the balance between time and frequency resolution.
% Hilbert Transform: Changing the filter parameters (e.g., filter order, cut-off frequencies) before applying the Hilbert transform can alter the frequency band of interest and potentially affect the temporal resolution of the extracted amplitude envelope.
% STFT: Modifying the window length and overlap in STFT can impact the time-frequency resolution, making the results appear more similar to wavelet-based analyses by optimizing the trade-off between temporal and spectral resolution.
