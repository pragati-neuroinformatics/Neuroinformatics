% Define the frequencies of interest
freqs_of_interest = linspace(2, 30, 8);

% Pre-allocate matrices for power and z-score normalized signals
power_at_freq_all = zeros(length(freqs_of_interest), size(hilbert_transformed_data, 2), size(hilbert_transformed_data, 3));
zScoreNormalizedSignals = zeros(size(power_at_freq_all));

% Loop through each frequency of interest
for f = 1:length(freqs_of_interest)
    % Find the index of the frequency of interest
    freq_idx = round(freqs_of_interest(f) / (EEG.srate / EEG.pnts));
    
    % Extract the power at the frequency of interest for each time point across all trials
    power_at_freq_all(f, :, :) = abs(hilbert_transformed_data(:, freq_idx, :)).^2;
    
    % Compute mean and standard deviation
    mu = mean(power_at_freq_all(f, :, :), 'all');
    sigma = std(power_at_freq_all(f, :, :), 0, 'all');
    
    % Apply z-score normalization
    zScoreNormalizedSignals(f, :, :) = (power_at_freq_all(f, :, :) - mu) / sigma;
end

% Plotting
figure;
for f = 1:length(freqs_of_interest)
    subplot(length(freqs_of_interest), 1, f)
    plot(EEG.times, squeeze(zScoreNormalizedSignals(f, :, :)));
    xlabel('Time (ms)')
    ylabel('Z-score normalized amplitude')
    title(['Z-score Normalized Trial-Averaged Signal at Electrode FCz - Frequency of Interest: ' num2str(freqs_of_interest(f)) ' Hz'])
end
