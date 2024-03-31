%1 a) Pick 3 electrodes and 5 time points and plot the distribution of phases 
% across trials
load sampleEEGdata
electrodes = [47, 57, 62]; %choosing FCz (47), Pz (57) and PO8 (62)
% center frequency
centerfreq = 12; % in Hz
times2plot = [ -200 150 360 600 800 ]; % in ms, from stimulus onset

% definte convolution parameters
n_wavelet_length     = EEG.pnts;
n_data_length        = EEG.pnts*EEG.trials;
n_convolution_length = n_wavelet_length+n_data_length-1;
n_conv_pow2   = pow2(nextpow2(n_convolution_length));

time    = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;
for k =1:3
    eegfft = fft(reshape(EEG.data(electrodes(k),:,:),1,[]),n_conv_pow2);
    
    % convolution
    eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
    eegcon = eegconv(1:n_convolution_length);
    eegconv = reshape(eegcon(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);
    
    % plot
    figure
    for subploti=1:5
        subplot(5,2,subploti*2 - 1)
        [junk,idx]=min(abs(EEG.times-times2plot(subploti)));
        polar([zeros(1,EEG.trials); angle(eegconv(idx,:))],[zeros(1,EEG.trials); ones(1,EEG.trials)])
        title([ 'Electrode number '  num2str(electrodes(k))  ' ITPC at ' num2str(times2plot(subploti)) ' ms = ' num2str(round(1000*abs(mean(exp(1i*angle(eegconv(idx,:))))))/1000) ])
        
        subplot(5,2,subploti*2)
        hist(angle(eegconv(idx,:)),20)
        set(gca,'xlim',[-pi pi]*1.1,'xtick',-pi:pi/2:pi,'xticklabel',{'-pi';'-pi/2';'0';'pi';'pi/2'},'ylim',[0 20])
    end
end
%%  Time-frequency plots of ITPC and decibel-baselined power for these electrodes,
% using the filter-Hilbert method.

load sampleEEGdata.mat
nyquist = EEG.srate / 2;
electrode_indx = [47 57 62];

%  frequencies of interest using logarithmic spacing
frequencies = logspace(log10(4), log10(40), 20);

s = logspace(log10(3), log10(10), length(frequencies)) ./ (2 * pi * frequencies);

% Defining baseline time period
baselinetime = [-300 -100];

baseidx = dsearchn(EEG.times', baselinetime(1)):dsearchn(EEG.times', baselinetime(2));

% Initialize matrices to store ITPC and power values for all electrodes
% (preallocation)
itpc_matrix_all = zeros(length(electrode_indx), length(frequencies), EEG.pnts);
power_matrix_all = zeros(length(electrode_indx), length(frequencies), EEG.pnts);

for elecInd = 1:3
    itpc_matrix = zeros(length(frequencies), EEG.pnts);
    power_matrix = zeros(length(frequencies), EEG.pnts);
    
    % Loop over each frequency
    for fi = 1:length(frequencies)
        % Create filter for current frequency
        lower_filter_bound = frequencies(fi) - 0.5 * s(fi);
        upper_filter_bound = frequencies(fi) + 0.5 * s(fi);
        transition_width = 0.2;
        filter_order = round(3 * (EEG.srate / lower_filter_bound));
        ffrequencies = [0 (1 - transition_width) * lower_filter_bound lower_filter_bound upper_filter_bound (1 + transition_width) * upper_filter_bound nyquist] / nyquist;
        idealresponse = [0 0 1 1 0 0];
        filterweights = firls(filter_order, ffrequencies, idealresponse);
             
        filtered_data = filtfilt(filterweights, 1, double(squeeze(EEG.data(electrode_indx(elecInd), :, :))));
        
        % Computing Hilbert transform
        hilbert_data = hilbert(filtered_data);
        
        % Compute ITPC and power for current frequency
        itpc_matrix(fi, :) = abs(mean(exp(1i * angle(hilbert_data)), 2));
        power_matrix(fi, :) = mean(abs(hilbert_data).^2, 2);
        power_matrix(fi, :) = 10*log10(power_matrix(fi, :) ./ mean(power_matrix(fi, baseidx), 2));
    end
    
    % Store ITPC and power matrices for the current electrode
    itpc_matrix_all(elecInd, :, :) = itpc_matrix(:, :);
    power_matrix_all(elecInd, :, :) = power_matrix(:,:);
end

% Plot ITPC and power matrices for each electrode
baselinetime = [-300 -100];
baseidx = dsearchn(EEG.times', baselinetime(1)):dsearchn(EEG.times', baselinetime(2));
for figInd = 1:length(electrode_indx)
    % Plotting ITPC
    figure;
    subplot(2, 1, 1)
    contourf(EEG.times, frequencies, squeeze(itpc_matrix_all(figInd, :, :)), 40, 'linecolor', 'none')
    set(gca, 'clim', [0.1 0.5], 'xlim', [-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('ITPC')
    
    % Plotting power
    subplot(2, 1, 2)
    contourf(EEG.times, frequencies, squeeze(power_matrix_all(figInd, :, :)), 40, 'linecolor', 'none')
    set(gca, 'clim', [-1 1], 'xlim', [-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('DB-power')
end

% The patterns of results from ITPC and power show differences rather than similarities. ITPC measures 
% the consistency of phase across trials, indicating temporal precision, while power reflects the 
% magnitude of oscillatory activity. Variations in phase coherence may not correspond to variations in power, 
% leading to distinct patterns between the two measures.

%%  compute wITPCz using reaction time as the trial-varying modulator. 

load sampleEEGdata.mat
centerfreq  = 6;
channel2use = 'po7';
times2save  = -200:50:1200;

% initialize matrix to store RTs
rts = zeros(1, EEG.trials); % Fixed size for storing RTs

for ei = 1:EEG.trials
    % find which event is time=0, and take the latency of the event thereafter.
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency) == 0);
    
    % use try-catch in case of no response
    try
        rts(ei) = EEG.epoch(ei).eventlatency{time0event + 1};
    catch
        rts(ei) = NaN; % Assign NaN if there's no response
    end
end

nyquist = EEG.srate / 2;
electrode_indx = [47 57 62];
frequencies = logspace(log10(4), log10(40), 20);
s = logspace(log10(3), log10(10), length(frequencies)) ./ (2 * pi * frequencies);

baselinetime = [-300 -100];
baseidx = dsearchn(EEG.times', baselinetime(1)):dsearchn(EEG.times', baselinetime(2));

witpc_matrix_all = zeros(length(electrode_indx), length(frequencies), EEG.pnts);
power_matrix_all = zeros(length(electrode_indx), length(frequencies), EEG.pnts);
for elecInd = 1:3
    witpc_matrix = zeros(length(frequencies), EEG.pnts);
    power_matrix = zeros(length(frequencies), EEG.pnts);
    for fi = 1:length(frequencies)
        % Create filter
        lower_filter_bound = frequencies(fi) - 0.5 * s(fi);
        upper_filter_bound = frequencies(fi) + 0.5 * s(fi);
        transition_width = 0.2;
        filter_order = round(3 * (EEG.srate / lower_filter_bound));
        ffrequencies = [0 (1 - transition_width) * lower_filter_bound lower_filter_bound upper_filter_bound (1 + transition_width) * upper_filter_bound nyquist] / nyquist;
        idealresponse = [0 0 1 1 0 0];
        filterweights = firls(filter_order, ffrequencies, idealresponse);
        
        % Apply filter to EEG data
        filtered_data = filtfilt(filterweights, 1, double(squeeze(EEG.data(electrode_indx(elecInd), :, :)))); % Assuming you're interested in the first electrode
        
        % Compute Hilbert transform
        hilbert_data = hilbert(filtered_data);
        
        % Compute ITPC with RTs
        witpc_matrix(fi, :) = abs(mean(repmat(rts, EEG.pnts, 1) .* exp(1i * angle(hilbert_data)), 2)); % Using repmat to match dimensions
        
        % Compute power
        power_matrix(fi, :) = mean(abs(hilbert_data).^2, 2);
        power_matrix(fi, :) = 10*log10(power_matrix(fi, :) ./ mean(power_matrix(fi, baseidx), 2));
    end
    witpc_matrix_all(elecInd, :, :) = witpc_matrix(:, :);
    power_matrix_all(elecInd, :, :) = power_matrix(:, :);
end

baselinetime = [-300 -100];
baseidx = dsearchn(EEG.times', baselinetime(1)):dsearchn(EEG.times', baselinetime(2));
for figInd = 1:length(electrode_indx)
    figure;
    subplot(2, 1, 1)
    contourf(EEG.times, frequencies, squeeze(witpc_matrix_all(figInd, :, :)), 40, 'linecolor', 'none')
    set(gca, 'clim', [min(witpc_matrix_all(:)), max(witpc_matrix_all(:))], 'xlim', [-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('WITPC')
    colorbar % Add colorbar to show the scale
    
    subplot(2, 1, 2)
    contourf(EEG.times, frequencies, squeeze(power_matrix_all(figInd, :, :)), 40, 'linecolor', 'none')
    set(gca, 'clim', [-1 1], 'xlim', [-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('DB-power ')
    colorbar % Add colorbar to show the scale
end




%% 

% Load sample EEG data
load sampleEEGdata.mat
center_freq = 6;
channel_to_use = 'po7';
times_to_save = -200:50:1200;

rts = zeros(1, EEG.trials); % Fixed size for storing RTs

for trial_index = 1:EEG.trials
    time0_event = find(cell2mat(EEG.epoch(trial_index).eventlatency) == 0);
    
    % Use try-catch in case of no response
    try
        rts(trial_index) = EEG.epoch(trial_index).eventlatency{time0_event + 1};
    catch
        rts(trial_index) = NaN; % Assign NaN if there's no response
    end
end

% Frequency parameters
nyquist = EEG.srate / 2;
electrode_indices = [47 57 62];
frequencies = logspace(log10(4), log10(40), 20);
s = logspace(log10(3), log10(10), length(frequencies)) ./ (2 * pi * frequencies);

% Baseline time
baseline_time = [-300 -100];
baseline_index = dsearchn(EEG.times', baseline_time(1)):dsearchn(EEG.times', baseline_time(2));

phase_matrix_all = zeros(length(electrode_indices), length(frequencies), EEG.pnts);
rt_matrix_all = zeros(length(electrode_indices), length(frequencies), EEG.pnts);

% Phase and reaction time for each electrode and frequency
for electrode_index = 1:length(electrode_indices)
    for frequency_index = 1:length(frequencies)
       
        lower_filter_bound = frequencies(frequency_index) - 0.5 * s(frequency_index);
        upper_filter_bound = frequencies(frequency_index) + 0.5 * s(frequency_index);
        transition_width = 0.2;
        filter_order = round(3 * (EEG.srate / lower_filter_bound));
        ffrequencies = [0 (1 - transition_width) * lower_filter_bound lower_filter_bound upper_filter_bound (1 + transition_width) * upper_filter_bound nyquist] / nyquist;
        ideal_response = [0 0 1 1 0 0];
        filter_weights = firls(filter_order, ffrequencies, ideal_response);
        
        filtered_data = filtfilt(filter_weights, 1, double(squeeze(EEG.data(electrode_indices(electrode_index), :, :)))); % Assuming you're interested in the first electrode
        
        hilbert_data = hilbert(filtered_data);
        
        phase_matrix_all(electrode_index, frequency_index, :) = mean(angle(hilbert_data)', 1);
        
        rt_matrix_all(electrode_index, frequency_index, :) = mean(repmat(rts, EEG.pnts, 1), 2); % Using repmat to match dimensions
    end
end



for electrode_index = 1:length(electrode_indices)
    figure;
    subplot(2, 1, 1)
    contourf(EEG.times, frequencies, squeeze(phase_matrix_all(electrode_index, :, :)), 40, 'linecolor', 'none')
    set(gca, 'clim', [-pi pi], 'xlim', [-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('Phase')
    colorbar % Show the scale
    
    subplot(2, 1, 2)
    plot(rts);
    xlabel('Trials')
    ylabel('Reaction time')
end



