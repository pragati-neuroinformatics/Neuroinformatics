% Question 1
load ('power.mat') % loading the power.matfile
Original_power_matrix = load ('power.mat');
Power_mat = Original_power_matrix.Power;

transposePower_mat = Power_mat'; %transposed the original power matrix
disp transposePower_mat

% To Evaluate if the observed power is significantly different across the
% two conditions 

% Power data
condition1 = transposePower_mat(:, 1); 
condition2 = transposePower_mat(:, 2); 

% Perform paired t-test (As we want to find out the difference between
% same subjects spectral power data for 2 conditions successful and unsuccessful) 
[h,p] = ttest(condition1,condition2)
% Findings: The returned value of h=1 indicates that ttest reject the null hypothesis
% at the default 5% signficant level (and the p value is less than 0.05 which means that 
% there is a signficant difference between the two conditions (successfull
% and unsuccessful) of that test subject. 
%%
% Question 2 
t= [0:0.001:4]; %(Time-period);
f = [2 8 12 25]; %(Frequrencies in Hz)
A = [0.5 1 1.5 2];
phase = [pi/1 pi/2 pi/3 pi/4];
figure;
for i = 1:4 
    x(i,:) = A(i)*sin(2*pi*f(i)*t+phase(i));
    subplot (7,1,i)
    plot (x(i,:))
    xlim([0 4000])
end 
%Mean extraction of all the 4 frequencies (mixed wave plot)
y = mean(x,1);
subplot (7,1,5)
plot (y)
xlim([0 4000])

% To recover the information from the mix-wave plot
Y = fft(y);
Z = abs (Y);
figure;
plot(Z);
plot([0:1:4000]*(1000/4000), Z); %to scale it down to the original version of frequencies  
%%
%Question-3
t= [0:0.001:4]; %(Time-period);
f = [2 8 12 25]; %(Frequrencies in Hz)
A = [0.5 1 1.5 2];
phase = [pi/1 pi/2 pi/3 pi/4];
noise_amplitude = [0.1 0.2 0.3 0.4]; %to add noise in the data
figure;
for i = 1:4 
    sw = A(i)*sin(2*pi*f(i)*t+phase(i));
    wn = noise_amplitude (i) * randn(size(t));
    sw(:, :) = sw(:,:) + wn(:,:);
    M = mean(sw,1);
    FR = fft(M);
    FRE = abs(FR);
    subplot(4,3,3*i-2)
    plot(sw)
    xlim([0 4001])
    subplot (4,3,3*i - 1)
    plot(FRE);
    plot([0:1:4000]*(1000/4000),FRE); 
    ni = abs(FR) < 200;
    FR(ni) = 0;
    filt = ifft(FR);
    subplot(4,3,3*i)
    plot(filt)
    xlim([0 4001])
 end 



