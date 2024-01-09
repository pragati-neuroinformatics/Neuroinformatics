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

% Perform independent t-test
[h, p, ci, stats] = ttest2(condition1, condition2);
%The p-value is 0.0605. 

% Question 2 
Fs= 2;
a = sin(2*pi*60*t) 
