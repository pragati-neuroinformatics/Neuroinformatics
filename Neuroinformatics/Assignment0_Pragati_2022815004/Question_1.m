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
disp = stats
%The p-value is 0.0605 and  h= 0 whichindicates that 
% ttest2 does not reject the null hypothesis at the default 5% significance level.

% Question 2 
t= [0:0.001:4]; %(Time-period);
f = [2 8 12 25]; %(Frequrencies in Hz)
A = [0.5 1 1.5 2];
phase = [pi/1 pi/2 pi/3 pi/4];
figure;
for i = 1:4 
    x(i,:) = A(i)*sin(2*pi*f(i)*t+phase(i));
    subplot (5,1,i)
    plot (x(i,:))
    xlim([0 4000])
end 
y = mean(x,1);
subplot (5,1,5)
plot (y)
xlim([0 4000])
%




