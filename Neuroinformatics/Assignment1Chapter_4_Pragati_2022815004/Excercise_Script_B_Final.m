% Exercise for Script B
%% Question6
%filepath= "C:\Users\Pragati\git\Neuroinformatics\Assignment1Chapter_4_Pragati_2022815004\amsterdam.bmp" 
Amsterdam_img = imread('C:\Users\Pragati\git\Neuroinformatics\Assignment1Chapter_4_Pragati_2022815004\amsterdam.bmp');
imshow(Amsterdam_img);
title('Amsterdam');
%% Question7 (to plot a thick red line from "Nieuwmarkt" to "Station Amsterdam Centraal"
Amsterdam_img = imread('C:\Users\Pragati\git\Neuroinformatics\Assignment1Chapter_4_Pragati_2022815004\amsterdam.bmp');
imshow(Amsterdam_img);
title('Amsterdam');
grid on;
grid minor;
axis on; 
% Nieuwmarkt_coordinates = [380,350]; Coordinates was accessible after
% axis was on. 
% Station_Amsterdam_Central = [375,75]; 
hold on;  % Allowing multiple plots/edits on the same figure
line([380 350],[375 75], 'color','r', 'LineWidth', 2);  % Plot the red line on the amsterdam map
plot(380, 350, 'ro', 'MarkerSize', 8); % placing the marker on point Nieuwmarkt  
plot(375, 75, 'ro', 'MarkerSize', 8); % placing the marker on point Station Amsterdam Central
hold off
%% Question8 (Plot a magenta star over the Waterlooplein metro station) 
Amsterdam_img = imread('C:\Users\Pragati\git\Neuroinformatics\Assignment1Chapter_4_Pragati_2022815004\amsterdam.bmp');
imshow(Amsterdam_img);
title('Amsterdam');
grid on;
axis on;
hold on; 
plot(390,492, 'p', 'MarkerSize', 20, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
hold off;
%% Question9 (Identifying the maximum value on each color dimension and plotted the circle using plot
% circle function and picking the pixel randomly for plotting 
Amsterdam_img = imread('C:\Users\Pragati\git\Neuroinformatics\Assignment1Chapter_4_Pragati_2022815004\amsterdam.bmp');
imshow(Amsterdam_img);
title('Amsterdam');

hold on;

% Find the maximum value on each color dimension
maxRed = max(max(Amsterdam_img(:,:,1)));
maxGreen = max(max(Amsterdam_img(:,:,2)));
maxBlue = max(max(Amsterdam_img(:,:,3)));

disp =maxRed;
disp =maxGreen;
disp =maxBlue;
[rowRed, colRed] = find(Amsterdam_img(:,:,1) == maxRed);

% Randomizing the selection of index x for maxRed
xRed = randi(numel(rowRed));
maxRedCoords = [rowRed(xRed), colRed(xRed)];

% Find the indices of the pixels with maximum green value
[rowGreen, colGreen] = find(Amsterdam_img(:,:,2) == maxGreen);

% Randomize the selection of index x for maxGreen
xGreen = randi(numel(rowGreen));
maxGreenCoords = [rowGreen(xGreen), colGreen(xGreen)];

% Find the indices of the pixels with maximum blue value
[rowBlue, colBlue] = find(Amsterdam_img(:,:,3) == maxBlue);

% Randomize the selection of index x for maxBlue
xBlue = randi(numel(rowBlue));
maxBlueCoords = [rowBlue(xBlue), colBlue(xBlue)];

plotCircle(maxRedCoords(1), maxRedCoords(2), 10, [maxRed,0,0])
hold on;
plotCircle(maxGreenCoords(1), maxGreenCoords(2), 10, [0,maxGreen,0])
hold on;

plotCircle(maxBlueCoords(1), maxBlueCoords(2), 10, [0,0,maxBlue])
hold off;
