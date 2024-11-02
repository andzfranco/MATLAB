set(0, 'DefaultAxesFontSize', 14); % Default font size for axes
set(0, 'DefaultTextFontSize', 14);  % Default font size for text (titles, labels)
set(0, 'DefaultLegendFontSize', 14); % Default font size for legends
%% READING IMAGES
%set folder name
folder = '[0] night';
filePattern = fullfile(folder, '*.CR2'); 
rawfiles = dir(filePattern);
disp(length(rawfiles));

N = numel(rawfiles);
Ss_array = 1:N;
Im_array = cell(1, N);

% %reading each image and assigning images to a cell
for k = 1:N
  %assigning file name
  baseFileName = rawfiles(k).name;
  fullFileName = fullfile(folder, baseFileName);
  disp(fullFileName)
  image = imread(fullFileName);
  imageArray = imresize(image, size(rgb2gray(image)));

  %storing the image to a cell, it's a bit different from an array
  Im_array{k} = imageArray;
  info = imfinfo(fullFileName); 

  digicam = info.DigitalCamera;
  shutterSpeed = digicam.ExposureTime; 
  Ss_array(k) = log(shutterSpeed);
end
%% OBTAINING POINTS

firstImage = rgb2gray(Im_array{1}); 

figure;
imshow(firstImage);
title('Select points');
[xPoints, yPoints] = ginput(8); 
points = [xPoints, yPoints];
%%
numPoints = size(points, 1);
grayValuesMatrix = zeros(numPoints, N);
redValuesMatrix = zeros(numPoints, N);
greenValuesMatrix = zeros(numPoints, N);
blueValuesMatrix = zeros(numPoints, N);
%% LOOP OVER IMAGES AND PLOT SHUTTER SPEED VS GRAY VALUES
for k = 1:N
    image = Im_array{k};
  
    %for grayscale and RGB
    grayImage = im2gray(image);
    redChannel = image(:, :, 1);
    greenChannel = image(:, :, 2);
    blueChannel = image(:, :, 3);

    %store values in the matrices
    for p = 1:numPoints
        grayValuesMatrix(p, k) = grayImage(round(points(p, 2)), round(points(p, 1)));
        redValuesMatrix(p, k) = redChannel(round(points(p, 2)), round(points(p, 1)));
        greenValuesMatrix(p, k) = greenChannel(round(points(p, 2)), round(points(p, 1)));
        blueValuesMatrix(p, k) = blueChannel(round(points(p, 2)), round(points(p, 1)));
    end

    %for visualization
    figure;
    imshow(grayImage);
    hold on;
    plot(points(:, 1), points(:, 2), 'ro', 'MarkerSize', 10, 'linewidth', 3);
    hold off
    title(['Image ', num2str(k), ' with selected points']);
end

%% PLOT SHUTTER SPEED VS GRAY LEVELS FOR EACH POINT
figure;
hold on;
for p = 1:numPoints
    scatter(grayValuesMatrix(p, :), Ss_array, 'filled'); 
end
hold off;
ylabel('g(Zij) = ln(delta t)');
xlabel('Gray pixel value, Zij');
title('g(Zij) vs Gray pixel value for selected points');
%legend('(1637, 431)', '(1991, 581)', '(1979, 1439)', '(3653, 533)', '(4661, 1283)', '(4361, 2075)', '(1355, 2369)', '(575, 1913)'); 
legend('(665, 617)', '(2063, 593)', '(3647, 527)', '(2753, 1907)', '(3971, 1997)', '(4121, 3047)', '(2171, 3137)', '(791, 3071)');  
grid on;

% figure;
% hold on;
% for p = 1:numPoints
%     scatter(redValuesMatrix(p, :), Ss_array, 'filled'); 
%     scatter(blueValuesMatrix(p, :), Ss_array, 'filled'); 
%     scatter(greenValuesMatrix(p, :), Ss_array, 'filled'); 
% end
% hold off;
% ylabel('Shutter Speed');
% xlabel('Gray Pixel Value');
% title('Shutter Speed vs RGB Level for Selected Points');
% legend('Point 1', 'Point 2', 'Point 3', 'Point 4', 'Point 5');  
% grid on;

 %% PREPARING VARIABLES FOR GSOLVE
% %Weighting function
%w = @(z) double(z <= 127) .* z + double(z > 127) .* (255 - z);
w = @ weighting_function;

%Smoothing
lambda = 10;  

[g, lE] = gsolve(grayValuesMatrix, Ss_array, lambda, w);
[g_red, lE_red] = gsolve(redValuesMatrix, Ss_array, lambda, w);
[g_green, lE_green] = gsolve(greenValuesMatrix, Ss_array, lambda, w);
[g_blue, lE_blue] = gsolve(blueValuesMatrix, Ss_array, lambda, w);
%%
%grayscale plot
figure;
hold on;
scatter(0:255, g);
hold off;
xlabel('Gray Pixel Value');
ylabel('g(Zij) = lnEi + ln(delta tj)');
title('Response Function (Gray)');
grid on;

%RGB combined plot
figure;
hold on;
scatter(0:255, g_red, 'r');
scatter(0:255, g_green, 'g');
scatter(0:255, g_blue, 'b');
hold off;
xlabel('Pixel Value');
ylabel('g(Zij) = lnEi + ln(delta tj)');
title('Response Function for RGB Channels');
legend('Red Channel', 'Green Channel', 'Blue Channel');
grid on;

%indiv RGB
figure;
subplot(3, 1, 1); 
hold on;
scatter(0:255, g_red, 'r');
hold off;
xlabel('Pixel Value');
ylabel('g(Zij)');
title('Red Channel');
grid on;

subplot(3, 1, 2);
hold on;
scatter(0:255, g_green, 'g');
hold off;
xlabel('Pixel Value');
ylabel('g(Zij)');
title('Green Channel');
grid on;

subplot(3, 1, 3); 
hold on;
scatter(0:255, g_blue, 'b');
hold off;
xlabel('Pixel Value');
ylabel('g(Zij)');
title('Blue Channel');
grid on;

 
%%
%initialize the radiance maps
[height, width] = size(Im_array{1}(:,:,1)); 
numImages = numel(Im_array);

ln_radiance_red = zeros(height, width);  
ln_radiance_green = zeros(height, width);  
ln_radiance_blue = zeros(height, width);  
ln_radiance = zeros(height, width);

%each pixel
for i = 1:height
    for j = 1:width
        weightNume = 0;
        weightDeno = 0;
        weightedSum_red = 0;
        weightedSum_green = 0;
        weightedSum_blue = 0;

        weightTotal_Red = 0;
        weightTotal_Green = 0;
        weightTotal_Blue =0;
        
        for k = 1:numImages
            image = Im_array{k}(i, j);
            grayImage = im2gray(image);
            %pixel_value = Im_array{k}(i, j);
            pixel_value = grayImage;
            weight = weighting_function(pixel_value);
            weightNume = weightNume + (weight * (g(pixel_value+1) - Ss_array(k)));
            weightDeno = weightDeno + weight;
           
            % %pixel values
            redValue = Im_array{k}(i, j, 1);   
            greenValue = Im_array{k}(i, j, 2);
            blueValue = Im_array{k}(i, j, 3);  

            % %weights for each channel
            weight_red = weighting_function(redValue);
            weight_green = weighting_function(greenValue);
            weight_blue = weighting_function(blueValue);

            % % %accumulate weighted sum for each channel
            weightedSum_red = weightedSum_red + weight_red * (g_red(redValue + 1) - Ss_array(k));
            weightedSum_green = weightedSum_green + weight_green * (g_green(greenValue + 1) - Ss_array(k));
            weightedSum_blue = weightedSum_blue + weight_blue * (g_blue(blueValue + 1) - Ss_array(k));
            
            %accumulate total weight
            %weightTotal_RGB = weightTotal_RGB + weight_red + weight_green + weight_blue;
            weightTotal_Red = weightTotal_Red + weight_red;
            weightTotal_Green = weightTotal_Green + weight_green;
            weightTotal_Blue = weightTotal_Blue + weight_blue;
        end
       
        ln_radiance(i, j) = weightNume / weightDeno;
        ln_radiance_red(i, j) = weightedSum_red / weightTotal_Red;
        ln_radiance_green(i, j) = weightedSum_green / weightTotal_Green;
        ln_radiance_blue(i, j) = weightedSum_blue / weightTotal_Blue;
        
    end
end

%exponentiate to get actual radiance values
radiance_red = exp(ln_radiance_red);
radiance_green = exp(ln_radiance_green);
radiance_blue = exp(ln_radiance_blue);
radiance = exp(ln_radiance);

rgb_hdr_image = cat(3, radiance_red, radiance_green, radiance_blue);
gray_hdr_image = radiance;
%%
figure;
imshow(gray_hdr_image);
title('Radiance Map (Grayscale)');
%%
figure;
imshow(rgb_hdr_image);
title('Radiance Map (RGB)');


%%
grayscale = tonemapfarbman(gray_hdr_image, Saturation=0.8, Exposure=2);
rgb = tonemapfarbman(rgb_hdr_image, Saturation = 0.8, Exposure=2);
%%
figure;
imshow(grayscale);
title('Tone Map (Grayscale)');
%%
figure;
imshow(rgb);
title('Tone Map (RGB)');

%%
% %% Plot for different lambda
% lambda_list = [0, 20, 50, 100];
% marker_list = {'o', 's', '^', 'd'};
% g_lambda_list = [];
% figure;
% hold on;
% 
% for i = 1:length(lambda_list)
%     lambda = lambda_list(i);
%     w = @ weighting_function;
%     [g_lambda_list, lE_lambda_list] = gsolve(grayValuesMatrix, Ss_array, lambda, w);
%     scatter(0:255, g_lambda_list,'DisplayName', ['\lambda = ', num2str(lambda)], 'Marker', marker_list{i});
% 
% end
% 
% xlabel('Pixel Value');
% ylabel('g(Zij)');
% title('Response Curve for Different \lambda Values');
% legend show;
% grid on;
% hold off;