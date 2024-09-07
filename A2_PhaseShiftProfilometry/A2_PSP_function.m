function Zw = get_phase(obj_images, ref_images, x, y, l, w, freq, rotation_angle)
%% GET_PHASE FUNCTION FOR GETTING Zw
%% Import and convert object images to grayscale and rotate
    I1o = imrotate(im2gray(imread(obj_images{1})), rotation_angle);
    I2o = imrotate(im2gray(imread(obj_images{2})), rotation_angle);
    I3o = imrotate(im2gray(imread(obj_images{3})), rotation_angle);
    I4o = imrotate(im2gray(imread(obj_images{4})), rotation_angle);
%% Import and convert reference images to grayscale and rotate
    I1 = imrotate(im2gray(imread(ref_images{1})), rotation_angle);
    I2 = imrotate(im2gray(imread(ref_images{2})), rotation_angle);
    I3 = imrotate(im2gray(imread(ref_images{3})), rotation_angle);
    I4 = imrotate(im2gray(imread(ref_images{4})), rotation_angle);
%% Crop images
    rect = [x, y, l, w];
    I1o = imcrop(I1o, rect);
    I2o = imcrop(I2o, rect);
    I3o = imcrop(I3o, rect);
    I4o = imcrop(I4o, rect);
    I1 = imcrop(I1, rect);
    I2 = imcrop(I2, rect);
    I3 = imcrop(I3, rect);
    I4 = imcrop(I4, rect);

    figure(1); 
    subplot(2,2,1); imshow(I1o);
    subplot(2,2,2); imshow(I2o);
    subplot(2,2,3); imshow(I3o);
    subplot(2,2,4); imshow(I4o);

    % Convert images to double-precision float
    I1o = double(I1o);
    I2o = double(I2o);
    I3o = double(I3o);
    I4o = double(I4o);
    I1 = double(I1);
    I2 = double(I2);
    I3 = double(I3);
    I4 = double(I4);
%% Calculate wrapped phase for object and background
    wrapped_object = atan2((I4o-I2o), (I1o-I3o));
    wrapped_bg = atan2((I4-I2), (I1-I3));
    figure(2); imagesc(wrapped_object);
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    h=colorbar;
    set(h,'fontsize',14);
    set(gca,'FontSize',13)
    figure(3); imagesc(wrapped_bg);
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    h=colorbar;
    set(h,'fontsize',14);
    set(gca,'FontSize',13);
%% Unwrap phase 
    uo = unwrap_TIE_FFT_DCT_iter(wrapped_object);
    ub = unwrap_TIE_FFT_DCT_iter(wrapped_bg);
    figure(4);imagesc(uo);
    figure(5); imagesc(ub);
%% Plot the unwrapped phase as a 3D surface
    [X, Y] = meshgrid(1:size(uo, 2), 1:size(uo, 1));
    
    figure(6);
    mesh(X, Y, uo);
    %title('3D Plot of Unwrapped Phase');
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    zlabel('Unwrapped Phase');
    set(gca,'FontSize',13)
    colormap();
    h=colorbar;
    set(h,'fontsize',14);
%% Plot the unwrapped phase as a 3D surface
    [X, Y] = meshgrid(1:size(ub, 2), 1:size(ub, 1));
    figure(7);
    mesh(X, Y, ub);
    h=colorbar;
    set(h,'fontsize',14);
    %title('3D Plot of Reference Unwrapped Phase');
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    zlabel('Unwrapped Phase');
    set(gca,'FontSize',13)
%% Own reference background
%se = strel('square', 40);
%se = strel('disk', 40);
%background = imopen(uo, se);
%figure(7); imagesc(background);
%figure(8); mesh(background)
%saveas(gcf, 'figure8.png'); 
%colormap(gca, gray);
%% Phase to height
    L = 65;
    p = freq;
    d = 26;
    z = 1.2;
    phi = uo - ub;
    figure(9); mesh(phi);
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    zlabel('Phase');

    % Zoom parameters
    L_prime = L/z;
    d_prime = d/z;
    p_prime = p/z;
    
    % Calculate height map
    Zw = ((L_prime * p_prime * phi) ./ ((2 * pi * (d_prime)) + (p_prime * phi))) + 2;

%% Plot height map of Zw
    figure(10);
    
    subplot(1,2,1); 
    [X, Y] = meshgrid(1:size(Zw, 2), 1:size(Zw, 1));
    mesh(X, Y, Zw);
    %title('Height Map of Zw');
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    zlabel('Height (Zw)');
    set(gca,'FontSize',13)
    %zlim([-10 4]);

    subplot(1,2,2); imagesc(Zw)
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
    set(gca, 'YDir', 'normal');
    set(gca,'FontSize',13)
    h=colorbar;
    set(h,'fontsize',14);

%% Making it enhanced
% Calculate the 2D Fourier Transform of the original image
F = fftshift(fft2(Zw));

% Apply a Gaussian filter in the frequency domain
sigma = 40; 
[x, y] = meshgrid(-size(F, 2)/2:size(F, 2)/2-1, -size(F, 1)/2:size(F, 1)/2-1);
gaussian_filter = exp(-(x.^2 + y.^2) / (2*sigma^2));
filtered_F = F .* gaussian_filter;
filtered_Zw = ifft2(ifftshift(filtered_F));

figure(11);
imagesc(abs(filtered_Zw));
set(gca, 'YDir', 'normal');
%colormap(gray);
axis off;
title('Filtered Image');

% 3D surface plot of the filtered image
figure(12);
s_ = surf(abs(filtered_Zw()));
s_.EdgeColor = "none";
%colormap(gray);
title('Filtered 3D Surface Plot');
%% Getting center
threshold = max(Zw(:)) * 0.9; 
binaryImage = Zw > threshold;

% Find properties
stats = regionprops(binaryImage, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');

imagesc(Zw);
set(gca, 'YDir', 'normal');
xlabel('X-axis (px)');
ylabel('Y-axis (px)');
colorbar;

if ~isempty(stats)
    % Get the centroid of the largest region 
    [~, idx] = max([stats.MajorAxisLength]); % Find largest component
    center = stats(idx).Centroid;
    
    % Plot the detected center on the height map
    hold on;
    plot(center(1), center(2), 'r*', 'MarkerSize', 10); 
    disp(['Detected center at: (X: ', num2str(center(1)), ', Y: ', num2str(center(2)), ')']);
end
hold off;

disp(Zw(int(center(2)), int(center(1))))
end

%% ACTUAL CODE HERE
% Object images
obj_images_0 = {"0Ao.JPG", "0Bo.JPG", "0Co.JPG", "0Do.JPG"};
obj_images_1 = {"1Ao.JPG", "1Bo.JPG", "1Co.JPG", "1Do.JPG"};
obj_images_2 = {"2Ao.JPG", "2Bo.JPG", "2Co.JPG", "2Do.JPG"};
obj_images_3 = {"3Ao.JPG", "3Bo.JPG", "3Co.JPG", "3Do.JPG"};
obj_images_4 = {"4Ao.JPG", "4Bo.JPG", "4Co.JPG", "4Do.JPG"};
obj_images_5 = {"5Ao.JPG", "5Bo.JPG", "5Co.JPG", "5Do.JPG"};

% Reference images
ref_images_0 = {"0A.jpg", "0B.JPG", "0C.JPG", "0D.JPG"};
ref_images_1 = {"1A.jpg", "1B.JPG", "1C.JPG", "1D.JPG"};
ref_images_2 = {"2A.JPG", "2B.JPG", "2C.JPG", "2D.JPG"};
ref_images_3 = {"3A.JPG", "3B.JPG", "3C.JPG", "3D.JPG"};
ref_images_4 = {"4A.JPG", "4B.JPG", "4C.JPG", "4D.JPG"};
ref_images_5 = {"5A.JPG", "5B.JPG", "5C.JPG", "5D.JPG"};

%% Cropping images
x_0 = 2000; % X top left
y_0 = 650; % Y top left
l_0 = 550; % Width 
w_0 = 950; % Height

x_1 = 2200; % X top left
y_1 = 800; % Y top left
l_1 = 550; % Width 
w_1 = 950; % Height

x_2 = 2100; % X top left
y_2 = 1200; % Y top left
l_2 = 300;  % Width 
w_2 = 300;  % Height

x_3 = 2100; % X top left
y_3 = 1200; % Y top left
l_3 = 300;  % Width 
w_3 = 300;  % Height

x_4 = 2100; % X top left
y_4 = 1200; % Y top left
l_4 = 300;  % Width 
w_4 = 300;  % Height

x_5 = 2070; % X top left
y_5 = 1600; % Y top left
l_5 = 400;  % Width 
w_5 = 300;  % Height

%% Use function
%Zw_0 = get_phase(obj_images_0, ref_images_0, x_0, y_0, l_0, w_0, 3, 0);
%Zw_1 = get_phase(obj_images_1, ref_images_1, x_1, y_1, l_1, w_1, 3, 0);
%Zw_2 = get_phase(obj_images_2, ref_images_2, x_2, y_2, l_2, w_2, 50, 0);
Zw_3 = get_phase(obj_images_3, ref_images_3, x_3, y_3, l_3, w_3, 4, 0);
%Zw_4 = get_phase(obj_images_4, ref_images_4, x_4, y_4, l_4, w_4, 50, 0);
%Zw_5 = get_phase(obj_images_5, ref_images_5, x_5, y_5, l_5, w_5, 50, 0);
