
%% Get surface normal function
function n = get_surfacenormal(V, I)
 g = inv(V'*V)*V'*I;
 n = normalize(g, 'norm'); %Euclidean norm
end

%% Get partial derivatives function

function [df_dx, df_dy] = get_df(n, imagesize)
%Reshape
nx = reshape(n(1, :), imagesize);
ny = reshape(n(2, :), imagesize);
nz = reshape(n(3, :), imagesize);

%Get partial derivatives
df_dx = -nx./nz;
df_dy = -ny./nz;

end
%% Line integral function
function f = get_lineintegral(df_dx, df_dy)
xcumsum = cumsum(df_dx, 2);
ycumsum = cumsum(df_dy, 1);
f = xcumsum + ycumsum;
end

%% Using photos.mat
%Load and extract images
images = load('Photos.mat');
I1 = images.I1; I2 = images.I2; I3 = images.I3; I4 = images.I4;
figure(1); 
subplot(2, 2, 1); imagesc(I1);
subplot(2, 2, 2); imagesc(I2);
subplot(2, 2, 3); imagesc(I3);
subplot(2, 2, 4); imagesc(I4);
colormap(gray);

n_imagesize = size(I1);

%Get V matrix
V1 = [0.085832, 0.17365, 0.98106];
V2 = [0.085832, -0.17365, 0.98106];
V3 = [0.17365, 0, 0.98481];
V4 = [0.16318, -0.34202, 0.92542];
V = [V1; V2; V3; V4];
disp(size(V));

%Get I Matrix
I = [I1(:)'; I2(:)'; I3(:)'; I4(:)'];
disp(size(I));

%Computation
n_photos = get_surfacenormal(V, I);
[n_dfdx, n_dfdy] = get_df(n_photos, n_imagesize);
n_f = get_lineintegral(n_dfdx, n_dfdy);

% Plots
figure(2);
imagesc(n_photos);
title('surface normal')
colormap(gray);

figure(3); 
mesh(n_f);
title('3D')

%% Using 20 cat images
filefolder = fullfile('A3_Turtle', 'Image_*.png'); 
imagefiles = dir(filefolder); 
nfiles = length(imagefiles); 
disp(nfiles);
images = cell(1, nfiles);

I_data = [];
for ii = 1:nfiles
   currentfilename = fullfile('A3_Turtle', imagefiles(ii).name); 
   currentimage = imread(currentfilename);  
   reshaped = currentimage(:)';
   I_data = double([I_data; reshaped]);
   images{ii} = currentimage;  
   
end

data_imagesize = size(images{1});
V_data = load('turtle_light_directions.txt')';
disp(size(V_data));
disp(size(I_data));

%Computation
n_data = get_surfacenormal(V_data, I_data);
[data_dfdx, data_dfdy] = get_df(n_data, data_imagesize);
data_f = get_lineintegral(data_dfdx, data_dfdy);

% % Plots
figure(3);
imagesc(n_data);
title('surface normal')
colormap(gray);

figure(3); 
mesh(data_f);
title('3D');
colormap('turbo');

% Frankot Chellappa
fc = frankotChellappa(data_dfdx, data_dfdy);
figure(4); mesh(fc); colormap("gray");


