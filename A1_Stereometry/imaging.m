folder = '10sided';
filePattern = fullfile(folder, '*.JPG'); 
rawfiles = dir(filePattern);
%n = 4; %four sided
n = 10; %six sided
%%
R = imread(fullfile(folder, rawfiles(1).name)); 
rect = [600, 1500, 2000, 2000];
R = imcrop(R, rect);
I = (im2gray(R));
imshow(I);
enableDefaultInteractivity(gca);

hold on;

%%
coord_a = zeros(n,2);
for i=1:n
    enableDefaultInteractivity(gca);
    [xi, yi, b] = ginput(1);
    coord_a(i,:) = [xi, yi];
    plot(coord_a(:, 1), coord_a(:,2), 'g+');
    zoom out;
end
hold off
exportgraphics(gca, 'a_corners.jpg');

%%
R = imread(fullfile(folder, rawfiles(2).name)); 
R = imcrop(R, rect);
I = (im2gray(R));
imshow(I);
enableDefaultInteractivity(gca);
hold on;
coord_b = zeros(n,2);

for i=1:n
    enableDefaultInteractivity(gca);
    [xi, yi, b] = ginput(1);
    coord_b(i,:) = [xi, yi];
    plot(coord_b(:, 1), coord_b(:,2), 'g+');
    zoom out;
end
hold off
exportgraphics(gca, 'b_corners.jpg');
%%
f = 26;
b = -10;

x_a = coord_a(:, 1);
x_b = coord_b(:, 1);
x = x_a - x_b;
z = [];
for i = 1:n
    p = (f*b)/x(i,1);
    z(i, 1)=p;
end

%% plot
y = coord_a(:,2);
xa = round(x_a);
ya = round(y);
x = linspace(min(xa), max(xa), 100);
y = linspace(min(ya), max(ya), 100);
[X, Y] = meshgrid(x, y);

Z = griddata(xa, ya, z, X, Y, 'linear');

figure;
% subplot(1,2,1);
mesh(X,Y,Z)
xlabel('X-axis (px)');
ylabel('Y-axis (px)');
zlabel('Depth');
set(gca,'FontSize',14);
h = colorbar;
set(h, 'fontsize', 14);
azimuth = 8; elevation = -70; % to make viewing angle similar to all trials
view([azimuth elevation])

% 
% subplot(1,2,2);
% mesh(X,Y,Z)
% xlabel('X-axis (px)');
% ylabel('Y-axis (px)');
% zlabel('Depth');
% set(gca,'FontSize',14);
% h = colorbar;
% set(h, 'fontsize', 14);
% ylabel(h, 'Depth', 'FontSize',14, 'Rotation',270);
% azimuth = 0; elevation = -90; % to make viewing angle similar to all trials
% view([azimuth elevation])
