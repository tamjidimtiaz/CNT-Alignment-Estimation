clc
clear all
close all
tic;

location = 'F:\carbon nano\Jacques Doumani data\Optical Image\P1\50x\*.bmp'; %  folder in which the images exist
     
ds = imageDatastore(location);         %  Creates a datastore for all images in your folder
noofportion = zeros(1,63);
count = 1;


% for 100 blocks we will calculate the power of fourier spectrum in 63 orientational direction. the orientation angle will be represented by radian

power = zeros(100,63);

% for each block we will get an individual order parameter so defining a
% vector to preserve each of the order parameter
S = zeros(1,100);

% Analysing each of the blocks preserved in a particular folder defined by
% the 'location' variable.
while hasdata(ds)

%process input image
src = read(ds);

% if rgb convert it to gray scale version
if numel(size(src))==3
    src = rgb2gray(src);
end

% figure(),imshow(src,[]);

% Image decomposition

lambda=0.8;
I = src;
I_s = decomposition_function(uint8(I),lambda,4,2);
I_s=double(I_s);% structural component
I_t=double(I)-double(I_s);%textural component
I = I_t;

[rows, columns, numberOfColorChannels] = size(I);

% figure,imshow(I_s,[]);
% figure,imshow(I_t,[]);

% make the image a squared one so that we can apply the radial summation operation uniformly

width = rows;
halfWidth = floor(width /2);
col1 = floor(rows/2) - halfWidth;
col2 = col1 + width;
row1 = 1;
row2 = columns;
if row1 < 1
    row1 = 1;
end
if row2 > rows
    row2 = rows;
end
if col1 < 1
    col1 = 1;
end
if col2 > columns
    col2 = columns;
end
I = I(row1:row2, col1:col2);



midx=floor(size(I,1)/2);
midy=floor(size(I,2)/2);
radius=floor(length(I)/2);
[xgrid, ygrid] = meshgrid(1:size(I,2), 1:size(I,1));
mask = ((xgrid-midx).^2 + (ygrid-midy).^2) <= radius.^2;
values = I(mask);
xCircle = xgrid(mask);
yCircle = ygrid(mask);
imgfpcirc=zeros(size(I));
for i=1:length(xCircle)
    imgfpcirc(xCircle(i),yCircle(i))=I(xCircle(i),yCircle(i));
end



level = 0.1;
a_bw = imgfpcirc;
a_bw(a_bw==-1)=0;
a_bw(a_bw~=0)=255;
L = a_bw;
level = 0.1;
a_bw = imbinarize(imgfpcirc,level);
Icomp = imcomplement(a_bw);
[B,L] = bwboundaries(a_bw);
% imshow(label2rgb(L, @jet, [.5 .5 .5]));
hold on
L(L~=1) = 255;
L  = imcomplement(L);

% Fourier Subtraction 

F2 = fft2(double(L));
shF = fftshift(F2);
shF(imgfpcirc==0)=0;
shF=shF/norm(shF);
imgfpcirc2 = abs(shF);

F1 = fft2(double(imgfpcirc));
shF = fftshift(F1);
shF(imgfpcirc==0)=0;
shF=shF/norm(shF);
imgfpcirc1 = abs(shF);


imgfpcirc = abs(imgfpcirc1-imgfpcirc2);
imgfpcirc = imrotate(imgfpcirc,-90);

% figure(), imshow(imgfpcirc,[]);

[N M] = size(imgfpcirc);

dimDiff = abs(N-M);
dimMax = max(N,M);


halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)
%% Compute radially average power spectrum
[X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho);
theta = round(theta,1);
theta(theta>0) = pi-theta(theta>0);
theta(theta<0) = -theta(theta<0);
theta(theta==pi) = 0;




% 
thetauni=unique(theta);


i = cell(length(thetauni), 1);
for r = 0:length(thetauni)-1
     i{r + 1} = find(theta == thetauni(r+1));
end

lengthofthetas = cellfun('length',i);
minthetas = min(lengthofthetas);
Pf = zeros(length(thetauni),1);
for r = 0:length(thetauni)-1
     Pf(r + 1,1) = sum(maxk(imgfpcirc(i{r+1}), minthetas));
     
end


x=thetauni;



Pf = medfilt1(Pf,5);
Pf = smooth(Pf);

y=(Pf-min(Pf))/(max(Pf)-min(Pf));
y(y<0.2)=0;
if thetauni(y==max(y))==0
    maxtheta=0;
else
    maxtheta = find(thetauni==randsample((thetauni(y==max(y))),1));
end


power(count,:) = y;

num = 0;
den = 0;
c=0;


% for j=1:length(thetauni)
%     c = cos((thetauni(j)))^2;
%     num = num + y(j)*c*sin(thetauni(j));
%     den = den + (y(j)*sin(thetauni(j)));
% end
% S(count) = abs((3*(num/den)-1)/2);
for j=1:length(thetauni)
    c = cos((thetauni(j)))^2;
    num = num + (y(j)*c);
    den = den + (y(j));
end
S(count) = abs((2*(num/den)-1));

% xi = linspace(min(x), max(x), 180);  
% y_1 = interp1(x, y, xi, 'spline', 'extrap');
% y_1 = (y_1-min(y_1))/(max(y_1)-min(y_1));
% y_1 = circshift(y_1,1-maxtheta);
% y_1(xi==max(xi))=y_1(xi==min(xi));
% x_1 = (xi*180)/pi;
% fig = figure('visible',false);
% plot(x_1, y_1, 'linewidth', 3);
% frame = getframe(fig);
% im = frame2im(frame);
% [A,map] = rgb2ind(im,256);
% imwrite(A, map, fullfile('F:\carbon nano\Jacques Doumani data\Optical Image\P4\P4_10x\distribution', [num2str(count), '.bmp']));


count = count+1



end

yfinal = sum(power)/100;

xi = linspace(min(x), max(x), 180);                     % Evenly-Spaced Interpolation Vector
y = interp1(x, yfinal, xi, 'spline', 'extrap');
yfinal = (yfinal-min(yfinal))/(max(yfinal)-min(yfinal));
maxtheta=find(thetauni==thetauni(yfinal==max(yfinal)));

y = (y-min(y))/(max(y)-min(y));
%%%% Please uncomment the below section if you want to calculate the S value
yfinal = circshift(yfinal,1-maxtheta);
%%%%
xi(xi==max(xi))=pi;
xi(xi==min(xi))=0;
y(xi==max(xi))=y(xi==min(xi));


yfinal(yfinal<0.2)=0;

num = 0;
den = 0;
c=0;

% for j=1:length(thetauni)
%     c = cos((thetauni(j)))^2;
%     num = num + (yfinal(j)*c*sin(thetauni(j)));
%     den = den + (yfinal(j)*sin(thetauni(j)));
% end
% Sfinal = abs((3*(num/den)-1)/2);
for j=1:length(thetauni)
    c = cos((thetauni(j)))^2;
    num = num + (yfinal(j)*c);
    den = den + (yfinal(j));
end
x = (xi*180)/pi;
Sfinal = abs((2*(num/den)-1));
figure(), plot(x, y, 'linewidth', 3);
ylabel('Normalized Power(AU)','Fontsize',30);
xlabel('orientation(Degree)','Fontsize',30);
set(gca,'fontsize',30);
xlim([0 177])
ylim([0 1])

toc;

