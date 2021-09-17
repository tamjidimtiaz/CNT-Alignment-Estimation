% to split the image into 10X10 blocks

currentimage = imcrop(currentimage,[1 1 886 1024]);
currentimage = imresize(currentimage,[880 1020]);

Height=102; width=88;
% currentimage = imread('F:\carbon nano\project\partially aligned\new\partially_aligned_new.png');
currentimage = imread('F:\carbon nano\Jacques Doumani data\Bad\P5\5um_-0.1254_-5.1254_cropped.jpg');

% Laplacian Filtering 
sigma = 0.4;
alpha = 0.5;
beta = 1.5;
currentimage = locallapfilt(currentimage, sigma, alpha, beta);
% figure(), imshow(currentimage,[]);

[c,r]=size( currentimage );
count=0;

for i=1:Height:r
    for j=1:width:c
        if(i<r && j<c)
            count=count+1;
            Image=imcrop( currentimage ,[i j Height-1 width-1]);
            imwrite(Image, fullfile('F:\carbon nano\Jacques Doumani data\Bad\P5\5um_1\', [num2str(count), '.bmp']));
        end
    end
end
