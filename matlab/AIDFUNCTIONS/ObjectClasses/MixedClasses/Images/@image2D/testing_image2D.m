%TEsting image2D
%% loading image

close all;clear all; clear classes;
cd C:\Work\images;
im = imread('Peter.jpg');
%figure;imshow(im);


obj = image2D('Image',im);
figure;imshow(obj);

%%


obj = image3D;
obj.Image = im;

newim = captureImage(obj);