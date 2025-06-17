%%% TOC presentation
clc; close all;
clear;

l=80;
p=[380,310];
a0=readTIF('E:\\data\\25.3.19\\nanobem24\\1\\deconvs_compares.tif');  

a(:,1:l)=a0(p(1)-l:p(1)+l,p(2)-l:p(2)-1,1);
a(:,1+l:2*l)=a0(p(1)-l:p(1)+l,p(2)-l:p(2)-1,5);

imshow(a,[]); colormap('parula');
%%% TOC presentation