%%% fig 1 presentation
clc; close all;
clear;
warning off;

m=1; n=2;
p=[20,1;20,40;35,1;35,40];
colorss={'k','r'};
a0=readTIF('deconvs_compares.tif');
ps=[301,443;301,453];
l=20;
xx=1:2*l+1; xx=xx*0.058;
idx=[1,5];
[lines,roi]=getROI(a0,ps,20,idx);

for i=1:2
    subplot(2,2,2*i-1); imshow(roi(:,:,i),[]); colorbar; clim([-0.3,0.1]);
    subplot(2,2,2*i); hold on;
    plot(xx,lines(:,2*i-1)); 
    mins=min(lines(:,2*i-1));
    fs1=sprintf('%.3f',mins);
    text(0.7,mins,fs1);
    plot(xx,lines(:,2*i)); ylim([-0.3,0.1]);
    hold off;
end
%%% fig 1 presentation

function [lines,roi]=getROI(a,p,l,idx)
    x1=p(1,1)-l:p(1,1)+l;
    y1=p(1,2)-l:p(1,2)+l;
    x2=p(2,1)-l:p(2,1)+l;
    for i=1:max(size(idx))
        roi(:,:,i)=a(y1,x1,idx(i));
        lines(:,2*i-1)=a(p(1,2),x1,idx(i));
        lines(:,2*i)=a(p(2,2),x2,idx(i));
    end
end
