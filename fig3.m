clc;close all;
clear;

% tic;
% y01=readTIF('1.tif'); 
% y0=y01(:,:,1:100);
% [psf,otf]=getPsf(520,0.9,0.058,0.9,0.9,size(y0,1));
% 
% names={'Raw.tif','RL.tif','Wie.tif','RL_Wie.tif','RL_Wie_filter.tif'};
% y=compares(y0,psf,otf,names);
% 
% toc;

%%%%% 从去背景后的图像获取颗粒对比度，先RVT + 强度阈值获取粗定位，后ROI高斯拟合获取contrast、localization
%%%%% The particle contrast was obtained from the image after background removal, the RVT intensity 
%%%%% threshold was obtained for coarse positioning, and then the ROI Gaussian fitting was obtained for 
%%%%% contrast and localization
a0=readTIF('deconvs_compares.tif');
rr=1:50;
t1=0.02; t2=0.05 ;r1=5;
means=1;
tic;
p=ContrastAndLocalizationFrom2DGaussian1(a0,r1,t1,t2,means,rr);
toc;

tic;
nn=0; t=1;
for i=1:size(p{1},1)
    contrasts(1)=p{1}(i,1);
    for j=2:5
        for k=1:min(size(p{j},1),size(p{1},1))
            diss(k)=hypot(p{1}(i,2)-p{j}(k,2),p{1}(i,3)-p{j}(k,3));
        end
        dis(j)=min(diss);
        idx=find(dis(j)==diss);
        contrasts(j)=p{j}(idx(1),1);

        clear diss;
    end  
    if max(dis)<t
        nn=nn+1;
        contrastss(nn,:)=contrasts;
    end
end
toc;

%%% fig 3
%%% shown enhanced contrast upon raw contrast
tic;
xx=1:nn;
figure; hold on;
colors={'b','r'};
for i=2:5
    scatter(contrastss(:,1),contrastss(:,i),'filled');
    tb=table(contrastss(:,1),contrastss(:,i));
    model=fitlm(tb);
    r2(i-1)=model.Rsquared.Adjusted;
    linearCoefficients = polyfit(contrastss(:,1), contrastss(:,i), 1);
    xFit = 0:0.01:0.8;
    % Get the estimated values with 

   
    yFit = polyval(linearCoefficients, xFit);
    % Plot the fit
    hold on;
    plot(xFit, yFit);
    fs1=sprintf('y=%.3f x + %.3f',linearCoefficients(1),linearCoefficients(2));
    text(0.05, 0.5+0.05*i, fs1);
end
xlim([0 0.4]); ylim([0 1]);
legend('Wiener','Wiener fit','RL','RL fit','Wie RL','Wie RL fit','Wie RL filter','Wie RL filter fit');
xlabel('Raw contrast');  ylabel('Deconvoluted contrast');

%%% comparison of contrast distribution
names1={'Raw','Wiener','RL'};
figure; t=10;
nbins=[6,10,10,12,12];
for i=1:5
    subplot(1,5,i); h=histfit(contrastss(:,i),nbins(i));
    % h(1).FaceColor = colorss1(i,:); h(2).Color ='k';
    pd = fitdist(contrastss(:,i),'Normal');
    fs1=sprintf('%.3f',pd.mu);
    text(0.1+pd.mu,20,fs1);
    xlabel('Contrast');  ylabel('Counts'); %title(names1{i});
    xlim([0 1]); ylim([0 120]); set(gca,'ytick',0:40:120);
end
toc;
%%% comparison of contrast distribution
%%% fig 3


%%% nearest points at a line y= a*x+b to points a0
function a=getNearLinePoints(l,a0)
    a(:,1)=a0(:,1)*0.5*l(1)+a0(:,2)*0.5+l(2)*0.5;
    a(:,2)=a0(:,1)*0.5+a0(:,2)*0.5/l(1)-0.5*l(2)/l(1);
end
function corrs=getCorrLinear(a,a0)
    a=getNearLinePoints(l,a0);
    corrs=corr(a,a0);
end

function a=getThresholded(a,r)
    for i=1:size(a,3)
        a1=a(:,:,i);
        t=max(max(a1))*r;
        a1(a1<t)=0;
        a(:,:,i)=a1;
    end
end
%%% locate particles in a frame at pixel precision, a0 the raw image, a the
%%% RVT processed image, r the block size of particle, t the threshold
%%% intehsity for particle identification in normalised RVT image
function [ps,rois]=curseLocalization(a,r,t,a0)
    a=normalMaxMin(a);
    [h,w]=size(a);
    maxs=max(max(a));
    idx=0;
    a1=a;
    % imshow(a1); colormap('turbo');
    while maxs>t
        
        [y0,x0]=find(a==maxs);      
        y=max(1,y0(1)-r):min(y0(1)+r,h);
        x=max(1,x0(1)-r):min(x0(1)+r,w);
        if y0>r&&y0<h-r&&x0>r&&x0<w-r
            idx=idx+1;
            ps(idx,:)=[y0,x0];
            rois(:,:,idx)=a0(y,x);

            % a1(y(1),x)=1; a1(max(y),x)=1; a1(y,x(1))=1; a1(y,max(x))=1; 
            % imshow(a1); colormap('turbo');
        end
        a(y,x)=0;
        maxs=max(max(a));

        % a1(y(1),x)=1; a1(max(y),x)=1; a1(y,x(1))=1; a1(y,max(x))=1; 
        % imshow(a1); colormap('turbo');
    end
end

function k=Gaussian2DfitLs5(a0,r,t,means)
    [h,w,nn]=size(a0); 
    x=1:w; y=1:h;
    [X Y]=meshgrid(x,y);
    x1=w/2;
    x2=h/2;
    p0=[1,x1,x2,x1,x2,1];
    options = optimoptions('fmincon','MaxIter',100,'Display','off');
    
    for i=1:nn
        a=a0(:,:,i);
        GaussianLoss = @(p)sum(sum(abs(a-p(6)-p(1)*exp(-0.5*(X-p(2)).*(X-p(2))./(p(4)*p(4))-0.5*(Y-p(3)).*(Y-p(3))./(p(5)*p(5))))));
        k(i,:)=fmincon(GaussianLoss,p0,[],[],[],[],[0.00001,x1-r,x2-r,0,0,means-t],[100,x1+r,x2+r,x1,x2,means+t],[],options);
    end
end

function k=ContrastAndLocalizationFrom2DGaussian(a,r,t1,t2,means,a0)
    [ps,rois]=curseLocalization(a,r,t1,means-a0);
    k=Gaussian2DfitLs5(rois,r,t2,means);
    k(:,2)=k(:,2)+ps(:,2)-r;
    k(:,3)=k(:,3)+ps(:,1)-r;
end

function k=ContrastAndLocalizationFrom2DGaussian1(a0,r,t1,t2,means,rr)
    a=RVT(a0,rr);
    a=getThresholded(a,t1);
    for i=1:size(a,3)
        [ps,rois]=curseLocalization(a(:,:,i),r,t1,means-a0(:,:,i));
        k{i}=Gaussian2DfitLs5(rois,r,t2,means);
        k{i}(:,2)=k{i}(:,2)+ps(:,2)-r;
        k{i}(:,3)=k{i}(:,3)+ps(:,1)-r;
    end
end

function k=ContrastAndLocalizationFrom2DGaussian2(a,r,t1,t2,means,a0)
    for i=1:size(a,3)
        [ps,rois]=curseLocalization(a(:,:,i),r,t1,means-a0(:,:,i));
        k{i}=Gaussian2DfitLs5(rois,r,t2,means);
        k{i}(:,2)=k{i}(:,2)+ps(:,2)-r;
        k{i}(:,3)=k{i}(:,3)+ps(:,1)-r;
    end
end

function a=compares(a0,psf,otf,names)
    a{1}=removeMedian(a0)-1;
    a{2}=removeMedian(RLdeconv(a0,psf,otf,5))-1;
    a0=WienerFilter(a0,otf);
    a{3}=removeMedian(a0)-1;
    a{4}=removeMedian(RLdeconv(a0,psf,otf,5))-1;
    a{5}=OTFfilter(a{4},0.55,0.058,520,0.9);
    for i=1:size(a,2)
        writeTIF(names{i},a{i});
        aa(:,:,i)=a{i}(:,:,1);
    end
    writeTIF('deconvs_compares.tif',aa);
end

function a=removeMedian(a0)
    med=getMedian(a0);
    for i=1:size(a0,3)
        a(:,:,i)=a0(:,:,i)./med;
    end
end

function a1=getSmoothed(a,r)
    [h,w]=size(a);
    
    h0=floor(h*0.5);
    w0=floor(w*0.5);
    l=hypot(h0,w0)*r;
    for i=1:h
        for j=1:w
            dis(i,j)=hypot(i-h0,j-w0);
        end
    end
    a1=a;
    a1(dis>l)=0;
end

function a1=OTFfilter(a,r,micronsPerPixel,lambda,NA)
    [h,w,nn]=size(a);
    cutoff=1000/(0.5*lambda/NA);
    cyclesPerMicron=1/(h*micronsPerPixel);
    otfr=r*cutoff/cyclesPerMicron;
    y=1:h;
    x=1:w;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-0.5*h,x-0.5*w);
    Mask1=zeros(h,w);
    Mask1(rad<otfr)=1;
    for i=1:nn
        a1(:,:,i)=real(ifft2(fftshift(fftshift(fft2(a(:,:,i))).*Mask1)));
    end
end


