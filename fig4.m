clc;close all; warning off;
% clear;

% %%% 1. identify and fit particles, this may take hours
% rr=1:50;
% t1=0.02; t2=0.05 ;r1=5;
% means=1;
% 
% pics={'Raw.tif','Wie.tif','RL.tif','RL_Wie.tif','RL_Wie_filter.tif'};
% 
% tic;
% nframes=10;
% for i=1:size(pics,2)
%     a0=readsTIF(pics{i},300);
%     tic;
%     p{i,:}=ContrastAndLocalizationFrom2DGaussian1(a0,r1,t1,t2,means,rr);
% end
% toc;
% %%% 1. identify and fit particles, this may take hours

%%% 2. load identified particle
% p0=struct2cell(load('tracked_particles.mat'));
% p=p0{1};
%%% 2. load identified particle

%%% 3. track particles in different frames
% tracks=trackParticles(p);
%%% 3. track particles in different frames

%%% 4. calculate nornalised localization error of tracked particles 
errs=getLocalizationError2(tracks,30,58);
%%% 4. calculate nornalised localization error of tracked particles 



function tracks=trackParticles(p)
    %%% the threshold particle distance between frames to determine whether 
    %%% the particles were the same particle in different frames 
    t=1.2;
    %%% frame numbers
    mm=size(p{1},2);
    %%% stack numbers:raw, Wiener, RL, Wiener_RL, SF
    nn=size(p,1);
    %%% tracks in each stack
    tracks=cell(1,nn);
    tic;
    for i=1:nn
        %%% in each stack (i), the number j particle in the 1st frame of
        %%% each stack (i) p{i}{1}(j,:) was set to be the start particle in
        %%% each track
        jj=0;
        %%% particle number (size(xx,1)) in the 1st frame of each stack (i)
        nparticles=size(p{i}{1},1);
        for j=1:nparticles
            x0=p{i}{1}(j,2);
            y0=p{i}{1}(j,3);
            ps(1,1:6)=p{i}{1}(j,:);
            ps(1,7:10)=0;
            %%% for frame 
            for k=2:mm
                %%% number of particles in the number k frame of the i stack
                nparticlesf=size(p{i}{k},1);
                dis=zeros(nparticlesf,1);
                for l=1:nparticlesf
                    dis(l)=hypot(x0-p{i}{k}(l,2),y0-p{i}{k}(l,3));
                end
                diss=min(dis);
                if diss<t
                    idx=find(dis==diss);
                    ps(k,1:6)=p{i}{k}(idx(1),:);
                    ps(k,7)=diss;
                    ps(k,8)=ps(k,2)-p{i}{1}(j,2);
                    ps(k,9)=ps(k,3)-p{i}{1}(j,3);
                    ps(k,10)=hypot(ps(k,8),ps(k,9));
                    x0=ps(k,2);
                    y0=ps(k,3);
                end
            end
            if size(ps,1)>20
                jj=jj+1;
                tracks{i}{jj}=ps;
            end
        end
    end
    toc;
end

function errs=getLocalizationErrors(tracks,steps)
    nn=size(tracks,2);
    errs=cell(nn,1);
    for i=1:nn
        mm=size(tracks{i},2);
        err=zeros(mm,1);
        for j=1:mm
            err(j)=tracks{i}{j}(steps,8);
        end
        mean_err=mean(err);
        
    end
end

function errs=getLocalizationError(tracks,steps)
    nn=size(tracks,2);
    errs=cell(nn,1);
    for i=1:nn
        mm=size(tracks{i},2);
        errs{i}=cell(mm,1);
        vectors=zeros(mm,2); %%% vectors(j)
        tic;
        for j=1:mm
            vectors(j,:)=[tracks{i}{j}(1,2)-tracks{i}{j}(1+steps,2),tracks{i}{j}(1,3)-tracks{i}{j}(1+steps,3)];
        end
        toc;
        mean_vec=mean(vectors);
        errs{i}=vectors-mean_vec;
    end
end

function errs=getLocalizationError2(tracks,steps,pixels)
    nn=size(tracks,2);
    errs=cell(nn,1);
    for i=1:nn
        mm=size(tracks{i},2);
        errs{i}=cell(mm,1);
        vectors=zeros(mm,2); 
        % tic;
        for j=1:mm
            vectors(j,1:2)=[tracks{i}{j}(steps,8),tracks{i}{j}(steps,9)];
        end
        % toc;
        mean_vec=mean(vectors);
        rad=hypot(mean_vec(1),mean_vec(2));
        errs{i}=pixels*(vectors-mean_vec)/rad;
        maxs(i)=max(max(abs(errs{i})));
    end
    ct=[0.5,0.8,1];
    getRatioCircle(errs,ct);
end

function getRatioCircle(ps,ct)
    nn=max(size(ps));
    lims=70;
    nbins=[10,10,10,10,10];
    for i=1:nn
        mm=max(size(ps{i}));
        dis=zeros(1,mm);
        for j=1:mm
            dis(j)=hypot(ps{i}(j,2),ps{i}(j,1));
        end
        diss=sort(dis);
        cts=ceil(diss(ceil(ct*mm)));
        ctss{i}=cts;
        subplot(2,nn,i); hold on;
        scatter(ps{i}(:,1),ps{i}(:,2),10,'filled'); xlim([-lims,lims]); ylim([-lims,lims]);
        set(gca,'ytick',-60:30:60); set(gca,'xtick',-60:30:60); axis equal;
        for k=1:max(size(cts))
            rectangle('Position', [-cts(k), -cts(k), 2*cts(k), 2*cts(k)], 'Curvature', [1 1], 'EdgeColor', 'r');
        end
        hold off;
        subplot(2,nn,i+nn); h=histfit(dis',nbins(i));
        % h(1).FaceColor = colorss1(i,:); h(2).Color ='k';
        pd = fitdist(dis','Normal');
        fs1=sprintf('%.3f',pd.mu);
        text(0.1+pd.mu,20,fs1);
        xlim([0,70]);
        xlabel('Absolute localization error');  ylabel('Counts'); %title(names1{i});
    end
end

%%% sub pixel localization throuth threshold restrained RVT and 2D Gaussian
%%% fitting, a0 the raw image, r the roi size, t1 the RVT threshold, t2 the
%%% vary range of background intensity of roi image, rr the RVT calculation
%%% range
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

%%% 2D Gaussian fit to approximate the sub-pixel localization and peak intensity
%%% a0 the roi image, background intensity were subtracted to around 0,  t
%%% the background intensity vary range from 0, r the size of roi
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