clc; close all;
clear;

%%% gig 5 2. load processed images and present the ROI  and the signal and
%%% background cut-line intensity
a0=readTIF('compares_deconvs_filters.tif'); 
a0=a0(397:460,391:454,:);
figure; m=5; n=5;
max1=max(max(max(a0))); min1=min(min(min(a0))); l=max1-min1;
max1=max1+0.1*l; min1=min1-0.1*l;
for i=1:m-2
    max1=max(max(max(a0(:,:,n*i-n+1:n*i)))); min1=min(min(min(a0(:,:,n*i-n+1:n*i)))); l=max1-min1;
    max1=max1+0.1*l; min1=min1-0.1*l;
    for j=1:n
        subplot(m,n,i*n+j-n); 
        imshow(a0(:,:,i*n+j-n),[]); colorbar; 
        clim([min1,max1]);
    end
end

p=[44,14;44,45;17,14;17,45];
ylims=[-0.25,0.12;-0.1,0.05];
colors={'k','g','r'};
for i=1:2
    for j=1:n
        subplot(m,n,i*n+2*n+j); hold on;
        for k=1:3
            a(:,k)=getCrossLine(a0(:,:,n*k-n+j),p(i*2-1:2*i,:));
        end
        maxs=max(max(a)); mins=min(min(a)); l=maxs-mins;
        maxs=maxs+0.1*l; mins=mins-0.1*l;
        xx1=1:max(size(a)); xx1=xx1*0.058;
        for k=1:3
           plot(xx1,a(:,k),colors{k}); 
           
        end
        ylim(ylims(i,:)); hold off;
    end
end
%%% gig 5 2. load processed images and present the ROI  and the signal and
%%% background cut-line intensity

% %%% fig 5 1. excute filters and deconvolutions and save the processed images for the
% %%% presentation, this miy take hours
% tic;
% y01=readTIF('1.tif'); 
% y0=y01(:,:,1:100);
% [psf,otf]=getPsf(520,0.9,0.058,0.9,0.9,size(y0,1));
% 
% noise_type =  'gw';
% noise_var = 0.001; % Noise variance
% seed = 0; % seed for pseudorandom noise realization
% 
% names={'Raw','RL','Wiener','RL_Wie','RL_Wie_filter'};
% filters={'.tif','_NLM.tif','_BM3D.tif'};
% y=compares(y0,psf,otf,noise_type,noise_var,seed,names,filters);
% toc;
% %%% fig 5 1. excute filters and deconvolutions and save the processed images for the
% %%% presentation, this miy take hours


function a=compares(a0,psf,otf,noise_type,noise_var,seed,names,filters)
    for i=1:size(a0,3)
        a1(:,:,i)=imnlmfilt(a0(:,:,i));
    end
    a2=double(getBM3D(noise_type, noise_var, seed, a0));

    toc;
    for i=1:size(filters,2)
        switch i
            case 1
                a00=a0;
            case 2
                a00=a1;
            otherwise
                a00=a2;
        end
        a{i*5-4}=removeMedian(a00)-1;
        a{i*5-3}=removeMedian(RLdeconv(a00,psf,otf,5))-1;
        a00=WienerFilter(a00,otf);
        a{i*5-2}=removeMedian(a00)-1;
        a{i*5-1}=removeMedian(RLdeconv(a00,psf,otf,5))-1;
        a{i*5}=OTFfilter(a{i*5-1},0.45,0.058,520,0.9);
    end

    for i=1:size(filters,2)
        for j=1:size(names,2)
            fs1=[names{j},filters{i}];
            writeTIF(fs1,a{5*i+j-5});
            aa(:,:,5*i+j-5)=a{5*i+j-5}(:,:,1);
        end
    end
    writeTIF('compares_deconvs_filters.tif',aa);
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

%% BM3D functions download from github, which were attached here in one matlab script
function a = getBM3D(noise_type, noise_var, seed, a0)
    [noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, seed, size(a0));
    a = BM3D(a0, PSD);
end
function [noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, realization, sz)
    % randn('seed',realization);
    randn();
    
    % Get pre-specified kernel
    kernel = getExperimentKernel(noise_type, noise_var, sz);
    
    % Create noisy image
    half_kernel = ceil(size(kernel) ./ 2);
    if(numel(sz) == 3 && numel(half_kernel) == 2)
        half_kernel(3) = 0;
    end
    
    % Crop edges
    noise = convn(randn(sz + 2 * half_kernel), kernel(end:-1:1,end:-1:1, :), 'same');
    noise = noise(1+half_kernel(1):end-half_kernel(1), 1+half_kernel(2):end-half_kernel(2), :);
    
    PSD = abs(fft2(kernel, sz(1), sz(2))).^2 * sz(1) * sz(2);
end
function kernel = getExperimentKernel(noiseType, noiseVar, sz)
    % Case gw / g0
    kernel = ones(1);
    
    noiseTypes = {'gw', 'g0', 'g1', 'g2', 'g3', 'g4', 'g1w', 'g2w', 'g3w', 'g4w'};
    
    found = false;
    for nt = noiseTypes
        if strcmp(noiseType, nt)
            found = true;
        end
    end
    if ~found
        disp('Error: Unknown noise type!')
        return;
    end
    
    if (~strcmp(noiseType, 'g4') && ~strcmp(noiseType, 'g4w')) || ~exist('sz', 'var')
        % Crop this size of kernel when generating,
        % unless pink noise, in which
        % Case we want to use the full image size
        sz = [101, 101];
       
    end
    
    % Sizes for meshgrids
    sz2 = -(1 - mod(sz, 2)) * 1 + floor(sz/2);
    sz1 = floor(sz/2);
    [uu, vv] = meshgrid(-sz1(1):sz2(1), -sz1(2):sz2(2)); 
    alpha = 0.8;
    
    switch (noiseType(1:2))
        case 'g1'
            % Horizontal line
            kernel = 16 - abs((1:31)-16);
            
        case 'g2'
            % Circular repeating pattern
            scale = 1;
            dist = (uu).^2 + (vv).^2;
            kernel = cos(sqrt(dist) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);
            
        case 'g3' 
            % Diagonal line pattern kernel
            scale = 1;
            kernel = cos((uu + vv) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);    
            
        case 'g4'
            % Pink noise
            dist = (uu).^2 + (vv).^2;
            n = sz(1)*sz(2);
            spec3 = sqrt((sqrt(n)*1e-2)./(sqrt(dist) +  sqrt(n)*1e-2));
            kernel = fftshift(ifft2(ifftshift(spec3)));
    end    
    
    % -- Noise with additional white component --
    if numel(noiseType) == 3 && noiseType(3) == 'w'
        kernel = kernel / norm(kernel(:));
        kalpha = sqrt((1 - alpha) + (alpha) * abs(fft2(kernel, sz(1), sz(2))).^2);
        kernel = fftshift(ifft2(kalpha));
    end
    
    % Correct variance
    kernel = kernel / norm(kernel(:)) * sqrt(noiseVar);
end
function [y_est, blocks] = BM3D(z, sigma_psd, profile, stage_arg, blockmatches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM3D is an algorithm for attenuation of additive spatially correlated
%  stationary (aka colored) Gaussian noise in grayscale and multichannel images.
%
%
%  FUNCTION INTERFACE:
%
%  y_est = BM3D(z, sigma_psd, profile)
%
%  INPUT ARGUMENTS:
%
%  -- required --
%
%         'z' : noisy image (M x N or M x N x C double array, intensities in range [0,1])
%               For multichannel images, block matching is performed on the first channel.
%  'sigma_psd' : noise power spectral density (M x N double nonnegative array)
%               OR
%               noise STD
%               OR
%               either of these in multichannel:
%               (M x N x C PSDs or 1 x C STDs)
%
% -- optional --
%
%   'profile' : 'np' --> Normal Profile (default)
%               'refilter' --> Apply refiltering
%               OR
%               a BM3DProfile object specifying the parameters
%               some other premade profiles also included from the previous versions
%               in BM3DProfile.m
%
%   'stage_arg' : Determines whether to perform hard-thresholding or wiener filtering.
%                 either BM3DProfile.HARD_THRESHOLDING, BM3DProfile.ALL_STAGES or an estimate
%                  of the noise-free image.
%                    - BM3DProfile.ALL_STAGES: Perform both.
%                    - BM3DProfile.HARD_THRESHOLDING: Perform hard-thresholding only.
%                    - ndarray, size of z: Perform Wiener Filtering with stage_arg as pilot.
%
%   'blockmatches' : Tuple {HT, Wiener}, with either value either:
%                      - false : Do not save blockmatches for phase
%                      (default)
%                      - true : Save blockmatches for phase
%                      - Pre-computed block-matching array returned by a
%                      previous call with [true]
%  OUTPUT:
%      'y_est'  denoised image  (M x N double array)
%      'y_est', {'blocks_ht', 'blocks_wie'} denoised image, plus HT and
%          Wiener blockmatches, if any storeBM values are set to True
%          (or [0] for missing block array, if only one calculated)
%
%
%  BASIC SIMULATION EXAMPLES:
%
%     Case 1)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise
%        generated as convution of AWGN against with kernel 'k'
%
%      k=[-1;2;-1]*[1 4 1]/100;   % e.g., a diagonal kernel
%      z=y+imfilter(randn(size(y)),k(end:-1:1,end:-1:1),'circular');
%
%      % define 'sigma_psd' from the kernel 'k'
%
%      sigma_psd=abs(fft2(k,size(z,1),size(z,2))).^2*numel(z);
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigma_psd);
%
%
%     Case 2)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise
%      % generated as convution of AWGN against with kernel 'k'
%      [x2, x1]=meshgrid(ceil(-size(y,2)/2):ceil(size(y,2)/2)-1,ceil(-size(y,1)/2):ceil(size(y,1)/2)-1)
%      sigma_psd=ifftshift(exp(-((x1/size(y,1)).^2+(x2/size(y,2)).^2)*10))*numel(y)/100;
%      z=y+real(ifft2(fft2(randn(size(y))).*sqrt(sigma_psd)/sqrt(numel(y))));
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigma_psd);
%
%     Case 3) If 'sigma_psd' is a singleton, this value is taken as sigma and
%             it is assumed that the noise is white variance sigma^2.
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive white Gaussian noise with variance sigma^2
%      sigma=0.1;
%      z=y+sigma*randn(size(y));
%
%      y_est = BM3D(z, sigma);
%
%      % or, equivalently,
%      sigma_psd = ones(size(z))*sigma^2*numel(z)
%      y_est = BM3D(z, sigma_psd)
%
%
%      Case 4)   MULTICHANNEL PROCESSING
%
%      y_est = BM3D(cat(3, z1, z2, z3), sigma_psd, 'np'); 
%
%      Multiple PSDs are optionally handled in the same way.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2019 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only
% be used for nonprofit noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes
% is prohibited.
%
% AUTHORS:
%     Y. Mäkinen, L. Azzari, K. Dabov, A. Foi
%     email: ymir.makinen@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('profile','var')
        profile         = 'np'; %% default profile
    end
    
    if isa(profile, 'string') || isa(profile, 'char')
        profile = BM3DProfile(profile);
    end
    
    if ~exist('stage_arg','var')
        stage_arg = profile.ALL_STAGES;  % By default, do both HT and Wie
    end
    
    if min(size(z, 1), size(z, 2)) < profile.N1 || min(size(z, 1), size(z, 2)) < profile.N1_wiener
        disp('Error: Image cannot be smaller than block size!')
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pro = convertToBM4DProfile(profile);
    
    channel_count = size(z, 3);
    dont_do_bm = ~exist('blockmatches','var') || numel(blockmatches) ~= 2 || ...
                    (blockmatches{1} == false && blockmatches{2} == false);
    
    if channel_count > 1  % Channel dimension should be 4
        z = permute(z, [1, 2, 4, 3]);
        if ndims(sigma_psd) == 3
            sigma_psd = permute(sigma_psd, [1, 2, 4, 3]);
        end
    
        if ~dont_do_bm
            disp('Warning: block-match data supplied with multichannel BM3D will be discarded. Call the function separately for each channel.');
        end
        [y_est] = BM4D_multichannel(z, sigma_psd, pro, stage_arg);
        y_est = squeeze(y_est);
        return
    end
    
    if dont_do_bm
        [y_est] = BM4D(z, sigma_psd, pro, stage_arg);
    else
        [y_est, blocks] = BM4D(z, sigma_psd, pro, stage_arg, blockmatches);
    end

end

function pro = convertToBM4DProfile(pro_in)
    pro = BM4DProfile('BM3D');

    pro.filter_strength = pro_in.filter_strength;

    pro.print_info = pro_in.print_info;

    pro.transform_3D_HT_name = pro_in.transform_2D_HT_name;
    pro.transform_3D_Wiener_name = pro_in.transform_2D_Wiener_name;
    pro.transform_NL_name = pro_in.transform_3rd_dim_name;

    pro.Nf = [pro_in.Nf, pro_in.Nf, 1];
    pro.Kin = pro_in.Kin;

    pro.denoise_residual = pro_in.denoise_residual;
    pro.residual_thr = pro_in.residual_thr;
    pro.max_pad_size = pro_in.max_pad_size;

    pro.gamma = pro_in.gamma;

    pro.N1 = [pro_in.N1, pro_in.N1, 1];
    pro.Nstep = [pro_in.Nstep, pro_in.Nstep, 1];

    pro.N2 = pro_in.N2;
    pro.Ns = [floor(pro_in.Ns / 2), floor(pro_in.Ns / 2), 0];
    pro.tau_match = pro_in.tau_match * pro_in.N1^2 / 255^2;

    pro.lambda_thr = pro_in.lambda_thr3D;
    pro.mu2 = pro_in.mu2;

    pro.lambda_thr_re = pro_in.lambda_thr3D_re;
    pro.mu2_re = pro_in.mu2_re;
    pro.beta = pro_in.beta;

    pro.N1_wiener = [pro_in.N1_wiener, pro_in.N1_wiener, 1];
    pro.Nstep_wiener = [pro_in.Nstep_wiener, pro_in.Nstep_wiener, 1];

    pro.N2_wiener = pro_in.N2_wiener;
    pro.Ns_wiener = [floor(pro_in.Ns_wiener / 2), floor(pro_in.Ns_wiener / 2), 0];
    pro.tau_match_wiener = pro_in.tau_match_wiener * pro_in.N1_wiener^2 / 255^2;
    pro.beta_wiener = pro_in.beta_wiener;
    pro.decLevel = pro_in.decLevel;

    pro.set_sharpen(pro_in.sharpen_alpha);
    pro.num_threads = pro_in.num_threads;

end
