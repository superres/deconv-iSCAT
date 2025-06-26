clc; close all;
clear;

a0=readTIF('compares_deconvs_filters.tif'); 
a0=a0(101:300,101:300,:);

titles={'Raw','RL','Wiener','RL Wiener','SF'};
colorss={'k','b','c','m','r'};
lines={'-','--'};

p=[49,160;84,96;69,80;137,70];

y0=floor(size(a0,1)*0.5);
x0=size(a0,2);
figure;
m=3; n=5;
for i=1:n
    a1=getFFT(a0(:,:,i));
    a2=getCrossLineProfile(a1);
    subplot(m,n,i); imshow(a0(:,:,i),[]);  title(titles{i}); colorbar; clim([-0.3,0.2]);
    subplot(m,n,i+n); imshow(a1,[]); colorbar; clim([0,6]);
    subplot(m,n,i+2*n); plot(a2(1,:)*0.058,a2(2,:),'Color',colorss{i}); ylim([0,6]);
end

% figure;
% for i=1:2
%     subplot(1,2,i); hold on;
%     for j=1:3
%         a1=getCrossLine(a0(:,:,3*j-2),p(2*i-1:2*i,:));
%         xx1=1:max(size(a1)); xx1=xx1*0.058;
%         plot(xx1',a1,'Color',colorss{j},'LineStyle',lines{i});
%         mina(j)=min(a1);
%         maxa(j)=max(a1);
%         ranges(j)=max(a1)-min(a1);
%     end
%     legend('Raw','Wiener','RL');
%     hold off;
% end

figure;
for i=1:n
    subplot(1,n,i); hold on;
    for j=1:2
        a1=getCrossLine(a0(:,:,3*i-2),p(2*j-1:2*j,:));
        xx1=1:max(size(a1)); xx1=xx1*0.058;
        plot(xx1',a1,'Color',colorss{i}); 
        xlim([0,4]); ylim([-0.35,0.1]);
    end
    legend('Raw','RL','Wiener','RL Wiener','SF');
    hold off;
end
