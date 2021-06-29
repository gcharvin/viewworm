function readworm(filename,frames)

if nargin==0
filename='worm_01_AL.tif';
end

[im map]=imread(filename,1);

%im=im(:,1:784);


s=imfinfo(filename);

if nargin==1
%mov=mov(:,:,:,frames); % cut 2 frames
end

%mov=uint8(zeros(size(im,1),size(im,2),1,50));
%length(s)

init=1;

% cropping image
 
hfig=figure;
imshow(im,[]);
 
title('Please double click on ROI when done !');
hax=gca;
windo=[100 100 50 50];
roi = imrect(hax,windo);
poly = wait(roi);
close;
 
minex=uint16(round(min(poly(1))));
maxex=uint16(round(max(poly(1)+poly(3)-1)));
 
miney=uint16(round(min(poly(2))));
maxey=uint16(round(max(poly(2)+poly(4)-1)));
 
%mov=mov(miney:maxey,minex:maxex,:,:);

mov=uint8(zeros(maxey-miney+1,maxex-minex+1,1,length(s)));
raw=uint16(zeros(maxey-miney+1,maxex-minex+1,1,length(s)));


imref=imread(filename,init);
%imref=imref(:,1:784);
%imref=imgradient(imref);
%figure, imshow(imref,[]);

%imrefraw=imref;
imrefraw=imref(miney:maxey,minex:maxex);

maxeref=65535*stretchlim(imref,[0.01 0.999999]);
%maxeref=65535*stretchlim(imref,[0 0.2]);

imref=imadjust(imref,[maxeref(1)/65536 maxeref(2)/65535],[0 1]);
imref8 = uint8(imref/256);

imref8=imref(miney:maxey,minex:maxex);
%imref8=imref;

cc=1;

  mov(:,:,1,cc)=imref8;
  raw(:,:,1,cc)=imrefraw;

cc=2;
%h=figure;
%hi=imshowpair(h,imref8,moved);


for i=init+1:init+size(mov,4)-1
   
  im=imread(filename,i);
  
  raw16=im;
  im=imadjust(im,[maxeref(1)/65536 maxeref(2)/65536],[]);
  
  img8 = uint8(im/256);

 c = normxcorr2(imrefraw,raw16);
%  size(imrefraw2)
%  size(raw16)
%  size(c)
%  minex
%  miney
 
% c(size(img8,1), size(img8,2))=c(size(img8,1)-1, size(img8,2)); % remove
 %white dot in the middle of image
  
 %figure, imshow(imref8,[]); 
 %figure, imshow(img8,[]);
 
 %figure, imshowpair(imref8,img8);

 %figure, imshow(c,[])
 
 %cbw=logical(zeros(size(c)));
 %cbw(round(size(cbw,1)/2-size(img8,1)/4):round(size(cbw,1)/2+size(img8,1)/4),round(size(cbw,2)/2-size(img8,2)/2):round(size(cbw,2)/2+size(img8,2)/2))=1;
 
% figure, imshow(cbw,[]);

 %c(~cbw)=0;
 %figure, imshow(c,[]);
 %return;
 
% [max_num,max_idx] = max(c(:));

%[row col]=ind2sub(size(c),max_idx);



%pos(1)=col;
%pos(2)=row;

[ypeak, xpeak] = find(c==max(c(:)));
%xpeak, ypeak

yoffset = ypeak-size(imrefraw,1);
xoffset = xpeak-size(imrefraw,2);


mov(:,:,1,cc)=img8(yoffset+1:yoffset+maxey-miney+1,xoffset+1:xoffset+maxex-minex+1);
raw(:,:,cc)=raw16(yoffset+1:yoffset+maxey-miney+1,xoffset+1:xoffset+maxex-minex+1);
  
%figure, imshow(raw(:,:,cc),[]);
%max_num
%pos


 %figure, surf(c), shading flat
 
 %figure, imshow(t,[]);
 
  %thr=0.28; % threshold for detected peaks

  %BW = im2bw(c,thr);

  %pp = regionprops(BW,'centroid');
  %pos = round(cat(1, pp.Centroid));

  %figure, imshow(BW,[]);
  
%   optimizer.InitialRadius = 0.01;
% optimizer.Epsilon = 1.5e-4;
% optimizer.GrowthFactor = 1.05;
% optimizer.MaximumIterations = 1000;

%size(img8)
  
  %moved=imregister(img8,imref8,'translation',optimizer,metric);
  %pos
  
 % -(pos(1)-size(img8,2))
 %-(pos(2)-size(img8,1))
  

  
%  moved=circshift(img8,-(pos(1)-size(img8,2)),2);
%  moved=circshift(moved,-(pos(2)-size(img8,1)),1);
  
%  rawmoved=circshift(raw16,-(pos(1)-size(raw16,2)),2);
%  rawmoved=circshift(rawmoved,-(pos(2)-size(raw16,1)),1);
  
%-(pos(1)-size(img8,2))
%-(pos(2)-size(img8,1))

 % figure, imshowpair(imref8,moved);
  
  %return;
  %pause(0.2);
  
 % figure, imshow(moved,[]);
  
 % return;
  %pause
  %close
  
  
  
  
  %mov(:,:,1,cc)=img8;
  
  %if mod(cc,50)==0
  %imref8=moved;
  %end
  
 % tic;
%   c = normxcorr2(imref8,img8);
%  % toc;
%  
%  thr=0.5; % threshold for detected peaks
% BW = im2bw(c,thr);
% 
% pp = regionprops(BW,'centroid');
% pos = round(cat(1, pp.Centroid))
% 
% 
%figure, imshowpair(moved,imref8)
% figure, surf(c), shading flat
%pause
%close

if mod(cc,50)==0
  fprintf('\n');
end

fprintf('.');
cc=cc+1;
end

%size(mov)
%size(raw)
%minex,maxex,miney,maxey

%aa=maxey+1-miney

mov=mov(1:maxey-miney+1,1:maxex-minex+1,:,:);
raw=raw(1:maxey-miney+1,1:maxex-minex+1,:);
%size(mov)
%size(raw)

fprintf('\n');

[pth fle ext ]=fileparts(filename);
% export registered data

save([fle '-registered.mat'],'raw'); % this image is raw data no filtering

% export registered movie 



%cmap=colormap(gray(256));
warning off all
v = VideoWriter([fle '-registered.mp4'],'MPEG-4');


%v.Colormap = cmap;

open(v)

writeVideo(v,mov);
close(v)
warning on all



