function viewworm(registeredmoviestr)


%findobj('Tag',['Trap' num2str(obj.id)])

if numel(findobj('Tag',registeredmoviestr)) % handle exists already
    h=findobj('Tag',registeredmoviestr);
    
    %h.Children(5).Tag
    
    hp=findobj(h,'Type','Axes');
    
%     hp(1)=findobj(h,'Tag','Axe1');
%     hp(2)=findobj(h,'Tag','Axe2');
%     hp(3)=findobj(h,'Tag','Axe3');
%     hp(4)=findobj(h,'Tag','Axe4');

% hp=h.Children
% 
% if strcmp(h.Children(1).Units,'normalized')
% hp(1:4)=h.Children(1:4);
% else
% hp(1:4)=h.Children(5:8);
% end
    
    him=h.UserData.him;
    updatedisplay(h,him,hp);
else


him=[];

t=[];
t.frame=1;
t.gfp=[];
t.pixtree=[];

h=figure('Tag',registeredmoviestr,'UserData',t);

im=buildimage(h); % returns a structure with all images to be displayed

hp(1)=subplot(2,2,1);
set(hp(1),'Tag','Axe1');

%aaa=hp(1).Tag
axis equal square

him.pixtraining=imshow(im.pixtraining);
hold on
him.pixtraining_raw=imshow(im.pixtraining_raw) ; % overlay raw gfp image to paint
him.pixtraining_raw.AlphaData=0.7;
hold off;
% 
 title(['Frame:' num2str(t.frame) ' - Int:' num2str(h.UserData.intensity)]);
% 
hp(2)=subplot(2,2,2);
set(hp(2),'Tag','Axe2');
axis equal square
him.pixclassif=imshow(im.pixclassif);
% title(['Pixel Classification - Trap: ' num2str(obj.id) '- Frame:' num2str(obj.frame)]);
% 
hp(3)=subplot(2,2,3);
set(hp(3),'Tag','Axe3');
axis equal square


him.trackclassif=imshow(im.trackclassif,[]);
% 
% %hold on
% %him.trackclassif=imshow(im.trackclassif) ; % overlay raw gfp image to paint
% %him.trackclassif.AlphaData=0.1;
% %hold off;
% 
% title(['Tracking Classification - Trap: ' num2str(obj.id) '- Frame:' num2str(obj.frame)]);
% 
 hp(4)=subplot(2,2,4);
 set(hp(4),'Tag','Axe4');
 xlabel('Time (frame)')
 ylabel('Fluo (A.U.)');
 set(hp(4),'FontSize',20);
 
% axis equal square
% him.overlay=imshow(im.overlay);
% title(['Raw image - Trap: ' num2str(obj.id) '- Frame:' num2str(obj.frame) ' - Int:' num2str(obj.intensity)]);
% 
% %himage=imshow(display(obj,hp));

h.Position(3)=1000;
h.Position(4)=1000;

h.UserData.him=him;

obj=[];
h.KeyPressFcn={@changeframe,h,him,hp};

%paint(him.pixtraining,h,hp,obj); % launches the function for pixel training

btnPaint = uicontrol('Style', 'togglebutton', 'String', 'Pixel Train',...
        'Position', [20 20 50 20],...
        'Callback', {@paint,him.pixtraining,h,hp,h.UserData}) ;
    
btnClassify = uicontrol('Style', 'pushbutton', 'String', 'Classify pixels',...
        'Position', [80 20 80 20],...
        'Callback', {@classify,h,him,hp}) ;
    
btnSave = uicontrol('Style', 'pushbutton', 'String', 'Clear training',...
        'Position', [180 20 80 20],...
        'Callback', {@cleartraining,h,him,hp}) ;

btnSave = uicontrol('Style', 'pushbutton', 'String', 'Save classifier',...
        'Position', [280 20 80 20],...
        'Callback', {@saveclassifier,h,registeredmoviestr}) ;
    
btnTrainObjects = uicontrol('Style', 'togglebutton', 'String', 'Select neurons',...
        'Position', [400 20 80 20],...
        'Callback', {@trainobjects,h,hp,obj,him}) ;
    
btnSetFrame = uicontrol('Style', 'edit', 'String', num2str(h.UserData.frame),...
        'Position', [630 20 50 20],...
        'Callback', {@setframe,h,him,hp},'Tag','frametext') ;
    
% btnSetDiv = uicontrol('Style', 'edit', 'String', 'No division',...
%         'Position', [450 20 80 20],...
%         'Callback', {},'Tag','divtext') ;
    
% btnTrainObjects2 = uicontrol('Style', 'pushbutton', 'String', 'Classify objects',...
%         'Position', [320 20 80 20],...
%         'Callback', {@classify,obj,him,hp}) ;

btnTrack = uicontrol('Style', 'pushbutton', 'String', 'Track',...
        'Position', [500 20 80 20],...
        'Callback', {@track,h}) ;
    
btnPlot = uicontrol('Style', 'pushbutton', 'String', 'Plot',...
        'Position', [700 20 80 20],...
        'Callback', {@plotresults,h,hp}) ;
    
btnSave = uicontrol('Style', 'pushbutton', 'String', 'Save',...
        'Position', [800 20 80 20],...
        'Callback', {@saveFigure,h,registeredmoviestr}) ;
    
btnMovie = uicontrol('Style', 'pushbutton', 'String', 'Movie',...
        'Position', [900 20 80 20],...
        'Callback', {@movie,h,registeredmoviestr}) ;
    
end
end

function cleartraining(handle, event, h,him,hp) 

h.UserData.train=zeros(size(h.UserData.train));

updatedisplay(h,him,hp);
end

function saveclassifier(handle, event, h,str) 

tree=h.UserData.tree;
save([str '-classifier.mat'],'tree');
end


function plotresults(handle, event, h,hp)
% this function assigns a training flag to all objects in the middle of the
% image


obj=h.UserData;

fluob=[];
frb=[];
fluor=[];
frr=[];

reverseStr='';
ce=1;

for i=1:size(obj.traintrack,4)
   
    tmp=obj.traintrack(:,:,3,i)>0;
    
    
    if sum(tmp(:))>0
     %   'ok'
       frb=[frb i];
       im=obj.gfp(:,:,i);
       fluob=[fluob sum(im(tmp))];
    end
    
    tmp=obj.traintrack(:,:,1,i)>0;
    
    if sum(tmp(:))>0
       frr=[frr i];
       im=obj.gfp(:,:,i);
       fluor=[fluor sum(im(tmp))];
    end
    
     msg = sprintf('%d / %d Frames tracked', ce , size(obj.traintrack,4) ); %Don't forget this semicolon
%     msg=[msg ' for trap ' obj.id];
     
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
      ce=ce+1;
      
end

axes(hp(4));

cla
plot(frb,fluob,'Color','b'); hold on 
plot(frr,fluor,'Color','r');

curframe=h.UserData.frame;

pb=find(frb==curframe);
pr=find(frr==curframe);

hb=[];
hr=[];

if numel(pb)
hb=plot(frb(pb),fluob(pb),'Color','b','Marker','.','MarkerSize',25); hold on 
end

if numel(pr)
hr=plot(frr(pr),fluor(pr),'Color','r','Marker','.','MarkerSize',25); hold on 
end

    
h.UserData.plot=[];
h.UserData.plot.fluob=fluob;
h.UserData.plot.fluor=fluor;
h.UserData.plot.frb=frb;
h.UserData.plot.frr=frr;
h.UserData.plot.hb=hb;
h.UserData.plot.hr=hr;

end

function track(handle, event, h)
% this function assigns a training flag to all objects in the middle of the
% image


obj=h.UserData;

curframe=h.UserData.frame;

frames=curframe+1:size(obj.traintrack,4);

nred=bwlabel(obj.traintrack(:,:,1,curframe)>0);

pred=regionprops(nred,'Centroid');

xr=0;
yr=0;
xb=0;
yb=0;

if numel(pred)~=0
   xr=pred.Centroid(1);
   yr=pred.Centroid(2);
end

nblue=bwlabel(obj.traintrack(:,:,3,curframe)>0);

pblue=regionprops(nblue,'Centroid');
if numel(pblue)~=0
   xb=pblue.Centroid(1);
   yb=pblue.Centroid(2);
end
    

newxr=xr;
newyr=yr;
newxb=xb;
newyb=yb;

reverseStr='';
ce=1;
for i=frames
    
   % size(obj.traintrack),i
    n=bwlabel(obj.traintrack(:,:,2,i)>0 | obj.traintrack(:,:,1,i)>0 |  obj.traintrack(:,:,3,i)>0);
    p=regionprops(n,'Centroid');
    
    %p.Centroid
    %figure, imshow(n,[])
    
    %return
    
     msg = sprintf('%d / %d Frames tracked', ce , numel(frames) ); %Don't forget this semicolon
%     msg=[msg ' for trap ' obj.id];
     
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
      ce=ce+1;
    
    if numel(p)==0
        continue
    end
    
    dred=[];
    for j=1:numel(p)
       % j
        dred(j)=sqrt((p(j).Centroid(1)-newxr)^2+(p(j).Centroid(2)-newyr)^2);
    end
    
   %d
  
    [dmin ix]=min(dred);
    
    nc=n==ix;
    
    obj.traintrack(:,:,1,i)=255*uint8(nc);
    tmp=obj.traintrack(:,:,2,i); tmp(nc)=0; obj.traintrack(:,:,2,i)=tmp;
    tmp=obj.traintrack(:,:,3,i); tmp(nc)=0; obj.traintrack(:,:,3,i)=tmp;
    
    pred=regionprops(nc,'Centroid');
    newxr=pred(1).Centroid(1);
    newyr=pred(1).Centroid(2);
    %
    
     dblue=[];
    for j=1:numel(p)
       % j
        dblue(j)=sqrt((p(j).Centroid(1)-newxb)^2+(p(j).Centroid(2)-newyb)^2);
    end
    
   %d
  
    [dmin ix]=min(dblue);
    nc=n==ix;
    
    obj.traintrack(:,:,3,i)=255*uint8(nc);
    tmp=obj.traintrack(:,:,2,i); tmp(nc)=0; obj.traintrack(:,:,2,i)=tmp;
    tmp=obj.traintrack(:,:,1,i); tmp(nc)=0; obj.traintrack(:,:,1,i)=tmp;
    
    pblue=regionprops(nc,'Centroid');
    newxb=pblue(1).Centroid(1);
    newyb=pblue(1).Centroid(2);
    
    %if mod(ce-1,50)==0
    
    %end
      
end

fprintf('\n');

h.UserData=obj;
end

function saveFigure(handle,event,h,str)

prompt = 'What is the offset (T0=?) ? ';
x = input(prompt);

fprintf('Saving....\n');

savefig(h,[str '-project.fig']);

hout=figure; 

if isfield(h.UserData,'plot')
frb=h.UserData.plot.frb;
fluob=h.UserData.plot.fluob;
frr=h.UserData.plot.frr;
fluor=h.UserData.plot.fluor;
end

plot(frb,fluob,'Color','b','LineWidth',2); hold on 
plot(frr,fluor,'Color','r','LineWidth',2);
xlabel('Time (frame)')
ylabel('Fluo (A.U.)');
set(gca,'FontSize',20);


title(['T0=' num2str(x)]);

savefig(hout,[str '-fluo.fig']);
saveas(hout,[str '-fluo.pdf']);
close

fprintf('Done !\n');

end

function movie(handle,event,h,str)

% export video
v = VideoWriter([str '-track.mp4'],'MPEG-4');
%v.Colormap = cmap;

bw=h.UserData.traintrack;

bwout=uint8(zeros(size(bw)));

gfp=h.UserData.gfp;
% 
% 
 meangfp=0.95*double(mean(gfp(:)));
 maxgfp=double(meangfp+h.UserData.intensity*(max(gfp(:))-meangfp));

 
%limi=stretchlim(rawgfp,[0.1 obj.intensity]);
%limi(2)=max(limi(2),1.2*limi(1));
gfpout=uint8(zeros(size(gfp)));

for i=1:size(gfp,3)
rawgfp=gfp(:,:,i);
gfpout(:,:,i) = imadjust(rawgfp,[meangfp/65535 maxgfp/65535],[0 1])/256;
end

%figure, imshow(gfpout(:,:,1),[])
%gfpout=255*gfpout;

gfpout=cat(4,gfpout,gfpout,gfpout);

gfpout=permute(gfpout,[1 2 4 3]);



for i=1:size(bw,4)
   tmp=bw(:,:,1,i)>0;
   %tmp=bwperim(tmp);
   bwout(:,:,1,i)=255*tmp;
   
   tmp=bw(:,:,3,i)>0;
   %tmp=bwperim(tmp);
   bwout(:,:,3,i)=255*tmp; 
end

%size(gfpout),size(bwout)
bwout=cat(1,gfpout,bwout);

open(v)
writeVideo(v,bwout);
close(v)


end

function classify(handle,event,h,him,hp)
%h,him,him
  %frames=h.UserData.frame;
  frames=1:size(h.UserData.gfp,3);
  
  %out=pixclassify(h,frames);
  out=pixclassify(h);
  
 % size(h.UserData.classi),size(out.classi)
  h.UserData.classi(:,:,:,frames)=out.classi(:,:,:,frames);
  h.UserData.traintrack(:,:,:,frames)=out.traintrack(:,:,:,frames);
  
  %figure,imshow(h.UserData.traintrack(:,:,:,h.UserData.frame),[]);
  
  updatedisplay(h,him,hp);
  % updates display
end

function im=buildimage(h)

% outputs a structure containing all displayed images
im=[];

frame=h.UserData.frame;

%obj.gfp

if numel(h.UserData.gfp)==0
    load([h.Tag '-registered.mat']);
    h.UserData.gfp=raw;
    
    h.UserData.intensity=0.7;
    h.UserData.train=zeros(size(h.UserData.gfp,1),size(h.UserData.gfp,2),3,size(h.UserData.gfp,3));
    h.UserData.traintrack=zeros(size(h.UserData.gfp,1),size(h.UserData.gfp,2),3,size(h.UserData.gfp,3));
    h.UserData.classi=zeros(size(h.UserData.gfp,1),size(h.UserData.gfp,2),3,size(h.UserData.gfp,3));
   % size(h.UserData.classi)
end


 totgfp=h.UserData.gfp;
% 
% 
 meangfp=0.95*double(mean(totgfp(:)));
 maxgfp=double(meangfp+h.UserData.intensity*(max(totgfp(:))-meangfp));

%limi=stretchlim(rawgfp,[0.1 obj.intensity]);
%limi(2)=max(limi(2),1.2*limi(1));

rawgfp=h.UserData.gfp(:,:,frame);

tracktrain=h.UserData.traintrack(:,:,:,frame);

temp = imadjust(rawgfp,[meangfp/65535 maxgfp/65535],[0 1]);

%temp = imadjust(rawgfp,[limi(1) limi(2)],[0 1]);


%temp=rawgfp;

%figure, imshow(temp,[]);

imrgb=cat(3,temp,temp,temp);

%figure, imshow(imrgb,[])
%size(h.UserData.train)
impaint=h.UserData.train(:,:,:,frame);

%impaint=cat(4,impaint,impaint);

im.pixtraining=impaint; % display RGB image of raw data to paint
im.pixtraining_raw=temp;

%  pix classification results
%size(h.UserData.classi)

%frame

im.pixclassif=h.UserData.classi(:,:,:,frame);

%figure, imshow(h.UserData.classi(:,:,:,frame),[]);

%  tracking classification results

im.trackclassif=tracktrain;

%  overlay

%rawphc = imadjust(rawphc,[meanphc/65535 maxphc/65535],[0 1]);

%imphc= uint16(cat(3,rawphc,rawphc,rawphc));
%imgfp=uint16(zeros(size(imphc)));
%imgfp(:,:,2)=temp;
%im.overlay=imgfp+imphc;
end


function setframe(handle,event,h,him,hp)

frame=str2num(handle.String);


if frame<size(h.UserData.gfp,3) & frame > 0
    
    h.UserData.frame=frame;
    updatedisplay(h,him,hp)
    
   % hl=findobj('Tag',['Trajline' obj.id]);
   % if numel(hl)>0
   % hl.XData=[obj.frame obj.frame];
   % end
end
end


function changeframe(handle,event,h,him,hp)

ok=0;

% if strcmp(event.Key,'uparrow')
% val=str2num(handle.Tag(5:end));
% han=findobj(0,'tag','movi')
% han.trap(val-1).view;
% delete(handle);
% end


if strcmp(event.Key,'rightarrow')
    if h.UserData.frame+1>size(h.UserData.gfp,3)
    return;
    end

    h.UserData.frame=h.UserData.frame+1;
    frame=h.UserData.frame+1;
    ok=1;
end

if strcmp(event.Key,'leftarrow')
    if h.UserData.frame-1<1
    return;
    end

    h.UserData.frame=h.UserData.frame-1;
    frame=h.UserData.frame-1;
    ok=1;
end

if strcmp(event.Key,'uparrow')
    h.UserData.intensity=max(0.01,h.UserData.intensity-0.01);
    ok=1;
end

if strcmp(event.Key,'downarrow')
    h.UserData.intensity=min(1,h.UserData.intensity+0.01);
  %  'ok'
    ok=1;
end

%  if strcmp(event.Key,'r') || strcmp(event.Key,'d' ) % reject divisions
%     if h.UserData.frame>1
%        % 'ok'
%         if obj.div.raw(obj.frame-1)==1 % putative division
%           %  'ok'
%             if obj.div.reject(obj.frame-1)==0 % frame is not rejected
%              %   'ok'
%              
%              if strcmp(event.Key,'r')
%                 obj.div.reject(obj.frame-1)=1;
%              else
%                
%                 obj.div.reject(obj.frame-1)=2;
%               %  aa=obj.div.reject(obj.frame-1)
%              end
%                 
% %                 hreject=findobj('Tag',['Rejectplot' num2str(obj.id)]);
% %                 hraw=findobj('Tag',['Rawplot' num2str(obj.id)]);
% %                 hreject.XData=[hreject.XData obj.frame-1];
% %                 hreject.YData=[hreject.YData hraw.YData(obj.frame-1)];
%             else
%                 obj.div.reject(obj.frame-1)=0;
%                 
%                 %hreject=findobj('Tag',['Rejectplot' num2str(obj.id)]);
%                 %hraw=findobj('Tag',['Rawplot' num2str(obj.id)]);
%                 
% %                 pix=find(hreject.XData==obj.frame-1);
% %                 
% %                 hreject.XData=hreject.XData( setxor(1:length(hreject.XData),pix));
% %                 hreject.YData=hreject.YData( setxor(1:length(hreject.YData),pix));
%             end
%             
%             hr=findobj('Tag',['Traj' obj.id]);
%             if numel(hr)>0
%             delete(hr);
%             end
%             
%             obj.traj;
%              h=findobj('Tag',['Trap' obj.id]);
%             figure(h);
%         end
%     end
%     ok=1;
%  end
    


if ok==1

updatedisplay(h,him,hp)

%  hl=findobj('Tag',['Trajline' obj.id]);
%     if numel(hl)>0
%     hl.XData=[obj.frame obj.frame];
%     end  
end
end

function updatedisplay(h,him,hp)

im=buildimage(h);

%him.overlay.CData=im.overlay;
him.pixtraining.CData=im.pixtraining;
him.pixtraining_raw.CData=im.pixtraining_raw;
him.pixclassif.CData=im.pixclassif;
him.trackclassif.CData=im.trackclassif;



 title(hp(1),['Frame:' num2str(h.UserData.frame) ' - Int:' num2str(h.UserData.intensity)]);
% title(hp(2),['Pixel Classification - Trap: ' obj.id '- Frame:' num2str(obj.frame)]);
% title(hp(3),['Tracking Classification - Trap: ' obj.id '- Frame:' num2str(obj.frame)]);
% title(hp(4),['Raw image - Trap: ' obj.id '- Frame:' num2str(obj.frame) ' - Int:' num2str(obj.intensity)]);

%h=findobj('Tag',['Trap' obj.id]);
    
t=findobj(h,'Tag','frametext');

t.String=num2str(h.UserData.frame);

% t=findobj(h,'Tag','divtext');
% t.String='No divison';

% if numel(obj.div.raw)>0
% if obj.frame>1
% if obj.div.raw(obj.frame-1)==1
%   t.String='Division ?';
% end
% 
% if obj.div.classi(obj.frame-1)==1
%   t.String='Division';
% end
% 
% if obj.div.reject(obj.frame-1)==1
%   t.String='Rejected division';
% end
% 
% if obj.div.reject(obj.frame-1)==2
%   t.String='Dead division';
% end

curframe=h.UserData.frame;

if isfield(h.UserData,'plot')
    
frb=h.UserData.plot.frb;
frr=h.UserData.plot.frr;
fluob=h.UserData.plot.fluob;
fluor=h.UserData.plot.fluor;

pb=find(frb==curframe);
pr=find(frr==curframe);

h.UserData.plot.hb.XData=frb(pb);
h.UserData.plot.hb.YData=fluob(pb);
h.UserData.plot.hr.XData=frr(pr);
h.UserData.plot.hr.YData=fluor(pr);
end



end



function paint(handle,event,imobj,figh,hp,obj)
%global tempimg

   if handle.Value==1
        
        set(figh,'WindowButtonDownFcn',@wbdcb);

        ah = hp(1); %axes('SortMethod','childorder');
       % axis ([1 10 1 10])
       %title('Click and drag')
   else
       figh.WindowButtonDownFcn='';
       figh.Pointer = 'arrow';
       figh.WindowButtonMotionFcn = '';
       figh.WindowButtonUpFcn = '';
       figure(figh); % set focus
   end
        
        function wbdcb(src,callbackdata)
     seltype = src.SelectionType;
   
     %obj
     obj.tree=[];
     
     if strcmp(seltype,'normal')
        src.Pointer = 'circle';
        cp = ah.CurrentPoint;
       
        xinit = cp(1,1);
        yinit = cp(1,2);
        
       % hl = line('XData',xinit,'YData',yinit,...
       % 'Marker','p','color','b');
        src.WindowButtonMotionFcn = {@wbmcb,2};
        src.WindowButtonUpFcn = @wbucb;
     

     end    
     if strcmp(seltype,'alt')
        src.Pointer = 'circle';
        cp = ah.CurrentPoint;
        xinit = cp(1,1);
        yinit = cp(1,2);
       % hl = line('XData',xinit,'YData',yinit,...
       % 'Marker','p','color','b');
        src.WindowButtonMotionFcn = {@wbmcb,1};
        src.WindowButtonUpFcn = @wbucb;

     end    
     
 
        function wbmcb(src,event,colortype)
           cp = ah.CurrentPoint;
           % For R2014a and earlier: 
           % cp = get(ah,'CurrentPoint');
           
           
           xdat = [xinit,cp(1,1)];
           ydat = [yinit,cp(1,2)];
           
           
           % enlarge pixel size
           
           %hl.XData = xdat;
           %hl.YData = ydat;
           % For R2014a and earlier: 
           % set(hl,'XData',xdat);
           % set(hl,'YData',ydat);
           
           % interpolate results
           
           finalX=xdat;
           finalY=ydat;

           in=finalX<=size(obj.train,2) & finalY<=size(obj.train,1) & finalX>0 & finalY>0;
           
           finalX=finalX(in);
           finalY=finalY(in);
           
           if numel(finalX)>=0
           
           imtemp=imobj.CData;
           
          % size(imtemp)
          % int32(finalY)
          % int32(finalX)
           %colortype*ones(1,length(finalX));
           
           linearInd = sub2ind(size(imtemp), int32(finalY), int32(finalX),colortype*ones(1,length(finalX)));
           imobj.CData(linearInd)=255;
           
           swi=3-colortype;
           linearInd = sub2ind(size(imtemp), int32(finalY), int32(finalX),swi*ones(1,length(finalX)));
           imobj.CData(linearInd)=0;

           obj.train(:,:,:,obj.frame)=imobj.CData;
          % hl.XData = xdat;
          % hl.YData = ydat;
          
         % figure, imshow(figh.UserData.train(:,:,:,figh.UserData.frame),[]);
          
           drawnow
           end
        end
   
        function wbucb(src,callbackdata)
           last_seltype = src.SelectionType;
           % For R2014a and earlier: 
           % last_seltype = get(src,'SelectionType');
           %if strcmp(last_seltype,'alt')
              src.Pointer = 'arrow';
              src.WindowButtonMotionFcn = '';
              src.WindowButtonUpFcn = '';
              % For R2014a and earlier:
              % set(src,'Pointer','arrow');
              % set(src,'WindowButtonMotionFcn','');
              % set(src,'WindowButtonUpFcn','');
              
              figh.UserData.train(:,:,:,figh.UserData.frame)=imobj.CData;
             %  figure, imshow(imobj.CData,[])
          % else
          %    return
           %end
        end
        end
end

function trainobjects(handle,event,figh,hp,obj,him)

   if handle.Value==1
        
        set(figh,'WindowButtonDownFcn',@wbdcbobjects);

       % 'ok'
        ah = hp(3); %axes('SortMethod','childorder');
       % axis ([1 10 1 10])
       %title('Click and drag')
   else
       figh.WindowButtonDownFcn='';
       figh.Pointer = 'arrow';
       figh.WindowButtonMotionFcn = '';
       figh.WindowButtonUpFcn = '';
       figure(figh); % set focus
   end
        
        function wbdcbobjects(src,callbackdata)
     seltype = src.SelectionType;
   
     
       % src.Pointer = 'circle';
        cp = ah.CurrentPoint;
       
        xinit = round(cp(1,1));
        yinit = round(cp(1,2));
        
        if xinit>0 && yinit>0 && xinit<=size(him.pixclassif.CData,2) && yinit<=size(him.pixclassif.CData,1)
            
        
        im=him.pixclassif.CData(:,:,2);
        im(im>0)=1;
        im=bwlabel(im);
        
        %figure, imshow(im,[])
        
        val=im(yinit,xinit);
        
        obj2=figh.UserData;
        
        if val>0 % make sure a real object is selescted
        
            if strcmp(seltype,'normal')
                
          tmp3=him.trackclassif.CData(:,:,3);
          tmp1=him.trackclassif.CData(:,:,1);
          
         if sum(tmp3(im==val))==0   %this object was not clicked
                 
        him.trackclassif.CData(:,:,3)=255*(im==val); % select it
        him.trackclassif.CData(:,:,2)=128*(im & ~(im==val) & ~tmp1); % deselect green channel
        him.trackclassif.CData(:,:,1)=255*(tmp1 & ~(im==val)); %deselect red channel
        
        obj2.traintrack(:,:,3,obj2.frame)=him.trackclassif.CData(:,:,3);
        obj2.traintrack(:,:,1,obj2.frame)=him.trackclassif.CData(:,:,1);
        obj2.traintrack(:,:,2,obj2.frame)=him.trackclassif.CData(:,:,2);
        
         else
        him.trackclassif.CData(:,:,3)=0*((im==val)); % deselect it
        him.trackclassif.CData(:,:,2)=128*(im & ~tmp1);
       % him.trackclassif.CData(:,:,2)=128*(im & ~(im==val)); % deselect green channel
       % him.trackclassif.CData(:,:,1)=255*(tmp1 & ~(im==val)); %deselect red channel
        
        obj2.traintrack(:,:,3,obj2.frame)=him.trackclassif.CData(:,:,3);
        obj2.traintrack(:,:,1,obj2.frame)=him.trackclassif.CData(:,:,1);
        obj2.traintrack(:,:,2,obj2.frame)=him.trackclassif.CData(:,:,2);
         end
        
         set(figh,'UserData',obj2);
            end
            
             if strcmp(seltype,'alt')
          tmp3=him.trackclassif.CData(:,:,3);
          tmp1=him.trackclassif.CData(:,:,1);
          
         if sum(tmp1(im==val))==0   %this object was not clicked
                 
        him.trackclassif.CData(:,:,1)=255*(im==val); % select it
        him.trackclassif.CData(:,:,2)=128*(im & ~(im==val) & ~tmp3); % deselect green channel
        him.trackclassif.CData(:,:,3)=255*(tmp3 & ~(im==val)); %deselect red channel
        
        obj2.traintrack(:,:,3,obj2.frame)=him.trackclassif.CData(:,:,3);
        obj2.traintrack(:,:,1,obj2.frame)=him.trackclassif.CData(:,:,1);
        obj2.traintrack(:,:,2,obj2.frame)=him.trackclassif.CData(:,:,2);
        
         else
        him.trackclassif.CData(:,:,1)=0*((im==val)); % deselect it
        him.trackclassif.CData(:,:,2)=128*(im & ~tmp3);
       % him.trackclassif.CData(:,:,2)=128*(im & ~(im==val)); % deselect green channel
       % him.trackclassif.CData(:,:,1)=255*(tmp1 & ~(im==val)); %deselect red channel
        
        obj2.traintrack(:,:,3,obj2.frame)=him.trackclassif.CData(:,:,3);
        obj2.traintrack(:,:,1,obj2.frame)=him.trackclassif.CData(:,:,1);
        obj2.traintrack(:,:,2,obj2.frame)=him.trackclassif.CData(:,:,2);
         end
        
         set(figh,'UserData',obj2);
        
        
            end
            
        end
       % hl = line('XData',xinit,'YData',yinit,...
       % 'Marker','p','color','b');
       % src.WindowButtonMotionFcn = {@wbmcb,2};
       % src.WindowButtonUpFcn = @wbucb;
        
     

     end    
%      if strcmp(seltype,'alt')
%         src.Pointer = 'circle';
%         cp = ah.CurrentPoint;
%         xinit = cp(1,1);
%         yinit = cp(1,2);
%        % hl = line('XData',xinit,'YData',yinit,...
%        % 'Marker','p','color','b');
%         src.WindowButtonMotionFcn = {@wbmcb,1};
%         src.WindowButtonUpFcn = @wbucb;
% 
%      end    
     
 
   
        end
end