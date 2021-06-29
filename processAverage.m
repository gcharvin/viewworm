function processAverage(low,high,offset,outputfilename)

%  low=low(1:2);
%  high=high(1:2);
%  offset=offset(1:2);

% binning frames
binning=2;

nframes = 80 ; % number of frames over which data are averaged to do normalization


% first, lets do basic synchronization 
synchro1=[];
synchro1.high=high; 
synchro1.low=low;
synchro1.lowpoolx=[];
synchro1.lowpooly=[];
synchro1.highpoolx=[];
synchro1.highpooly=[];

for i=1:numel(offset)
synchro1.high(i).data(1,:)=synchro1.high(i).data(1,:)-offset(i);
synchro1.low(i).data(1,:)=synchro1.low(i).data(1,:)-offset(i);

synchro1.highpoolx=[synchro1.highpoolx synchro1.high(i).data(1,:)];
synchro1.highpooly=[synchro1.highpooly synchro1.high(i).data(2,:)];

synchro1.lowpoolx =[synchro1.lowpoolx synchro1.low(i).data(1,:)];
synchro1.lowpooly =[synchro1.lowpooly synchro1.low(i).data(2,:)];
end

[synchro1.highpoolx ix]=sort(synchro1.highpoolx);
synchro1.highpooly=synchro1.highpooly(ix);
[synchro1.lowpoolx ix]=sort(synchro1.lowpoolx);
synchro1.lowpooly=synchro1.lowpooly(ix);

mine=min(min(synchro1.highpoolx),min(synchro1.lowpoolx));
maxe=max(max(synchro1.highpoolx),max(synchro1.lowpoolx));

edges=mine:binning:maxe;

highid=discretize(synchro1.highpoolx,edges);
lowid=discretize(synchro1.lowpoolx,edges);

higharrmean=zeros(1,length(edges)-1);
higharrstd=zeros(1,length(edges)-1);
lowarrmean=zeros(1,length(edges)-1);
lowarrstd=zeros(1,length(edges)-1);


for i=1:length(higharrmean)
    higharrmean(i)=mean(synchro1.highpooly(highid==i));
    higharrstd(i)=std(synchro1.highpooly(highid==i))./sqrt(sum(highid==i));
    
    lowarrmean(i)=mean(synchro1.lowpooly(lowid==i));
    lowarrstd(i)=std(synchro1.lowpooly(lowid==i))./sqrt(sum(lowid==i));
end

x=0.5*(edges(1:end-1)+edges(2:end));

figure('Color','w'), plot(synchro1.highpoolx,synchro1.highpooly,'Color','r','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,higharrmean,higharrstd,{'k-','LineWidth',2}); hold on 
xlabel('Time (frames)');
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);
a=gca;
yli=a.YLim;
line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');

figure('Color','w'), plot(synchro1.lowpoolx,synchro1.lowpooly,'Color','b','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,lowarrmean,lowarrstd,{'k-','LineWidth',2}); hold on 
xlabel('Time (frames)');
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);
a=gca;
yli=a.YLim;
line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');


% now sychronize cells, normalize cells to their initial level 

synchro2=[];
synchro2.high=high; 
synchro2.low=low;
synchro2.lowpoolx=[];
synchro2.lowpooly=[];
synchro2.highpoolx=[];
synchro2.highpooly=[];

for i=1:numel(offset) % synchronization 
synchro2.high(i).data(1,:)=synchro2.high(i).data(1,:)-offset(i);

synchro2.high(i).data(2,:)=synchro2.high(i).data(2,:)./mean(synchro2.high(i).data(2,1:nframes));

synchro2.low(i).data(1,:)=synchro2.low(i).data(1,:)-offset(i);

synchro2.low(i).data(2,:)=synchro2.low(i).data(2,:)./mean(synchro2.low(i).data(2,1:nframes));

end

for i=1:numel(offset) 
synchro2.highpoolx=[synchro2.highpoolx synchro2.high(i).data(1,:)];
synchro2.highpooly=[synchro2.highpooly synchro2.high(i).data(2,:)];

synchro2.lowpoolx =[synchro2.lowpoolx synchro2.low(i).data(1,:)];
synchro2.lowpooly =[synchro2.lowpooly synchro2.low(i).data(2,:)];
end

[synchro2.highpoolx ix]=sort(synchro2.highpoolx);
synchro2.highpooly=synchro2.highpooly(ix);
[synchro2.lowpoolx ix]=sort(synchro2.lowpoolx);
synchro2.lowpooly=synchro2.lowpooly(ix);


mine=min(min(synchro2.highpoolx),min(synchro2.lowpoolx));
maxe=max(max(synchro2.highpoolx),max(synchro2.lowpoolx));

edges=mine:binning:maxe;

highid=discretize(synchro2.highpoolx,edges);
lowid=discretize(synchro2.lowpoolx,edges);

higharrmean=zeros(1,length(edges)-1);
higharrstd=zeros(1,length(edges)-1);
lowarrmean=zeros(1,length(edges)-1);
lowarrstd=zeros(1,length(edges)-1);


for i=1:length(higharrmean)
    higharrmean(i)=mean(synchro2.highpooly(highid==i));
    higharrstd(i)=std(synchro2.highpooly(highid==i))./sqrt(sum(highid==i));
    
    lowarrmean(i)=mean(synchro2.lowpooly(lowid==i));
    lowarrstd(i)=std(synchro2.lowpooly(lowid==i))./sqrt(sum(lowid==i));
end

x=0.5*(edges(1:end-1)+edges(2:end));

figure('Position',[100 100 500 700],'Color','w')
subplot(2,1,1);

plot(synchro2.highpoolx,synchro2.highpooly,'Color','r','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,higharrmean,higharrstd,{'k-','LineWidth',2}); hold on 
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);

a=gca;
yli=a.YLim;

line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
line([0 100], [max(higharrmean) max(higharrmean)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
text(100, max(higharrmean), num2str(max(higharrmean),3),'FontSize',24);

ylim([0 3]);

subplot(2,1,2);

plot(synchro2.lowpoolx,synchro2.lowpooly,'Color','b','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,lowarrmean,lowarrstd,{'k-','LineWidth',2}); hold on 
xlabel('Time (frames)');
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);

a=gca;
yli=a.YLim;
line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
line([0 100], [max(lowarrmean) max(lowarrmean)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
text(100, max(lowarrmean), num2str(max(lowarrmean),3),'FontSize',24);


ylim([0 3]);


% now sychronize cells, normalize cells to their initial level 



synchro3=[];
synchro3.high=high; 
synchro3.low=low;
synchro3.lowpoolx=[];
synchro3.lowpooly=[];
synchro3.highpoolx=[];
synchro3.highpooly=[];

for i=1:numel(offset) % synchronization 
synchro3.high(i).data(1,:)=synchro3.high(i).data(1,:)-offset(i);

temp=mean(synchro3.high(i).data(2,1:nframes));

synchro3.high(i).data(2,:)=synchro3.high(i).data(2,:)./mean(synchro3.high(i).data(2,1:nframes));

synchro3.low(i).data(1,:)=synchro3.low(i).data(1,:)-offset(i);

synchro3.low(i).data(2,:)=synchro3.low(i).data(2,:)./temp;
end

for i=1:numel(offset) 
synchro3.highpoolx=[synchro3.highpoolx synchro3.high(i).data(1,:)];
synchro3.highpooly=[synchro3.highpooly synchro3.high(i).data(2,:)];

synchro3.lowpoolx =[synchro3.lowpoolx synchro3.low(i).data(1,:)];
synchro3.lowpooly =[synchro3.lowpooly synchro3.low(i).data(2,:)];
end

[synchro3.highpoolx ix]=sort(synchro3.highpoolx);
synchro3.highpooly=synchro3.highpooly(ix);
[synchro3.lowpoolx ix]=sort(synchro3.lowpoolx);
synchro3.lowpooly=synchro3.lowpooly(ix);


mine=min(min(synchro3.highpoolx),min(synchro3.lowpoolx));
maxe=max(max(synchro3.highpoolx),max(synchro3.lowpoolx));

edges=mine:binning:maxe;

highid=discretize(synchro3.highpoolx,edges);
lowid=discretize(synchro3.lowpoolx,edges);

higharrmean=zeros(1,length(edges)-1);
higharrstd=zeros(1,length(edges)-1);
lowarrmean=zeros(1,length(edges)-1);
lowarrstd=zeros(1,length(edges)-1);


for i=1:length(higharrmean)
    higharrmean(i)=mean(synchro3.highpooly(highid==i));
    higharrstd(i)=std(synchro3.highpooly(highid==i))./sqrt(sum(highid==i));
    
    lowarrmean(i)=mean(synchro3.lowpooly(lowid==i));
    lowarrstd(i)=std(synchro3.lowpooly(lowid==i))./sqrt(sum(lowid==i));
end

x=0.5*(edges(1:end-1)+edges(2:end));

figure('Position',[100 100 500 700],'Color','w')
subplot(2,1,1);

plot(synchro3.highpoolx,synchro3.highpooly,'Color','r','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,higharrmean,higharrstd,{'k-','LineWidth',2}); hold on 
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);

a=gca;
yli=a.YLim;


line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
line([0 100], [max(higharrmean(x>0)) max(higharrmean(x>0))],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
text(100, max(higharrmean(x>0)), num2str(max(higharrmean(x>0)),3),'FontSize',24);

ylim([0 3]);

subplot(2,1,2);

plot(synchro3.lowpoolx,synchro3.lowpooly,'Color','b','LineStyle','none','Marker','.'); hold on;
shadedErrorBar(x,lowarrmean,lowarrstd,{'k-','LineWidth',2}); hold on 
xlabel('Time (frames)');
ylabel('Average fluo level (A.U.)')

set(gca,'FontSize',24);

a=gca;
yli=a.YLim;
line([0 0], [yli(1) yli(2)],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
line([0 100], [max(lowarrmean(x>0)) max(lowarrmean(x>0))],'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
text(100, max(lowarrmean(x>0)), num2str(max(lowarrmean(x>0)),3),'FontSize',24);


%ylim([0 3]);

% 
% shadedErrorBar(Xp,arrmeanhigh,arrstdhigh,{'b-o','markerfacecolor','b'});
% 
% xlabel('Time (frames)');
% ylabel('Average fluo level (A.U.)')
% 
% xlim([-100 200])
% ylim([0 120])


% ok, now synchronize all the traces
% 
% maxe=0;
% for i=1:numel(offset)
%     maxe=max(size(low(i).data,2),maxe);
% end
% 
% arrlow=-ones(length(low),2*maxe);
% arrhigh=-ones(length(high),2*maxe);
% 
% for i=1:length(low)
%    
%     pixmax=size(low(i).data,2);
%     
%    % size(maxe+1-offset(i):maxe+pixmax-offset(i))
%     
%     arrlow(i,maxe+1-offset(i):maxe+pixmax-offset(i))=low(i).data(2,:);
%     arrhigh(i,maxe+1-offset(i):maxe+pixmax-offset(i))=high(i).data(2,:);
%     
% end

%arrlow
%now average curves and propose binning !! 

% arrmeanlow=[];
% arrmeanhigh=[];
% 
% arrstdlow=[];
% arrstdhigh=[];
% 
% 
% edges=linspace(0,2*maxe,2*maxe/binning);
% 
% x=discretize(1:2*maxe,edges); 
% 
% tim=1:2*maxe;
% 
% %Xp=zeros(1,max(x));
% Xp=[];
% 
% cc=1;
% 
% for j=1:max(x)
%     pix=find(x==j);
%    % j
%     %aa=mean(tim(pix))
%     tmp=arrlow(:,pix);
%     
%     ok=0;
%     for k=pix
%        if sum(arrlow(:,k)>0)>1
%            ok=1;
%        end
%     end
%     
%     tmp=tmp(tmp>0);
%     
%     if ok==1 % 2 objects necessary for averaging
%         
%         Xp(cc)=mean(tim(pix));
%         
%         arrmeanlow(cc)=mean(tmp);
%         arrstdlow(cc)=std(tmp)./sqrt(length(tmp));
%         
%         cc=cc+1;
%     end
%     
%     %arrmeanlow
% end
% 
% Xp=[];
% 
% cc=1;
% 
% for j=1:max(x)
%     pix=find(x==j);
%    % j
%     %aa=mean(tim(pix))
%     tmp=arrhigh(:,pix);
%     
%     ok=0;
%     for k=pix
%        if sum(arrhigh(:,k)>0)>1
%            ok=1;
%        end
%     end
%     
%     tmp=tmp(tmp>0);
%     
%     if ok==1 % 2 objects necessary for averaging
%         
%         Xp(cc)=mean(tim(pix));
%         
%         arrmeanhigh(cc)=mean(tmp);
%         arrstdhigh(cc)=std(tmp)./sqrt(length(tmp));
%         
%         cc=cc+1;
%     end
%     
%     %arrmeanlow
% end
% 
% 
% Xp=Xp-maxe;
% 
% %Xp,arrmeanlow,arrstdlow
% Xp2=Xp(1:length(arrmeanlow));



% figure, shadedErrorBar(Xp2,arrmeanlow,arrstdlow,{'r-o','markerfacecolor','r'}); hold on 
% 
% shadedErrorBar(Xp,arrmeanhigh,arrstdhigh,{'b-o','markerfacecolor','b'});
% 
% xlabel('Time (frames)');
% ylabel('Average fluo level (A.U.)')
% 
% xlim([-100 200])
% ylim([0 120])
% 
% set(gca,'FontSize',24);
% 
% saveas(gcf,[foldername '.pdf']);
% savefig(gcf,[foldername '.fig']);
% 
