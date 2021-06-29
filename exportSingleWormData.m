function exportSingleWormData(low,high,offset,outputfilename)

% binning frames
binning=1; % KEEP BINNING TO & HERE §§§§

nframes = 80 ; % number of frames over which data are averaged to do normalization

interva=[0 50]; % temporal window compute the max of the curve , units are 

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

% for i=1:numel(offset) 
% synchro2.highpoolx=[synchro2.highpoolx synchro2.high(i).data(1,:)];
% synchro2.highpooly=[synchro2.highpooly synchro2.high(i).data(2,:)];
% 
% synchro2.lowpoolx =[synchro2.lowpoolx synchro2.low(i).data(1,:)];
% synchro2.lowpooly =[synchro2.lowpooly synchro2.low(i).data(2,:)];
% end
% 
% [synchro2.highpoolx ix]=sort(synchro2.highpoolx);
% synchro2.highpooly=synchro2.highpooly(ix);
% [synchro2.lowpoolx ix]=sort(synchro2.lowpoolx);
% synchro2.lowpooly=synchro2.lowpooly(ix);
% 
% 
% mine=min(min(synchro2.highpoolx),min(synchro2.lowpoolx));
% maxe=max(max(synchro2.highpoolx),max(synchro2.lowpoolx));
% 
% edges=mine:binning:maxe;
% x=0.5*(edges(1:end-1)+edges(2:end));
% 
% highid=discretize(synchro2.highpoolx,edges);
% lowid=discretize(synchro2.lowpoolx,edges);
% 
% higharrmean=zeros(1,length(edges)-1);
% higharrstd=zeros(1,length(edges)-1);
% lowarrmean=zeros(1,length(edges)-1);
% lowarrstd=zeros(1,length(edges)-1);
% 
% 
% for i=1:length(higharrmean)
%     higharrmean(i)=mean(synchro2.highpooly(highid==i));
%     higharrstd(i)=std(synchro2.highpooly(highid==i))./sqrt(sum(highid==i));
%     
%     lowarrmean(i)=mean(synchro2.lowpooly(lowid==i));
%     lowarrstd(i)=std(synchro2.lowpooly(lowid==i))./sqrt(sum(lowid==i));
% end

% calculate max around the time 0 for each single curve

wormhigh=[];
wormlow=[];

for i=1:numel(offset) % synchronization 
pix=synchro2.high(i).data(1,:);
pix=pix>=interva(1) & pix<interva(2);

tmp=synchro2.high(i).data(2,:);
[m ix]=max(tmp(pix));
wormhigh=[wormhigh m];

pix=synchro2.low(i).data(1,:);
pix=pix>=interva(1) & pix<interva(2);

tmp=synchro2.low(i).data(2,:);
[m ix]=max(tmp(pix));
wormlow=[wormlow m];

end

%pix=x>=0 & x<50;


% now extract stitstics for 

disp('Exporting single worm max data to xls file')




% [m ix]= max(higharrmean(pix));
% ix=ix+find(pix==1,1,'first')
% highid
% wormhigh=synchro2.highpooly(highid==ix)
% 
% [m ix]= max(lowarrmean(pix));
% ix=ix+find(pix==1,1,'first')
% lowid
% wormlow=synchro2.highpooly(lowid==ix)

T=table(wormhigh',wormlow','VariableNames',{'High','Low'})

writetable(T,[outputfilename '.xls'],'Sheet',1);


