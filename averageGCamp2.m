function [low high offset]=averageGCamp(foldername)
% goes over all cuvres in folder to extract neuron data

if ~isfolder(foldername) % tests if folder is a right foldername. If not, return...
    return;
end

list=dir([foldername '/*.fig']); % lists all .fig files in folder


cc=1;

binning=2;


low=[] ; % low neuron 
high=[] ; % low neuron 

low.data=[];
high.data=[];

offset=[];

cc=1; % counter 




for i=1:numel(list) % loops on .fig files in folder 
    
   
    openfig(fullfile(foldername,list(i).name)); % open figure
    
    pause(0.2);
    
    hfig=gca; % get the handle of the axis on current figure; 
    t=hfig.Title.String;
    
    if strcmp(class(t),'cell') % check the class of object and modify it if necessary
        t=t{1};
    end
    
    ix=strfind(t,'T0='); % finds the offset time in object
    
    offset(cc)=str2num(t(ix+3:end));
    
    for j=1:numel(hfig.Children) % look for the two fluo plots
       
        if hfig.Children(j).Color(3)==1 % this is the red plot (RGB color mode) --> low neuron 
            low(cc).data=[hfig.Children(j).XData ; hfig.Children(j).YData]; % gets the X and Y coordinates of fluo points
        end
        if hfig.Children(j).Color(1)==1 % this is the blue plot (RGB color mode) . --> high neuron
            high(cc).data=[hfig.Children(j).XData ; hfig.Children(j).YData]; % gets the X and Y coordinates of fluo points
        end
        
        
    end
    
   % ww=hfig.Children(j).YData(end)
    close % close current figure before loading another one
    cc=cc+1;
end
