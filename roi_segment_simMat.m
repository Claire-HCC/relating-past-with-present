%% this script computes segment x segment pattern similarity matrix
% set parameters;
loc='cluster';
set_parameters;
timeUnit='segment' ;
self_other='selfother';
froidir='shen';
srm='SRM';
if strmatch(srm,'SRM');
    rnames=dir([expdir '/fMRI/timeseries/' timeUnit '/roi/' froidir '/*_pauseAudResid_zscore_srm.mat']);
else
    rnames=dir([expdir '/fMRI/timeseries/' timeUnit '/roi/' froidir '/*_pauseAudResid_zscore.mat']);
end

% get the names of ROIs
rnames={rnames.name}';
rnames=strrep(rnames,'.mat','');
roi_labels=regexp(rnames,'[0-9]*','match');
roi_labels=str2num(cell2mat([roi_labels{:}]'));

% create empty segment x segment x subject pattern similarity matrix
sims=zeros(length(rnames),45*45,25);

% loop over ROIs
for ri=1:length(rnames);
    rname=rnames{ri};
    fr =  sprintf('%s/fMRI/timeseries/%s/roi/%s/%s.mat',expdir,timeUnit,froidir,rname);
    load(fr,'gdata');
    
    % normalize the time series across time
    gdata_z=zscore(gdata,0,2);
    
    if ~isempty(gdata);
        corr_all=zeros(45,45,25);
        
        % loop over subjects
        for si = 1:25;
            subj=subjects{si};
            
            othersi=~ismember(1:25,si);
            others=mean(gdata_z(:,:,othersi),3);
            self=gdata_z(:,:,si);
            
            % inter- or intra-subject analysis
            if strcmp(self_other,'selfother');
                sim(:,:,si)=corr(self,others);
                
            elseif strcmp(self_other,'selfself');
                sim(:,:,si)=corr(self,self);
            end;
            
            sim_temp=sim(:,:,si);
            sims(ri,:,si)=sim_temp(:);
        end
    end
    
    % fisher's z transformation
    sim=0.5*log((1+sim)./(1-sim));
    
    save([expdir '/fMRI/simMat/roi/' self_other '/' froidir '/' timeUnit '/simMatZ_' rname '.mat' ], 'sim','roi_labels','srm');
end





