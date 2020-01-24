tic
clearvars
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
    % rnames=dir([expdir '/fMRI/timeseries/' timeUnit '/roi/' froidir '/zscore_*.mat']);
end

rnames={rnames.name}';
rnames=strrep(rnames,'.mat','');

sims=zeros(length(rnames),45*45,25);
for ri=1:length(rnames);
    rname=rnames{ri};
    fr =  sprintf('%s/fMRI/timeseries/%s/roi/%s/%s.mat',expdir,timeUnit,froidir,rname);
    load(fr,'gdata');
    
    gdata_z=zscore(gdata,0,2);
    
    if ~isempty(gdata);
        corr_all=zeros(45,45,25);
        
        for si = 1:25;
            subj=subjects{si};
            
            othersi=~ismember(1:25,si);
            others=mean(gdata_z(:,:,othersi),3);
            self=gdata_z(:,:,si);
            
            if strcmp(self_other,'selfother');
                sim(:,:,si)=corr(self,others);
            elseif strcmp(self_other,'selfself');
                sim(:,:,si)=corr(self,self);
            elseif strcmp(self_other,'authorAudience');
                sim(:,:,si)=corr(author,self);
            end;
            
            sim_temp=sim(:,:,si);
            sims(ri,:,si)=sim_temp(:);
        end
    end
    sim=0.5*log((1+sim)./(1-sim));
    
    save([expdir '/fMRI/simMat/roi/' self_other '/' froidir '/' timeUnit '/' srm '/simMatZ_' rname '.mat' ], 'sim','roi_labels','srm');
end

sim=sims;
% save([expdir '/fMRI/simMat/roi/' self_other '/' froidir '/' timeUnit '/' srm '/simMatZ_allROIs.mat' ], 'sim','roi_labels','srm');

toc
beep




