%% this script runs group t-test for the storyline & storyline x time effect 
% set parameters;
loc='cluster';
set_parameters;
self_others={'selfother','selfself'};
froidir='shen';
srms={'SRM','noSRM'};

% set contrasts
connames={'storyline','storylineXtime_2timeBins'};
load([expdir '/scripts_claire/segment_simMat_contrasts.mat'],'xnames','hp_array');

% inter or intra-subject analyses
for sf=1;
    self_other=self_others{sf};
    
    % SRM or not
    for sm=1;
        srm=srms{sm};
        
        % extract ROI names
        rnames=dir(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/simMatZ_roi*zscore_srm.mat',expdir,self_other,froidir));
        rnames={rnames(:).name}';
        rnames=strrep(cellstr(rnames),'.mat','');
        roi_labels=regexp(strrep(rnames,'simMat_zscore_roi_',''),'[0-9]*','match');
        roi_labels=str2num(cell2mat([roi_labels{:}]'));
        
        if  ismember(srm,{'noSRM'});
            rnames=strrep(rnames,'_srm','');
        end
        
        % loop over contrasts
        for ci=1:length(connames);
            
            % get the contrast 
            conname=connames{ci};
            con=zeros(length(rnames),25);
            contrast=hp_array(:,find(strcmp(xnames,conname)));
            
            % apply the contrast to each ROI
            for ri=1:length(rnames);
                roi_label=roi_labels(ri);
                rname=rnames{ri};
                
                load(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/%s.mat',expdir,self_other,froidir,rname),'sim');
                
                % loop over subjects
                for si=1:25;
                    sim_temp=sim(:,:,si);
                    con(ri,si)=nansum(sim_temp(:).*contrast(:));
                end
            end
        end
        
        % group t-test
        [~,p]=ttest(con',0,'tail','right');
        m=mean(con,2);
        
        % FWE correction
        sig_fwe=(p<(0.05/length(roi_labels)));
        
        % threshold storyline x time effect within regions showing significant storyline effect
        if ~ismember(conname,{'storyline'});
            sig_fwe_storyline=load(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_stats.mat',expdir,self_other,froidir,srm,'storyline'),'sig_fwe');
            sig_fwe_storyline=sig_fwe_storyline.sig_fwe;
            
            sig_fwe_withinSotrylineFwe= zeros(length(p),1);
            sig_fwe_withinSotrylineFwe(sig_fwe_storyline==1)=(p(sig_fwe_storyline==1)<(0.05/sum(sig_fwe_storyline)));
            
            save(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_stats.mat',expdir,self_other,froidir,srm,conname),'p','rnames','con','m','sig_fwe','roi_labels','sig_fwe_withinSotrylineFwe');
            
        else
            save(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_stats.mat',expdir,self_other,froidir,srm,conname),'p','rnames','con','m','sig_fwe','roi_labels');
        end
        
    end
end

