clear all

loc='myPC';
set_parameters;

self_other='selfother';
modelName='A12B12C';
froidir='shen';
srm='SRM';
bnames={{'A1A1','A2A2','B1B1','B2B2','A1A2','B1B2',   'A1B1','A2B2','A1B2','A2B1'},...
    {'A2A2','B2B2','A1B1',   'A1A1','B1B1','A2B2'}};
conrasts={[1/(6)*ones(1,6)  -1/(4)*ones(1,4)],...
    [1/2 1/2 1 -1/2 -1/2 -1]};
connames={'storyline','storylineXtime'};

load([expdir '/fMRI/simMat_orig/wholeBrain/' self_other '/segment/subjects_glm/' modelName '/hp.mat'],'xnames','hp_array');

if strmatch(srm,'SRM');
    rnames=dir(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/simMatZ_roi*zscore_srm.mat',expdir,self_other,froidir));
else
    rnames=dir(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/simMatZ_roi*zscore.mat',expdir,self_other,froidir));
end

rnames={rnames(:).name}';
rnames=strrep(cellstr(rnames),'.mat','');
roi_labels=regexp(strrep(rnames,'simMat_zscore_roi_',''),'[0-9]*','match');
roi_labels=str2num(cell2mat([roi_labels{:}]'));

for ci=1:length(bnames);
    conname=connames{ci};
    con=zeros(length(roi_labels),25);
    p=nan(length(roi_labels),1);
    
    for ri=1:length(rnames);
        roi_label=roi_labels(ri);
        rname=rnames{ri};
        
        % 45 segments x 45 segments x 25 subjects pattern similarity matrix computed with leave-one-subject-out-method
        load(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/%s.mat',expdir,self_other,froidir,rname),'sim');
        
        for subj=1:25;
            sim_subj=sim(:,:,subj);
            
            for bi=1:length(bnames{ci});
                bname=bnames{ci}{bi};
                b=hp_array(:,find(strcmp(xnames,bname)));
                b=reshape(b,45,45);
                
                con(ri,subj)= con(ri,subj)+ conrasts{ci}(bi)*mean(sim_subj(b==1));
            end
        end
    end
    
    [~,p]=ttest(con',0,'tail','right');
    m=mean(con,2);
    
    sig_fdr=fdr0(p,0.05);
    sig_fwe=(p<(0.05/length(roi_labels)));
    
    % storyline x time effect, within regions showing storyline effect
    if strcmp(conname,'storylineXtime');
        sig_fwe_storyline=load(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s.mat',expdir,self_other,froidir,srm,'storyline'),'sig_fwe');
        sig_fwe_storyline=sig_fwe_storyline.sig_fwe;
        sig_fwe_withinSotryline= zeros(length(p),1);
        sig_fwe_withinSotryline(sig_fwe_storyline==1)=(p(sig_fwe_storyline==1)<(0.05/sum(sig_fwe_storyline)));
        
        save(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s.mat',expdir,self_other,froidir,srm,conname),'p','rnames','con','m','sig_fdr','sig_fwe','roi_labels','sig_fwe_withinSotryline');
        nii=roiTable2wholeBrainNii_shen([roi_labels(sig_fwe_withinSotryline==1 ) m(sig_fwe_withinSotryline==1)]);
        save_nii(nii,sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_fwe_withinSotryline.nii',expdir,self_other,froidir,srm,conname));
        
    else
        save(sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s.mat',expdir,self_other,froidir,srm,conname),'p','rnames','con','m','sig_fdr','sig_fwe','roi_labels');
    end
    
    nii=roiTable2wholeBrainNii_shen([roi_labels m]);
    save_nii(nii,sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s.nii',expdir,self_other,froidir,srm,conname))
    
    nii=roiTable2wholeBrainNii_shen([roi_labels(sig_fdr==1) m(sig_fdr==1)]);
    save_nii(nii,sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_fdr.nii',expdir,self_other,froidir,srm,conname))
    
    nii=roiTable2wholeBrainNii_shen([roi_labels(sig_fwe==1 ) m(sig_fwe==1 )]);
    save_nii(nii,sprintf('%s/fMRI/simMat/roi/%s/%s/segment/ttestSubj/%s/%s_fwe.nii',expdir,self_other,froidir,srm,conname))
end

