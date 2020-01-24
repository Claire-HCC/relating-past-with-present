% motif used to be called as memory plant (mp)....
clear all
close all
tic
% loc='cluster';
set_parameters;
froidir='shen';
self_other='selfother';
srm='noSRM';
load(sprintf('%s/fMRI/simMat/roi/%s/%s/mp/%s/mp_data.mat',expdir,self_other,froidir,srm),'mp_data');

mp_d=sortrows(mp_data,{'subj','mp_token_id_C','mp_token_id_AB'});
mp_d(:,{'mp_C','mp_AB'})=[];
mp_d(:,'mp_id_AB_shuffled')=table(0);

rnames=regexp(mp_d.Properties.VariableNames,'.*_pauseAudResid_zscore','match');
roi_coli=find(~cellfun(@isempty,rnames));
rnames=[rnames{roi_coli}]';
roi_labels=regexp(rnames,'[0-9]*','match');
roi_labels=str2num(cell2mat([roi_labels{:}]'));

% remove similarity between unrelated mp
mp_d(mp_d.cu==0,:)=[];
mp_d_CtokenXsubj=grpstats([mp_d(:,{'subj','mp_token_id_C','cu'}) mp_d(:,roi_coli)],{'subj','mp_token_id_C','cu'});

mp_cu2=grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu==2,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu==2))';
mp_cu1=grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu==1,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu==1))';

iters=10000;
for iter=1:iters;
    
    for ab=[13 23];
        % get all the C mp_token_id in one story line
        token_id_C=table2array(unique(mp_d(mp_d.mp_ab_C==ab & mp_d.subj==1,'mp_token_id_C')));
        
        % get all the mp_id_AB in the same story line
        ids_AB=table2array(mp_d(mp_d.mp_token_id_C==token_id_C(1,1) & mp_d.subj==1,'mp_id_AB'));
        
        % shuffle the mp_id_AB in the same way for all subjects X mp_token_id_C
        shuff1 = randperm(length(ids_AB));
        mp_d.mp_id_AB_shuffled(mp_d.mp_ab_C==ab)=repmat(ids_AB(shuff1),25*length(token_id_C),1);
    end
    
    mp_d.cu_shuffled(:,1)=1;
    mp_d.cu_shuffled(mp_d.mp_id_C==mp_d.mp_id_AB_shuffled,1)=2;
    
    mp_d_CtokenXsubj=grpstats([mp_d(:,{'subj','mp_token_id_C','cu_shuffled'}) mp_d(:,roi_coli)],{'subj','mp_token_id_C','cu_shuffled'});
    
    null_cu2(:,iter)=mean(grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu_shuffled==2,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu_shuffled==2)));
    null_cu1(:,iter)=mean(grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu_shuffled==1,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu_shuffled==1)));
    
end

toc
mpEffect=mp_cu2-mp_cu1;
mpEffect_m=mean(mpEffect,2);
null_mpEffect=null_cu2-null_cu1;

save([expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_stats.mat'],'mp_cu2','mp_cu1','null_cu2','null_cu1','iters','srm','mpEffect','rnames');

nii=roiTable2wholeBrainNii_shen([roi_labels mpEffect_m]);
save_nii(nii,[expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_simz.nii']);

p=sum(min(null_mpEffect)<mpEffect_m,2)/iters;
sig_mp_qmax_neg=p<0.05;
nii=roiTable2wholeBrainNii_shen([roi_labels(sig_mp_qmax_neg) mpEffect_m(sig_mp_qmax_neg)]);
save_nii(nii,[expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_simz_qmax05_neg.nii']);

p=sum(max(null_mpEffect)>mpEffect_m,2)/iters;
sig_mp_qmax=p<0.05;
nii=roiTable2wholeBrainNii_shen([roi_labels(sig_mp_qmax) mpEffect_m(sig_mp_qmax)]);
save_nii(nii,[expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_simz_qmax05_pos.nii']);

p=sum(null_mpEffect>mpEffect_m,2)/iters;
sig_mp_qfwe=p<(0.05/length(roi_labels));
nii=roiTable2wholeBrainNii_shen([roi_labels(sig_mp_qfwe) mpEffect_m(sig_mp_qfwe)]);
save_nii(nii,[expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_simz_qfwe05_pos.nii']);

p=sum(null_mpEffect>mpEffect_m,2)/iters;
sig_mp_qfdr=(fdr0(p,0.05)==1);
nii=roiTable2wholeBrainNii_shen([roi_labels(sig_mp_qfdr) mpEffect_m(sig_mp_qfdr)]);
save_nii(nii,[expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_simz_qfdr05_pos.nii']);

save([expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_stats2.mat'],'mpEffect','mp_cu1','mp_cu2','null_mpEffect','null_cu1','null_cu2','roi_labels','sig_mp_qmax','sig_mp_qfwe','sig_mp_qfdr','iters','srm');

beep
