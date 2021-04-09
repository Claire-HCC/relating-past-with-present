%% this script compare the pattern similarity between C-AB moments with shared motif or only shared storyline (but different motifs) 

% set parameters;
loc='cluster';
set_parameters;
froidir='shen'
self_other='selfother';
srm='SRM';

% load the long table listing all C x AB motif tokens for each subject
load(sprintf('%s/fMRI/simMat/roi/%s/%s/mp/%s/mp_data.mat',expdir,self_other,froidir,srm),'mp_data');
mp_d=sortrows(mp_data,{'subj','mp_token_id_C','mp_token_id_AB'});
mp_d(:,{'mp_C','mp_AB'})=[];
mp_d(:,'mp_id_AB_shuffled')=table(0);

% list the names of all ROIs
rnames=regexp(mp_d.Properties.VariableNames,'.*_pauseAudResid_zscore','match');
roi_coli=find(~cellfun(@isempty,rnames));
rnames=[rnames{roi_coli}]';
roi_labels=regexp(rnames,'[0-9]*','match');
roi_labels=str2num(cell2mat([roi_labels{:}]'));

% remove unrelated pairs
mp_d(mp_d.cu==0,:)=[];
% average C motif tokens of the same motif type
mp_d_CtokenXsubj=grpstats([mp_d(:,{'subj','mp_token_id_C','cu'}) mp_d(:,roi_coli)],{'subj','mp_token_id_C','cu'});

% compute the averaged pattern similarities between C-AB motif pairs with shared motif and pairs sharing only the storyline (but with different motifs)
mp_cu2=grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu==2,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu==2))';
mp_cu1=grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu==1,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu==1))';

% null distribution generated by shuffling motif labels
iters=10000;
for iter=1:iters;
    
    % shuffle A and B motifs respectively
    for ab=[13 23];
        % list all the motif tokens in part C
        token_id_C=table2array(unique(mp_d(mp_d.mp_ab_C==ab & mp_d.subj==1,'mp_token_id_C')));
        
        % list all the motif types in storyline A or B
        ids_AB=table2array(mp_d(mp_d.mp_token_id_C==token_id_C(1,1) & mp_d.subj==1,'mp_id_AB'));
        
        % shuffle A or B sotryline motif types in the same way for all subjects, and all motif tokens in C
        shuff1 = randperm(length(ids_AB));
        mp_d.mp_id_AB_shuffled(mp_d.mp_ab_C==ab)=repmat(ids_AB(shuff1),25*length(token_id_C),1);
    end
    
    % categrize the shuffleded motif pairs
    % 2=shared motif; 1=shared storyline (different motif)
    mp_d.cu_shuffled(:,1)=1;
    mp_d.cu_shuffled(mp_d.mp_id_C==mp_d.mp_id_AB_shuffled,1)=2;
    
    % average pattern similarities by pair category 
    % 2=shared motif; 1=shared storyline (different motifs)
    mp_d_CtokenXsubj=grpstats([mp_d(:,{'subj','mp_token_id_C','cu_shuffled'}) mp_d(:,roi_coli)],{'subj','mp_token_id_C','cu_shuffled'});
    
    null_cu2(:,iter)=mean(grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu_shuffled==2,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu_shuffled==2)));
    null_cu1(:,iter)=mean(grpstats(table2array(mp_d_CtokenXsubj(mp_d_CtokenXsubj.cu_shuffled==1,5:end)),mp_d_CtokenXsubj.subj(mp_d_CtokenXsubj.cu_shuffled==1)));
end

% compute the real motif effect
mpEffect=mp_cu2-mp_cu1;
mpEffect_m=mean(mpEffect,2);

% compute the null motif effect
null_mpEffect=null_cu2-null_cu1;

% p-value and FWE correction 
p=sum(max(null_mpEffect)>mpEffect_m,2)/iters;
sig_mp_qmax=p<0.05;

save([expdir 'fMRI/simMat/roi/' self_other '/' froidir '/mp/' srm '/mp_stats.mat'],'mpEffect','mp_cu1','mp_cu2','null_mpEffect','null_cu1','null_cu2','roi_labels','sig_mp_qmax','iters','srm');

