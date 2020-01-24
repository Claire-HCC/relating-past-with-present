tic
% motif used to be called as memory plant (mp)....

clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
type='Pearson';
fragmentScanN=5;
self_other='selfother';
srm='SRM';
froidir='shen';

if strmatch(srm,'SRM');
    rnames=dir([expdir '/fMRI/timeseries/' timeUnit '/roi/' froidir '/*_pauseAudResid_zscore_srm.mat']);
else
    rnames=dir([expdir '/fMRI/timeseries/' timeUnit '/roi/' froidir '/roi*_pauseAudResid_zscore.mat']);
end
rnames={rnames.name};
rnames=strrep(rnames,'.mat','');

load([dir_texts '/segment_time.mat']);
mp_table=readtable([dir_texts '/hypothesis/mp_token_table_used.csv']);
mpN=size(mp_table,1);

% create a table with all the pairs of c mp & ab mp from the same storyline for each subject
mp_c_inds=find(mp_table.abc==3 & ismember(mp_table.mp_ab,[13 23]));
n=1;
for si=1:25;
    
    for mp_c_i=1:length(mp_c_inds);
        mp_c_ind=mp_c_inds(mp_c_i);
        mp_ab_inds=find(ismember(mp_table.abc,[1 2]) );
        
        % create table to save the data
        nrow=length(mp_ab_inds);
        mp_data(n:(n+nrow-1),{'mp_ab_C','mp_C','mp_id_C','mp_token_id_C','abc_C','emotion_C','Character_C','MemRepetition_C'})=repmat(mp_table(mp_c_ind,{'mp_ab','mp','mp_id','mp_token_id','abc','EMOTIONALWEIGHTOFELEMENT','CharacterMain','MemRepetition'}),nrow,1);
        mp_data(n:(n+nrow-1),{'mp_ab_AB','mp_AB','mp_id_AB','mp_token_id_AB','abc_AB','emotion_AB','Character_AB','MemRepetition_AB'})=mp_table( mp_ab_inds,{'mp_ab','mp','mp_id','mp_token_id','abc','EMOTIONALWEIGHTOFELEMENT','CharacterMain','MemRepetition'});
        mp_data(n:(n+nrow-1),'subj')=table(si);
        
        n=n+nrow;
    end
end

mp_data.cu(:,1)=0;
mp_data.cu(mp_data.mp_ab_C==mp_data.mp_ab_AB)=1;
mp_data.cu(mp_data.mp_id_C==mp_data.mp_id_AB)=2;

mp_data.Character_cu(:,1)=double(0);
mp_data.Character_cu(ismember(mp_data.Character_AB,{'Clara','Gary'})==ismember(mp_data.Character_C,{'Clara','Gary'}))=1;
mp_data.Character_cu(ismember(mp_data.Character_AB,{'Margaret','Sasha','Steven','Alexander'})==ismember(mp_data.Character_C,{'Margaret','Sasha','Steven','Alexander'}))=1;
mp_data.Character_cu(strcmp(mp_data.Character_C,mp_data.Character_AB))=2;

% fill in the fisher's z transformed r in the table
for ri=1:length(rnames);
    rname=rnames{ri};
    
    fr =  sprintf('%s/fMRI/timeseries/%s/roi/%s/%s.mat',expdir,timeUnit,froidir,rname);
    load(fr,'gdata');
    gdata_z=zeros(size(gdata));
    gdata_z(:,segmentv_inTr>0,:)=zscore(gdata(:,segmentv_inTr>0,:),0,2);
    
    for i=1:size(mp_data,1);
        
        mp_tr_C=round(table2array(mp_table(mp_table.mp_token_id==mp_data.mp_token_id_C(i),'mp_start_tr')));
        mp_tr_C=(mp_tr_C+1):(mp_tr_C+fragmentScanN);
        mp_tr_C(ismember(mp_tr_C,find(segmentv_inTr==0)))=[];
        mp_data.mp_tr_C(i)=mp_tr_C(1);
        
        mp_tr_AB=round(table2array(mp_table(mp_table.mp_token_id==mp_data.mp_token_id_AB(i),'mp_start_tr')));
        mp_tr_AB=(mp_tr_AB+1):(mp_tr_AB+fragmentScanN);
        mp_tr_AB(ismember(mp_tr_AB,find(segmentv_inTr==0)))=[];
        mp_data.mp_tr_AB(i)=mp_tr_AB(1);
        
        si=mp_data.subj(i);
        othersi=1:25;
        othersi=othersi(~ismember(othersi,si));
        
        self_C=mean(gdata_z(:,mp_tr_C,si),2);
        self_AB=mean(gdata_z(:,mp_tr_AB,si),2);
        others_C=mean(mean(gdata_z(:,mp_tr_C,othersi),2),3);
        others_AB=mean(mean(gdata_z(:,mp_tr_AB,othersi),2),3);
        
        if strcmp(self_other,'selfother');
            sim1=corr(self_C,others_AB);
            sim2=corr(self_AB,others_C);
            sim_z=(0.5*log((1+sim1)./(1-sim1))+0.5*log((1+sim2)./(1-sim2)))/2;
        elseif strcmp(self_other,'selfself');
            sim=corr(self_AB,self_C);
            sim_z=0.5*log((1+sim)./(1-sim));
        end
        
        mp_data.(rname)(i)=sim_z';
        
    end
end
% save(sprintf('%s/fMRI/simMat/roi/%s/%s/mp/%s/mp_data.mat',expdir,self_other,froidir,srm),'mp_data','srm','-v7.3');
toc
beep;