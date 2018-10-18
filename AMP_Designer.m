% First of all
% load the AMP Forest model, as the evolutionary landscape
% Forest 100 contains 100 different decision trees;
% it's also optional to use AMP_Forest_50.mat to load a smaller model to
% accelearte the designing process.
load AMP_Forest_50.mat
% load the reference antimicrobial peptide sequences where the AMP Forests
% were built on
load AMP_reference.mat

%% define the parameters that will be input into the AMP Evolver to
% generatea the mutant library, which should be setup by users
generation = 8;% we run 100 in silico mutation-selection evolutionary cycles
mutation_rate = 100;% for each in silico mutagenesis method, we randomly repeat 50 times
selection_strengh = 10;% how many individuals will survive into the next round
parent = "AAAAAAAAAARAAAAAAAAA";% one example

[AMP_optimized,S0] = AMP_Evolver(parent,mdl,AMP_ref,generation,mutation_rate,selection_strengh);

%% All the functions
function [AMP_result,top_Score] = AMP_Evolver(parent,mdl,AMP_ref,generation,replicates,threshold)
% AMP_Evolver properties:
% parent - the original sequence where we start to mutate and evolve,cell format
% mdl - what kind of model we will use the score the peptide sequences
% AMP_ref - the reference for distance calcualtion
% generation - how many generation the algorithm will iterate
% threshold - how many top candidates will be retained for the next generation
% replicates - how many mutant will be randomly created in the mutagenesis function 
top_seq = parent;
top_Score = cell(generation,1);
top_Score(1) = {scorepredict(mdl,parent,AMP_ref)};
for i = 2:generation
    tic
    mutant = AMP_mutagenerator(top_seq,replicates);
    % estimate the arbitrary score of the bacterial survival under a
    % certain peptide sequences
    Score = scorepredict(mdl,mutant,AMP_ref);
    % sort all the mutants, and select ones that have the lowest survival
    Score = [cell2mat(top_Score(i-1));Score];
    mutant = [cellstr(top_seq);mutant];
    [Score,index] = sort(Score);
    top_Score(i) = {Score(1:threshold)};
    top_seq = string(mutant(index(1:threshold)));
    % store the scores
    toc
end
AMP_result = top_seq;
% calculate the predicted efficiency of parent sequence
end

%%
function [AMP_result,top_Score] = AMP_Evolvercluster(parent,mdl,AMP_ref,generation,replicates,threshold)
% AMP_Evolver properties:
% parent - the original sequence where we start to mutate and evolve,cell format
% mdl - what kind of model we will use the score the peptide sequences
% AMP_ref - the reference for distance calcualtion
% generation - how many generation the algorithm will iterate
% threshold - how many top candidates will be retained for the next generation
% replicates - how many mutant will be randomly created in the mutagenesis function 
top_seq = parent;
top_Score = cell(generation,1);
top_Score(1) = {scorepredict(mdl,parent,AMP_ref)};
for i = 2:generation
    tic
    mutant = AMP_clustermutator(top_seq,replicates);
    % estimate the arbitrary score of the bacterial survival under a
    % certain peptide sequences
    Score = scorepredict(mdl,mutant,AMP_ref);
    % sort all the mutants, and select ones that have the lowest survival
    Score = [cell2mat(top_Score(i-1));Score];
    mutant = [cellstr(top_seq);mutant];
    [Score,index] = sort(Score);
    top_Score(i) = {Score(1:threshold)};
    top_seq = string(mutant(index(1:threshold)));
    % store the scores
    toc
end
AMP_result = top_seq;
% calculate the predicted efficiency of parent sequence
end

%%
function mut = AMP_randmutagenesis(AMP_seq,replicates)
% define 20 commonly-used amino acids, and special groups:
AA = 'ARNDCQEGHILKMFPSTWYV';
number_AMP = numel(AMP_seq);
mut = {};
for k = 1:number_AMP
    seq = AMP_seq(k);
length_AMP = numel(seq);
mutant = cell(replicates,1);
for i = 1:replicates
   temp = char(seq);
   mutation_number = randi(length_AMP);
   for j = 1:mutation_number
   temp(randi(length_AMP)) = AA(randi(20));
   end
   mutant(i) = {temp};
end
mut = [mut;mutant];
end

end

%%
function mutant = AMP_clustermutator(AMP_seq,replicates)
Positive = 'RK';
Negative = 'ED';

AMP_seq = cellstr(AMP_seq);
Index_RK = cellfun(@(x) [strfind(x,'R'),strfind(x,'K')],AMP_seq,'UniformOutput',false);
Count_RK = cellfun(@sum,Index_RK);

AMP_charge_x = {};AMP_charge = {};
for i = 1:numel(AMP_seq)
    for j = 1:Count_RK % j define how many positive aa will be chosen as the center of the mutation
        %for k = 1:258 % replicates
        for k = 1:replicates
                distance = round(2*randn(j,1));% randomly choose the distance between the mutation and positive aa in wt
                Index_charge = randsample(Index_RK{i},j) + distance'; % for each chosen positive aa, mutate one aa close to it.
                Index_charge(find(Index_charge <1)) = 1;
                Index_charge(find(Index_charge > length(AMP_seq{i})))= length(AMP_seq{i}); % cut off
                AMP_temp = AMP_seq{i};
                AMP_temp(Index_charge) = 'x'; % replace aa by x
                AMP_charge_x = [AMP_charge_x;AMP_temp];
        end
    end
end
AMP_charge_x = unique(AMP_charge_x);
% apply different replace rules
for i = 1:2
    AMP_charge = [AMP_charge;strrep(AMP_charge_x,'x',Positive(i))];
end
for i = 1:2
    AMP_charge = [AMP_charge;strrep(AMP_charge_x,'x',Negative(i))];
end
AMP_charge = unique(AMP_charge);
AMP_hydrophobic = {};
AMP_hydrophobic = [AMP_hydrophobic;strrep(AMP_charge_x,'x','L')];
AMP_hydrophobic = unique(AMP_hydrophobic);
%% Find the unique list of all the mutants along with wild-type.
mutant = [cellstr(AMP_seq);AMP_charge;AMP_hydrophobic];

end

%%
% estimate the arbitrary score of the bacterial survival under a certain
% peptide
% mutant - in cell array; cell elements must be character vectors
function Score = scorepredict(mdl,mutant,ref)
    D = seqdist_ref(mutant,ref);
    Score = 2.1419*predict(mdl,D');
end
function D = seqdist_ref(x,ref)
D = zeros(numel(ref),numel(x));
for i = 1:numel(ref)
    parfor j = 1:numel(x)
        D(i,j) = seqpdist([ref(i),x(j)]);
    end
end
end