% Trial select for BST
function [trials] = BST_TS_9(bst,varargin)
% Generates a unique list of trials corresponding to requested
% stimulus parameters. When 9 STDs have been played before a DEV
%
% These values should be as follows: 'Block',[1 3 5],'Lev1',[0 5
% 10],'Frq1',[1000 2000],... All
%
% When called without inputs, function returns a list of all
% trials.

% AH Updated 230601 to allow text searching too. With the try catch at
% line 26
% 

Parameter=varargin(1:2:end);
Value=varargin(2:2:end);
for j=1:length(Parameter)
    In.(Parameter{j})=Value{j};
end

EpocNamesL=fieldnames(In);
findList=[];

for i = 1:length(EpocNamesL)
    try
        findList=[findList; find(ismember(bst.Epocs.Values.(EpocNamesL{i}),In.(EpocNamesL{i})) ) ];
    catch
        findList=[findList; find(strcmp([bst.Epocs.Values.(EpocNamesL{i})], In.(EpocNamesL{i})) )]; % single line engine
    end
   
end
findCount=hist(findList,1:1:bst.NTrials);
trials_Complete=find(findCount==length(EpocNamesL));

n=2;
trials(1,1) = trials_Complete(1,1);

for ii = 2:length(trials_Complete)
if (trials_Complete(1,ii-1) - trials_Complete(1,ii)) < -8
   trials(1,n) = trials_Complete(ii);
   n = n +1;
end
end

if isempty(trials)
    trials=nan;
end

end