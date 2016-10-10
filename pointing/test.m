clear all;
close all
load fisheriris
xdata = meas(1:70,3:4);
group = species(1:70);
svmStruct = svmtrain(xdata,group, 'kernel_function', 'linear','ShowPlot',true);
succeed =0;
fail =0;
for i =70:length(species)
[C,F] = svmclassifydist(svmStruct, meas(i,3:4) );
    if ismember(C,species(i)) 
   % fprintf('succeed \n')
        succeed = succeed + 1;
    else
  %  fprintf('fail \n')
        fail = fail + 1;
    end
end

%% SVM Multiclass Example
% SVM is inherently one vs one classification.
% This is an example of how to implement multiclassification using the
% one vs all approach.
TrainingSet=[ 1 10;2 20;3 30;4 40;5 50;6 66;3 30;4.1 42];
TestSet=    [3 34; 1 14; 2.2 25; 6.2 63];
GroupTrain=[1;1;2;2;3;3;2;2];
results = multisvm(TrainingSet, GroupTrain, TestSet);
disp('multi class problem');
disp(results);

% N = size(true_features,1);
% results = multisvm( true_features(1:N/2,:), true_label(1:N/2), true_features(N/2+1:end,:) );
