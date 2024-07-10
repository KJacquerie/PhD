clear all 
close all 
clc

%%
% load('MNIST_[0, 1, 2]_n=[100, 100, 100]_crop=4_downsample=2_binarize=True_onehot=True.mat')

% NB_neurons = size(data,2);
% NB_data = size(data,1);
% dataset = zeros(NB_neurons, NB_data/2);
% testset = zeros(NB_neurons, NB_data/2);
% 
% dataset(:,1:50) =data(1:50,:)'; 
% dataset(:,51:100) =data(101:150,:)'; 
% dataset(:,101:150) =data(201:250,:)'; 
% 
% testset(:,1:50) =data(51:100,:)'; 
% testset(:,51:100) =data(151:200,:)'; 
% testset(:,101:150) =data(251:300,:)'; 

load('MNIST_[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]_n=[100, 100, 100, 100, 100, 100, 100, 100, 100, 100]_crop=3_downsample=1_binarize=True_onehot=True.mat');

%load('MNIST_[0, 1, 2, 3, 4, 5]_n=[100, 100, 100, 100, 100, 100]_crop=3_downsample=1_binarize=True_onehot=True.mat');
%load('MNIST_[0, 1, 2, 3, 4]_n=[100, 100, 100, 100, 100]_crop=3_downsample=1_binarize=True_onehot=True.mat');


%%
NB_neurons = size(data,2);
NB_data = size(data,1);
dataset = zeros(NB_neurons, NB_data/2);
testset = zeros(NB_neurons, NB_data/2);

%%
dataset(:,1:50) =data(1:50,:)'; 
dataset(:,51:100) =data(101:150,:)'; 
dataset(:,101:150) =data(201:250,:)'; 
dataset(:,151:200) =data(301:350,:)'; 
dataset(:,201:250) =data(401:450,:)'; 
dataset(:,251:300) =data(501:550,:)'; 
dataset(:,301:350) =data(601:650,:)'; 
dataset(:,351:400) =data(701:750,:)'; 
dataset(:,401:450) =data(801:850,:)'; 
dataset(:,451:500) =data(901:950,:)'; 


testset(:,1:50) =data(51:100,:)'; 
testset(:,51:100) =data(151:200,:)'; 
testset(:,101:150) =data(251:300,:)'; 
testset(:,151:200) =data(351:400,:)'; 
testset(:,201:250) =data(451:500,:)'; 
testset(:,251:300) =data(551:600,:)'; 
testset(:,301:350) =data(651:700,:)'; 
testset(:,351:400) =data(751:800,:)'; 
testset(:,401:450) =data(851:900,:)'; 
testset(:,451:500) =data(951:1000,:)'; 


%%

mm = randi([1 500],1,13);


for idx = 1:1:length(mm)
   figure(idx)
    %subplot(length(mm),1,idx)
     X = reshape(dataset(:,mm(idx)), [sqrt(NB_neurons),sqrt(NB_neurons)]); 
     
     imshow(X')
end

% build pixel ON 
for idx=1:1:size(dataset,2)
    X = reshape(dataset(:,idx), [sqrt(NB_neurons),sqrt(NB_neurons)]); 
    pixelON(idx) = length(find(X>0.5)); 
end
for idx=1:1:size(dataset,2)
    X = reshape(testset(:,idx), [sqrt(NB_neurons),sqrt(NB_neurons)]); 
    pixelON_test(idx) = length(find(X>0.5)); 
end
%%%%%
% before save dataset clean the file in your folder otherwise it stacked at
% the end

%% SAVE NEW DATA SET -- erase the previous one & uncomment

%dlmwrite('dataset_MNIST.dat', dataset,'-append','delimiter',' ','roffset',1)
%dlmwrite('testset_MNIST.dat', dataset,'-append','delimiter',' ','roffset',1)

%%
%idx_wl = randi(NB_neurons*5, 1, 100);
%dlmwrite('idx_wl.dat', round(idx_wl),'-append','delimiter',' ','roffset',1)

 
