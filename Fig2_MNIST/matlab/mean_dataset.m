clc
clear all 
close all 

%%


load('dataset_MNIST.dat');
load('testset_MNIST.dat'); 

%%
NB_pixels = size(dataset_MNIST,1); 
N = 50; 
NB_dgt = size(dataset_MNIST,2)/N; 
NB_neurons = sqrt(NB_pixels);
mean_data = zeros(NB_pixels, NB_dgt); 
mean_testset = zeros(NB_pixels, NB_dgt); 

%%
F = 1/16*[-1 2 -1; 2 4 2; -1 2 -1]; 

for idx=1:1:NB_dgt
    mean_data(:,idx) = mean(dataset_MNIST(:,(idx-1)*N+1:N*idx),2); 
    X = reshape(mean_data(:,idx), [NB_neurons,NB_neurons]);
    Y = conv2(X,F);
    figure(idx)
    subplot(1,2,1)
    imagesc(X')
    subplot(1,2,2)
    imagesc(Y')
 
end

%%

dlmwrite('mean_MNIST.dat', mean_data,'-append','delimiter',' ','roffset',1)
%dlmwrite('filt_MNIST.dat', dataset,'-append','delimiter',' ','roffset',1)

%%
% 
% close all
% 
% for idx=1:1:50
% figure(idx)
% imagesc(reshape(dataset_MNIST(:,100+idx), [22,22]))
% end
% 
% figure(100)
% subplot(3,1,1)
% imagesc(reshape(neurons_freq(5,2:end-10), [22,22]))
% subplot(3,1,2)
% imagesc(reshape(neurons_freq(25,2:end-10), [22,22]))
% subplot(3,1,3)
% imagesc(reshape(neurons_freq(45,2:end-10), [22,22]))


%%

mean_MNIST_presented = zeros(484, 10); 

for idx=1:1:NB_dgt
    mean_MNIST_presented(:,idx) = mean(dataset_MNIST(:,idx_deM(idx,:)),2);
end

dlmwrite('mean_MNIST_presented.dat', mean_MNIST_presented,'-append','delimiter',' ','roffset',1)


%%

for idx=1:1:NB_dgt
    
    X = reshape(mean_data(:,idx), [NB_neurons,NB_neurons]);
    Y = reshape(mean_MNIST_presented(:,idx), [NB_neurons,NB_neurons]);
    figure(idx)
    subplot(1,2,1)
    imagesc(X')
    subplot(1,2,2)
    imagesc(Y')
 
end

