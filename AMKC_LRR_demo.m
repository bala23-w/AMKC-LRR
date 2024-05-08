clear
clc
warning off;

path = '.\';
addpath(genpath(path));
dataName = 'proteinFold';
%% caltech101_nTrain5_48
load([path,'datasets\',dataName,'_Kmatrix.mat'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
qnorm = 2;
% 创建存放结果的文件夹
result_dir = fullfile(pwd, '..', [dataName, '_AMKC_LRR']);
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

tic

lambdaset8 = 2.^[-10:1:10];
for tau = 0.1:0.05:1
    for il =1:length(lambdaset8)
        for arfa = 10.^[-5:1:5]
            tic
            mkkm_rm_result_file = fullfile(result_dir, [num2str(tau),'_', num2str(lambdaset8(il)), '_', num2str(arfa) '_result_AMKC_LRR.mat']);
            [H_normalized8,gamma81,obj81,Z] = myregmultikernelclustering(KH,numclass,lambdaset8(il),tau,arfa);
            eval = myNMIACC(H_normalized8,Y,numclass);
            save(mkkm_rm_result_file, 'eval');
            toc;
        end
    end
end