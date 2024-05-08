function [H_normalized,gamma,obj,Z] = myregmultikernelclustering(K,cluster_count,lambda,tau,alpha)

nbkernel = size(K,3);
gamma0 = ones(nbkernel,1)/nbkernel;
flag = 1;
iter = 0;
sample_num = size(K,1);
K_mu  = mycombFun(K,gamma0);
numSel = round(tau*sample_num);
NS = genarateNeighborhood(K_mu,numSel);
%%--Calculate Neighborhood--%%
Ai_sum = zeros(sample_num);
for i =1:sample_num
	Ai_sum(NS(:,i),NS(:,i)) = Ai_sum(NS(:,i),NS(:,i))+1;
end
% Ai_sum = Ai_sum./sample_num;
% alpha = 1;
while flag
    iter = iter + 1;
    KC  = mycombFun(K,gamma0.^2);
    KC = Ai_sum.*KC;
    %% 求 H
    H = mykernelkmeans(KC,cluster_count);
    %% 求 Z
    tic
    resultMatrix = calK(K,tau);  % 计算K
    X=cell(1,1);
    X{1}=resultMatrix;
    Z = Z_demo(X,gamma0,alpha,lambda);
    toc
    obj(iter)  = calObj(H,K,Z{1},gamma0,lambda,Ai_sum);  
    %% 求γ
    tic
    [gamma]= updatekernelweightsV2(H,K,Z{1},lambda,Ai_sum);  % 求伽马
    gamma0 = gamma;
    toc
    %%
    if iter>2 && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4)
        flag =0;
    end
    if iter>100
        flag =0;
    end
end

H_normalized = H./ repmat(sqrt(sum(H.^2, 2)), 1,cluster_count);


% 产生样本近邻
function indx_0 = genarateNeighborhood(KC,tau,DIRECTION)
if nargin<3
    DIRECTION = 'descend';
end
num = size(KC,1);
KC0 = KC - 10^8*eye(num);
[val,indx] = sort(KC0,DIRECTION);
indx_0 = indx(1:tau,:);  % 选择矩阵 indx 中从第一行到第 tau 行的所有行，以及所有列。
end
end