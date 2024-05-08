function resultMatrix = calK(KH,tau)

nbkernel = size(KH,3);
gamma0 = ones(nbkernel,1)/nbkernel;
sample_num = size(KH,1);


K_mu  = mycombFun(KH,gamma0);
numSel = round(tau*sample_num);
NS = genarateNeighborhood(K_mu,numSel);
%%--Calculate Neighborhood--%%
Ai_sum = zeros(sample_num);
for i =1:sample_num
	Ai_sum(NS(:,i),NS(:,i)) = Ai_sum(NS(:,i),NS(:,i))+1;
end
Ai_sum = Ai_sum./sample_num;



numker = size(KH,3);

rows = numel(Ai_sum); 


% 根据计算结果初始化resultMatrix
resultMatrix = zeros(rows, numker);

% 遍历第三维度
for k = 1:numker
    slice = KH(:, :, k); % 获取当前切片
    modifiedSlice = Ai_sum .* slice; % 将切片逐元素乘以矩阵,选取邻域元素
    columnVector = modifiedSlice(:); % 将修改后的切片按列展开成一列
    resultMatrix(:, k) = columnVector; % 将列向量放入结果矩阵的对应列
end




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