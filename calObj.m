function obj = calObj(T,K,H0,gamma0,lambda,Ai_sum)  % T——H; K——KH; H0——M

nbkernel = size(K,3);
num = size(K,1);
f = zeros(nbkernel,1);
for p = 1 : nbkernel
    f(p) = trace(T'*(Ai_sum.*K(:,:,p))*T);
    % f(p) = (trace(K(:,:,p))-trace(T'*K(:,:,p)*T));
end
    
obj = (1/2)*gamma0'*(lambda*H0+2*diag(f))*gamma0;
