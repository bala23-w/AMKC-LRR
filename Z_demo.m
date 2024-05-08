function [Z] =Z_demo(X,F,alpha,lll)
% X is union of data about cell
% for example: X=cell(3,1),X{1}=200*300;X{2}=800*300;X{3}=600*300;
%   Detailed explanation goes here
c=length(F);
A=size(X);[~,N]=size(X{1});%N is the number of samples;
K=A(1);                     %K is number of data;
Z=cell(K,1);lambda=cell(K,1);
J=Z;S=lambda;L=Z;

V=cell(K,1);t=0;
 
 U=cell(K,1);I=eye(N,N);h=ones(N,1);
 for i=1:K  
       J{i}=zeros(N,N);    
       Z{i}=zeros(N,N);
       S{i}=zeros(N,N);
       lambda{i}=zeros(N,N);
       U{i}=zeros(N,N);
       V{i}=zeros(N,c);                 
 end
 
W=ones(K,1)*(1/K);
H=ones(K,1);R1=ones(K,1);R2=ones(K,1);
mu=0.05;
mumax=10^10;le=1.2;delta=10^(-6);
while ( max(R1)>=delta || max(R2)>=delta)
    M1=zeros(N,N);
    N1=zeros(N,c);
    for i=1:K
        
       J{i}= my_SVD(Z{i}+(1/mu)*lambda{i},alpha/mu);

        S{i}=(mu*I+2*W(i)*X{i}'*X{i})\(mu*Z{i}+U{i}+2*W(i)*X{i}'*X{i});
         
        V{i}=zeros(N,c);

        partial2=(F+V{i})*(F+V{i})'*(lll/2);
        A=(1/2)*(J{i}+S{i}-(lambda{i}+U{i}+partial2)/mu);
        A_bar=A-diag(diag(A));
        Z{i}=max((A_bar+A_bar')/2,0);
        H(i)=norm(X{i}-X{i}*S{i},'fro')^2;
        L{i}=Z{i};
        M1=L{i}+M1;
        N1=L{i}*V{i}+N1;
        
    end
% [F, eigenvalue] = eigs(M1,c);
% [~,b] = sort(diag(eigenvalue),'ascend');%good
% eigenvalue = eigenvalue(:,b);
% F = F(:,b);
% F= F(:,1:c);
    for i=1:K
        lambda{i}=lambda{i}+mu*(Z{i}-J{i});
        U{i}=U{i}+mu*(Z{i}-S{i});
        R1(i)=norm(Z{i}-J{i},inf);
        R2(i)=norm(Z{i}-S{i},inf);
    end
    mu=min(le*mu,mumax);
    t=t+1;
    if t>300
        break
    end
end  


end


        
        
        
        
        
        
        
        
        
        
        
    