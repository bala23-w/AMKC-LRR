function D= my_SVD(Y,a)
% Da(Y)=argmin(1/2*norm(X-Y)+a*trace norm(X)),make X =Da(Y)
% Da(Y) is a operation,Y's single value decomposition UBV
% B's diag element bi,Da=U*diag(max(bi-a,0))*V
[I1,I2]=size(Y);
[U,B,V]=svd(Y);b=diag(B);m=size(b,1);
   a=a*ones(m,1);
   c=max(b-a,0);
   C=spdiags(c,0,I1,I2);
   D=U*C*V';
