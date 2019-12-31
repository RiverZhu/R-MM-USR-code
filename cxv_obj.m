function [X_est,y]=cxv_obj(a,s,N)
n = 2;
C = nan(n+1,n+1,N);
for i = 1:N
    C(:,:,i) = [eye(2),-a(:,i);-a(:,i)',a(:,i)'*a(:,i)];
end



cvx_begin
   variable X(n+1,n+1) semidefinite ;
   variable G(N+1,N+1)  semidefinite ;
   minimize trace(G)-2*G(N+1,1:N)*s ;
   subject to
      for i = 1:N
          G(i,i)==trace(C(:,:,i)*X);
      end
     G(N+1,N+1)==1;
     X(n+1,n+1)==1;
cvx_end
[V_x,D_x] = eig(X);
[d_x,ind_d] = sort(diag(D_x),'descend');
X_est_extend = sqrt(d_x(1))*V_x(:,ind_d(1));
X_est = X_est_extend(1:end-1)/X_est_extend(end);
d = sqrt(sum((a-repmat(X_est,1,N)).^2));
d = d';
y=sum((s-d).^2);
