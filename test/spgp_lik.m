% Edward Snelson's SPGP implementation in Octave/Matlab

function [fw,dfw] = spgp_lik(w,y,x,n,del)

% spgp_lik_3: neg. log likelihood for SPGP and gradients with respect to
% pseudo-inputs and hyperparameters. Gaussian covariance with one
% lengthscale per dimension.
%
% y -- training targets (N x 1)
% x -- training inputs (N x dim)
% n -- number of pseudo-inputs
% w -- parameters = [reshape(xb,n*dim,1); hyp]
%
%      hyp -- hyperparameters (dim+2 x 1)
%      hyp(1:dim) = log( b )
%      hyp(dim+1) = log( c )
%      hyp(dim+2) = log( sig )
%
%      where cov = c * exp[-0.5 * sum_d b_d * (x_d - x'_d)^2]
%                       + sig * delta(x,x')
%
%      xb -- pseudo-inputs (n*dim)
%
% del -- OPTIONAL jitter (default 1e-6)
%
% fw -- likelihood
% dfw -- OPTIONAL gradients
%
% Edward Snelson (2006)

if nargin < 5; del = 1e-6; end % default jitter

[N,dim] = size(x); xb = reshape(w(1:end-dim-2),n,dim);
b = exp(w(end-dim-1:end-2)); c = exp(w(end-1)); sig = exp(w(end));

xb = xb.*repmat(sqrt(b)',n,1);
x = x.*repmat(sqrt(b)',N,1);

Q = xb*xb';
Q = repmat(diag(Q),1,n) + repmat(diag(Q)',n,1) - 2*Q;
Q = c*exp(-0.5*Q) + del*eye(n);

K = -2*xb*x' + repmat(sum(x.*x,2)',n,1) + repmat(sum(xb.*xb,2),1,N);
K = c*exp(-0.5*K);

L = chol(Q)';
V = L\K;
ep = 1 + (c-sum(V.^2)')/sig;
K = K./repmat(sqrt(ep)',n,1);
V = V./repmat(sqrt(ep)',n,1); y = y./sqrt(ep);
Lm = chol(sig*eye(n) + V*V')';
invLmV = Lm\V;
bet = invLmV*y;

% Likelihood
fw = sum(log(diag(Lm))) + (N-n)/2*log(sig) + ...
      (y'*y - bet'*bet)/2/sig + sum(log(ep))/2 + 0.5*N*log(2*pi);

% OPTIONAL derivatives
if nargout > 1

% precomputations
Lt = L*Lm;
B1 = Lt'\(invLmV);
b1 = Lt'\bet;
invLV = L'\V;
invL = inv(L); invQ = invL'*invL; clear invL
invLt = inv(Lt); invA = invLt'*invLt; clear invLt
mu = ((Lm'\bet)'*V)';
sumVsq = sum(V.^2)'; clear V
bigsum = y.*(bet'*invLmV)'/sig - sum(invLmV.*invLmV)'/2 - (y.^2+mu.^2)/2/sig ...
         + 0.5;
TT = invLV*(invLV'.*repmat(bigsum,1,n));

% pseudo inputs and lengthscales
for i = 1:dim
% dnnQ = (repmat(xb(:,i),1,n)-repmat(xb(:,i)',n,1)).*Q;
% dNnK = (repmat(x(:,i)',n,1)-repmat(xb(:,i),1,N)).*K;
dnnQ = dist(xb(:,i),xb(:,i)).*Q;
dNnK = dist(-xb(:,i),-x(:,i)).*K;

epdot = -2/sig*dNnK.*invLV; epPmod = -sum(epdot)';

dfxb(:,i) = - b1.*(dNnK*(y-mu)/sig + dnnQ*b1) ...
    + sum((invQ - invA*sig).*dnnQ,2) ...
    + epdot*bigsum - 2/sig*sum(dnnQ.*TT,2);

dfb(i,1) = (((y-mu)'.*(b1'*dNnK))/sig ...
           + (epPmod.*bigsum)')*x(:,i);

dNnK = dNnK.*B1; % overwrite dNnK
dfxb(:,i) = dfxb(:,i) + sum(dNnK,2);
dfb(i,1) = dfb(i,1) - sum(dNnK,1)*x(:,i);

dfxb(:,i) = dfxb(:,i)*sqrt(b(i));

dfb(i,1) = dfb(i,1)/sqrt(b(i));
dfb(i,1) = dfb(i,1) + dfxb(:,i)'*xb(:,i)/b(i);
dfb(i,1) = dfb(i,1)*sqrt(b(i))/2;
end

% size
epc = (c./ep - sumVsq - del*sum((invLV).^2)')/sig;

dfc = (n + del*trace(invQ-sig*invA) ...
     - sig*sum(sum(invA.*Q')))/2 ...
    - mu'*(y-mu)/sig + b1'*(Q-del*eye(n))*b1/2 ...
      + epc'*bigsum;

% noise
dfsig = sum(bigsum./ep);

dfw = [reshape(dfxb,n*dim,1);dfb;dfc;dfsig];

end
