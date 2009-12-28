% Octave script for testing Gaussian process regression results
%
% Copyright (C) 2009-  Markus Mottl
% email: markus.mottl@gmail.com
% WWW:   http://www.ocaml.info

format long

global log_sf2;

load data/inputs
load data/targets
load data/inducing_points
load data/sigma2
load data/log_ell
load data/log_sf2

global epsilon = 1e-6;
sigma = sqrt(sigma2);
global log_sf2 = log_sf2;
global log_sf2_e = log_sf2 + epsilon;
global inv_ell2 = exp(-2 * log_ell);
global inv_ell2_e = exp(-2*(log_ell + epsilon));
log_inv_ell2 = log(inv_ell2);
[dim, N] = size(inputs);
[dim, M] = size(inducing_points);

function res = eval_rbf2(r2, a, b)
  res = exp(a + -0.5 * b * r2);
end

function res = kf(x, y, a, b)
  [dim, n1] = size(x);
  n2 = size(y, 2);
  r2 = repmat(sum(x' .* x', 2), 1, n2) - 2 * x' * y + repmat(sum(y' .* y', 2)', n1, 1);
  res = eval_rbf2(r2, a, b);
  [dim, N] = size(res);
  if (dim == N)
    jitter = 1e-6;
    res = res + jitter*eye(N);
  endif
  res = res';
end

function res = k(x, y)
  global log_sf2 inv_ell2;
  res = kf(x, y, log_sf2, inv_ell2);
end

function res = k_e(x, y)
  global log_sf2 log_sf2_e inv_ell2 inv_ell2_e epsilon;
  res = kf(x, y, log_sf2, inv_ell2_e);
end

function res = kf_diag(x, a, b)
  r2 = zeros(size(x, 2), 1);
  res = eval_rbf2(r2, a, b);
end

function res = k_diag(x)
  global log_sf2 inv_ell2;
  res = kf_diag(x, log_sf2, inv_ell2);
end

function res = k_diag_e(x)
  global log_sf2 log_sf2_e inv_ell2 inv_ell2_e epsilon;
  res = kf_diag(x, log_sf2, inv_ell2_e);
end


%%%%%%%%%%%%%%%%%%%% Covariance matrices %%%%%%%%%%%%%%%%%%%%

Km = k(inducing_points, inducing_points);
Km_e = k_e(inducing_points, inducing_points);
dKm = (Km_e - Km) / epsilon;

Knm = k(inducing_points, inputs);
Knm_e = k_e(inducing_points, inputs);
dKnm = (Knm_e - Knm) / epsilon;

Kn_diag = k_diag(inputs);
Kn_e_diag = k_diag_e(inputs);
dKn_diag = (Kn_e_diag - Kn_diag) / epsilon;


%%%%%%%%%%%%%%%%%%%% Main definitions %%%%%%%%%%%%%%%%%%%%

y = targets;

cholKm = chol(Km);
V = Knm / cholKm;

r = Kn_diag - sum(V .^ 2, 2);
s = r + sigma2;
is = ones(size(s, 1), 1) ./ s;
is_2 = sqrt(is);

inv_lam_sigma = repmat(is_2, 1, size(Knm, 2));

Knm_ = inv_lam_sigma .* Knm;

[Q, R] = qr([Knm_; chol(Km)], 1);
SF = diag(sign(diag(R)));
Q = Q(1:N,1:end)*SF;
R = SF*R;
S = inv_lam_sigma .* Q / R';

B = Km + Knm_' * Knm_;

%%%%%%%%%%%%%%%%%%%% Standard %%%%%%%%%%%%%%%%%%%%

%%%%%% Log evidence

l1 = ...
  -0.5*(...
    2*sum(log(diag(R))) - 2*sum(log(diag(cholKm))) + sum(log(s)) ...
    + N * log(2*pi))

y_ = is_2 .* y;
t = S'*y;
u = y_ - Q*(Q'*y_);
l2 = -0.5*(u'*y_)

l = l1 + l2


%%%%%% Log evidence derivative

T = inv(Km) - inv(B);

U = V / cholKm';

v1 = is .* (ones(size(Q, 1), 1) - sum(Q .^ 2, 2));
U1 = repmat(sqrt(v1), 1, size(U, 2)) .* U;
W1 = T - U1'*U1;
X1 = S - repmat(v1, 1, size(U, 2)) .* U;

dl1 = -0.5*(v1' * dKn_diag - trace(W1'*dKm)) - trace(X1'*dKnm)

w = is_2 .* u;
v2 = w .* w;
U2 = repmat(w, 1, size(U, 2)) .* U;
W2 = t*t' - U2'*U2;
X2 = w*t' - repmat(v2, 1, size(U, 2)) .* U;

dl2 = 0.5*(v2' * dKn_diag - trace(W2'*dKm)) + trace(X2'*dKnm)

dl = dl1 + dl2


%%%%%% Log evidence derivative wrt. noise

dls1 = -0.5*sum(v1)
dls2 = 0.5*sum(v2)
dls = dls1 + dls2


%%%%%%%%%%%%%%%%%%%% Variational %%%%%%%%%%%%%%%%%%%%

%%%%%% Log evidence

vl1 = l1 + -0.5*is'*r
vl = vl1 + l2


%%%%%% Log evidence derivative

vv1 = is .* (2*ones(size(Q, 1), 1) - is .* r - sum(Q .* 2, 2));
vU1 = repmat(sqrt(vv1), 1, size(U, 2)) .* U;
vW1 = T - vU1'*vU1;
vX1 = S - repmat(vv1, 1, size(U, 2)) .* U;

vdl1 = -0.5*(vv1' * dKn_diag - trace(vW1'*dKm)) - trace(vX1'*dKnm)
vdl = vdl1 + dl2


%%%%%% Log evidence derivative wrt. noise

vdls1 = -0.5*(sum(vv1) - sum(is))
vdls = vdls1 + dls2


%%%%%%%%%%%%%%%%%%%%%%%%% Ed Snelson's stuff %%%%%%%%%%%%%%%%%%%%%%%%%

hyp = [log_inv_ell2; log_sf2; log(sigma2)];
ew = [reshape(inducing_points', M*dim, 1); hyp];
[eds_neg_log_likelihood, dfw] = spgp_lik(ew, y, inputs', M);
eds_evidence = -eds_neg_log_likelihood
eds_dlog_ell = -(-dfw(end - 2) * 2)
eds_dlog_sf2 = -dfw(end - 1)
eds_dsigma2 = -dfw(end) / sigma2
