format long

global log_sf2;

load data/inputs
load data/targets
load data/inducing_inputs
load data/sigma2
load data/log_ell
load data/log_sf2

global epsilon = 1e-6;
sigma = sqrt(sigma2);
global log_sf2 = log_sf2;
global inv_ell2 = exp(-2 * log_ell);
global inv_ell2_e = exp(-2*(log_ell + epsilon));
log_inv_ell2 = log(inv_ell2);
[dim, N] = size(inputs);
[dim, M] = size(inducing_inputs);

function res = eval_rbf2(r2, a, b)
  res = exp(a + -0.5 * b * r2);
end

function res = kf(x, y, a, b)
  [dim, n1] = size(x);
  n2 = size(y, 2);
  repmat(sum(y' .* y', 2), 1, n1);
  r2 = repmat(sum(x' .* x', 2), 1, n2) - 2 * x' * y + repmat(sum(y' .* y', 2)', n1, 1);
  res = eval_rbf2(r2, a, b);
  [dim, N] = size(res);
  if (dim == N)
    jitter = 1e-6;
    res = res + jitter*eye(N);
  endif
end

function res = k(x, y)
  global log_sf2 inv_ell2;
  res = kf(x, y, log_sf2, inv_ell2);
end

function res = k_e(x, y)
  global log_sf2 inv_ell2 inv_ell2_e epsilon;
  res = kf(x, y, log_sf2, inv_ell2_e);
end


%%%%%% Covariance matrices

Km = k(inducing_inputs, inducing_inputs);

Km_e = k_e(inducing_inputs, inducing_inputs);
dKm = (Km_e - Km) / epsilon;

Kmn = k(inducing_inputs, inputs);
Kmn_e = k_e(inducing_inputs, inputs);
dKmn = (Kmn_e - Kmn) / epsilon;

Kn = k(inputs, inputs);
Kn_e = k_e(inputs, inputs);
dKn = (Kn_e - Kn) / epsilon;


%%%%%% Main computations

y = targets;

cholKm = chol(Km);
V = cholKm' \ Kmn;  % with Qn
Qn = V' * V;

lam = diag(diag(Kn - Qn));
lam_sigma2 = lam + sigma2 * eye(N);
inv_lam_sigma2 = inv(lam_sigma2);

Kmn_ = Kmn * sqrt(inv_lam_sigma2);
Kmn__ = Kmn * inv_lam_sigma2;

B = Km + Kmn_ * Kmn_';
cholB = chol(B);

S = inv(Km) - inv(B);
T = cholB' \ Kmn__;
U = cholB \ T;
W = cholKm \ V;  % with Qn

r = diag(lam);
s = diag(lam_sigma2);
is = diag(inv_lam_sigma2);
dsh = 0.5 * diag(dKn) + diag(W'*(0.5 * dKm*W - dKmn));


%%%%%% Log evidence

t = U*y;
u = (Kmn' * t - y) .* is;

l1 = 0.5*(log(det(Km)) - log(det(B)) - sum(log(s)) - N * log(2*pi))
l2 = 0.5*(u'*y)
l = l1 + l2
vl1 = l1 - 0.5*is'*r
vl = vl1 + l2


%%%%%% Log evidence derivative

v = diag(T' * T) - is;
w = u .* u;

dl1c = 0.5*trace(S'*dKm) - trace(U'*dKmn);
vv = v - is + is.*(is.*r);

dl1 = v' * dsh + dl1c
dl2 = w' * dsh - t'*(dKm*0.5*t + dKmn*u)
dl = dl1 + dl2
vdl1 = vv' * dsh + dl1c
vdl = vdl1 + dl2

%%%%%% Log evidence derivative wrt. noise

dls1 = 0.5*(sum(v))
dls2 = 0.5*(sum(w))
dls = dls1 + dls2
vdls1 = 0.5*(sum(vv) + sum(is))
vdls = vdls1 + dls2


%%%%%% Ed's stuff

hyp = [log_inv_ell2; log_sf2; log(sigma2)];
ew = [reshape(inducing_inputs', M*dim, 1); hyp];
[eds_neg_log_likelihood, dfw] = spgp_lik(ew, y, inputs', M);
eds_evidence = -eds_neg_log_likelihood
eds_dlog_ell = -(-dfw(end - 2) * 2)
eds_dlog_sf2 = -dfw(end - 1)
eds_dsigma2 = -dfw(end) / sigma2
