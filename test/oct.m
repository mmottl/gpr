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
global log_sf2_e = log_sf2 + epsilon;
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
  global log_sf2 log_sf2_e inv_ell2 inv_ell2_e epsilon;
  res = kf(x, y, log_sf2_e, inv_ell2);
end


%%%%%%%%%%%%%%%%%%%% Covariance matrices %%%%%%%%%%%%%%%%%%%%

Km = k(inducing_inputs, inducing_inputs);

Km_e = k_e(inducing_inputs, inducing_inputs);
dKm = (Km_e - Km) / epsilon;

Kmn = k(inducing_inputs, inputs);
Kmn_e = k_e(inducing_inputs, inputs);
dKmn = (Kmn_e - Kmn) / epsilon;

Kn = k(inputs, inputs);
Kn_e = k_e(inputs, inputs);
dKn = (Kn_e - Kn) / epsilon;


%%%%%%%%%%%%%%%%%%%% Main definitions %%%%%%%%%%%%%%%%%%%%

y = targets;

cholKm = chol(Km);
V = cholKm' \ Kmn;
Qn = V' * V;

lam = diag(diag(Kn - Qn));
lam_sigma2 = lam + sigma2 * eye(N);
inv_lam_sigma2 = inv(lam_sigma2);

Kmn_ = Kmn * sqrt(inv_lam_sigma2);

B = Km + Kmn_ * Kmn_';
cholB = chol(B);

r = diag(lam);
s = diag(lam_sigma2);
is = diag(inv_lam_sigma2);


%%%%%%%%%%%%%%%%%%%% Standard %%%%%%%%%%%%%%%%%%%%

%%%%%% Log evidence

l1 = -0.5*(log(det(B)) - log(det(Km)) + sum(log(s)) + N * log(2*pi))

R = cholB' \ Kmn;
S = cholB \ R * inv_lam_sigma2;
t = S*y;
u = is .* (y - Kmn' * t);
l2 = -0.5*(u'*y)

l = l1 + l2


%%%%%% Log evidence derivative

T = inv(Km) - inv(B);

U = cholKm \ V;

v1 = is .* (is .* (s - diag(R' * R)));
U1 = U*diag(sqrt(v1));
W1 = T - U1*U1';
X1 = S - U*diag(v1);

dl1 = -0.5*(v1' * diag(dKn) - trace(W1'*dKm)) - trace(X1'*dKmn)

v2 = u .* u;
U2 = U*diag(u);
W2 = t*t' - U2*U2';
X2 = t*u' - U*diag(v2);

dl2 = 0.5*(v2' * diag(dKn) - trace(W2'*dKm)) + trace(X2'*dKmn)

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

vv1 = is .* (is .* (s + s - r - diag(R' * R)));
vU1 = U*diag(sqrt(vv1));
vW1 = T - vU1*vU1';
vX1 = S - U*diag(vv1);

vdl1 = -0.5*(vv1' * diag(dKn) - trace(vW1'*dKm)) - trace(vX1'*dKmn)
vdl = vdl1 + dl2


%%%%%% Log evidence derivative wrt. noise

vdls1 = -0.5*(sum(vv1) - sum(is))
vdls = vdls1 + dls2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ed's stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyp = [log_inv_ell2; log_sf2; log(sigma2)];
ew = [reshape(inducing_inputs', M*dim, 1); hyp];
[eds_neg_log_likelihood, dfw] = spgp_lik(ew, y, inputs', M);
eds_evidence = -eds_neg_log_likelihood
eds_dlog_ell = -(-dfw(end - 2) * 2)
eds_dlog_sf2 = -dfw(end - 1)
eds_dsigma2 = -dfw(end) / sigma2
