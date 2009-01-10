load data/inputs
load data/targets
load data/inducing_inputs
load data/sigma2

sigma = sqrt(sigma2);
[dim, N] = size(inputs);
[dim, M] = size(inducing_inputs);

function res = eval_rbf2(r2, a, b)
  res = exp(a + b * r2);
end

function res = kf(x, y, a, b)
  [dim, n1] = size(x);
  n2 = size(y, 2);
  repmat(sum(y' .* y', 2), 1, n1);
  r2 = repmat(sum(x' .* x', 2), 1, n2) - 2 * x' * y + repmat(sum(y' .* y', 2)', n1, 1);
  res = eval_rbf2(r2, a, b);
end

function res = k(x, y)
  a = -0.5;
  b = -1.0;
  res = kf(x, y, a, b);
end

function res = k_e(x, y)
  a = -0.5;
  b = -1.0;
  epsilon = 10e-6;
  res = kf(x, y, a + epsilon, b);
end

km = k(inducing_inputs, inducing_inputs);
km_e = k_e(inducing_inputs, inducing_inputs);
dkm = (km_e - km) / 10e-6;

kmn = k(inducing_inputs, inputs);
kmn_e = k_e(inducing_inputs, inputs);
dkmn = (kmn_e - kmn) / 10e-6;

kn = k(inputs, inputs);
kn_e = k_e(inputs, inputs);
dkn = (kn_e - kn) / 10e-6;

jitter = 10e-6;
km = km + jitter*eye(M);
kn = kn + jitter*eye(N);

km_chol = chol(km);

qn = kmn' * inv(km) * kmn;
lam = diag(diag(kn - qn));
lam_sigma2 = lam + sigma2 * eye(N);
inv_lam_sigma2 = inv(lam_sigma2);
inv_lam_sigma = sqrt(inv_lam_sigma2);
kmn_ = kmn * inv_lam_sigma;
y = targets;
y_ = inv_lam_sigma * y;
log_det_lam_sigma2 = log(det(lam_sigma2));
b = km + kmn * inv_lam_sigma2 * kmn';
b_chol = chol(b);
kmn_y_ = kmn_ * y_;

%

dlam__ = diag(inv_lam_sigma2 * diag(dkn - 2*dkmn'*inv(km)*kmn + kmn'*inv(km)*dkm*inv(km)*kmn));
dkmn_trace = trace(inv_lam_sigma2 * kmn' * inv(b)*(2*dkmn - kmn * dlam__));
dlam__trace = trace(dlam__);
dkm_trace = trace(((inv(b) - inv(km)) * dkm));
trace_sum = dkmn_trace + dlam__trace + dkm_trace

y__=inv_lam_sigma2*y;
x=inv(b)*kmn*y__;
z=kmn'*x;

dlam_factor = ((2*z - y).*y__ - inv_lam_sigma2*z.*z)
dkmn_factor = 2*(inv_lam_sigma2*z - y__)
dkm_nll = x'*dkm*x
dkmn_nll = x'*dkmn*dkmn_factor
dlam_nll = dlam_factor'*diag(dlam__)
dl2 = 0.5*(dkm_nll + dkmn_nll + dlam_nll);


%

log_det_b = log(det(b));
log_det_km = log(det(km));
log_det_lam_sigma2 = log(det(lam_sigma2));

model_nll = (log_det_b - log_det_km + log_det_lam_sigma2 + N * log(2*pi)) / 2
trained_nll = (y' * inv(qn + lam_sigma2) * y) / 2;

nll = (model_nll + trained_nll);
evidence = - nll;

model_nll_dsigma2 = trace(inv(qn + lam_sigma2)) / 2;
model_evidence_dsigma2 = - model_nll_dsigma2;

% Trained
evidence

% Ed's stuff
hyp = [-log(0.5); -1 / 2; log(sigma2)];
w = [reshape(inducing_inputs', M*dim, 1); hyp];
[eds_neg_log_likelihood, dfw] = spgp_lik(w, y, inputs', M);
eds_evidence = -eds_neg_log_likelihood
eds_dsigma2 = -dfw(end) / sigma2
