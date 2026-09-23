function R = run_phase_a(btc_id, varargin)
% RUN_PHASE_A  Fit one BTC with Volponi's SMIMfit through an adapter that reproduces
% her Lead_SMIM.m workflow, and export everything the R validation needs.
%
%   R = run_phase_a(btc_id, 'name', value, ...)
%
% Options (defaults):
%   'data_dir'     private/data                      where btc_corrected.csv / fit_windows.csv live
%   'out_root'     private/validation/matlab_out     one sub-folder per run is created
%   'objective'    'volponi'                         'volponi' (her objFunctionSMIM) | 'log' (objFunctionSMIM_log)
%   't_end_mode'   'ref'                             'ref' | 'lod' | 'loq' | 'record' | numeric seconds
%   'seed'         1000                              rng(seed,'twister') before the fit
%   'instrumented' false                             also run lsqnonlin from every start point separately
%   'display'      'iter'                            lsqnonlin Display (runFit hard-codes 'iter'; only
%                                                    affects the instrumented pass)
%
% Workflow (mirrors Lead_SMIM.m; line numbers refer to her file):
%   1. keep samples with t >= 0                                    (lines 68-73, IndexRelease)
%   2. Conc(Conc < 0) = 0                                          (line 57)
%   3. Conc(t < HitTime) = 0                                       (lines 76-77)
%   4. drop samples with t > t_end                                 (lines 85-88, CutTime)
%   5. [ccNorm, mrf, Qdg] = cNorm(t, Conc, 'c', 'cMass', M, 'Q', Q) (lines 93-97)
%   6. params_guess = (upper + lower)/2                            (line 100)
%   7. SMIMfit(t, ccNorm, 'c', 'L', L, 't_end', t(end), ...)       (line 103)
%   8. covariance / SE                                             (lines 107-140)
% Decisions applied here (PLAN.md §10): 4 free parameters, log10(t1) and log10(t2) fixed
% at log10(1.5) and log10(14000) by setting lb == ub; SE uses only the 4 free columns.
%
% Nothing in Volponi's folder is modified. Her folder is put on the path read-only.

p = inputParser;
addParameter(p, 'data_dir', '');
addParameter(p, 'out_root', '');
addParameter(p, 'objective', 'volponi', @(s) any(strcmp(s, {'volponi', 'log'})));
addParameter(p, 't_end_mode', 'ref');
addParameter(p, 'seed', 1000, @isnumeric);
addParameter(p, 'instrumented', false, @islogical);
addParameter(p, 'display', 'iter');
addParameter(p, 'smoke', false, @islogical);   % NOT A FIT: 2-iteration lsqnonlin to exercise the export path
parse(p, varargin{:});
o = p.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
her_dir   = fullfile(repo_root, 'inputs', 'SabrinaVolponi-SMIMfit-058f87d');
addpath(her_dir);                       % her core, read-only
addpath(fileparts(mfilename('fullpath')));
if isempty(o.data_dir), o.data_dir = fullfile(repo_root, 'private', 'data'); end
if isempty(o.out_root), o.out_root = fullfile(repo_root, 'private', 'validation', 'matlab_out'); end

%% ---- load ---------------------------------------------------------------------------
B = load_btc(btc_id, o.data_dir);
switch class(o.t_end_mode)
    case 'char'
        switch o.t_end_mode
            case 'ref',    t_end = B.t_end_ref_s;
            case 'lod',    t_end = B.t_end_lod_s;
            case 'loq',    t_end = B.t_end_loq_s;
            case 'record', t_end = B.t_end_record_s;
            otherwise, error('unknown t_end_mode %s', o.t_end_mode);
        end
        if isnan(t_end), t_end = B.t_end_record_s; end   % rule did not bind -> whole record
        tag = o.t_end_mode;
    otherwise
        t_end = o.t_end_mode; tag = sprintf('t%g', t_end);
end
run_name = sprintf('%s__%s__%s__seed%d', btc_id, o.objective, tag, o.seed);
if o.smoke, run_name = [run_name '__SMOKE_NOT_A_FIT']; end
out_dir  = fullfile(o.out_root, run_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
diary(fullfile(out_dir, 'matlab_log.txt')); diary on
fprintf('=== run_phase_a %s  (%s)\n', run_name, datestr(now));

%% ---- Lead_SMIM.m cleaning, reproduced step by step ---------------------------------
t    = B.t_s;
Conc = B.C_ppb;
% step 1: IndexRelease (Lead line 68-73): drop everything before release
keep = t >= 0;  t = t(keep);  Conc = Conc(keep);
n_after_release = numel(t);
% step 2: negatives -> 0 (Lead line 57)
n_neg_clamped = sum(Conc < 0);
Conc(Conc < 0) = 0;
% step 3: pre-hit -> 0 (Lead line 76-77)
n_prehit_zeroed = sum(t < B.hit_time_s & Conc > 0);
Conc(t < B.hit_time_s) = 0;
% step 4: cut (Lead line 85-88)
cut = t > t_end;  t(cut) = [];  Conc(cut) = [];
% step 5: normalise (Lead line 93-97); cMass in ug, Q in L/s, C in ug/L
[ccNorm, massRecoveryFraction, Qdg] = cNorm(t, Conc, 'c', 'cMass', B.mass_ug, 'Q', B.Q_L_s);

%% ---- bounds / guess (Lead lines 14, 17, 100) with t1, t2 fixed (decision Q1) --------
params_upper = [1.5 2.3 0.5 2 0.2 12];
params_lower = [0.05 2E-10 0.0001 0  0 0.5];
logT1 = log10(1.5);  logT2 = log10(14000);
params_upper(5) = logT1;  params_lower(5) = logT1;
params_upper(6) = logT2;  params_lower(6) = logT2;
params_guess = (params_upper + params_lower) / 2;
free = [1 2 3 4];  pfree = numel(free);

%% ---- start points: regenerate exactly what runFit.m:51-52 will draw -----------------
rng(o.seed, 'twister');
points = lhsdesignbnd(60, 6, params_lower, params_upper);
points(61, :) = params_guess;
write17(fullfile(out_dir, 'starts.csv'), {'v','D','Lambda','beta','logT1','logT2'}, points);

%% ---- fit --------------------------------------------------------------------------------
tfit = tic;
if o.smoke
    % exercise every downstream step with a deliberately truncated single lsqnonlin run
    tt = t; cc = ccNorm;
    if tt(1) == 0, tt = tt(2:end); cc = cc(2:end); end
    idx = find(cc <= 0); tt(idx) = []; cc(idx) = [];
    data = struct('tobs', tt, 'cobs', cc, 'L', B.L_m, 'injectDuration', B.inject_duration_s);
    obj_function = @objFunctionSMIM; if strcmp(o.objective, 'log'), obj_function = @objFunctionSMIM_log; end
    f = @(fp) obj_function(fp, 1, data, @TPLmodel);
    opts = optimoptions(@lsqnonlin, 'FiniteDifferenceType', 'central', 'MaxIterations', 2, 'Display', 'off');
    [x, resnorm, residual, ~, ~, lambda, jacobian] = lsqnonlin(f, params_guess, params_lower, params_upper, opts);
    M = struct('params_fit', x, 'resnorm', resnorm, 'resid', residual, 'jacobian', jacobian, 'lambda', lambda, ...
               'ccfit', TPLmodel(x, data), 'tcfit', tt, 'cobs', cc);
    o.objective = [o.objective '_SMOKE'];
end
switch o.objective
    case 'volponi'
        rng(o.seed, 'twister');
        M = SMIMfit(t, ccNorm, 'c', 'L', B.L_m, 't_end', t(end), 'model_type', 'TPL', ...
                    'params_guess', params_guess, 'params_upper', params_upper, 'params_lower', params_lower);
        obj_function = @objFunctionSMIM;
    case 'log'
        % SMIMfit.m hard-wires objFunctionSMIM through createModel.m, so for the log
        % objective we reproduce SMIMfit.m lines 119-142 (data preparation) and call her
        % runFit.m directly with the new objective handle. Her core is untouched.
        tt = t; cc = ccNorm;
        idx = find(tt > t(end)); tt(idx) = []; cc(idx) = [];
        if tt(1) == 0, tt = tt(2:end); cc = cc(2:end); end
        idx = find(cc <= 0); tt(idx) = []; cc(idx) = [];
        data.tobs = tt; data.cobs = cc; data.L = B.L_m; data.injectDuration = B.inject_duration_s;
        rng(o.seed, 'twister');
        M = runFit('c', @objFunctionSMIM_log, data, @TPLmodel, params_guess, params_upper, params_lower);
        M.tcfit = tt; M.cobs = cc;
        obj_function = @objFunctionSMIM_log;
end
fit_seconds = toc(tfit);

%% ---- post-processing ---------------------------------------------------------------------
data = struct('tobs', M.tcfit, 'cobs', M.cobs, 'L', B.L_m, 'injectDuration', B.inject_duration_s);
n = numel(M.cobs);
J = full(M.jacobian);
% Lead_SMIM.m lines 110-123, with the two fixed columns removed (they are identically 0)
Jf = J(:, free);
SSWR = sum((M.cobs - M.ccfit).^2);          % her formula: unweighted, ccfit NOT renormalised
CEV_volponi = SSWR / (n - pfree);
Cov_volponi = CEV_volponi * inv(Jf.' * Jf);
SE_volponi = nan(1, 6); SE_volponi(free) = sqrt(diag(Cov_volponi))';
% internally consistent alternative (diagnostic only)
CEV_consistent = M.resnorm / (n - pfree);
Cov_consistent = CEV_consistent * inv(Jf.' * Jf);
SE_consistent = nan(1, 6); SE_consistent(free) = sqrt(diag(Cov_consistent))';

ccfit_norm = M.ccfit / trapz(M.tcfit, M.ccfit);   % what the objective actually compares
r2_lin = 1 - sum((M.cobs - ccfit_norm).^2) / sum((M.cobs - mean(M.cobs)).^2);
lo = log(M.cobs); lf = log(max(ccfit_norm, realmin));
r2_log = 1 - sum((lo - lf).^2) / sum((lo - mean(lo)).^2);
wmae = sum(M.resid);
bound_hit = (abs(M.params_fit(free) - params_lower(free)) < 1e-8 * max(1, abs(params_lower(free)))) | ...
            (abs(M.params_fit(free) - params_upper(free)) < 1e-8 * max(1, abs(params_upper(free))));

%% ---- exports ------------------------------------------------------------------------------
write17(fullfile(out_dir, 'params_fit.csv'), {'v','D','Lambda','beta','logT1','logT2'}, M.params_fit(:)');
write17(fullfile(out_dir, 'se_volponi.csv'),    {'v','D','Lambda','beta','logT1','logT2'}, SE_volponi);
write17(fullfile(out_dir, 'se_consistent.csv'), {'v','D','Lambda','beta','logT1','logT2'}, SE_consistent);
write17(fullfile(out_dir, 'curves.csv'), {'tobs','cobs','ccfit_raw','ccfit_norm','resid'}, ...
        [M.tcfit(:) M.cobs(:) M.ccfit(:) ccfit_norm(:) M.resid(:)]);
write17(fullfile(out_dir, 'jacobian.csv'), {'v','D','Lambda','beta','logT1','logT2'}, J);
write17(fullfile(out_dir, 'cleaned_input.csv'), {'t','Conc_ppb','ccNorm'}, [t(:) Conc(:) ccNorm(:)]);
write17(fullfile(out_dir, 'bounds.csv'), {'v','D','Lambda','beta','logT1','logT2'}, [params_lower; params_upper; params_guess]);

S = struct();
S.btc_id = char(B.btc_id); S.objective = o.objective; S.t_end_mode = tag; S.t_end_s = t_end;
S.seed = o.seed; S.L_m = B.L_m; S.Q_L_s = B.Q_L_s; S.mass_ug = B.mass_ug; S.hit_time_s = B.hit_time_s;
S.n_after_release = n_after_release; S.n_neg_clamped = n_neg_clamped; S.n_prehit_zeroed = n_prehit_zeroed;
S.n_in_window = numel(t); S.n_fitted = n; S.n_dropped_le0_in_fit = numel(t) - 1 - n;
S.mass_recovery_fraction = massRecoveryFraction; S.Q_dilution_gauging = Qdg;
S.resnorm = M.resnorm; S.wmae = wmae; S.r2_lin = r2_lin; S.r2_log = r2_log;
S.bound_hit_free = bound_hit; S.fit_seconds = fit_seconds;
S.matlab_version = version; S.date = datestr(now, 31);
v = ver; S.toolboxes = strjoin(arrayfun(@(k) sprintf('%s %s', v(k).Name, v(k).Version), 1:numel(v), 'UniformOutput', false), '; ');
fid = fopen(fullfile(out_dir, 'summary.json'), 'w'); fprintf(fid, '%s', jsonencode(S, 'PrettyPrint', true)); fclose(fid);

%% ---- instrumented pass (optional): every start point separately ------------------------
if o.instrumented
    fprintf('--- instrumented pass: %d starts\n', size(points, 1));
    maxiter = 10000; tol = 1e-14;                                     % runFit.m:37-39
    f = @(fp) obj_function(fp, 1, data, @TPLmodel);                    % runFit.m:42
    opts = optimoptions(@lsqnonlin, 'FunctionTolerance', tol, 'FiniteDifferenceType', 'central', ...
        'MaxFunctionEvaluations', maxiter, 'OptimalityTolerance', tol, 'StepTolerance', tol, ...
        'MaxIterations', maxiter, 'Display', o.display);               % runFit.m:44-45
    ns = size(points, 1);
    PS = nan(ns, 6 + 5);
    for k = 1:ns
        tk = tic;
        [xk, rnk, ~, efk, outk] = lsqnonlin(f, points(k, :), params_lower, params_upper, opts);
        PS(k, :) = [xk(:)' rnk efk outk.iterations outk.funcCount toc(tk)];
        fprintf('start %2d: resnorm %.10g exitflag %d iter %d feval %d (%.0f s)\n', k, rnk, efk, outk.iterations, outk.funcCount, PS(k, end));
    end
    write17(fullfile(out_dir, 'per_start.csv'), {'v','D','Lambda','beta','logT1','logT2','resnorm','exitflag','iterations','funcCount','seconds'}, PS);
    % MultiStart selection rule: lowest Fval among exitflag > 0, else lowest overall
    ok = PS(:, 8) > 0;
    if any(ok), cand = find(ok); else, cand = (1:ns)'; end
    [~, j] = min(PS(cand, 7)); best = cand(j);
    [xr, rnr] = lsqnonlin(f, PS(best, 1:6), params_lower, params_upper, opts);   % runFit.m:61 refine
    write17(fullfile(out_dir, 'instrumented_refined.csv'), {'v','D','Lambda','beta','logT1','logT2','resnorm','best_start'}, [xr(:)' rnr best]);
    fprintf('instrumented best start %d -> refined resnorm %.10g vs SMIMfit %.10g (max |dparam| = %.3g)\n', ...
            best, rnr, M.resnorm, max(abs(xr(:)' - M.params_fit(:)')));
end

%% ---- quick-look plot -------------------------------------------------------------------------
fig = figure('Visible', 'off', 'Position', [0 0 1100 420]);
subplot(1, 2, 1); plot(M.tcfit, M.cobs, 'o', 'MarkerSize', 3); hold on; plot(M.tcfit, ccfit_norm, '-', 'LineWidth', 2);
xlabel('t (s)'); ylabel('normalised C'); title(strrep(run_name, '_', ' '), 'FontSize', 8); legend('obs', 'fit');
subplot(1, 2, 2); loglog(M.tcfit, M.cobs, 'o', 'MarkerSize', 3); hold on; loglog(M.tcfit, ccfit_norm, '-', 'LineWidth', 2);
xlabel('t (s)'); ylabel('normalised C'); title(sprintf('R^2 lin %.4f  log %.4f', r2_lin, r2_log));
saveas(fig, fullfile(out_dir, 'fit.png')); close(fig);

fprintf('=== done %s: params %s resnorm %.10g fit %.0f s\n', run_name, mat2str(M.params_fit, 8), M.resnorm, fit_seconds);
diary off
R = struct('M', M, 'summary', S, 'points', points, 'out_dir', out_dir);
end

function write17(path, header, X)
% CSV with 17 significant digits so R reads bit-identical doubles.
fid = fopen(path, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fmt = [strjoin(repmat({'%.17g'}, 1, size(X, 2)), ','), '\n'];
fprintf(fid, fmt, X.');
fclose(fid);
end
