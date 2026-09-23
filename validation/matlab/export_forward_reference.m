function export_forward_reference(btc_ids, data_dir, out_dir, n_design, seed)
% EXPORT_FORWARD_REFERENCE  Evaluate Volponi's forward model (TPLmodel) and both
% objectives on a fixed parameter design for each BTC, for equivalence items 1-2.
%
%   export_forward_reference(btc_ids, data_dir, out_dir, n_design, seed)
%
% For every BTC: the fitted time vector is prepared exactly as SMIMfit.m does (t >= 0,
% clamp, pre-hit zero, cut at t_end_ref, cNorm, drop t == 0, drop C <= 0), then for each
% design row the following are written with 17 significant digits:
%   <id>_design.csv     parameter rows (v, D, Lambda, beta, logT1, logT2)
%   <id>_forward.csv    columns = design rows, rows = tobs: TPLmodel(params, data)
%   <id>_obj.csv        per design row: sum(f.^2) for objFunctionSMIM and objFunctionSMIM_log,
%                       plus the full residual vectors in <id>_resid_volponi.csv / _log.csv
%   <id>_tobs.csv       tobs, cobs
% The design is LHS over the Lead_SMIM bounds for the 4 free parameters (t1, t2 fixed),
% plus deliberate edge rows: beta = 0, 1, 2, 1 +/- 1e-6, the bounds mid-point, and
% v/D/Lambda at their bounds.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'inputs', 'SabrinaVolponi-SMIMfit-058f87d'));
addpath(fileparts(mfilename('fullpath')));
if nargin < 4, n_design = 40; end
if nargin < 5, seed = 4242; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

lb = [0.05 2E-10 0.0001 0  log10(1.5) log10(14000)];
ub = [1.5 2.3 0.5 2 log10(1.5) log10(14000)];
mid = (lb + ub) / 2;
rng(seed, 'twister');
X = lhsdesignbnd(n_design, 6, lb, ub);
edge = [mid; mid; mid; mid; mid; lb; ub];
edge(1, 4) = 0; edge(2, 4) = 1; edge(3, 4) = 2; edge(4, 4) = 1 - 1e-6; edge(5, 4) = 1 + 1e-6;
edge(6, 4) = 0.5; edge(7, 4) = 1.5;   % bounds rows but with an interior beta
X = [X; edge];

for i = 1:numel(btc_ids)
    id = btc_ids{i};
    B = load_btc(id, data_dir);
    t = B.t_s; C = B.C_ppb;
    keep = t >= 0; t = t(keep); C = C(keep);
    C(C < 0) = 0; C(t < B.hit_time_s) = 0;
    cut = t > B.t_end_ref_s; t(cut) = []; C(cut) = [];
    cn = cNorm(t, C, 'c', 'cMass', B.mass_ug, 'Q', B.Q_L_s);
    if t(1) == 0, t = t(2:end); cn = cn(2:end); end
    idx = cn <= 0; t(idx) = []; cn(idx) = [];
    data = struct('tobs', t, 'cobs', cn, 'L', B.L_m, 'injectDuration', B.inject_duration_s);

    F = nan(numel(t), size(X, 1)); Rv = F; Rl = F; obj = nan(size(X, 1), 2);
    for k = 1:size(X, 1)
        F(:, k) = TPLmodel(X(k, :), data);
        fv = objFunctionSMIM(X(k, :), 1, data, @TPLmodel);     Rv(:, k) = fv(:); obj(k, 1) = sum(fv.^2);
        fl = objFunctionSMIM_log(X(k, :), 1, data, @TPLmodel); Rl(:, k) = fl(:); obj(k, 2) = sum(fl.^2);
    end
    w17(fullfile(out_dir, [id '_design.csv']), {'v','D','Lambda','beta','logT1','logT2'}, X);
    w17(fullfile(out_dir, [id '_tobs.csv']), {'tobs','cobs'}, [t(:) cn(:)]);
    w17(fullfile(out_dir, [id '_forward.csv']), cellstr("d" + string(1:size(X,1))), F);
    w17(fullfile(out_dir, [id '_resid_volponi.csv']), cellstr("d" + string(1:size(X,1))), Rv);
    w17(fullfile(out_dir, [id '_resid_log.csv']), cellstr("d" + string(1:size(X,1))), Rl);
    w17(fullfile(out_dir, [id '_obj.csv']), {'ss_volponi','ss_log'}, obj);
    fprintf('%s: n=%d design=%d done\n', id, numel(t), size(X, 1));
end
end

function w17(path, header, X)
fid = fopen(path, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fprintf(fid, [strjoin(repmat({'%.17g'}, 1, size(X, 2)), ','), '\n'], X.');
fclose(fid);
end
