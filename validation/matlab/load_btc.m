function B = load_btc(btc_id, data_dir)
% LOAD_BTC  Read one corrected breakthrough curve and its fit-window metadata.
%
%   B = load_btc(btc_id, data_dir)
%
% Reads <data_dir>/btc_corrected.csv and <data_dir>/fit_windows.csv, which are written
% by the single preprocessing script (private/preprocess.R for the private run;
% private/make_example_data.R for the public A-D example). Gain and background
% corrections are ALREADY applied in those files; nothing is corrected here.
%
% Returns a struct with
%   .t_s, .C_ppb       all samples (including pre-release t < 0), in time order
%   .L_m, .Q_L_s, .mass_ug, .inject_duration_s
%   .hit_time_s, .t_end_ref_s, .t_end_lod_s, .t_end_loq_s, .t_end_record_s, .fit_t_end_s
%   .sigma_ppb, .baseline_ppb
%
% This is the adapter equivalent of Lead_SMIM.m lines 40-52 (xlsread of the template).

d = readtable(fullfile(data_dir, 'btc_corrected.csv'), 'TextType', 'string', ...
              'VariableNamingRule', 'preserve');
w = readtable(fullfile(data_dir, 'fit_windows.csv'), 'TextType', 'string', ...
              'VariableNamingRule', 'preserve');

rows = d.btc_id == string(btc_id);
if ~any(rows), error('load_btc:notFound', 'No rows for btc_id "%s"', btc_id); end
d = sortrows(d(rows, :), 't_s');
wr = w(w.btc_id == string(btc_id), :);
if height(wr) ~= 1, error('load_btc:window', 'Expected 1 fit_windows row for "%s", got %d', btc_id, height(wr)); end

B.btc_id            = string(btc_id);
B.t_s               = d.t_s(:);
B.C_ppb             = d.C_ppb(:);
B.L_m               = d.reach_length_m(1);
B.Q_L_s             = d.discharge_L_s(1);
B.mass_ug           = d.tracer_mass_mg(1) * 1000;   % mg -> ug so that Q [L/s] * C [ug/L] * t [s] matches
B.inject_duration_s = d.inject_duration_s(1);
B.baseline_ppb      = d.baseline_ppb(1);
B.hit_time_s        = wr.hit_time_s;
B.t_end_ref_s       = wr.t_end_ref_s;
B.t_end_lod_s       = wr.t_end_lod_s;
B.t_end_loq_s       = wr.t_end_loq_s;
B.t_end_record_s    = wr.t_end_record_s;
B.fit_t_end_s       = wr.fit_t_end_s;
B.sigma_ppb         = wr.sigma_ppb;
end
