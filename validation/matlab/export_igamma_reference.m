% export_igamma_reference.m
% Reference values of MATLAB's igamma(a, z) (Symbolic Math Toolbox, upper incomplete
% gamma, non-regularised) for negative real order a = -beta and complex z, on the kind of
% grid LapPsiFuncTPL.m uses. Consumed by tests/testthat/test-igamma.R to validate the
% R implementation igamma_complex() to < 1e-12 relative error.
%
% z = tau + T1*s with tau = T1/T2, s = gamma + i*pi*k/T (de Hoog), for the time scales
% of our BTCs (T = 2*max(t) with max t in 10^1 .. 10^4) plus a few generic points.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(repo_root, 'validation', 'matlab_out', 'igamma_reference.csv');

betas = [0 1e-9 1e-6 1e-3 0.1 0.27646 0.5 0.78 0.9 0.99 0.999999 1 1.000001 1.01 1.3 1.5 1.9 1.999999 2];
T1 = 1.5; T2 = 14000; tau = T1 / T2;
rows = [];
for Tmax = [19.9 99 999 2820 9999 13500]
    T = 2 * Tmax; g = -log(1e-10) / (2 * T);
    s = g + 1i * pi * (0:50)' / T;
    z = tau + T1 * s;
    for b = betas
        y = igamma(-b, z);
        rows = [rows; repmat(b, numel(z), 1) real(z) imag(z) real(y) imag(y)]; %#ok<AGROW>
    end
end
% generic points (also |z| ~ 1 and larger, and real z)
zg = [0.001; 0.01; 0.1; 0.5; 1; 2; 5; 0.001+0.05i; 0.1+1i; 1+1i; 2-3i; 0.05i; 0.5i; 5i];
for b = betas
    y = igamma(-b, zg);
    rows = [rows; repmat(b, numel(zg), 1) real(zg) imag(zg) real(y) imag(y)]; %#ok<AGROW>
end
fid = fopen(out, 'w');
fprintf(fid, 'beta,re_z,im_z,re_G,im_G\n');
fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g\n', rows.');
fclose(fid);
fprintf('wrote %d rows to %s\n', size(rows, 1), out);
