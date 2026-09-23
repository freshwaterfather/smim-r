% repro_her_example.m
% Phase A, step 1: run Sabrina Volponi's Lead_SMIM.m workflow, unmodified in substance,
% on her own example sheet ('example_Release1' in TemplateExample.xlsx), so that the
% MATLAB environment (Optimization, Global Optimization, Statistics, Symbolic) is proven
% before anything is adapted to our data.
%
% Differences from Lead_SMIM.m (all logged in validation/ADAPTATION_LOG.md):
%   1. No `cd(fileparts(matlab.desktop.editor.getActiveFilename))` (fails in -batch);
%      her repo is added to the path read-only instead.
%   2. Only the 'example_Release1' sheet is processed (Lead loops over all sheets,
%      which would error on the 'Instructions' and 'Template' sheets).
%   3. rng(seed,'twister') is set before SMIMfit so the run is reproducible, and the
%      identical start-point matrix is regenerated and exported.
%   4. Outputs go to an explicit output folder, never her repo.
%   5. Plot section is replaced by a PNG export (her .fig save needs a display).
% Lines 40-104 and 107-157 of Lead_SMIM.m are otherwise copied verbatim.
%
% Usage (from repo root):  matlab -batch "run('validation/matlab/repro_her_example.m')"

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
her_dir   = fullfile(repo_root, 'inputs', 'SabrinaVolponi-SMIMfit-058f87d');
out_dir   = fullfile(repo_root, 'validation', 'matlab_out', 'her_example');
addpath(her_dir);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
seed = 20260923;

%% Part 1: User definitions (Lead_SMIM.m lines 11-17, verbatim values)
    fileName = fullfile(her_dir, 'TemplateExample.xlsx');
    params_upper = [1.5 2.3 0.5 2 0.2 12];
    params_lower = [0.05 2E-10 0.0001 0  0 0.5];

%% Part 2: Run the SMIM (Lead_SMIM.m lines 27-104, sheet restricted)
    SiteDescriptions = "example_Release1";
    p = 6;
    SMIMSum = zeros(numel(SiteDescriptions),15);
    ScaledData = cell(numel(SiteDescriptions), 2);

for i=1:numel(SiteDescriptions)
    %Load in excel file data
        SensorLoc = xlsread(fileName, SiteDescriptions(i), 'F4:F4');
        Conc = xlsread(fileName, SiteDescriptions(i), 'K:K');
        ReleaseTime = xlsread(fileName, SiteDescriptions(i), 'F3:F3');
        HitTime = xlsread(fileName, SiteDescriptions(i), 'F6:F6');
        RunTime= xlsread(fileName, SiteDescriptions(i), 'B:B');
        CutTime = xlsread(fileName, SiteDescriptions(i), 'F9:F9');
        % [adaptation 6] Lead_SMIM.m:46 reads Day from column C, which is empty in her
        % example sheet (dates are text in column A). Derive day-of-month from column A.
        [~, txtA] = xlsread(fileName, SiteDescriptions(i), 'A:A');
        dA = datetime(txtA, 'InputFormat', 'M/d/yyyy');
        Day = day(dA(~isnat(dA)));
        DayStart = xlsread(fileName, SiteDescriptions(i), 'F7:F7');
        DayStop = xlsread(fileName, SiteDescriptions(i), 'F10:F10');
        Q = xlsread(fileName, SiteDescriptions(i),'F12:F12');
        TimeStep = xlsread(fileName, SiteDescriptions(i),'F13:F13');
        TracerMass = xlsread(fileName, SiteDescriptions(i), 'F5:F5');

    %Clean data
            %Remove negative values
                Conc(find(Conc < 0)) = 0;
            %Make sure vectors are the same size
                ConcSize = 1:length(Conc);
                RunTime = RunTime(ConcSize);
                Day = Day(ConcSize);
            % Set the first concentration measurement to the release start.
                IndexRelease = find((RunTime < ReleaseTime) & (Day == DayStart));
                Conc(IndexRelease) = [];
                    RunTime(IndexRelease) = [];
                    Day(IndexRelease) = [];
          % Set concentrations before the hit time to zero.
                IndexHit = find((RunTime < HitTime)&(Day == DayStart));
                Conc(IndexHit) = 0;
          % Rewrite the time vector such that the first time is zero
                t = (0:TimeStep:((TimeStep*length(Conc))-1))';
          % Remove data past the cutoff time.
                IndexCut = find(((RunTime > CutTime)& (Day == DayStop)) | (Day > DayStop));
                Conc(IndexCut) = [];
                t(IndexCut) = [];
                Day(IndexCut) = [];

    %Running the SMIM
          % Normalize the data
            if ~isempty(Q)
                [ccNorm, massRecoveryFraction, Qdg] = cNorm(t, Conc, 'c', 'cMass', TracerMass, 'Q', Q);
            else
                 [ccNorm, massRecoveryFraction, Qdg] = cNorm(t, Conc, 'c', 'cMass', TracerMass);
            end
          %Initial guess for optimization
          params_guess  = [(params_upper + params_lower)/2];

          % [adaptation 3] reproducible start points; export the identical matrix
            rng(seed, 'twister');
            points = lhsdesignbnd(60, 6, params_lower, params_upper);
            points(61, :) = params_guess;
            writematrix(points, fullfile(out_dir, 'starts.csv'));
            rng(seed, 'twister');

         % Run the model and save output
            tfit = tic;
            ModCon =SMIMfit(t, ccNorm, 'c', 'L', SensorLoc, 't_end', t(length(t)),'model_type', 'TPL', 'params_guess', params_guess, 'params_upper', params_upper, 'params_lower', params_lower)
            fit_seconds = toc(tfit);
            save(fullfile(out_dir, append(SiteDescriptions(i), '.mat')));

   % Calculate Covariance matrix (Lead_SMIM.m lines 107-123)
               Jacobian = full(ModCon.jacobian);
               SSWR = sum((ModCon.cobs - ModCon.ccfit).^2);
               n = numel(ModCon.cobs);
               CEV = SSWR/(n-p);
               Covariance = CEV*inv(Jacobian.'*Jacobian);

  % Save Model Data (Lead_SMIM.m lines 127-157)
            SMIMSum(i,1) = ModCon.params_fit(1); %U
            SMIMSum(i,2) = ModCon.params_fit(2); % D
            SMIMSum(i,3) = ModCon.params_fit(3); %Lambda
            SMIMSum(i,4) = ModCon.params_fit(4); %Beta
            SMIMSum(i,5) = ModCon.params_fit(5); %logT1
            SMIMSum(i,6) = ModCon.params_fit(6); %logT2
            SMIMSum(i, 7) = sqrt(Covariance(1,1));
            SMIMSum(i, 8) = sqrt(Covariance(2,2));
            SMIMSum(i, 9) = sqrt(Covariance(3,3));
            SMIMSum(i, 10) = sqrt(Covariance(4,4));
            SMIMSum(i,11) = sqrt(Covariance(5,5));
            SMIMSum(i,12) = sqrt(Covariance(6,6));
            SMIMSum(i,13) = sum(ModCon.resid);
            SMIMSum(i,14) = massRecoveryFraction;
            SMIMSum(i,15) = Qdg;

            SMIMTable = table(SiteDescriptions', ...
                SMIMSum(:,1), SMIMSum(:,2), SMIMSum(:,3), SMIMSum(:,4), SMIMSum(:,5), ...
                SMIMSum(:,6), SMIMSum(:,7), SMIMSum(:,8), SMIMSum(:,9), SMIMSum(:,10), SMIMSum(:,11), ...
                SMIMSum(:,12), SMIMSum(:,13), SMIMSum(:,14), SMIMSum(:,15));
            SMIMTable.Properties.VariableNames(1:16) = {'SiteDescription', 'U', 'D','Lambda', ...
                'Beta', 'logT1', 'logT2', 'SE_U', 'SE_D', 'SE_Lambda', 'SE_B', 'SE_logT1', 'SE_logT2', ...
                'WMAE', 'MassRecoveryFractionQAv', 'QEstimated'};
            writetable(SMIMTable, fullfile(out_dir, 'SMIMResults_repro.csv'))

  % Extra exports for comparison with her SMIMResults.xls / example_Release1.mat
            writematrix([ModCon.tcfit(:) ModCon.cobs(:) ModCon.ccfit(:) ModCon.resid(:)], fullfile(out_dir, 'curves.csv'));
            writematrix(ModCon.params_fit(:)', fullfile(out_dir, 'params_fit.csv'));
            fid = fopen(fullfile(out_dir, 'run_info.txt'), 'w');
            fprintf(fid, 'seed=%d\r\nresnorm=%.17g\r\nfit_seconds=%.1f\r\nn_obs=%d\r\nmatlab=%s\r\n', seed, ModCon.resnorm, fit_seconds, n, version);
            fclose(fid);

 %  Plot (Lead_SMIM.m lines 162-216, PNG instead of .fig)
                hitIndex = find(Conc, 1, 'first');
                concMass = Conc; tMass = t;
                concMass(1:(hitIndex-1)) = []; tMass(1:(hitIndex-1)) = [];
            [predPeak, predPeakIndex] = max(ModCon.ccfit);
            [obsPeak, obsPeakIndex] = max(concMass);
            scalefactor = (obsPeak/predPeak);
            scaledPred = ModCon.ccfit .* scalefactor;
            NlogFig = figure('Visible', 'off', 'Position', [0 0 1100 450]);
            subplot(1,2,1); plot(ModCon.tcfit, scaledPred, 'linewidth', 3); hold on
            plot(tMass, concMass, 'o', 'markersize', 4); xlabel('Time (s)'); ylabel('Tracer RWT Concentration [mg/L]'); legend('SMIM Predicted','Observed')
            subplot(1,2,2); plot(ModCon.tcfit, scaledPred, 'linewidth', 3); hold on
            plot(tMass, concMass, 'o', 'markersize', 4); set(gca, 'YScale', 'log'); xlabel('Time (s)'); ylabel('Tracer RWT Concentration [mg/L]');
            saveas(NlogFig, fullfile(out_dir, append(SiteDescriptions(i), '.png')));
            close(NlogFig)
end
fprintf('\nDONE her example: fit_seconds=%.1f\n', fit_seconds);
