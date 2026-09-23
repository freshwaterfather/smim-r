%===============================================================================
% objFunctionSMIM_log  --  NEW FILE (not part of Volponi's SMIMfit)
%
% Log-residual objective for lsqnonlin, with the same signature and data handling as
% Volponi's objFunctionSMIM.m so that runFit.m can use it unchanged:
%
%   f_i = log(cobs_i) - log(cfit_i),   cfit normalised to unit area over tobs by trapz
%
% Conventions copied from objFunctionSMIM.m:
%   - K_mass is unused (always 1) and kept only for the interface
%   - c_fit is divided by trapz(data.tobs, c_fit)
%   - any NaN in c_fit  -> residual vector filled with 1e16
% Addition:
%   - any c_fit <= 0 after normalisation (log undefined) -> the same 1e16 penalty vector
% cobs is strictly positive by construction (SMIMfit drops C <= 0 before fitting).
%===============================================================================

function f = objFunctionSMIM_log(params, K_mass, data, pdf_function)

cobs = data.cobs;
c_fit = K_mass .* pdf_function(params, data);
if (sum(isnan(c_fit)) ~= 0)
    f = repelem(1e16, numel(cobs));
    return
end

scalefactor = trapz(data.tobs, c_fit);
c_fit = c_fit / scalefactor;

if any(c_fit <= 0)
    f = repelem(1e16, numel(cobs));
    return
end

f = log(cobs) - log(c_fit);
end
