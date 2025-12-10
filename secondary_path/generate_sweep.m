function [sweepSig, invFilter, info] = generate_sweep(cfg)
% generate_sweep
% 产生 ESS (指数正弦扫频) 及其逆滤波器 (Farina 方法)
% 必要字段: fs,T,f1,f2,padLeading,padTrailing,amplitude
fs = cfg.fs; T = cfg.T; f1 = cfg.f1; f2 = cfg.f2;
N  = round(T*fs);
t  = (0:N-1)'/fs;

L = log(f2/f1);
K = T / L;
phase = 2*pi*f1*K*(exp(t*L/T)-1);
sweepCore = sin(phase);

sweepCore = sweepCore / max(abs(sweepCore)+1e-12) * cfg.amplitude;

padLead  = zeros(round(cfg.padLeading*fs),1);
padTrail = zeros(round(cfg.padTrailing*fs),1);
sweepSig = [padLead; sweepCore; padTrail];

% 逆滤波器 (简化 Farina)
w = exp(t*L/T);
invCore = flipud(sweepCore) .* w;
invCore = invCore / max(abs(invCore)+1e-12);
invFilter = [invCore; zeros(length(padLead)+length(padTrail),1)];

info = struct();
info.Ncore = N;
info.fullLength = length(sweepSig);
info.f1 = f1; info.f2 = f2;
info.padLeadingSamples  = numel(padLead);
info.padTrailingSamples = numel(padTrail);
info.amplitude = cfg.amplitude;
info.fs = fs;
end