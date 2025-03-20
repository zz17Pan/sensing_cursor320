function tx_signal = generate_fmcw_signal(params)
% GENERATE_FMCW_SIGNAL 生成FMCW发射信号
%   tx_signal = GENERATE_FMCW_SIGNAL(params) 根据参数生成FMCW信号
%   FMCW信号：s_tx(t) = A * exp(j2π(f_c*t + 0.5*μ*t^2))

% 计算每个chirp的时间点
t_chirp = (0:params.N_samples-1) / params.fs;

% 计算FMCW信号的相位
phase = 2 * pi * (params.fc * t_chirp + 0.5 * params.mu * t_chirp.^2);

% 生成复信号
tx_chirp = exp(1j * phase);

% 对所有chirp重复
tx_signal = repmat(tx_chirp, params.N_chirps, 1);

end 