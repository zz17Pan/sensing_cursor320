function kf = kalman_filter_update(kf, z, params)
% KALMAN_FILTER_UPDATE 使用扩展卡尔曼滤波(EKF)更新目标状态估计
%   kf = KALMAN_FILTER_UPDATE(kf, z, params) 根据观测值z更新卡尔曼滤波器状态
%   z为观测向量 [R; theta; phi] (角度需要用弧度)

fprintf('  卡尔曼滤波更新 - 输入观测: R=%.2fm, θ=%.2f°, φ=%.2f°\n', ...
        z(1), z(2)*180/pi, z(3)*180/pi);

% ==== 预测步骤 ====
% 1. 状态预测
x_pred = kf.F * kf.x;

% 位置部分的预测值 (笛卡尔坐标)
px_pred = x_pred(1);
py_pred = x_pred(4);
pz_pred = x_pred(7);

% 从笛卡尔预测位置计算球坐标
r_pred = sqrt(px_pred^2 + py_pred^2 + pz_pred^2);
theta_pred = atan2(px_pred, py_pred);
phi_pred = atan2(pz_pred, sqrt(px_pred^2 + py_pred^2));

fprintf('  预测状态 (球坐标): R=%.2fm, θ=%.2f°, φ=%.2f°\n', ...
        r_pred, theta_pred*180/pi, phi_pred*180/pi);

% 2. 协方差预测
P_pred = kf.F * kf.P * kf.F' + kf.Q;

% ==== 更新步骤 (使用扩展卡尔曼滤波处理非线性观测) ====
% 3. 计算当前状态下的观测预测
z_pred = state_to_measurement(x_pred);

% 4. 计算观测矩阵的雅可比矩阵 (线性化)
H = compute_observation_jacobian(x_pred);

% 5. 计算卡尔曼增益
S = H * P_pred * H' + kf.R;
% 确保S是正定矩阵 (数值稳定性)
S = (S + S') / 2;  % 确保对称
[V, D] = eig(S);   % 特征分解
d = diag(D);
d(d < 1e-6) = 1e-6; % 限制最小特征值
S_stable = V * diag(d) * V';

K = P_pred * H' / S_stable;

% 6. 状态更新
innovation = z - z_pred;

% 角度环绕处理 (确保角度差在-pi到pi之间)
innovation(2) = wrapToPi(innovation(2));
innovation(3) = wrapToPi(innovation(3));

fprintf('  创新向量: ΔR=%.2fm, Δθ=%.2f°, Δφ=%.2f°\n', ...
        innovation(1), innovation(2)*180/pi, innovation(3)*180/pi);

% 限制创新幅度 (防止偏离过大)
max_innovation = [5.0; pi/6; pi/6]; % 最大允许创新: 5m, 30度, 30度
for i = 1:3
    if abs(innovation(i)) > max_innovation(i)
        innovation(i) = sign(innovation(i)) * max_innovation(i);
        fprintf('  创新值 %d 被限制到 %.2f\n', i, innovation(i));
    end
end

x_update = x_pred + K * innovation;

% 7. 协方差更新
P_update = (eye(size(P_pred)) - K * H) * P_pred;
% 确保协方差矩阵对称正定
P_update = (P_update + P_update') / 2;

% 更新卡尔曼滤波器状态
kf.x = x_update;
kf.P = P_update;
kf.H = H; % 保存当前的雅可比矩阵

% 更新后状态的球坐标
px_update = x_update(1);
py_update = x_update(4);
pz_update = x_update(7);
r_update = sqrt(px_update^2 + py_update^2 + pz_update^2);
theta_update = atan2(px_update, py_update);
phi_update = atan2(pz_update, sqrt(px_update^2 + py_update^2));

fprintf('  更新后状态 (球坐标): R=%.2fm, θ=%.2f°, φ=%.2f°\n', ...
        r_update, theta_update*180/pi, phi_update*180/pi);
fprintf('  速度估计: vx=%.2f, vy=%.2f, vz=%.2f m/s\n', ...
        x_update(2), x_update(5), x_update(8));

end

function z_pred = state_to_measurement(x)
% 将状态向量转换为观测向量 [r; theta; phi]
px = x(1);  % x位置
py = x(4);  % y位置
pz = x(7);  % z位置

% 计算距离
r = sqrt(px^2 + py^2 + pz^2);

% 计算方位角
theta = atan2(px, py);

% 计算俯仰角
rho_h = sqrt(px^2 + py^2);
phi = atan2(pz, rho_h);

% 组合观测向量
z_pred = [r; theta; phi];
end

function H = compute_observation_jacobian(x)
% 计算观测方程对状态的雅可比矩阵
px = x(1);  % x位置
py = x(4);  % y位置
pz = x(7);  % z位置

% 计算中间值
r = sqrt(px^2 + py^2 + pz^2);
rho_h = sqrt(px^2 + py^2);

% 初始化雅可比矩阵
H = zeros(3, 9);

% 距离对状态的导数
if r > 1e-6  % 避免除零
    H(1, 1) = px / r;     % dr/dx
    H(1, 4) = py / r;     % dr/dy
    H(1, 7) = pz / r;     % dr/dz
else
    H(1, 1) = 1;
    H(1, 4) = 0;
    H(1, 7) = 0;
end

% 方位角对状态的导数
if rho_h > 1e-6  % 避免除零
    H(2, 1) = py / (px^2 + py^2);      % dtheta/dx
    H(2, 4) = -px / (px^2 + py^2);     % dtheta/dy
    % dtheta/dz = 0
else
    H(2, 1) = 0;
    H(2, 4) = 0;
end

% 俯仰角对状态的导数
if r > 1e-6  % 避免除零
    H(3, 1) = -px * pz / (r^2 * rho_h);  % dphi/dx
    H(3, 4) = -py * pz / (r^2 * rho_h);  % dphi/dy
    H(3, 7) = rho_h / r^2;              % dphi/dz
else
    H(3, 1) = 0;
    H(3, 4) = 0;
    H(3, 7) = 1;
end

end 