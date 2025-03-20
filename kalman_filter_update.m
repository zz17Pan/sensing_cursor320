function kf = kalman_filter_update(kf, z, params)
% KALMAN_FILTER_UPDATE 更新卡尔曼滤波器状态
%   kf = KALMAN_FILTER_UPDATE(kf, z, params) 使用最新观测向量更新卡尔曼滤波器
%   z = [R; theta; phi] - 观测向量（距离，方位角，俯仰角，弧度制）

% 显示输入观测，确保单位和格式正确
fprintf('卡尔曼滤波器更新:\n');
fprintf('  观测: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    z(1), z(2)*180/pi, z(3)*180/pi);

% 保存上一个状态以备调试
x_prev = kf.x;

% ===== 预测步骤 =====
% 1. 状态预测
x_pred = kf.F * kf.x;

% 2. 协方差预测
P_pred = kf.F * kf.P * kf.F' + kf.Q;

% 确保协方差矩阵是对称的
P_pred = (P_pred + P_pred')/2;

% ===== 计算当前状态下的预测位置和角度 =====
x = x_pred(1); % x位置
y = x_pred(4); % y位置
z_pos = x_pred(7); % z位置
vx = x_pred(2); % x速度
vy = x_pred(5); % y速度
vz = x_pred(8); % z速度

% 计算预测的距离
r_pred = sqrt(x^2 + y^2 + z_pos^2);

% 计算预测的方位角，标准定义：逆时针，相对x轴
theta_pred = atan2(y, x);

% 计算预测的俯仰角
phi_pred = atan2(z_pos, sqrt(x^2 + y^2));

% 报告预测值
fprintf('  状态预测: 位置=[%.2f, %.2f, %.2f]m, 速度=[%.2f, %.2f, %.2f]m/s\n', ...
    x, y, z_pos, vx, vy, vz);
fprintf('  预测观测: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    r_pred, theta_pred*180/pi, phi_pred*180/pi);

% ===== 计算雅可比矩阵 =====
% 观测方程对状态的偏导数
dr_dx = x / r_pred;
dr_dy = y / r_pred;
dr_dz = z_pos / r_pred;

dtheta_dx = -y / (x^2 + y^2);
dtheta_dy = x / (x^2 + y^2);
dtheta_dz = 0;

dphi_dx = -x * z_pos / (r_pred^2 * sqrt(x^2 + y^2));
dphi_dy = -y * z_pos / (r_pred^2 * sqrt(x^2 + y^2));
dphi_dz = sqrt(x^2 + y^2) / r_pred^2;

% 构建雅可比矩阵 (3x9)
H = zeros(3, 9);
H(1, 1) = dr_dx;     H(1, 4) = dr_dy;     H(1, 7) = dr_dz;
H(2, 1) = dtheta_dx; H(2, 4) = dtheta_dy; H(2, 7) = dtheta_dz;
H(3, 1) = dphi_dx;   H(3, 4) = dphi_dy;   H(3, 7) = dphi_dz;

% ===== 更新步骤 =====
% 0. 构建观测向量和预测观测
z_pred = [r_pred; theta_pred; phi_pred];

% 1. 计算创新向量（测量残差）
innovation = z - z_pred;

% 1.1 处理角度环绕（确保角度差在-pi和pi之间）
innovation(2) = wrapToPi(innovation(2));
innovation(3) = wrapToPi(innovation(3));

% 1.2 检查创新向量是否有异常值，限制创新大小
max_innovation = [10.0; 20*pi/180; 20*pi/180]; % 最大允许创新量：距离10m，角度20度
innovation = min(max(innovation, -max_innovation), max_innovation);

fprintf('  创新向量: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    innovation(1), innovation(2)*180/pi, innovation(3)*180/pi);

% 2. 计算创新协方差
S = H * P_pred * H' + kf.R;

% 2.1 确保S是对称正定的
S = (S + S')/2;
[V, D] = eig(S);
D = max(D, 1e-6 * eye(size(D))); % 确保最小特征值大于零
S = V * D * V';

% 3. 计算卡尔曼增益
K = P_pred * H' / S; % 直接使用矩阵除法，更加稳定

% 4. 更新状态
kf.x = x_pred + K * innovation;

% 5. 更新状态协方差矩阵 (使用Josheph稳定形式)
I = eye(size(kf.P));
kf.P = (I - K * H) * P_pred * (I - K * H)' + K * kf.R * K';

% 5.1 确保P是对称正定的
kf.P = (kf.P + kf.P')/2;
[V, D] = eig(kf.P);
D = max(D, 1e-6 * eye(size(D))); % 确保最小特征值大于零
kf.P = V * D * V';

% ===== 报告更新后的状态 =====
x = kf.x(1); 
y = kf.x(4);
z_pos = kf.x(7);
vx = kf.x(2);
vy = kf.x(5);
vz = kf.x(8);
ax = kf.x(3);
ay = kf.x(6);
az = kf.x(9);

fprintf('  更新后状态: 位置=[%.2f, %.2f, %.2f]m, 速度=[%.2f, %.2f, %.2f]m/s\n', ...
    x, y, z_pos, vx, vy, vz);
fprintf('  加速度: [%.2f, %.2f, %.2f]m/s²\n', ax, ay, az);
fprintf('  速度比较: 更新前=[%.2f, %.2f, %.2f]m/s, 更新后=[%.2f, %.2f, %.2f]m/s\n', ...
    x_prev(2), x_prev(5), x_prev(8), vx, vy, vz);

% 计算更新后的球坐标
r_est = sqrt(x^2 + y^2 + z_pos^2);
theta_est = atan2(y, x);
phi_est = atan2(z_pos, sqrt(x^2 + y^2));

fprintf('  对应球坐标: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    r_est, theta_est*180/pi, phi_est*180/pi);

end

function angle = wrapToPi(angle)
% 将角度包装到-pi到pi的范围内
    angle = mod(angle + pi, 2*pi) - pi;
end 
