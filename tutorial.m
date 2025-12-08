% =========================================================================
% LJ-SUITE-2D: LOW-Q INVESTIGATION (Large Box, No Binning)
% =========================================================================
clear; close all; clc;

% =========================================================================
% 1. LABORATORY PARAMETERS
% =========================================================================

make_movie   = false;
vis_interval = 10;   % Update less often because N is large

% --- System State ---
rho    = 0.30;       
N_side = 30;         % <--- CHANGE: 1600 particles (Big Box = Low q)
T      = 0.8;        

% --- Dynamics ---
D0     = 0.05;       
dt     = 0.002;      
t_max  = 50.0;       % Shorter time (sim is heavier)

% --- Potential ---
sigma   = 1.0;       
epsilon = 5;
r_cut   = 3 * sigma; 

% --- Analysis ---
q_max   = 15.0;       % <--- CHANGE: Focus only on Low q
corr_t  = 2.0; 

% =========================================================================
% 2. SYSTEM INITIALIZATION
% =========================================================================
fprintf('--- INITIALIZING LARGE SYSTEM ---\n');

N = N_side^2;
L = sqrt(N/rho);
steps = ceil(t_max / dt);
n_corr_steps = ceil(corr_t / dt);

fprintf('Particles: %d (Heavier Sim)\n', N);
fprintf('Box Size L: %.2f\n', L);
fprintf('Min Accessible q (2pi/L): %.3f\n', 2*pi/L);

% Grid Start
[x, y] = meshgrid(linspace(0.5, L-0.5, N_side), linspace(0.5, L-0.5, N_side));
pos = [x(:), y(:)];
pos = pos + 0.02 * randn(N,2); 
pos = mod(pos, L);

pos_unwrapped = pos; 

% Storage (Less frequent saving to save RAM)
n_save = 50;
save_interval = floor(steps / n_save);
if save_interval < 1, save_interval = 1; end
traj_wrapped   = zeros(N, 2, n_save);

% =========================================================================
% 3. VISUALIZATION SETUP
% =========================================================================
h_fig = figure('Color','k', 'Name', 'Large Scale Sim', 'Position', [50 300 600 600]);
ax = axes('Parent', h_fig);
axis equal; axis([0 L 0 L]); box on; hold on;
title('Running Large System...');

% Only draw dots for speed (rectangles are too slow for N=1600)
h_dots = plot(pos(:,1), pos(:,2), 'b.', 'MarkerSize', 8);

% =========================================================================
% 4. SIMULATION LOOP
% =========================================================================
fprintf('--- RUNNING SIMULATION ---\n');
frame_idx = 1;

for t = 1:steps
    % --- Optimized Force Calculation ---
    % For N=1600, full matrix is heavy. But Matlab handles it OK.
    
    dx = pos(:,1) - pos(:,1)'; 
    dx = dx - L * round(dx/L);
    
    dy = pos(:,2) - pos(:,2)'; 
    dy = dy - L * round(dy/L);
    
    r2 = dx.^2 + dy.^2;
    r2(1:N+1:end) = Inf; 
    
    mask = (r2 < r_cut^2);
    
    % Force Calculation
    r2_c = r2(mask);
    r2_c(r2_c < 0.25) = 0.25; % Floor
    
    inv_r2 = (sigma^2) ./ r2_c;
    inv_r6 = inv_r2.^3;
    f_div_r = (24 * epsilon ./ r2_c) .* (2*(inv_r6.^2) - inv_r6);
    
    fx = zeros(N,N); fy = zeros(N,N);
    fx(mask) = f_div_r .* dx(mask);
    fy(mask) = f_div_r .* dy(mask);
    
    forces = [sum(fx, 2), sum(fy, 2)];
    
    % Cap
    forces(forces > 5000) = 5000;
    forces(forces < -5000) = -5000;
    
    % --- Update ---
    noise = sqrt(2*D0*dt) * randn(N,2);
    displacement = (D0/T)*forces*dt + noise;
    
    pos_new = pos + displacement;
    pos = mod(pos_new, L);
    
    if mod(t, save_interval) == 0 && frame_idx <= n_save
        traj_wrapped(:,:,frame_idx)   = pos;
        frame_idx = frame_idx + 1;
    end
    
    if mod(t, vis_interval) == 0
        set(h_dots, 'XData', pos(:,1), 'YData', pos(:,2));
        drawnow limitrate;
    end
end

% =========================================================================
% 5. LOW-Q ANALYSIS (NO BINNING)
% =========================================================================
fprintf('\n--- ANALYZING DISCRETE MODES ---\n');

% 1. Generate ALL discrete wavevectors
dq = 2*pi/L;
n_max = floor(q_max / dq);

q_results = []; % Will store [q_mag, S(q), Deff]

% Loop over integer indices (n_x, n_y) directly
for nx = -n_max:n_max
    for ny = -n_max:n_max
        if nx==0 && ny==0, continue; end
        
        qvec = [nx, ny] * dq;
        qmag = norm(qvec);
        
        if qmag > q_max, continue; end
        
        % Calculate S(q) and Deff for this SPECIFIC mode
        rho_t = zeros(1, n_save);
        
        % Vectorized phase sum
        qx = qvec(1); qy = qvec(2);
        term = exp(-1i * (qx*traj_wrapped(:,1,:) + qy*traj_wrapped(:,2,:))); 
        rho_t = squeeze(sum(term, 1))'; 
        
        % S(q)
        S_val = mean(abs(rho_t).^2) / N;
        
        % F(q,t) by looking at the autocorrelation of the fourier
        % coefficients
        n_frames_corr = floor(corr_t / (dt * save_interval));
        d_rho = rho_t - mean(rho_t);
        ac = xcorr(d_rho, n_frames_corr, 'unbiased');
        center = length(ac)/2 + 0.5;
        F_qt = real(ac(center : center + n_frames_corr - 1));
        F_qt = F_qt / F_qt(1);
        
        % find characteristic lifetime of the F(q,t) to extract Gamma
        idx_decay = find(F_qt < 0.367, 1);
        if isempty(idx_decay), idx_decay = length(F_qt); end
        
        t_lag = (0:n_frames_corr-1) * (dt * save_interval);
        Gamma = 1 / t_lag(idx_decay);

        % calculate Deff from gamma
        Deff_val = Gamma / (qmag^2);
        
        % Store Result
        q_results = [q_results; qmag, S_val, Deff_val];
    end
end

% =========================================================================
% 6. PLOTTING THE DISCRETE LOW-Q PHYSICS
% =========================================================================
figure('Color','k', 'Name', 'Low Q Investigation', 'Position', [100 100 1000 500]);

% Plot S(q)
subplot(1,2,1);
plot(q_results(:,1), q_results(:,2), 'bo', 'MarkerSize', 4);
hold on;
xline(2*pi/L, 'r--', 'Min q (2\pi/L)');
xlabel('q'); ylabel('S(q)');
title('Discrete S(q) at Low q');
subtitle('Note the gap between 0 and the red line');
grid on;

% Plot Deff
subplot(1,2,2);
plot(q_results(:,1), q_results(:,3), 'mo', 'MarkerSize', 4);
hold on;
yline(D0, 'g--', 'D0');
xline(2*pi/L, 'r--', 'Min q');
xlabel('q'); ylabel('D_{eff}');
yscale log
title('Collective Diffusion at Low q');
subtitle('Scattered points due to low statistics of single modes');
grid on;
ylim([0 max(q_results(:,3))*1.2]);