% =========================================================================
% LJ-SUITE-2D: ENSEMBLE ANISOTROPY (The "Noise Killer")
% =========================================================================
clear; close all; clc;

% --- 1. SETUP ---
num_runs = 20;          % <--- THE KEY: Average over 20 independent universes
N_side   = 20;          % 400 particles
rho      = 0.7;         
T        = 0.8;           
D0       = 0.05;       
dt       = 0.002;      
t_max    = 40.0;       

% Potential
sigma = 1.0; epsilon = 1.0; r_cut = 2.5*sigma; 

% Analysis Settings
n_projections = 360;    % 1 degree resolution
steps = ceil(t_max/dt);
rp=sigma/2;
dr=sqrt(4*D0*dt);
dr_save_interval=0.1*2*rp;
dt_save_interval=dr_save_interval^2/(4*D0);
save_interval=floor(dt_save_interval/dt);
nsave=ceil(steps/save_interval);
time_axis = (0:n_save-1) * dt * save_interval;

% Initialize Accumulators
% We will sum the MSD(theta, t) over all runs here
Sum_MSD_Proj = zeros(n_projections, n_save);

% Pre-calculate projection vectors
angles = linspace(0, 2*pi, n_projections);
unit_vecs = [cos(angles)', sin(angles)'];

% Indices for [1,0] and [1,1]
[~, idx_10] = min(abs(angles - 0));
[~, idx_11] = min(abs(angles - pi/4));

% --- 2. ENSEMBLE LOOP ---
fprintf('Starting Ensemble Average over %d runs...\n', num_runs);
total_timer = tic;

for run = 1:num_runs
    fprintf('  Run %d / %d ... ', run, num_runs);
    
    % --- A. Run Simulation ---
    N = N_side^2; L = sqrt(N/rho);
    [x,y] = meshgrid(linspace(0.5,L-0.5,N_side), linspace(0.5,L-0.5,N_side));
    pos = [x(:),y(:)] + 0.01*randn(N,2); pos = mod(pos,L);
    pos_unwrapped = pos;
    
    traj_unwrapped = zeros(N, 2, n_save);
    
    for t = 1:steps
        dx = pos(:,1)-pos(:,1)'; dx = dx-L*round(dx/L);
        dy = pos(:,2)-pos(:,2)'; dy = dy-L*round(dy/L);
        r2 = dx.^2+dy.^2; r2(1:N+1:end)=Inf;
        mask = r2 < r_cut^2;
        r2_c = r2(mask); r2_c(r2_c<0.25)=0.25;
        inv_r2 = 1./r2_c; inv_r6 = inv_r2.^3;
        f_mag = (24*epsilon./r2_c).*(2*inv_r6.^2 - inv_r6);
        f_mag(f_mag>5000)=5000;
        fx=zeros(N,N); fy=zeros(N,N);
        fx(mask)=f_mag.*dx(mask); fy(mask)=f_mag.*dy(mask);
        
        disp = (D0/T)*[sum(fx,2), sum(fy,2)]*dt + sqrt(2*D0*dt)*randn(N,2);
        pos = pos + disp;
        pos_unwrapped = pos_unwrapped + disp;
        pos = mod(pos, L);
        
        if mod(t,save_interval)==0 
            traj_unwrapped(:,:,t/save_interval)=pos_unwrapped; 
        end
    end
    
    % --- B. Project and Accumulate ---
    % Calculate Dr for this run
    Dr = traj_unwrapped - traj_unwrapped(:,:,1); % (N, 2, Frames)
    
    % Vectorized Projection
    % Reshape Dr to (N*Frames, 2) for speed
    Dr_flat = reshape(permute(Dr, [1 3 2]), [], 2);
    
    for k = 1:n_projections
        u = unit_vecs(k, :);
        % Project: (Dr . u)^2
        proj_sq = (Dr_flat * u').^2;
        % Reshape back to (N, Frames) and average over N
        proj_sq_mat = reshape(proj_sq, N, n_save);
        msd_t = mean(proj_sq_mat, 1);
        
        % Add to accumulator
        Sum_MSD_Proj(k, :) = Sum_MSD_Proj(k, :) + msd_t;
    end
    
    fprintf('Done.\n');
end
toc(total_timer);

% --- 3. FINALIZE STATISTICS ---
Avg_MSD_Proj = Sum_MSD_Proj / num_runs;

% Calculate Mean over all angles and Residuals
Global_Mean_MSD = mean(Avg_MSD_Proj, 1);
Residuals = Avg_MSD_Proj - Global_Mean_MSD;

% =========================================================================
% 4. VISUALIZATION
% =========================================================================
figure('Color','w', 'Position', [100 100 1200 600]);

% Plot 1: Ensemble Residuals
subplot(1,2,1);
plot(time_axis, Residuals(idx_10, :), 'r-', 'LineWidth', 2); hold on;
plot(time_axis, Residuals(idx_11, :), 'b-', 'LineWidth', 2);
yline(0, 'k--');
title(sprintf('Ensemble Average (%d Runs)', num_runs));
legend('[1,0] (Box Axis)', '[1,1] (Diagonal)', 'Mean Ref');
xlabel('Time \tau'); ylabel('Average \Delta MSD');
grid on;
subtitle('If this stays flat at 0, PBCs have NO effect on D_s');

% Plot 2: Mobility Landscape (End of Run)
subplot(1,2,2);
final_residuals = Residuals(:, end);
percent_dev = 100 * final_residuals / Global_Mean_MSD(end);

polarplot(angles, percent_dev, 'k-', 'LineWidth', 1);
hold on;
fill(angles(percent_dev>0), percent_dev(percent_dev>0), 'r', 'FaceAlpha', 0.4);
fill(angles(percent_dev<0), percent_dev(percent_dev<0), 'b', 'FaceAlpha', 0.4);

title('Mobility Landscape (Ensemble Averaged)');
subtitle(sprintf('Deviation from Mean (%%). Run Time: %.1f', t_max));
thetaticks(0:45:315);