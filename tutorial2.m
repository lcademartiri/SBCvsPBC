% =========================================================================
% LJ-SUITE-2D: ANISOTROPY MAP (CLEANED via SYMMETRY)
% =========================================================================
clear; close all; clc;

% --- 1. SETUP (Longer Run for Better Low-Q Stats) ---
N_side = 40;            
rho    = 0.7;           
T      = 1.0;           
D0     = 0.05;       
dt     = 0.002;      
t_max  = 50.0;       % INCREASED: Run longer to catch slow modes
n_grid = 6; 

% Potential
sigma = 1.0; epsilon = 1.0; r_cut = 1.12*sigma; 

% --- 2. SIMULATION ---
fprintf('--- INITIALIZING LONG RUN (Anisotropy) ---\n');
N = N_side^2; L = sqrt(N/rho); steps = ceil(t_max/dt);
[x,y] = meshgrid(linspace(0.5,L-0.5,N_side), linspace(0.5,L-0.5,N_side));
pos = [x(:),y(:)] + 0.05*randn(N,2); pos = mod(pos,L);

% Save less often (we only need low frequency data for low q)
n_save = 2000; save_interval = floor(steps/n_save);
traj = zeros(N,2,n_save);

fprintf('Running Dynamics (N=%d)...\n', N);
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
    pos = pos + (D0/T)*[sum(fx,2), sum(fy,2)]*dt + sqrt(2*D0*dt)*randn(N,2);
    pos = mod(pos, L);
    if mod(t,save_interval)==0, traj(:,:,t/save_interval)=pos; end
    if mod(t,10)==0
        disp(100*t/steps)
    end
end

% --- 3. MAPPING (Raw Calculation) ---
fprintf('\n--- MAPPING & SYMMETRIZING ---\n');
dq = 2*pi/L;
range = -n_grid:n_grid;
raw_Deff = nan(length(range), length(range));
center = n_grid + 1;

for i = 1:length(range)
    nx = range(i);
    for j = 1:length(range)
        ny = range(j);
        if nx==0 && ny==0, continue; end
        
        qx = nx * dq; qy = ny * dq; qmag = sqrt(qx^2 + qy^2);
        
        % Correlation
        term = exp(-1i * (qx*traj(:,1,:) + qy*traj(:,2,:)));
        rho_t = squeeze(sum(term,1))';
        d_rho = rho_t - mean(rho_t);
        ac = xcorr(d_rho, floor(n_save/4), 'unbiased');
        center_idx = length(ac)/2 + 0.5;
        F_qt = real(ac(center_idx:end)); F_qt = F_qt / F_qt(1);
        
        idx_decay = find(F_qt < 0.367, 1);
        if isempty(idx_decay), idx_decay = length(F_qt); end
        decay_time = idx_decay * dt * save_interval;
        
        raw_Deff(j,i) = (1/decay_time) / qmag^2;
    end
end

% --- 4. SYMMETRIZATION (The "Noise Cleaner") ---
% We enforce the physical law that D(nx, ny) == D(ny, nx) == D(-nx, ny)
clean_Deff = raw_Deff; 

for i = 1:n_grid+1 % Loop over positive quadrant
    nx_idx = center + (i-1);
    for j = 1:n_grid+1
        ny_idx = center + (j-1);
        
        % Collect all 8 symmetric partners
        vals = [];
        % Indices for +/- nx and +/- ny
        ix_pos = center + (i-1); ix_neg = center - (i-1);
        iy_pos = center + (j-1); iy_neg = center - (j-1);
        
        indices = [ix_pos, iy_pos; ix_neg, iy_pos; ix_pos, iy_neg; ix_neg, iy_neg; ...
                   iy_pos, ix_pos; iy_neg, ix_pos; iy_pos, ix_neg; iy_neg, ix_neg];
               
        for k = 1:8
            val = raw_Deff(indices(k,2), indices(k,1)); % Swap Row/Col for matrix indexing
            if ~isnan(val), vals = [vals, val]; end
        end
        
        avg_val = mean(vals);
        
        % Write back to all 8 partners
        for k = 1:8
            clean_Deff(indices(k,2), indices(k,1)) = avg_val;
        end
    end
end

% --- 5. VISUALIZATION ---
figure('Color','k', 'Position', [100 100 1000 500]);

% Plot Raw (Noisy)
subplot(1,2,1);
pcolor(range, range, raw_Deff);
axis square; colorbar; shading flat;
title('Raw Data (Noisy)');
xlabel('n_x'); ylabel('n_y');

% Plot Symmetrized (Clean)
subplot(1,2,2);
h = pcolor(range, range, clean_Deff);
set(h, 'EdgeColor', 'k', 'LineWidth', 0.5); % Grid lines
axis square; colorbar;
title('Symmetrized Data (Clean)');
xlabel('n_x'); ylabel('n_y');
subtitle('Averaged over 8 symmetries');

% Force the color scale to highlight the center structure
caxis([0 prctile(clean_Deff(:), 95)]); 
colormap(gca, 'jet');

% Extract line profiles
D10 = clean_Deff(center+1:center+n_grid-2, center);      % nx = 1:n_grid, ny=0
D11 = diag(clean_Deff(center+1:center+n_grid-2, center+1:center+n_grid-2));

q10 = (1:n_grid-2) * dq;
q11 = (1:n_grid-2) * dq * sqrt(2);
figure
scatter(q10, D10, 30,'filled');
hold
scatter(q11, D11, 30,'filled');
legend('[10]','[11]');
xlabel('|q|'); ylabel('D_{eff}(q)');
xscale log
yscale log

