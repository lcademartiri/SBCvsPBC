
    % Uses Newtonian Dynamics (Momentum Conserving) to reveal 
    % the Lattice Anisotropy in collective relaxation.

    % --- Simulation Parameters ---
    N = 900;          % Larger system for better hydrodynamics
    rho = 0.8;        % Higher density (Liquid state, not gas)
    L = sqrt(N/rho);
    dt = 0.002;       % Smaller timestep for NVE stability
    
    % Equilibration vs Production
    n_equil = 5000;
    n_prod = 10000;
    
    % WCA Potential
    rc2 = (2^(1/6))^2;
    
    % --- Initialization ---
    fprintf('Initializing NVE MD (N=%d, L=%.2f)...\n', N, L);
    
    % Lattice Start
    n_side = ceil(sqrt(N));
    x = linspace(-L/2, L/2, n_side+1); x = x(1:end-1);
    [X, Y] = meshgrid(x, x);
    pos = [X(:), Y(:)]; pos = pos(1:N, :);
    pos_real = pos; % Unwrapped coordinates
    
    % Initial Velocities (T=1.0)
    vel = randn(N, 2);
    vel = vel - mean(vel); % Zero total momentum
    T_target = 1.0;
    vel = vel * sqrt(T_target / (0.5*mean(sum(vel.^2, 2))));
    
    acc = get_acc(pos, L, rc2);
    
    % --- 1. Equilibration (Berendsen-like scaling) ---
    fprintf('Equilibrating...\n');
    for step = 1:n_equil
        vel = vel + 0.5 * acc * dt;
        pos = pos + vel * dt;
        pos_real = pos_real + vel * dt;
        pos = mod(pos + L/2, L) - L/2;
        acc = get_acc(pos, L, rc2);
        vel = vel + 0.5 * acc * dt;
        
        % Velocity Rescaling (Simple Thermostat)
        if mod(step, 50) == 0
            T_curr = 0.5 * mean(sum(vel.^2, 2));
            scale = sqrt(T_target / T_curr);
            vel = vel * scale;
        end
    end
    
    % --- 2. Production (NVE - Momentum Conserving) ---
    fprintf('Running Production (Hydrodynamics On)...\n');
    
    n_frames = 250;
    save_freq = 40; % Save every 40 steps
    traj_real = zeros(N, 2, n_frames);
    
    frame = 1;
    for step = 1:n_prod
        % Standard Velocity Verlet (No Thermostat!)
        vel = vel + 0.5 * acc * dt;
        dr = vel * dt;
        pos = pos + dr;
        pos_real = pos_real + dr; % Track Unwrapped
        
        pos = mod(pos + L/2, L) - L/2;
        
        acc_new = get_acc(pos, L, rc2);
        vel = vel + 0.5 * acc_new * dt;
        acc = acc_new;
        
        if mod(step, save_freq) == 0 && frame <= n_frames
            traj_real(:,:,frame) = pos_real;
            frame = frame + 1;
        end
    end
    
    % --- Analysis: Windowed F(k,t) ---
    fprintf('Analyzing Anisotropy...\n');
    
    % We analyze at the first diffraction peak (structural) or low k (hydro)
    % Let's look at k = 2pi/L (The Box Mode)
    k_mag = 2*pi/L;       
    
    % Use a large circular window
    R_window = L/2 * 0.95; 
    
    angles = 0:10:360;    
    tau_list = zeros(size(angles));
    t_axis = (0:n_frames-1) * dt * save_freq;
    
    % Initial Config
    pos0 = traj_real(:,:,1);
    mask0 = sum(pos0.^2, 2) <= R_window^2;
    
    for i = 1:length(angles)
        theta = deg2rad(angles(i));
        k_vec = k_mag * [cos(theta), sin(theta)];
        
        F_t = zeros(1, n_frames);
        rho0 = sum(exp(1i * (pos0(mask0,:) * k_vec')));
        
        for f = 1:n_frames
            pos_t = traj_real(:,:,f);
            rho_t = sum(exp(1i * (pos_t(mask0,:) * k_vec')));
            F_t(f) = real(rho_t * conj(rho0));
        end
        
        % Normalize
        F_t = F_t / F_t(1);
        
        % Better Decay Extraction: Interpolate to find 1/e crossing
        % This handles the "Integer Frame" quantization issue
        below_idx = find(F_t < 1/exp(1), 1);
        
        if isempty(below_idx)
            tau_list(i) = t_axis(end); % Decay too slow
        elseif below_idx == 1
            tau_list(i) = 0; % Decay too fast
        else
            % Linear Interpolation for sub-frame accuracy
            y1 = F_t(below_idx-1); t1 = t_axis(below_idx-1);
            y2 = F_t(below_idx);   t2 = t_axis(below_idx);
            target = 1/exp(1);
            
            % Slope
            m = (y2 - y1) / (t2 - t1);
            tau_list(i) = t1 + (target - y1) / m;
        end
    end
    
    % --- Plotting ---
    figure('Position', [100, 100, 900, 500], 'Color', 'w');
    
    % Smooth the data slightly for visualization of the trend
    % (MD is noisy without ensemble averaging)
    tau_smooth = smooth(tau_list, 3);
    
    plot(angles, tau_list, 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4); hold on;
    plot(angles, tau_smooth, '-', 'LineWidth', 2.5, 'Color', [0 0.447 0.741]);
    
    title('\textbf{Dynamical Anisotropy in NVE Liquid}', 'Interpreter', 'latex', 'FontSize', 16);
    subtitle(['Momentum Conserving Dynamics + Circular Window'], 'Interpreter', 'latex');
    xlabel('Angle (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Relaxation Time $\tau(\theta)$', 'Interpreter', 'latex', 'FontSize', 14);
    
    grid on; box on; xlim([0 360]); xticks(0:45:360);
    
    % Annotations
    y_avg = mean(tau_smooth);
    y_scale = max(tau_smooth) - min(tau_smooth);
    ylim([y_avg - y_scale*1.5, y_avg + y_scale*1.5]);
    
    set(gcf, 'InvertHardcopy', 'off');


function acc = get_acc(pos, L, rc2)
    N = size(pos, 1);
    acc = zeros(N, 2);
    % Vectorized Force Calc
    for i = 1:N
        d = bsxfun(@minus, pos, pos(i,:));
        d = d - L*round(d/L); 
        r2 = sum(d.^2, 2);
        mask = r2 < rc2 & r2 > 0;
        if any(mask)
            r2_in = r2(mask);
            d_in = d(mask, :);
            inv_r2 = 1./r2_in;
            inv_r6 = inv_r2.^3;
            f_mag = 48 * (inv_r6.^2 - 0.5*inv_r6) .* inv_r2;
            acc(i,:) = acc(i,:) + sum(bsxfun(@times, f_mag, d_in), 1);
        end
    end
end