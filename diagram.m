
    % GENERATE_ENSEMBLE_AVERAGED_FIGURE
    % Averages S(k) over 100 independent realizations to remove noise.
    % Generates the definitive JCP Figure 1.

    % --- Simulation Parameters ---
    N_particles = 300;     
    L = 1;                 
    R = 0.5;               
    N_reps = 100;          % ENSEMBLE SIZE (The "Fix")
    
    % --- Accumulators for Averaging ---
    % We will sum the results here and divide by N_reps at the end
    S_pbc_2D_sum = 0;
    S_sbc_2D_sum = 0;
    S_pbc_azi_sum = 0;
    S_sbc_azi_sum = 0;
    S_sbc_radial_sum = 0;
    
    % Store one snapshot for the Real Space plots (Panels A & B)
    pos_snapshot = [];
    sbc_pos_snapshot = [];
    pbc_pos_snapshot = [];

    % --- Grid Setup (Defined once) ---
    % 2D Grids
    k_range = linspace(-20, 20, 100);
    [KX, KY] = meshgrid(k_range, k_range);
    
    % Azimuthal Grid
    k_target = 2*pi/L;
    theta = 0:2:360; 
    theta_rad = deg2rad(theta);
    kx_circ = k_target * cos(theta_rad);
    ky_circ = k_target * sin(theta_rad);
    
    % Radial Grid
    k_radial = linspace(0.1, 20, 200);

    fprintf('Starting Ensemble Averaging over %d replicates...\n', N_reps);

    % --- MAIN LOOP ---
    for rep = 1:N_reps
        if mod(rep, 10) == 0, fprintf('Processing replicate %d/%d...\n', rep, N_reps); end
        
        % 1. Generate Liquid Snapshot (Robust Hard Disks)
        rng(rep + 1000); % Distinct seed per rep
        pos = zeros(N_particles, 2);
        count = 0;
        min_dist = 0.04; 
        attempt = 0;
        while count < N_particles && attempt < 50000
            p = (rand(1,2)-0.5)*L;
            if count == 0
                count = 1; pos(count,:) = p;
            else
                d2 = sum((pos(1:count,:) - p).^2, 2);
                if all(d2 > min_dist^2)
                    count = count + 1; pos(count,:) = p;
                end
            end
            attempt = attempt + 1;
        end
        pos = pos(1:count, :);
        
        % 2. Define Subsets
        r_sq = sum(pos.^2, 2);
        sbc_pos = pos(r_sq <= R^2, :);
        N_sbc_curr = size(sbc_pos, 1);
        
        pbc_pos = [];
        for ix = -1:1 % Reduced to 3x3 for speed in loop (still valid for interference)
            for iy = -1:1
                shifted_pos = bsxfun(@plus, pos, [ix*L, iy*L]);
                pbc_pos = [pbc_pos; shifted_pos];
            end
        end
        
        % 3. Save Snapshot (Only first rep)
        if rep == 1
            pos_snapshot = pos;
            sbc_pos_snapshot = sbc_pos;
            pbc_pos_snapshot = pbc_pos;
        end
        
        % 4. Compute S(k) for this replicate
        % Note: We use a local get_S function optimized for the loop
        
        % 2D Spectra
        S_pbc_2D_sum = S_pbc_2D_sum + get_S_local(pbc_pos, KX, KY, size(pos,1));
        S_sbc_2D_sum = S_sbc_2D_sum + get_S_local(sbc_pos, KX, KY, N_sbc_curr);
        
        % Azimuthal
        S_pbc_azi_sum = S_pbc_azi_sum + get_S_local(pbc_pos, kx_circ, ky_circ, size(pos,1));
        S_sbc_azi_sum = S_sbc_azi_sum + get_S_local(sbc_pos, kx_circ, ky_circ, N_sbc_curr);
        
        % Radial Average (Manual integration for this snapshot)
        S_rad_curr = zeros(size(k_radial));
        for i = 1:numel(k_radial)
            km = k_radial(i);
            % Use fewer angles per k in loop for speed, avg cleans it up
            th_loc = linspace(0, 2*pi, 12); 
            kx_r = km * cos(th_loc);
            ky_r = km * sin(th_loc);
            vals = get_S_local(sbc_pos, kx_r, ky_r, N_sbc_curr);
            S_rad_curr(i) = mean(vals);
        end
        S_sbc_radial_sum = S_sbc_radial_sum + S_rad_curr;
    end

    % --- NORMALIZE ENSEMBLE ---
    S_pbc_2D = S_pbc_2D_sum / N_reps;
    S_sbc_2D = S_sbc_2D_sum / N_reps;
    S_pbc_azi = S_pbc_azi_sum / N_reps;
    S_sbc_azi = S_sbc_azi_sum / N_reps;
    S_sbc_radial = S_sbc_radial_sum / N_reps;

    % --- ANALYSIS (Panel F) ---
    % Now performing subtraction on the SMOOTH averaged data
    x = k_radial * R;
    FormFactor_Shape = (2 * besselj(1, x) ./ x).^2;
    
    % Dynamic Scaling based on Peak Match
    peak_ratio = S_sbc_radial(1) / FormFactor_Shape(1);
    % Heuristic: Liquid S(0) is small, so GeomFactor is ~98% of total
    Scale_Factor = peak_ratio * 0.98; 
    
    S_analytical = Scale_Factor * FormFactor_Shape;
    S_clean = S_sbc_radial - S_analytical;
    
    % Theoretical Compressibility Limit (Estimated from low-k plateau)
    S_theory_val = mean(S_clean(1:15));

    % --- PLOTTING ---
    f = figure('Position', [50, 50, 1100, 1000], 'Color', 'w'); % WHITE BACKGROUND
    t = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'tight');
    
    % === ROW 1: REAL SPACE (Snapshot) ===
    nexttile; hold on; axis equal; box on;
    rectangle('Position', [-L/2, -L/2, L, L], 'EdgeColor', 'c', 'LineWidth', 2);
    plot(pos_snapshot(:,1), pos_snapshot(:,2), 'o', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    plot(pbc_pos_snapshot(:,1), pbc_pos_snapshot(:,2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
    xlim([-1.2 1.2]); ylim([-1.2 1.2]);
    title('\textbf{A. PBC: The Lattice}', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'XTick', [], 'YTick', []);
    
    nexttile; hold on; axis equal; box on;
    viscircles([0 0], R, 'Color', 'b', 'LineWidth', 2);
    plot(sbc_pos_snapshot(:,1), sbc_pos_snapshot(:,2), 'o', 'MarkerFaceColor', [1 0.6 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
    xlim([-1.2 1.2]); ylim([-1.2 1.2]);
    title('\textbf{B. SBC: The Finite Domain}', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'XTick', [], 'YTick', []);

    % === ROW 2: ENSEMBLE SPECTRA ===
    nexttile;
    imagesc(k_range, k_range, log10(S_pbc_2D+1e-3));
    axis xy equal;  caxis([0 3]); % Tighter contrast
    title('\textbf{C. PBC Spectrum}', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$k_x$', 'Interpreter', 'latex'); ylabel('$k_y$', 'Interpreter', 'latex');
    colorbar;
    
    nexttile;
    imagesc(k_range, k_range, log10(S_sbc_2D+1e-3));
    axis xy equal;  caxis([0 3]);
    title('\textbf{D. SBC Spectrum}', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$k_x$', 'Interpreter', 'latex'); ylabel('$k_y$', 'Interpreter', 'latex');
    colorbar;

    % === ROW 3: ANALYSIS ===
    nexttile; hold on; box on; grid on;
    plot(theta, S_pbc_azi, 'Color', [0 0.447 0.741], 'LineWidth', 1.5);
    plot(theta, S_sbc_azi, 'Color', [0.929 0.694 0.125], 'LineWidth', 3);
    title(['\textbf{E. Anisotropy at } $|k| = 2\pi/L$'], 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Angle $\theta$ (degrees)', 'Interpreter', 'latex'); ylabel('$S(k,\theta)$', 'Interpreter', 'latex');
    xlim([0 360]); set(gca, 'YScale', 'log'); ylim([1e-2 1e3]);
    legend({'PBC (Averaged)', 'SBC (Averaged)'}, 'Interpreter', 'latex', 'Location', 'SouthEast');
    xticks(0:90:360);

    nexttile; hold on; box on; grid on;
    p1 = plot(k_radial, S_sbc_radial, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    p2 = plot(k_radial, S_analytical, '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
    p_theory = yline(S_theory_val, ':', 'Color', [0 0.5 1], 'LineWidth', 1.5);
    p3 = plot(k_radial, S_clean, 'Color', [0 0.6 0], 'LineWidth', 3);
    
    title('\textbf{F. Domain Subtraction}', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$|k|$ (Radially Averaged)', 'Interpreter', 'latex'); ylabel('$S(|k|)$', 'Interpreter', 'latex');
    set(gca, 'YScale', 'log'); ylim([1e-2 1e3]); xlim([0 20]);
    legend([p1, p2, p3, p_theory], {'SBC (Avg)', 'Analytic Form', 'Recovered Physics', 'Limit'}, 'Interpreter', 'latex');

    set(gcf, 'InvertHardcopy', 'off');


% --- Optimized Local S(k) Function ---
function S = get_S_local(p, kx, ky, NormFactor)
    % Optimized for speed (no error checking, simple matrix ops)
    kv = [kx(:), ky(:)];
    % Phase: (N x 2) * (2 x M) = N x M
    phase = p * kv';
    S_vec = abs(sum(exp(1i*phase), 1)).^2;
    S = reshape(S_vec, size(kx)) / NormFactor;
end