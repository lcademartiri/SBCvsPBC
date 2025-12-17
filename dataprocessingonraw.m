%% CLEANUP

clear all
close all
clc

%% DATA SELECTION

ic=26;
filenameseries='SBCvsPBC_%d.mat';
filenameseriesdata='SBCvsPBC_%d_DATA.mat';

%% FLAGS
plottingenabled=true;

%% FOLDERS

if exist('D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
end
if exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
end
    
toolbox_folder = '..\ARBD_toolbox';
utilities_folder =  fullfile('..', '..', 'Utilities');
utilities_folder = genpath(utilities_folder);
addpath(utilities_folder)
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)
cmaplibrary=loadColormaps();

%% LOAD FILE

filename=sprintf(filenameseries,ic);
load(filename,'EDGES')
p=POS;
clear POS PDF PDFT SSF
[N, dim, T_steps] = size(p);
cleaner=true(ic,1);
cleaner(ic,1)=false;
EDGES(cleaner,:)={[]};

%% POST LOAD

% --- 1. COM CORRECTION & UNWRAPPING ---
clear pt pu put dp sp spt spu sput
COMw = mean(p,1);   % wrapped COM
pt=p-COMw; % COM-subtracted wrapped coordinates
slices=(1:T_steps/10:T_steps)';
for i0=1:10
    sp{i0,1}=p(:,:,slices(i0):slices(i0)+T_steps/10-1);
    spt{i0,1}=pt(:,:,slices(i0):slices(i0)+T_steps/10-1);
end
if S.bc==2
    for irep=1:10
        spu{irep,1} = zeros(size(sp{irep,1})); % initialize unwrapped positions array
        spu{irep,1}(:,:,1) = sp{irep,1}(:,:,1); % starting positions
        for t = 2:T_steps/10
            dp = sp{irep,1}(:,:,t) - sp{irep,1}(:,:,t-1); % displacements
            dp(dp >  S.br) = dp(dp >  S.br) - 2*S.br; % MIC #1
            dp(dp < -S.br) = dp(dp < -S.br) + 2*S.br; % MIC #1
            spu{irep,1}(:,:,t) = spu{irep,1}(:,:,t-1) + dp;
        end
        COMu = (mean(spu{irep,1},1));
        sput{irep,1}=spu{irep,1}-COMu;
    end
elseif S.bc==3
    for irep=1:10
        spu{irep,1} = zeros(size(sp{irep,1})); % initialize unwrapped positions array
        spu{irep,1}(:,:,1) = sp{irep,1}(:,:,1); % starting positions
        A_mat = S.fcc.A; 
        invA_mat = S.fcc.invA;
        for t = 2:T_steps/10
            dp = sp{irep,1}(:,:,t) - sp{irep,1}(:,:,t-1); % A. Calculate raw Cartesian displacement
            dp_frac = dp * invA_mat; % B. Convert displacement to Fractional coordinates
            dp_frac = dp_frac - round(dp_frac); % C. Apply Minimum Image Convention in Fractional Space
            dp_mic = dp_frac * A_mat; % D. Convert back to Cartesian displacement
            spu{irep,1}(:,:,t) = spu{irep,1}(:,:,t-1) + dp_mic; % E. Accumulate
        end
        COMu = mean(spu{irep,1}, 1);
        sput{irep,1} = spu{irep,1} - COMu;
    end
end

%% WINDOWS

% K-WINDOW
k_fundamental = 2*pi/(2*S.br);
k_max=pi/S.rp;
k_mags=(k_fundamental:k_fundamental:k_max)';
k_mags=sort([k_mags;k_fundamental.*[0.5;sqrt(2);sqrt(3);pi]]);
nK = length(k_mags);

% AZIMUTH-WINDOW
az_deg = (0:10:170)';
az_rad = deg2rad(az_deg); 
nAz = length(az_rad);

% ELEVATION-WINDOW
el_deg = (-90:10:90)';
el_rad = deg2rad(el_deg); 
nEl = length(el_rad);

% DYNAMICS
max_lag = T_steps/100; 

% OUTPUT ARRAY INITIALIZATION
% output arrays containing the means across replicates
S_std = zeros(nAz, nEl, nK);  
S_mask = zeros(nAz, nEl, nK); 
Deff_std = zeros(nAz, nEl, nK);
Deff_mask = zeros(nAz, nEl, nK);

% Scalar Maps (Keep Replicates for Error Bars)
    % Dimensions: Theta x Phi x K x Replicates
SSTD      = zeros(nAz, nEl, nK, 10);
DEFFSTD   = zeros(nAz, nEl, nK, 10);
GAMMASTD  = zeros(nAz, nEl, nK, 10);

% Full Time Correlation Function (Average Only)
% Dimensions: Theta x Phi x K x Time
% Using 'single' to save 50% RAM
F_AVG_STD = zeros(nAz, nEl, nK, max_lag+1);

if S.bc == 2 || S.bc==3
    SMASK     = zeros(nAz, nEl, nK, 10);
    DEFFMASK  = zeros(nAz, nEl, nK, 10);
    GAMMAMASK = zeros(nAz, nEl, nK, 10);
    F_AVG_MASK= zeros(nAz, nEl, nK, max_lag+1);
end

%% MASK

% Mask is defined on the COM-corrected 'p'
if S.bc == 2 % PBC
    for irep=1:10
	    dist_sq = sum(spt{irep,1}.^2, 2); 
	    mask{irep,1} = squeeze(dist_sq < S.br^2); 
	    N_eff{irep,1} = mean(sum(mask{irep,1}, 1));
    end
elseif S.bc==3
    for irep=1:10
        box_heights = abs(dot(S.fcc.A, S.fcc.normals, 2));
        R_mask = min(box_heights) / 2; % The inscribed radius is half the smallest height
        dist_sq = sum(spt{irep,1}.^2, 2); % (N x 1 x T_seg)
        mask{irep,1} = squeeze(dist_sq < R_mask^2);
        N_eff_t = sum(mask{irep,1}, 1);
        N_eff{irep,1} = mean(N_eff_t);
        if N_eff{irep,1} < 1, N_eff{irep,1} = 1; end
    end
elseif S.bc==1 % SBC
    for irep=1:10
        mask{irep,1} = true(S.N, 1, T_steps/10);
        N_eff{irep,1} = S.N;
    end
end

%% SUPERLOOP
for irep=1:10
    % initialize replicate storage
    seg_S = zeros(nAz, nEl, nK);
    seg_D = zeros(nAz, nEl, nK);
    seg_G = zeros(nAz, nEl, nK);
    if S.bc==2 || S.bc==3
        seg_Sm = zeros(nAz, nEl, nK);
        seg_Dm = zeros(nAz, nEl, nK);
        seg_Gm = zeros(nAz, nEl, nK);
    end

    for az_idx = 1:nAz
        az = az_rad(az_idx);
        sin_th = sin(az); 
        cos_th = cos(az);
        
        for el_idx = 1:nEl
            el = el_rad(el_idx);
            nx = sin_th * cos(el); % Unit Direction Vector n_hat
            ny = sin_th * sin(el); % Unit Direction Vector n_hat
            nz = cos_th; % Unit Direction Vector n_hat
            if S.bc==3
                [nfcc,~]=FCCrotate([nx,ny,nz],[1,1,1]./norm([1,1,1]));
                nx=nfcc(1);
                ny=nfcc(2);
                nz=nfcc(3);
            end
            
            for k_idx = 1:nK
                k_val = k_mags(k_idx);
                qx = k_val * nx; % 3D q-vector
                qy = k_val * ny; % 3D q-vector
                qz = k_val * nz; % 3D q-vector
            
                % -- STATIC S(k) --
                % Use wrapped coords (pt)
                px = squeeze(spt{irep,1}(:,1,:)); 
                py = squeeze(spt{irep,1}(:,2,:)); 
                pz = squeeze(spt{irep,1}(:,3,:));
                
                % Dot product q*r
                phase = -(qx.*px + qy.*py + qz.*pz);
                E = exp(1i * phase);
                
                % Standard
                rho_stat = sum(E, 1);
                seg_S(az_idx, el_idx, k_idx) = mean(abs(rho_stat).^2) / N;
                
                if S.bc == 2 || S.bc==3
                    rho_mask = sum(mask{irep,1} .* E, 1);
                    seg_Sm(az_idx, el_idx, k_idx) = mean(abs(rho_mask).^2) / N_eff{irep,1};
                end

                % -- DYNAMICS Deff(k) --
                % Use unwrapped coords (put)
                if S.bc~=1
                    pux = squeeze(sput{irep,1}(:,1,:));
                    puy = squeeze(sput{irep,1}(:,2,:));
                    puz = squeeze(sput{irep,1}(:,3,:));
                else
                    pux = squeeze(spt{irep,1}(:,1,:));
                    puy = squeeze(spt{irep,1}(:,2,:));
                    puz = squeeze(spt{irep,1}(:,3,:));
                end
                
                phase_dyn = -(qx.*pux + qy.*puy + qz.*puz);
                E_dyn = exp(1i * phase_dyn);
                rho_dyn = sum(E_dyn, 1);
                [F_curve, ~] = compute_acf(rho_dyn, max_lag);
                [D_val, G_val] = fit_dynamics(F_curve, S.timestep*10, k_val);
                
                seg_D(az_idx, el_idx, k_idx) = D_val;
                seg_G(az_idx, el_idx, k_idx) = G_val;
                F_AVG_STD(az_idx, el_idx, k_idx, :) = F_AVG_STD(az_idx, el_idx, k_idx, :) + reshape((F_curve), 1, 1, 1, []);
                
                %seg_D(az_idx, el_idx, k_idx) = get_deff(rho_dyn, max_lag, S.timestep, k_val);
                if S.bc == 2 || S.bc==3
                    rho_dyn_m = sum(mask{irep,1} .* E_dyn, 1);
                    [F_curve_m, ~] = compute_acf(rho_dyn_m, max_lag);
                    [D_m, G_m] = fit_dynamics(F_curve_m, S.timestep*10, k_val);
                    seg_Dm(az_idx, el_idx, k_idx) = D_m;
                    seg_Gm(az_idx, el_idx, k_idx) = G_m;
                    F_AVG_MASK(az_idx, el_idx, k_idx, :) = F_AVG_MASK(az_idx, el_idx, k_idx, :) + reshape(single(F_curve_m), 1, 1, 1, []);
                    %seg_Dm(az_idx, el_idx, k_idx) = get_deff(rho_dyn_m, max_lag, S.timestep, k_val);
                end

                disp([irep,az_idx, el_idx, k_idx])
            end
        end
    end
    % Store segment
    FSTD(:,:,:,:,irep)=F_AVG_STD;
    SSTD(:,:,:,irep) = seg_S;
    DEFFSTD(:,:,:,irep) = seg_D;
    GAMMASTD(:,:,:,irep) = seg_G;
    if S.bc == 2 || S.bc==3
        FMASK(:,:,:,:,irep)=F_AVG_MASK;
        SMASK(:,:,:,irep) = seg_Sm;
        DEFFMASK(:,:,:,irep) = seg_Dm;
        GAMMAMASK(:,:,:,irep) = seg_Gm;
    end
end

% cleanup
clear sp* phase* pt pu* px py pz seg_* F_AVG_*

% test=mean(DEFFSTD,4);
% test0el=squeeze(test(:,10,:));
% az_full = [az_rad; az_rad + pi;2*pi];
% [AZ_GRID, K_GRID] = meshgrid(rad2deg(az_full),k_mags);
% test0el_full = [test0el; test0el;test0el(1,:)];
% test0el_full=test0el_full';
% OriginData=[AZ_GRID(:),K_GRID(:),test0el_full(:)];

%% --- SSF ANALYSIS --------------------------------------------------

load(filename,'SSF')

% --- 0. SETUP & BASICS ---
fprintf('Starting SSF Analysis...\n');
if ~isfield(SSF, 'nsnap') || SSF.nsnap < 1
    error('SSF Error: No snapshots found in structure.');
end

N_part = S.N; 
norm_fac = 1 / (SSF.nsnap * N_part);

% Define High-Contrast Colors for Black Background
c_raw   = [0.9 0.9 0.9];       % White/Gray (Raw Data)
c_stat  = [0.2 1.0 0.2];       % Neon Green (Static/Bragg/FormFactor)
c_fluct = [1.0 0.2 1.0];       % Magenta (Pure Fluctuations)
c_ref   = [0.4 0.8 1.0 0.5];   % Cyan transparent (Reference)

% --- 1. PROCESS CARTESIAN PROBES (Common to all) ---
% A. "Random Directions" (Off-Lattice / Continuum)
% In SBC: Used for subtraction method. In PBC: Shows spectral leakage.
k_rand = SSF.kmag_sampling;
S_rand_tot  = SSF.sampling.sum_abs2 * norm_fac;
S_rand_stat = (abs(SSF.sampling.sum_rho / SSF.nsnap).^2) / N_part;
S_rand_fluc = S_rand_tot - S_rand_stat;

% B. "Reciprocal Lattice" (Cubic/FCC Grid)
% In SBC: Meaningless. In PBC: The correct Eigenmodes.
k_lat = SSF.kmag_lattice;
S_lat_tot  = SSF.lattice.sum_abs2 * norm_fac;
S_lat_stat = (abs(SSF.lattice.sum_rho / SSF.nsnap).^2) / N_part;
S_lat_fluc = S_lat_tot - S_lat_stat;

% --- 2. BRANCH: SPHERICAL BOUNDARY CONDITIONS (SBC) ---
if S.bc == 1
    fprintf('Detected SBC. Computing Bessel Eigenmodes...\n');
    
    % --- BESSEL NORMALIZATION (Neumann + Volume + Degeneracy) ---
    R = S.br;
    V_sphere = (4/3)*pi*R^3;
    modes = SSF.sbc_modes.data; % [k, l, n]
    k_bes = modes(:,1);
    l_vals = modes(:,2);
    
    % Calculate Geometric Weights
    W_ln = zeros(size(k_bes));
    for i = 1:numel(k_bes)
        l = l_vals(i);
        x = k_bes(i) * R;
        if x < 1e-4
            W_ln(i) = 1.0;
        else
            val = sqrt(pi/(2*x)) * besselj(l + 0.5, x);
            % Neumann Factor: (1 - l(l+1)/x^2)
            geo_fac = (1 - (l*(l+1))/(x^2));
            W_ln(i) = ((R^3/2) * val^2 * geo_fac) / V_sphere;
        end
    end
    
    % Apply Normalizations
    raw_bes = SSF.sbc_modes.sum_power / SSF.nsnap;
    S_bes = (raw_bes ./ W_ln) / N_part;  % Geometric
    S_bes = S_bes ./ (2*l_vals + 1);     % Degeneracy Correction!

    % Separate Bulk (l>0) from Wall Layering (l=0)
    mask_bulk = l_vals > 0 & k_bes > 1e-4;
    mask_wall = l_vals == 0 & k_bes > 1e-4;

    % --- PLOT SBC ---
    figure('Color','k', 'Position', [50 50 1400 600]);
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % PANEL 1: CARTESIAN SUBTRACTION (Method A)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Raw
    plot(k_rand, S_rand_tot, 'o', 'Color', c_raw, 'MarkerSize', 2, ...
        'DisplayName', 'Raw (Includes Sphere Form Factor)');
    % Static Background
    plot(k_rand, S_rand_stat, '-', 'Color', c_stat, 'LineWidth', 1.5, ...
        'DisplayName', 'Static Profile (Form Factor)');
    % Result
    plot(k_rand, S_rand_fluc, '.', 'Color', c_fluct, 'MarkerSize', 10, ...
        'DisplayName', 'Cleaned S(k) (Subtraction)');
    
    title('Method A: Cartesian Basis (Subtraction)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 10]); xlim([0 max(k_rand)]); set(gca, 'YScale', 'log');

    % PANEL 2: BESSEL EIGENMODES (Method B)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Bulk Modes (Rainbow)
    scatter(k_bes(mask_bulk), S_bes(mask_bulk), 30, l_vals(mask_bulk), 'filled', ...
        'DisplayName', 'Bulk Modes (l > 0)');
    colormap(gca, 'turbo'); 
    cb = colorbar; cb.Label.String = 'Angular Momentum (l)'; 
    cb.Color = 'w'; cb.Label.Color = 'w';
    
    % Wall Modes (White Outliers)
    plot(k_bes(mask_wall), S_bes(mask_wall), 'wo', 'MarkerSize', 6, 'LineWidth', 1.5, ...
        'DisplayName', 'Wall Layering (l = 0)');
    
    % Overlay Method A Result for comparison
    plot(k_rand, smooth(S_rand_fluc, 5), '-', 'Color', c_ref, 'LineWidth', 2, ...
        'DisplayName', 'Ref: Cartesian Cleaned');
    
    title('Method B: Bessel Basis (Natural Separation)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    yline(1, 'w--');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w', 'Location', 'Southeast');
    ylim([0.01 10]); xlim([0 max(k_bes)]); set(gca, 'YScale', 'log');
    
    sgtitle('SBC Analysis: Recovering Bulk S(k) in Confinement', 'Color', 'w', 'FontSize', 14);


% --- 3. BRANCH: PERIODIC BOUNDARY CONDITIONS (PBC/FCC) ---
else
    fprintf('Detected PBC/FCC. Comparing Lattice vs Random Probes...\n');
    
    figure('Color','k', 'Position', [50 50 1400 600]);
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % PANEL 1: OFF-LATTICE (Showing Leakage)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Plot Random Vectors (Off-lattice)
    % These will look messy because they hit the sinc-function zeros of the Box
    plot(k_rand, S_rand_tot, 'o', 'Color', c_raw, 'MarkerSize', 3, ...
        'DisplayName', 'Off-Lattice (Random Directions)');
    
    % Overlay Lattice Points (To show where the peaks SHOULD be)
    plot(k_lat, S_lat_tot, 'r.', 'MarkerSize', 5, 'DisplayName', 'Lattice Points (Ref)');
    
    title('Probe A: Off-Lattice Vectors (Spectral Leakage)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 100]); xlim([0 max(k_rand)]); set(gca, 'YScale', 'log', 'XScale', 'log');
    
    % PANEL 2: RECIPROCAL LATTICE (Correct Eigenmodes)
    nexttile; hold on; box on;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.2);
    
    % Use Stems for discrete lattice modes
    % Raw Total
    stem(k_lat, S_lat_tot, 'Color', [0.5 0.5 0.5], 'Marker', 'none', 'LineWidth', 1, ...
        'DisplayName', 'Total (Raw)');
    
    % Static (Bragg Peaks - Crystal Order)
    if max(S_lat_stat) > 0.1
        stem(k_lat, S_lat_stat, 'Color', c_stat, 'Marker', 's', 'MarkerSize', 4, ...
            'LineStyle', 'none', 'DisplayName', 'Static (Bragg Order)');
    end
    
    % Fluctuations (Liquid Physics)
    stem(k_lat, S_lat_fluc, 'Color', c_fluct, 'Marker', 'o', 'MarkerSize', 4, ...
        'LineWidth', 1.5, 'LineStyle', 'none', 'DisplayName', 'Pure Fluctuations');
        
    title('Probe B: Reciprocal Lattice (Eigenmodes)', 'Color', 'w');
    xlabel('k', 'Color','w'); ylabel('S(k)', 'Color','w');
    yline(1, 'w--');
    grid on; legend('Color','k','TextColor','w','EdgeColor','w');
    ylim([0.01 100]); xlim([0 max(k_lat)]); set(gca, 'YScale', 'log', 'XScale', 'log');
    
    sgtitle('PBC Analysis: Importance of Lattice Commensurability', 'Color', 'w', 'FontSize', 14);
end

%% --- PDF ANALYSIS --------------------------------------------------

load(filename,'PDFT','PDF')
% --- general purpose data
N=S.N;
globV=S.bv;
dens=N/globV;
kbT=1.38e-23*298.15;
nosnap=size(PDFT,1); % number of snapshots in PDF
% ---
% --- sequences
for i0=1:nosnap
    AZ(:,i0)=PDFT{i0,1};
    EL(:,i0)=PDFT{i0,2};
    RHO(:,i0)=PDFT{i0,3};
    AZS(:,i0)=PDFT{i0,4};
    ELS(:,i0)=PDFT{i0,5};
end
% clear PDFT
% ---
% --- means
MEAN_AZ=mean(AZ,2);
MEAN_EL=mean(EL,2);
MEAN_RHO=mean(RHO,2);
MEAN_AZS=mean(AZS,2);
MEAN_ELS=mean(ELS,2);
% ---
% --- moving averages
MMEAN_AZ=movmean(AZ,100,2);
MMEAN_AZ=MMEAN_AZ./sum(MMEAN_AZ);
MMEAN_EL=movmean(EL,100,2);
MMEAN_EL=MMEAN_EL./sum(MMEAN_EL);
MMEAN_RHO=movmean(RHO,100,2);
MMEAN_RHO=MMEAN_RHO./sum(MMEAN_RHO);
MMEAN_AZS=movmean(AZS,100,2);
MMEAN_AZS=MMEAN_AZS./sum(MMEAN_AZS);
MMEAN_ELS=movmean(ELS,100,2);
MMEAN_ELS=MMEAN_ELS./sum(MMEAN_ELS);
% ---
% --- variance of moving averages
STD_MMEAN_AZ=std(MMEAN_AZ)';
STD_MMEAN_EL=std(MMEAN_EL)';
STD_MMEAN_RHO=std(MMEAN_RHO)';
STD_MMEAN_AZS=std(MMEAN_AZS)';
STD_MMEAN_ELS=std(MMEAN_ELS)';
% --------------------------------------------------------------

% --- VERIFY AND IDENTIFY THERMALIZATION -----------------------
fprintf('condition: %d - boundary condition: %d -  phi: %.3f - verify and identify thermalization\n', ic, S.bc, S.phi);
% --- find a good lognormal maximum using a coarse binning
edges_azs=linspace(log10(min(STD_MMEAN_AZS)),log10(max(STD_MMEAN_AZS)),100)';
[counts_azs,edges_azs]=histcounts(log10(STD_MMEAN_AZS),edges_azs);
hist_azs=[edges_azs(1:end-1,:),counts_azs'./sum(counts_azs)];
peakvalue=hist_azs(hist_azs(:,2)==max(hist_azs(:,2)),1);
peakvalue=peakvalue(1);
% ---
% --- fine log10 histogram for windowed lognormal fit
edges_azs=linspace(log10(min(STD_MMEAN_AZS)),log10(max(STD_MMEAN_AZS)),1000)';
[counts_azs,edges_azs]=histcounts(log10(STD_MMEAN_AZS),edges_azs);
hist_azs=[edges_azs(1:end-1,:),counts_azs'./sum(counts_azs)];
[peakid,edges_azs]=histcounts((peakvalue),edges_azs);
% ---
% --- find peak bin
peakid=[edges_azs(1:end-1,:),peakid'];
peakid=find(peakid(:,2));
peakA=hist_azs(peakid,2);
peakA=peakA(1);
% ---
% --- define gaussian fit options
ft = fittype(...
'A * ( eta/(1+((x-mu)/gamma)^2) + (1-eta)*exp(-(x-mu)^2/(2*sigma^2)) )', ...
'independent','x', ...
'coefficients',{'A','mu','sigma','gamma','eta'});

opts = fitoptions(ft);

opts.Lower = [ ...
    0.5*peakA, ...                 % A
    1.5*peakvalue, ...             % mu
    0.001*abs(peakvalue), ...       % sigma
    0.001*abs(peakvalue), ...       % gamma
    0 ];                           % eta

opts.Upper = [ ...
    1.5*peakA, ...                 % A
    0.8*peakvalue, ...             % mu
    1*abs(peakvalue), ...          % sigma
    1*abs(peakvalue), ...          % gamma
    1 ];                           % eta

opts.StartPoint = [ ...
    peakA, ...                     % A
    peakvalue, ...                 % mu
    0.2*abs(peakvalue), ...        % sigma
    0.2*abs(peakvalue), ...        % gamma
    0.5 ];                         % eta  (equal L/G mix initially)
% ---
% --- windowed fit
xwins=[10:10:1000]';
xall=hist_azs(:,1);
yall=hist_azs(:,2);
SStot = sum((yall - mean(yall)).^2);
R2best=-1e6;
bestfit=[];
for iwin=1:numel(xwins)
    xwin=xwins(iwin);
    bottom=peakid-xwin;
    if bottom<1
        bottom=1;
    end
    top=peakid+xwin;
    if top>size(hist_azs,1)
        top=size(hist_azs,1);
    end
    x=hist_azs(bottom:top,1);
    y=hist_azs(bottom:top,2);
    [cf, gof] = fit(x, y, ft, opts);
    yfit_all = cf.A * ( cf.eta ./ (1 + ((xall - cf.mu)/cf.gamma).^2) + ...
                    (1-cf.eta) .* exp(-(xall - cf.mu).^2/(2*cf.sigma^2)) );
    SSres = sum((yall - yfit_all).^2);        
    R2(iwin) = 1 - SSres/SStot;
    if R2(iwin)>R2best
        bestfit.cf=cf;
        bestfit.gof=gof;
        R2best=R2(iwin);
    end
end
if plottingenabled==1
    figure
    plot(xall,yall)
    hold
    yfit_all = bestfit.cf.A * ( bestfit.cf.eta ./ (1 + ((xall - bestfit.cf.mu)/bestfit.cf.gamma).^2) + ...
                            (1-bestfit.cf.eta) .* exp(-(xall - bestfit.cf.mu).^2/(2*bestfit.cf.sigma^2)) );
    plot(xall,yfit_all)
    hold
end
confints=confint(bestfit.cf);
thermalrange=[10^(confints(1,2)),10^(confints(2,2))];
% ---
% --- find thermalized data
for idtherm=1:nosnap
    if STD_MMEAN_AZS(idtherm)>thermalrange(1) && STD_MMEAN_AZS(idtherm)<thermalrange(2)
        break
    end
end
% ---
% --- cleanup
clear peak* opts R2* bottom top x* y* edges* gof cf de iwin i0 ft count_azs SSres SStot
% ---

% --- THERMALIZED DATA ----------------------------------------------
fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculate aggregated, fully thermalized data\n', ic, S.bc, S.phi);
% --- thermalized means
MEAN_AZ=mean(AZ(:,idtherm:end),2);
MEAN_AZ=MEAN_AZ./sum(MEAN_AZ);
MEAN_EL=mean(EL(:,idtherm:end),2);
MEAN_EL=MEAN_EL./sum(MEAN_EL);
MEAN_RHO=mean(RHO(:,idtherm:end),2);
MEAN_RHO=MEAN_RHO./sum(MEAN_RHO);
MEAN_AZS=mean(AZS(:,idtherm:end),2);
MEAN_AZS=MEAN_AZS./sum(MEAN_AZS);
MEAN_ELS=mean(ELS(:,idtherm:end),2);
MEAN_ELS=MEAN_ELS./sum(MEAN_ELS);
% ---
% --- plotting
if plottingenabled==1
    figure
    plot(MEAN_AZS)
    hold
    plot(MEAN_AZ)
    figure
    plot(MEAN_ELS)
    hold
    plot(MEAN_EL)
end
% ---
% -------------------------------------------------------------------

% -------------------------------------------------------------------
       
% --- determining the ideal gas distance distribution by Montecarlo ---
fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculating pair distribution function\n', ic, S.bc, S.phi);
if S.bc==1
    filepdfdenom = sprintf('PDFdenom_SBC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
elseif S.bc==2
    filepdfdenom = sprintf('PDFdenom_PBCc_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
elseif S.bc==3
    filepdfdenom = sprintf('PDFdenom_PBCFCC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
end
load(filepdfdenom,'gdenominator');
% ---
% --- determine the g
fprintf('condition: %d - boundary condition: %d -  phi: %.3f - average pairs numerator\n', ic, S.bc, S.phi);
gnumerator=sum(RHO(:,idtherm:end),2)./(size(RHO,2)-idtherm); % average number of pairs in that bin
g=gnumerator./gdenominator;
g=[PDF.pdfedges{3}(1:end-1),g];
% ---
% --- plot g
if plottingenabled==1
    figure
    plot(g(:,1),g(:,2))
    xlim([2*S.rp 2*S.br+2*S.rc]);
    ylim([0.5 max(g(:,2))*1.2]);
    yline(1)
    xline(2*S.br)
    xline(2*(S.br-S.rp))
end
% ---
% cleanup
clear PDFT STD_MMEAN* MMEAN* AZ EL RHO AZS ELS

%% --- WAITING TIME ANALYSIS --------------------------------------------------
       
% --- per particle ---
colls=EDGES{ic,1};
clear EDGES
idxswap=colls(:,3)<colls(:,2);
colls(idxswap,[2 3])=colls(idxswap,[3 2]);
colls=sortrows(colls,[2 1],"ascend");
ds=zeros(size(colls,1)-1,1);
q=1;
[~,b]=ismember((1:S.N)',colls(:,2));
for ib=1:size(b,1)-1
    xmin=b(ib);
    xmax=b(ib+1)-1;
    s=colls(xmin:xmax,1);
    ds(q:q+xmax-xmin-1,1)=diff(s);
    q=xmax+1;
    disp(ib)
end
ds(ds==0,:)=[];
wtppart.edges=(1:max(ds))';
[wtppart.counts,wtppart.edges]=histcounts(ds,wtppart.edges);
wtppart.dist=[wtppart.edges(1:end-1,1),wtppart.counts',wtppart.counts'./sum(wtppart.counts)];
wtppart.dist(wtppart.dist(:,2)==0,:)=[];
if plottingenabled==1
    figure
    plot(wtppart.dist(:,1),wtppart.dist(:,3))
    xscale log
    yscale log
end
% ----------------------
% --- per pair ---
ds=zeros(size(colls,1)-1,1);
q=1;
ncolls=size(colls,1);
colls=sortrows(colls,[2,3,1]);
qds=1;
while true
    coll1=colls(q,2);
    coll2=colls(q,3);
    qs=q+1;
    while qs>0
        if qs>ncolls
            qs=qs-1;
            s=colls(q:qs,1);
            ds(qds:qds+numel(s)-2,1)=diff(s);
            break
        end
        if colls(qs,3)==coll2
            if colls(qs,2)==coll1
                qs=qs+1;                   
            else
                qs=-qs;
            end
        else
            qs=-qs;
        end
    end
    if qs~=ncolls
        qs=-qs-1;
        s=colls(q:qs,1);
        ds(qds:qds+numel(s)-2,1)=diff(s);
        qds=qds+numel(s)-1;
    end
    q=qs+1;
    if q>ncolls
        break
    end

end
ds(ds==0,:)=[];
wtppair.edges=(1:max(ds))';
[wtppair.counts,wtppair.edges]=histcounts(ds,wtppair.edges);
wtppair.dist=[wtppair.edges(1:end-1,1),wtppair.counts',wtppair.counts'./sum(wtppair.counts)];
wtppair.dist(wtppair.dist(:,2)==0,:)=[];
if plottingenabled==1
    figure
    plot(wtppair.dist(:,1),wtppair.dist(:,3))
    xscale log
    yscale log
end
% ----------------------

% --- CONSECUTIVE COLLISION ANALYSIS ----
ds(ds~=1,:)=0;
ds=logical(ds);
d = diff([0, ds(:)', 0]);
startIndex = find(d == 1);
endIndex   = find(d == -1);
runLengths = endIndex - startIndex;
cc.edges=(1:max(runLengths))';
[cc.counts,cc.edges]=histcounts(runLengths,cc.edges);
cc.dist=[cc.edges(1:end-1,1),cc.counts',cc.counts'./sum(cc.counts)];
cc.dist(cc.dist(:,2)==0,:)=[];
if plottingenabled==1
    figure
    plot(cc.dist(:,1),cc.dist(:,3))
    xscale log
    yscale log
end
% --------------------------------------------------------------

%% COMPILING DATA

fprintf('condition: %d - boundary condition: %d -  phi: %.3f - collating data\n', ic, S.bc, S.phi);
ISO(ic).g=g;
ISO(ic).pdfedges=PDF.pdfedges;
ISO(ic).thermalized_mAZ=[PDF.pdfedges{1}(1:end-1),MEAN_AZ];
ISO(ic).thermalized_mAZS=[PDF.pdfedges{1}(1:end-1),MEAN_AZS];
ISO(ic).thermalized_mEL=[PDF.pdfedges{2}(1:end-1),MEAN_EL];
ISO(ic).thermalized_mELS=[PDF.pdfedges{2}(1:end-1),MEAN_ELS];
ISO(ic).thermalized_mRHO=[PDF.pdfedges{3}(1:end-1),MEAN_RHO];
ISO(ic).thermalizationtime=idtherm*S.kt;
ISO(ic).thermalizationfit=bestfit;
ISO(ic).S=S;
ISO(ic).P=P;
ISO(ic).V=V;
ISO(ic).C=C;
ISO(ic).F=FSTD;
ISO(ic).SSF=SSTD;
ISO(ic).DEFF=DEFFSTD;
ISO(ic).GAMMA=GAMMASTD;
if S.bc==2 || S.bc==3
    ISO(ic).Fm=FMASK;
    ISO(ic).SSFm=SMASK;
    ISO(ic).DEFFm=DEFFMASK;
    ISO(ic).GAMMAm=GAMMAMASK;
end
ISO(ic).SSF=SSF;
ISO(ic).PDF=PDF;
if S.bc==1
    ISO(ic).Sk_rand=[k_rand,S_rand_tot',S_rand_stat',S_rand_fluc'];
    ISO(ic).Sk_bes=[k_bes,S_bes,mask_bulk,mask_wall];
elseif S.bc==2 || S.bc==3
    ISO(ic).Sk_rand=[k_rand,S_rand_tot'];
    ISO(ic).Sk_lat=[k_lat,S_lat_tot',S_lat_stat',S_lat_fluc'];
end
ISO(ic).WT.ppart=wtppart;
ISO(ic).WT.ppair=wtppair;
ISO(ic).CC=cc;

filenamedata=sprintf(filenameseriesdata,ic);
save([output_folder,'\',filenamedata],'ISO','-v7.3');

%% HELPER FUNCTIONS

function D = get_deff(rho_t, max_lag, dt, k_val)
    if any(isnan(rho_t)), D=NaN; return; end
    rho_c = rho_t - mean(rho_t);
    acf = xcorr(rho_c, max_lag, 'biased');
    acf = acf(max_lag+1:end);
    
    if abs(acf(1)) < 1e-9, D=NaN; return; end
    
    F_norm = abs(acf) / abs(acf(1));
    time = (0:max_lag) * dt;
    idx = F_norm < 0.9 & F_norm > 0.2;
    
    if sum(idx) < 5
        D = NaN; 
    else
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        D = -p(1) / k_val^2;
    end
end

function [F_norm, lags_vec] = compute_acf(rho_t, max_lag)
    if any(isnan(rho_t))
        F_norm = nan(max_lag+1, 1); lags_vec = 0:max_lag; return; 
    end
    rho_c = rho_t - mean(rho_t);
    acf = xcorr(rho_c, max_lag, 'biased');
    acf = acf(max_lag+1:end);
    if abs(acf(1)) > 1e-9
        F_norm = abs(acf) / abs(acf(1));
    else
        F_norm = zeros(size(acf));
    end
    lags_vec = 0:max_lag;
end

function [D, Gamma] = fit_dynamics(F_norm, dt, k_val)
    time = (0:(length(F_norm)-1)) * dt;
    idx = F_norm < 0.9 & F_norm > 0.2;
    if sum(idx) < 5
        D = NaN; Gamma = NaN;
    else
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        Gamma = -p(1);
        D = Gamma / k_val^2;
    end
end