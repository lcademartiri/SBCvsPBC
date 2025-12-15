%% CLEANUP

clear all
close all
clc

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

load('SBCvsPBC_27.mat')
p=POS;
clear POS PDF PDFT SSF
[N, dim, T_steps] = size(p);

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
                [D_val, G_val] = fit_dynamics(F_curve, S.timestep, k_val);
                
                seg_D(az_idx, el_idx, k_idx) = D_val;
                seg_G(az_idx, el_idx, k_idx) = G_val;
                F_AVG_STD(az_idx, el_idx, k_idx, :) = F_AVG_STD(az_idx, el_idx, k_idx, :) + reshape((F_curve), 1, 1, 1, []);
                
                %seg_D(az_idx, el_idx, k_idx) = get_deff(rho_dyn, max_lag, S.timestep, k_val);
                if S.bc == 2 || S.bc==3
                    rho_dyn_m = sum(mask{irep,1} .* E_dyn, 1);
                    [F_curve_m, ~] = compute_acf(rho_dyn_m, max_lag);
                    [D_m, G_m] = fit_dynamics(F_curve_m, S.timestep, k_val);
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
    SSTD(:,:,:,irep) = seg_S;
    DEFFSTD(:,:,:,irep) = seg_D;
    if S.bc == 2 || S.bc==3
        SMASK(:,:,:,irep) = seg_Sm;
        DEFFMASK(:,:,:,irep) = seg_Dm;
    end
end

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