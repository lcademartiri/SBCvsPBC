%% CLEANUP

clear all
close all
clc

%% DATA SELECTION

ic=27;
filenameseries='SBCvsPBC_fixedPOT_%d.mat';
filenameseriesdata='SBCvsPBC_fixedPOT_%d_DATA.mat';

%% FLAGS
plottingenabled=true;
dodeff=true;
dovanhove=true;
dopdf=true;
dowt=true;


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
load(filename)
p=POS;
clear POS PDF PDFT SSF
[~, dim, T_steps] = size(p);
clearinghouse=true(ic,1);
clearinghouse(ic,1)=false;
EDGES(clearinghouse,:)={[]};
clear clearinghouse


%%  COM CORRECTION & UNWRAPPING ---
clear pt pu put dp sp spt spu sput
[spt,sput]=comCorrectionSpliceUnwrap(p,10,S);

% p     : whole     non-COM-corrected   wrapped     trajectories
% sp    : sliced    non-COM-corrected   wrapped     trajectories
% spt   : sliced    COM-corrected       wrapped     trajectories
% spu   : sliced    non-COM-corrected   unwrapped   trajectories
% sput  : sliced    COM-corrected       unwrapped   trajectories

%% MASK

[mask,N_eff]=sphericalMask(S,spt);


%% WAVEVECTORS

if dodeff
    % K-WINDOW
    k_fundamental = 2*pi/(2*S.br);
    k_max=pi/S.rp;
    k_mags=(k_fundamental:k_fundamental:k_max)';
    k_mags=sort([k_mags;k_fundamental.*[sqrt(2);sqrt(3);pi]]);
    nK = length(k_mags);
    
    % Fibonacci vectors
    az=fibonacci_sphere(600);
    [az,el,~]=cart2sph(az(:,1),az(:,2),az(:,3));
    azel=[az,el];
    azel(azel(:,1)<0,:)=[];

    % add equator vectors
    equator=linspace(0,pi,181)';
    equator(end,:)=[];
    equator(1,2)=0;
    azel_rad=[azel;equator];
    azel_deg=rad2deg(azel_rad);
    nAzel=size(azel_rad,1);
    clear az el equator azel
    
    % DYNAMICS
    max_lag = T_steps/100; 
end

%% STATICS

if S.bc==1
    [SSTD,~]=statics(S,spt,azel_rad,k_mags,0,0);
else
    [SSTD,SMASK]=statics(S,spt,azel_rad,k_mags,mask,N_eff);
end

    % test=mean(ISO(25).DEFF,3);
    % idxazel=azel_rad(:,2)==0;
    % tempaz=azel_rad(idxazel,1);
    % tempaz=sort([tempaz-pi;tempaz;pi]);
    % [AZ_GRID, K_GRID] = meshgrid(rad2deg(tempaz),k_mags);
    % test0el=test(idxazel,:);
    % test0el_full = [test0el; test0el;test0el(1,:)];
    % test0el_full=test0el_full';
    % OriginData=[AZ_GRID(:),K_GRID(:),test0el_full(:)];
%% DYNAMICS

if S.bc==1
    [FSTD,GAMMASTD,DEFFSTD,~,~,~]=dynamics(S,sput,azel_rad,k_mags,max_lag,0,0);
else
    [FSTD,GAMMASTD,DEFFSTD,FMASK,GAMMAMASK,DEFFMASK]=dynamics(S,sput,azel_rad,k_mags,max_lag,mask,N_eff);
end

%% --- VAN HOVE CORRELATION FUNCTION --------------------------------------------

if dovanhove

    % --- NON_GAUSSIAN PARAMETER
    
    % Setup
    % Assuming 'data' is your 1x10 cell array containing Nx3xT matrices
    % 'dt' is your timestep (e.g., 200 * tau)
    
    n_cells = length(sput);
    [alphaN, alphadim, n_steps] = size(sput{1});
    
    % Define maximum lag time (e.g., 10% of trajectory or until statistics get poor)
    % If n_steps is huge (10^6), restrict this to the relevant timescale (e.g., 5000 frames)
    max_tau = 100; 
    
    % Initialize accumulators for moments
    alphasum_r2 = zeros(max_tau, 1);
    alphasum_r4 = zeros(max_tau, 1);
    N_samples = zeros(max_tau, 1); % To track number of samples per lag
    
    % Calculation Loop
    fprintf('Computing Alpha2...\n');
    
    for alphac = 1:n_cells
        traj = sput{alphac}; % Get trajectory from cell (must be UNWRAPPED)
        
        % Loop over lag times tau
        % (We vectorize the calculation over all start times t0 for speed)
        for tau = 1:max_tau
            
            % 1. Compute displacements for this lag tau
            % delta_r will be size: [N, 3, (n_steps - tau)]
            % This effectively uses ALL valid starting times t0
            delta_r = traj(:, :, 1+tau:end) - traj(:, :, 1:end-tau);
            
            % 2. Compute Squared Displacement (r^2) per particle per window
            % Sum over dimension 2 (x, y, z)
            dr2 = sum(delta_r.^2, 2); 
            
            % Flatten to a vector (collapses N particles and all t0 windows)
            dr2_vec = dr2(:);
            
            % 3. Accumulate Moments
            alphasum_r2(tau) = alphasum_r2(tau) + sum(dr2_vec);
            alphasum_r4(tau) = alphasum_r4(tau) + sum(dr2_vec.^2);
            N_samples(tau) = N_samples(tau) + length(dr2_vec);
            disp([alphac,tau])
        end
    
        fprintf('Processed Cell %d/%d\n', alphac, n_cells);
    end
    
    % Final Compute
    % Average Moments FIRST
    M2 = alphasum_r2 ./ N_samples; % Mean Squared Displacement <r^2>
    M4 = alphasum_r4 ./ N_samples; % Mean Quartic Displacement <r^4>
    
    % Calculate Alpha2 (3D formula)
    % alpha2 = (3 * <r^4>) / (5 * <r^2>^2) - 1
    alpha2 = (3 .* M4) ./ (5 .* M2.^2) - 1;
    
    % Plotting
    time_axis = (1:max_tau); % Multiply by your actual dt if needed
    
    figure;
    plot(time_axis, alpha2, 'LineWidth', 2, 'Color', 'w'); % White line for black bg
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 2, 'FontSize', 14);
    xlabel('Time (steps)');
    ylabel('\alpha_2(t)');
    title('Non-Gaussian Parameter (Smoothed)');
    grid on;
    
    % Find the peak (Cage Breaking Time)
    [max_val, idx] = max(alpha2);
    t_star = time_axis(idx);
    xline(t_star, '--r', 'Cage Breaking t*', 'LineWidth', 2);
    
    % ---
    % Setup
    % Use the t_star you found (approx 22 steps based on your plot)
    target_tau = t_star; 
    
    % Accumulate all displacements dx, dy, dz for this specific lag
    all_dx = [];
    all_dy = [];
    all_dz = [];
    
    fprintf('Collecting displacements at tau = %d...\n', target_tau);
    
    for c = 1:n_cells
        traj = sput{c}; % Must be UNWRAPPED and COM-corrected
        
        % Calculate displacement for ALL windows of length target_tau
        % Vectorized for speed
        delta_r = traj(:, :, 1+target_tau:end) - traj(:, :, 1:end-target_tau);
        
        % Flatten and append
        all_dx = [all_dx; reshape(delta_r(:,1,:), [], 1)];
        all_dy = [all_dy; reshape(delta_r(:,2,:), [], 1)];
        all_dz = [all_dz; reshape(delta_r(:,3,:), [], 1)];
    end
    limit=max([max(abs(all_dx)),max(abs(all_dy)),max(abs(all_dz))]);
    
    % Plotting the 2D Heatmap (XY Plane)
    figure('Color', 'k');
    
    % 2D Histogram parameters
    bin_edges = linspace(-limit, limit, 500); % Adjust range based on your particle size (sigma)
    
    % Compute 2D histogram
    h = histogram2(all_dx, all_dy, bin_edges, bin_edges, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
    
    % Visual Styling for "The Kill Shot"
    colormap(hot); % 'hot' or 'parula' work well on black
    axis equal;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('\Delta x', 'Color', 'w', 'FontWeight', 'bold');
    ylabel('\Delta y', 'Color', 'w', 'FontWeight', 'bold');
    title(['Cage Breaking Geometry at t^* = ' num2str(target_tau)], 'Color', 'w', 'FontWeight', 'bold');
    colorbar('Color', 'w', 'LineWidth', 3, 'FontWeight', 'bold');
    set(gca, 'ColorScale', 'log');

    % Quantitative Cuts: Axis vs Diagonal
    % Define width of the slice (how thick is our "line")
    slice_width = bin_edges(2) - bin_edges(1);
    
    % 1. Cut along X-axis (where dy is near 0)
    % Find indices where dy is within +/- 1 bin of 0
    axis_mask = abs(all_dy) < slice_width;
    axis_vals = abs(all_dx(axis_mask));
    
    % 2. Cut along Diagonal (where |dx| approx |dy|)
    % Find indices where |abs(dx) - abs(dy)| is small
    diag_mask = abs(abs(all_dx) - abs(all_dy)) < slice_width;
    diag_vals = sqrt(all_dx(diag_mask).^2 + all_dy(diag_mask).^2); % Radial distance
    
    % Binning for 1D plot
    r_edges = linspace(0, limit, 50);
    [counts_axis, ~] = histcounts(axis_vals, r_edges, 'Normalization', 'pdf');
    [counts_diag, ~] = histcounts(diag_vals, r_edges, 'Normalization', 'pdf');
    r_centers = r_edges(1:end-1) + diff(r_edges)/2;
    
    % Plot Comparison
    figure('Color', 'k');
    semilogy(r_centers, counts_axis, 'r-', 'LineWidth', 3, 'DisplayName', 'Axis (0 deg)');
    hold on;
    semilogy(r_centers, counts_diag, 'b-', 'LineWidth', 3, 'DisplayName', 'Diagonal (45 deg)');
    legend('TextColor', 'w', 'Color', 'none', 'FontSize', 14, 'FontWeight', 'bold');
    
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 3, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Displacement Distance r', 'Color', 'w', 'FontWeight', 'bold');
    ylabel('Probability P(r)', 'Color', 'w', 'FontWeight', 'bold');
    title('Directional Probability Decay', 'Color', 'w', 'FontWeight', 'bold');
    grid on;
end

%% --- THERMALIZATION ANALYSIS --------------------------------------------------



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
    
    %% --- ANGULARLY RESOLVED PDF ANALYSIS --------------------------------------------------
    
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
    
    %% --- RADIAL PDF ANALYSIS --------------------------------------------------
           
    % --- determining the ideal gas distance distribution by Montecarlo ---
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - calculating pair distribution function\n', ic, S.bc, S.phi);
    if S.bc~=1
        if S.bc==1
            filepdfdenom = sprintf('PDFdenom_SBC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
        elseif S.bc==2
            filepdfdenom = sprintf('PDFdenom_PBCc_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
        elseif S.bc==3
            filepdfdenom = sprintf('PDFdenom_PBCFCC_%.0e_%.0e_%.0f.mat',S.rp,S.phi,S.N);
        end
        if exist(filepdfdenom,'file')
            load(filepdfdenom,'gdenominator');
        else
            gdenominator=PDFdenom(S,PDF,1e5,data_folder);
        end
    else
        shellcenters=PDF.pdfedges{3}(1:end-1)+0.5*(PDF.pdfedges{3}(2)-PDF.pdfedges{3}(1));
        shellvolumes=(4/3)*pi*PDF.pdfedges{3}(2:end).^3-(4/3)*pi*PDF.pdfedges{3}(1:end-1).^3;
        ndens=S.N/S.bv;
        gdenominator=ndens.*shellvolumes.*(S.N/2);
    end
    % ---
    % --- determine the g
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - average pairs numerator\n', ic, S.bc, S.phi);
    gnumerator=sum(RHO(:,idtherm:end),2)./(size(RHO,2)-idtherm); % average number of pairs in that bin
    g=gnumerator./gdenominator;
    g=[PDF.pdfedges{3}(1:end-1)+0.5*(PDF.pdfedges{3}(2:end)-PDF.pdfedges{3}(1:end-1)),g];
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

if dowt
       
    % --- per particle ---
    colls=EDGES{ic,1};
    clear EDGES
    colls(colls(:,1)==0,:)=[]; % eliminate filler
    idmap=uint32(zeros(2*S.N,2*S.N));
    q=1;
    for i1=1:2*S.N
        for i2=i1:2*S.N
            idmap(i1,i2)=q;
            idmap(i2,i1)=q;
            q=q+1;
        end
    end
    PPAIR=cell(q-1,1);
    clear i1 i2
    ncolls=size(colls,1);
    for i0=1:ncolls
        colls(i0,2)=idmap(colls(i0,2),colls(i0,3));
    end
    colls(:,3)=[];
    colls=sortrows(colls, [2 1]);
    idswitch=diff(colls(:,2))~=0;
    switches=find(idswitch);
    startend=[1,switches(1);switches+1,[switches(2:end,1);ncolls]];
    npairs=size(startend,1);
    for i0=1:npairs
        pairid=colls(startend(i0,1),2);
        PPAIR{pairid,1}=colls(startend(i0,1):startend(i0,2),1);
    end
    clear colls idswitch startend switches
    %
    ds=[];
    dspp=[];
    ccpp=[];
    qdspp=1;
    qccpp=1;
    for in=1:S.N
        cells=idmap(:,in);
        ds=diff(sort(vertcat(PPAIR{cells})));
        cells=idmap(1:S.N,in);
        ncells=numel(cells);
        dspp=nan(numel(ds),1);
        ccpp=nan(numel(ds),1);
        qdspp=1;
        qccpp=1;
        for icell=1:ncells
            dstemp=diff(PPAIR{cells(icell)});
            dspp(qdspp:qdspp-2+numel(PPAIR{cells(icell)}),1)=dstemp;
            dstemp=(dstemp==1);
            dstemp=regionprops(dstemp,'Area');
            dstemp=[dstemp.Area]';
            ccpp(qccpp:qccpp-1+numel(dstemp),1)=dstemp;
        end
        DSPP{in,1}=dspp;
        CCPP{in,1}=ccpp;
        DS{in,1}=ds;
        disp(in)
    end
    DSPP=vertcat(DSPP{:});
    DS=vertcat(DS{:});
    CCPP=vertcat(CCPP{:});
    
    perparte=(0:max(DS))';
    perpaire=(0:max(DSPP))';
    perpaircce=(0:max(CCPP))';
    clear dspp
    clear ds
    [perpartc,perparte]=histcounts(DS,perparte);
    perpart=double([perparte(1:end-1,1),perpartc']);
    perpart(:,3)=perpart(:,2)./sum(perpart(:,2));
    [perpairc,perpaire]=histcounts(DSPP,perpaire);
    perpair=double([perpaire(1:end-1,1),perpairc']);
    perpair(:,3)=perpair(:,2)./sum(perpair(:,2));
    [perpairccc,perpaircce]=histcounts(CCPP,perpaircce);
    perpaircc=double([perpaircce(1:end-1,1),perpairccc']);
    perpaircc(:,3)=perpaircc(:,2)./sum(perpaircc(:,2));
    perpart(perpart(:,2)==0,:)=[];
    perpair(perpair(:,2)==0,:)=[];
    perpaircc(perpaircc(:,2)==0,:)=[];
end

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
ISO(ic).WT.ppart=perpart;
ISO(ic).WT.ppair=perpair;
ISO(ic).CC=perpaircc;

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