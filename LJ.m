%% CLEANUP

clear all
close all
rng('shuffle')
warning('off', 'all');
% debugging=true

%% ADD PATHS TO STORAGE FOLDER

if exist('D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
else
    data_folder = 'G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
end
if exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
else
    output_folder='C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
end
    
    toolbox_folder = '..\ARBD_toolbox';
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)
% data_folder = 'C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\failure of Smoluchowskis collisional model in molecular systems\data\code\dt_rescaling_final';
% addpath(data_folder)

%% SERIES NAME

filenamecollisionseries='SBCvsPBC_%d.mat';
filenamecollisionseries_temp='SBCvsPBC_temp_%d.mat';

%% FIXED PHYSICAL PARAMETERS

C.T=298.15; %K (temperatura)
C.kb=1.38e-23; %J*K^-1 (costante di Boltzmann)
C.kbT=C.kb*C.T;
C.epsilon=1*C.kbT;
C.hydrationlayer=2.5e-10;

%% OPTIONAL SWITCHES

P.pdf=1;
P.ssf=1;
P.dens=0;
P.equipartition=0;
P.cluster=0;
P.exvol=0;
P.collkin=1;
P.correctionwindow=0; % window of the correction; 0 = no correction, 1e6 = S.rc, any other value X=X*S.rp;
P.convergencemode=1; % 1  by step numbers, 2 by coll numbers, 3 by coll rates convergence
P.pdfbins=[180,90]; % number of bins for azimuth, elevation
P.pdfthermints=10; % interval of steps between snapshots of pdfs to ascertain thermalization
P.rfstepsforeffdiffest=10000; % number of random walk steps that I need to collect to extract the alpha parameter

%% SIMULATION PARAMETERS

P.nodisp=1e7;  % size of displacement library
P.reps=1; % number of replicates
P.maxcoll=1e7; % minimum number of collisions to measure
P.maxsteps=1e6; % maximum number of steps to make when doing pdfs
P.kuhnmultiplier=200; % multiplier of taur that is used to either include or exclude persistence
P.kuhnmultiplierVACF=200; % multiplier of taur that is used to calculate autocorrelation
P.convergencewindows=[10000,1000]; % number of timesteps over which to describe the evolution of k in time; number of steps over which to evaluate convergence
P.evinterval=1e4; % ncoll interval between excluded volume estimations
P.noevpoints=1e5; % number of points generated to estimate accessible volume
P.ratio=1000; % maximum number of physical trajectory steps in each relaxation time to calculate the effective diffusivity. This is the number that ensures that the timescale is adaptive

%% VARIABLE SPACE and ASSUMPTIONS

C.solventproperties=solventlibrary();
C.phasesproperties=phaseslibrary();

V.solvent_list=[1]; % list of solvents according to the library above
V.phase_list=[2]; % list of phases according to the library above
V.r_range=[1e-08,1e-8,1,1]; % min, max, steps of particle radii, then log switch (1 if logscale)
V.N_range=[1000,1000,1,1]; % min, max, steps of particle numbers, then log switch (1 if logscale)
V.phi_range=[0.01,0.40,5,1]; % min, max, steps of volume fractions, then log switch (1 if logscale)
V.bc_range=[1,3,3]; % min, max, steps of boundary conditions (1:SBC; 2:PBC_cubic; 3:PBC_fcc; 4:BB)
V.bbm_range=[1,1,1]; % big box multiplier
V.LJpot_range=[1,2,2]; % LJ potential trigger (0 = HS; 1 = LJ; 2 = WCA)
V.T_range=[298.15,298.15,1]; % min, max, steps of temperatures

V.assumptions(1)=1; % exclude those where the Kuhn time is smaller than the trajectory step time
V.assumptions(2)=1; % exclude those where the particle is smaller than the solvent
V.assumptions(3)=1; % exclude those where the particle is lighter than the solvent

[V,CONDS]=simconditions_LJ(V,C,P);

%% EFFECTIVE DIFFUSIVITY CALCULATION

disp('collection and estimation of effective diffusivities for the conditions to be analyzed')

CONDS=effective_diffusivity(data_folder,CONDS,P,C);

%% GHOST PARTICLE MATRIX 

GPMAT=ghostparticlematrix();

%% SIMULATION EXECUTION

for ic=21
    
    if CONDS.alpha(ic,1)==0
        continue
    end
    
    % --- filename    
    filename=sprintf(filenamecollisionseries,ic);
    tempfilename=sprintf(filenamecollisionseries_temp,ic);
    %---

    % --- continue if filename exists
    if exist(filename,'file')>0 
       continue
    end
    % ---

    % --- initialize main variables from CONDS struct array
    S.timestep=V.c{ic,11};
    S.kbT=C.kbT;
    S.bc=CONDS.bc(ic,1);
    S.br=CONDS.br(ic,1);
    S.bv=CONDS.bv(ic,1);    
    S.rp=CONDS.rp(ic,1);
    S.N=CONDS.N(ic,1);
    S.phi=CONDS.phi(ic,1);
    S.solvent=CONDS.solvent(ic,1);
    S.phase=CONDS.phase(ic,1);
    S.alpha=CONDS.alpha(ic,1);
    S.esdiff=C.kbT/(6*pi*CONDS.eta(ic,1)*S.rp);
    S.stdx=sqrt(2*S.alpha*S.esdiff*CONDS.kt(ic,1));
    S.kt=CONDS.kt(ic,1);
    S.bbm=CONDS.bbm(ic,1);
    S.potential=CONDS.pot(ic);
    S.pot_corr=false;
    if S.potential>0
        S.pot_epsilon=C.epsilon;
        S.pot_sigma=2*S.rp+C.hydrationlayer;
        S.rc=3*S.pot_sigma;
        P.pdfsigmaplusthresh=2.76; % particle radius multipler within which to collect PDFsigma data
    else
        S.rc=2*S.rp;
        P.pdfsigmaplusthresh=3.2; % particle radius multipler within which to collect PDFsigma data
    end
    if P.correctionwindow==0
        S.correctionwindow=0;
    elseif P.correctionwindow==1e6
        S.correctionwindow=S.rc;
    else
        S.correctionwindow=P.correctionwindow*S.rp;
    end
    if S.bc==2 || S.bc==3
        [~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
        tokens = regexp(cmdout, '\d+', 'match');
        S.cacheSizeMB = max(str2double(tokens))/1024;
    end
    S.gtrig=S.rc;
    if S.bc==3 % THIS IF CONDITION IS VERIFIED AND CORRECT
        S.fcc.unitvecs=V.fcc_unitvecs;
        S.fcc.a1 = 2*S.br.*S.fcc.unitvecs(1,:);
        S.fcc.a2 = 2*S.br.*S.fcc.unitvecs(2,:);
        S.fcc.a3 = 2*S.br.*S.fcc.unitvecs(3,:);
        S.fcc.A=[S.fcc.a1; S.fcc.a2; S.fcc.a3]';
        S.fcc.invA = inv(S.fcc.A);
        S.fcc.reciprocal_vectors = S.fcc.invA';
        S.fcc.normals = [cross(S.fcc.a2,S.fcc.a3);
           cross(S.fcc.a3,S.fcc.a1);
           cross(S.fcc.a1,S.fcc.a2)];
        S.fcc.normals = S.fcc.normals ./ vecnorm(S.fcc.normals,2,2);
        S.fcc.frac_thresh = S.rc ./ abs(sum(S.fcc.A'.*S.fcc.normals,2))';
        [n1, n2, n3] = ndgrid(-1:1);
        S.fcc.shift_coeffs = [n1(:), n2(:), n3(:)];
        S.fcc.shift_coeffs(all(S.fcc.shift_coeffs==0,2),:) = []; % A robust way to remove [0,0,0]
        S.fcc.cart_shifts = S.fcc.shift_coeffs * S.fcc.A;
        clear n1 n2 n3
    elseif S.bc==2
        % for cubic PBC use diag(L,L,L) for MIC transforms
        L = 2*S.br;
        S.fcc.A = diag([L,L,L]);
        S.fcc.invA = diag([1/L,1/L,1/L]); 
    end
    clear nx ny nz dx dy dz SX SY SZ neighbor_offsets neighbor_linear valid lin all_ids remainder cx cy cz
    
    IDlist=double(linspace(1,S.N,S.N)'); % id list for all particles in the boundary to use when the position matrix 'p' is rebuilt at the beginning of every time step
    % ---

    % --- CALCULATION OF LJ FORCE -----
    if S.potential~=0
        clear H
        H=pot_force(S.potential,S.rc,30000,S.pot_sigma,S.pot_epsilon);
        H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
    end
    % ---

    % --- loop over replicates
    for irep=1:P.reps % loop over replicates 

        % ---- CALCULATION OF DISPLACEMENT LIBRARIES ----
        DISP=build_noise_library(S.stdx,P.nodisp);
        % ----

        % --- MONTECARLO CLAMP DETERMINATION ---
        if S.potential~=0
            S.pot_clamp=mcdClamp(P.nodisp,DISP,S.esdiff,S.timestep,H(H(:,1) >= 0.8 * S.pot_sigma, :),C.kbT);            
        end
        % ---

        % ---- INITIALIZING PDF, SSF, DENSITY, EQUIPARTITION ---- 
        if P.pdf==1
            PDF=pdf_initialization(S,P);
        end
        if P.ssf==1
            SSF=ssf_initialization(S);
        end
        % calculation of bins for density diagnostics must depend on the
        % shape of the boundary. In the case of SBC we look at the density
        % (mass and number) distribution as a function of distance from the
        % origin as well as azimuth and elevation. In the case of PBCs we
        % do the same (within the sphere circumscribed by the boundary) AND
        % we also look at the density distributions in slabs parallel to
        % the walls.
        if P.dens==1
            DCOMP = densityEngine('setup', [], S, P, [],data_folder);
            disp('density mapping setup - completed')
        end
        if P.equipartition==1
            EQUIP = equipartitionEngine('setup', [], S, P, []);
            disp('equipartition mapping setup - completed')
        end
        % ---

        % --- STARTING POSITIONS & DETERMINING SBC RADIAL DISPLACEMENT ASYMMETRY CORRECTION ---
        if S.bc==1
            if S.potential==1, potname='lj'; elseif S.potential==2, potname='wca'; else potname='hs'; end
            filestartingconfiguration = sprintf('START_SBC_%s_%.0e_%.0e_%.0f_%.1f_%.1e.mat',...
                        potname,S.rp,S.phi,S.N,S.pot_epsilon/S.kbT,S.pot_sigma);
            if exist(filestartingconfiguration,'file')
                load(filestartingconfiguration,'p','pgp')
            else
                [p,pgp]=sbc_setup_sgd_v9(S,PDF,[],data_folder);
            end
            disp('initial particle configuration - determined')
            if S.pot_corr,disp('radial displacement asymmetry - determined'); end
            % --- CREATE MASTER GHOSTLIST ----
            pgp0= struct('p',{},'ID',{},'active',{},'lastStep',{});
            pgp0(1).p=pgp;
            pgp0(1).ID=IDlist;
            pgp0(1).active=zeros(S.N,1); % ACTIVE flag
            pgp0(1).lastStep=zeros(S.N,1);% LAST STEP ACTIVE
            disp('master ghost list - created')
            % ---
        else
            p=startingPositions_lj(S);
        end
        % ---

        % ---- INITIALIZE SERVICE VARIABLES AND COUNTERS ----
        disp('initializing service variables and counters')
        qd=1; % displacement counter
        qs=1; % step counter
        qc=0; % collisions counter (only counts valid two-body collisions)
        qev=1; % counter of AV cell array data storage;
        kconverged=0;% flag to identify when the value of k has converged and when, therefore the sim can be concluded
        counterflag=0; % flag to tell me when to output the sim state - without it, the simstate gets broadcasted too many times
        DEGREES{ic,irep}=zeros(100,1); % initialize the DEGREES cell array that stores the histogram of the node degrees of the clusters formed in the simulation
        CLUSTERS{ic,irep}=[]; % initialize the CLUSTERS cell array that stores the numbers of nodes (col 2), the number of edges (col 3) and the completeness (col 4) of ALL collision clusters occurring and the step at which they occur (col 1) 
        EDGES{ic,irep}=uint32(zeros(1e6,3)); % initialize the EDGES cell array that stores the IDs of all colliders. This is essential when we are going to havce to look at the pair waiting times at low phi.
        qedges=1;
        AV{ic,irep}=[]; % accessible volume estimations;
        
        % ---- RUN UNTIL CONVERGENCE -----        
        disp('starting simulation')
        while kconverged==0
            
            % --- SBC GHOST ACTIVATION AND DEACTIVATION -----
            if S.bc==1
                prho=vecnorm(p,2,2); % calculate distances from origin of reals
                idxgp=prho>(S.br-S.rc); % mask of active ghosts
                pgp=pgp0.p(idxgp,:); % active ghosts positions
                ptemp=[p;pgp]; % concatenated list of reals and active ghosts
                S.Na=size(ptemp,1); % number of particles (reals + active ghosts)
            end
            % ---
            
            % ---- CALCULATION OF DISPLACEMENTS BY POTENTIALS --------------
            if S.potential~=0
                if S.bc==1 % SBC
                    % calculate the displacement components due to potentials for all particles (reals and active ghosts)
                    disppot=potential_displacements_v13(ptemp, S, H, H_interpolant, 0);
                    % extract potential displacements for active ghosts
                    disppotgp=disppot(S.N+1:end,:);
                    % extract potential displacements for reals
                    disppot=disppot(1:S.N,:);
                else % NON-SBC conditions (PBC)
                    % Under MIC we evaluate pairwise displacements using
                    % minimum-image convention on the real particle list.
                    disppot=potential_displacements_v13(p, S, H, H_interpolant, 0,S.cacheSizeMB);
                end
            end
            % --------------------------------------------------------------

            % ---- DISPLACE PARTICLE POSITIONS -----------------------------
            if S.bc~=1 % non-SBC conditions
                displacements=DISP(qd:qd+S.N-1,:); % extract displacements for reals
                qd=qd+S.N; % update displacement library counter
                if S.potential~=0 % if not HS add potential displacement component
                    displacements=displacements+disppot;
                end
                p=[p,IDlist,displacements,p+displacements];
                % ----  BIG BOX BOUNDARY INTERACTIONS -----
                if S.bc==4
                    % this function teleports particle that touch the hard
                    % boundary to avoid progressive condensation and therefore
                    % depletion from the core of the domain
                    p=bigboxTeleport(S,p); 
                end
            elseif S.bc==1 % SBC
                % extract brownian displacements for reals and ghosts
                    displacements=DISP(qd:qd+S.Na-1,:); % extract displacements for reals and ghosts
                    qd=qd+S.Na; % update displacement library counter by reals + active ghosts
                    dispgp=displacements(S.N+1:end,:); % extract brownian displacement for active ghosts
                    displacements=displacements(1:S.N,:); % extract brownian displacement for reals
                % add potential displacements, if necessary, to both reals and ghosts displacements
                    if S.potential~=0 
                        displacements=displacements+disppot;
                        dispgp=dispgp+disppotgp;
                    end
                % convert ghost displacements to tangential only
                    pgp_norm = vecnorm(pgp,2,2);
                    pgp_norm(pgp_norm==0) = eps;
                    pgp_dir = pgp ./ pgp_norm;
                    v_rad_component = sum(dispgp .* pgp_dir, 2); % extract radial component to the total ghost displacement          
                    v_tan = dispgp - (v_rad_component .* pgp_dir); % subtract radial component from the 3D displacement
                % apply displacements to active ghosts 
                    pgp2_temp = pgp + v_tan; % move ghosts tangentially creating tentative position
                % correct real displacements by ASYMCORR
                    if S.pot_corr
                        dispcorr=zeros(size(p,1),3);
                        % mask particles within 2*S.rp of the boundary
                        idxasym = prho > (S.br - S.correctionwindow);
                        if any(idxasym)
                            % interpolate correction value for each masked particle
                            delta_r = interp1(ASYMCORR(:,1), ASYMCORR(:,2), prho(idxasym), 'pchip', 0);
                            % compute unit radial vectors
                            rho_hat = p(idxasym,:) ./ prho(idxasym);
                            % apply correction only to the radial component of the displacement
                            dispcorr(idxasym,:) = dispcorr(idxasym,:) + delta_r .* rho_hat;
                        end
                        displacements=displacements+dispcorr;
                    end
                    p2=p+displacements; % move reals creating tentative position
                % reset tether distance of initial ghosts
                    p2rho = vecnorm(p2,2,2); % get tentative real position norms
                    pgp2_temp_norm = vecnorm(pgp2_temp,2,2); % get tentative ghost positions norms
                    pgp2_temp_norm(pgp2_temp_norm==0) = eps; % (guard for tiny norms)
                    pgp2_dir = pgp2_temp ./ pgp2_temp_norm; % get versors of tentative ghost positions
                    pgp2 = pgp2_dir .* (2*S.br-p2rho(idxgp,:)); % set distance of tentative ghosts to the tether distance from tentative real positions
                    pgp2disps = pgp2-pgp; % store displacements for eventual reset of positions in case of HS potential
                    pgp=[pgp,IDlist(idxgp,:),pgp2disps,pgp2];
                % update ghost master list
                    pgp0.active(:,:)=0;
                    pgp0.active(pgp(:,4),1)=1;
                    pgp0.lastStep(pgp(:,4),:)=qs;
                    p=[p,IDlist,displacements,p+displacements]; % the p matrix contains [positions at the beginning of step, particle ID, displacements, positions at the end of the step]
                    p=[p;pgp];
            end
            % --------------------------------------------------------------

            % ---- PROCESSING COLLISIONS
            if P.collkin==1
                dists=pdist(p(:,8:10))';
                idxcoll=sign(2*S.rp-dists)+1; % calculate all distances, subtract them from 2r, and pick the sign. If -1 then distances>2r, if +1 then distance<2r. Then add 1.
                if sum(idxcoll,'all')>0
    
                    % --- labeling individual particles to ensure only the right particles are considered
                    Nreal=size(p,1); % total number of particles (real + ghost)
                    p(:,11)=zeros(Nreal,1); % index to account for particles that are moved or not during a single step
                    if S.bc==1
                        p(:,12:13)=ones(Nreal,2); % index to label the particles that are inside the boundary either real or ghosts)
                        p(vecnorm(p(:,1:3),2,2)>(S.br+S.rp),12)=0; % col 12 becomes zero if initial positions are FULLY outside the box
                        p(vecnorm(p(:,8:10),2,2)>(S.br+S.rp),13)=0; % col 13 becomes zero if final positions are FULLY outside the box
                    elseif S.bc==2
                        p(:,12:13)=ones(Nreal,2); % index to label the particles that are inside the boundary either real or ghosts)
                        p(sum(abs(p(:,1:3))>(S.br+S.rc),2)>0,12)=0; % col 12 becomes zero if initial positions are FULLY outside the box
                        p(sum(abs(p(:,8:10))>(S.br+S.rc),2)>0,13)=0; % col 13 becomes zero if final positions are FULLY outside the box
                    elseif S.bc==3
                        p(:,12:13)=ones(Nreal,2); % index to label the particles that are inside the boundary either real or ghosts)
                        p(oob_fcc(p(:,1:3),S,-S.rc),12)=0; % col 12 becomes zero if initial positions are FULLY outside the box
                        p(oob_fcc(p(:,8:10),S,-S.rc),13)=0; % col 13 becomes zero if final positions are FULLY outside the box
                    end
                    % ---
                    
                    flagcoll=1; % flag that makes me exit the collision analysis loop when no more valid collisions are detected; 1 - collisions are present; 2 - collisions have been processed but some might remain; 3 - collisions
                    edges=[]; % initialize list of collision types
                    while flagcoll>0 % repeat collider identification and motion until no collision remains                    
                        % -- COLLISION CHECK AT THE BEGINNING OF EVERY LOOP except the first ----  
                        if flagcoll==2 && S.potential==0  % avoid rechecking if this is the first cycle of collision detangling                        
                            idxcoll=pdist(p(:,8:10))'-(2*S.rp)<0;
                            if S.bc~=4 && (sum(idxcoll)==0)
                                flagcoll=0;
                                counterflag=0;
                                break
                            elseif S.bc==4 && (sum(idxcoll)==0 || isempty(colliders)==1)
                                flagcoll=0;
                                counterflag=0;
                                break
                            end
                        elseif flagcoll==2 && S.potential~=0 % do not reprocess collisions in the presence of soft potentials
                            counterflag=0;
                            break
                        end
                        % ---
    
                        % --- IDENTIFY THE COLLIDERS
                        Nreal=size(p,1);
                        collpairs=find(idxcoll);
                        bin=ceil(-0.5*sqrt(8*nchoosek(Nreal,2)-8*collpairs+1)+Nreal-0.5);
                        coll1=bin;
                        binedges=0.5*(Nreal-bin)-0.5*(Nreal-bin).^2+nchoosek(Nreal,2);
                        coll2=Nreal-binedges+collpairs;
                        % ---
                        
                        % --- MAKE A LIST OF THE COLLIDERS (position of collider #1, position of collider #2, row# of collider #1 and #2), theN identify those who collide OUTSIDE the boundary, and then eliminate the midpoint coordinates from the collider matrix (for economy) ---
                        if S.bc~=4
                            % COL9=RESET BIT FOR COLLIDER 1 (1 means particle HAS BEEN RESET ALREADY)
                            % COL10=RESET BIT FOR COLLIDER 2 (1 means particle HAS BEEN RESET ALREADY)
                            % COL11=START-OF-STEP-OOB BIT FOR COLLIDER 1 (0 means particle WAS OOB at start of step)
                            % COL12=START-OF-STEP-OOB BIT FOR COLLIDER 2 (0 means particle WAS OOB at start of step)
                            % COL13=END-OF-STEP-OOB BIT FOR COLLIDER 1 (0 means particle IS OOB at end of step)
                            % COL14=END-OF-STEP-OOB BIT FOR COLLIDER 2 (0 means particle IS OOB at end of step)
                            colliders=[p(coll1,8:10),p(coll2,8:10),coll1,coll2,p(coll1,11),p(coll2,11),p(coll1,12),p(coll2,12),p(coll1,13),p(coll2,13)];
                            % identify all collisions among non-reals
                            collidersnonreal=colliders(:,7)>S.N & colliders(:,8)>S.N;
                            % if the first sweep through the collisions has
                            % been done and you find that the only collisions
                            % left are those among non-reals then proceed. This
                            % implies that some unresolved collisions will stay
                            % unresolved which is necessary. Remember that
                            % ghosts will never collide if their reals are not
                            % colliding but compensators will collide with
                            % ghosts and with each other. The compensators are
                            % repelling each other by HS at the beginning of
                            % cycle so that shouldn't be an issue. but the
                            % collisions between ghosts and compensators have
                            % to be ignored to avoid spooky action at a
                            % distance. In the case of soft potentials the
                            % compensators will still move away from ghosts.
                            if flagcoll==2 & sum(collidersnonreal)==size(colliders,1)
                                break
                            end
                        else
                            colliders=[p(coll1,8:10),p(coll2,8:10),coll1,coll2,p(coll1,11),p(coll2,11)];
                        end                    
                        % ---
    
                        % --- CALCULATE COLLISION MIDPOINTS 
                        if S.bc==1
                            mp=vecnorm(colliders(:,1:3)+0.5*(colliders(:,4:6)-colliders(:,1:3)),2,2); % midpoint calculation                     
                            idxmp=mp>S.br; % index matrix of collisions occuring outside the boundary
                        elseif S.bc==2 || S.bc==3
                            % Under PBC (cubic or FCC), midpoint OOB checks are invalid.
                            % MIC decides collisions entirely. Keep midpoint only for diagnostics.
                            tempvector = colliders(:,4:6) - colliders(:,1:3);
                            mp = colliders(:,1:3) + 0.5 * tempvector;  % midpoint (diagnostic only)
                            idxmp = false(size(mp,1),1);               % no OOB rejections under PBC
                        elseif S.bc==4
                            tempvector=(colliders(:,4:6)-colliders(:,1:3));
                            mp=colliders(:,1:3)+(tempvector./norm(tempvector).*(0.5*norm(tempvector))); % midpoint calculation
                            idxmp=sum(abs(mp)>(S.br/S.bbm),2)>0; % index matrix of collisions occuring outside the boundary - which are excluded from consideration only in the final count.
                        end
                        % IDS{qs,1}=[IDS{qs,1};colliders(:,7:8)];
                        if sum(sum(colliders(~idxmp,9:10),2)==2)>0 % edge case check to see if sim finds collisions between particles that already have been moved back
                            disp('houston, we have a problem')
                            pause
                        end
                        % if the first sweep through the collisions has been
                        % done and all the collisions left are outside the
                        % boundary (even if they involve reals), then you
                        % proceed. I am not sure if this is unavoidable but it
                        % seems it is. I guess it weakens the compensation by
                        % compensators in the case of HS potentials but it is
                        % an edge case to avoid breaking the overlap resolution
                        if flagcoll==2 & isempty(colliders(~idxmp,:))==1
                            flagcoll=0;
                            counterflag=0;
                            break
                        end
                        % ---
    
                        % --- MOVE BACK COLLIDERS --- code check OK
                        if isempty(colliders)==0
                            % --- store info about the OVERLAPPERS whether in HS or LJ
                            ids=colliders(:,7:8); % row #s of the two colliders (columns) in each two-body collisions (rows)                        
                            edges=[edges;ids(~idxmp,:)]; % add to the iterative list of the row #s of the two colliders (columns) in each two-body collisions (rows), UNLESS the collisions happened outside the boundary. Therefore 'edges' works as a list of the valid collisions that accumulate during a time-step
                            % --- move back all particles (real and associated ghost) that have collided if on HS 
                            if S.potential==0
                                if ~isempty(ids)
                                    if S.bc~=4
                                        ids=nonzeros(unique(reshape(ids,[numel(ids) 1]))); % reshape as column vector
                                        ids=unique(nonzeros(p(ids,4))); % use the row #s of all colliders to find out all the particle IDS of the colliders.  
                                        ids=find(ismember(p(:,4),ids)); % use the particle IDS to find all particles (i.e., the row #s), real or ghost whose ID matches that of a collider.
                                    end
                                    p(ids,8:10)=p(ids,1:3); % move those back
                                    p(ids,11)=1;% update the list of moved particles
                                    if S.bc~=4
                                        p(ids,13)=p(ids,12); % the columns indicating out of boundary must be updated as well.
                                    end
                                end
                            end
                            flagcoll=2;
                            % ---
                        end
                        % ---
                    end
                    % ---
    
                    % --- MANY-BODY COLLISIONS ANALYSIS and COLLISION DATA STORAGE ---
                    % this part takes place after all collisions have been processed and no more collisions are found.
                    % The 'edges' variable contains the pairs of particles that have collided during this timestep
                    if isempty(edges)==0
                        tempedges=[ones(size(edges,1),1).*qs,edges];
                        qc=qc+size(tempedges,1);
                        if qc>size(EDGES{ic,irep},1)
                            EDGES{ic,irep}(qc+1e6,1)=0;
                        end
                        EDGES{ic,irep}(qc-size(tempedges,1)+1:qc,:)=tempedges;
                        if P.cluster==1
                            [degrees,gq,twobodycolls]=graphAnalysis(edges,qs);
                            DEGREES{ic,irep}(1:length(degrees),1)=DEGREES{ic,irep}(1:length(degrees),:)+degrees; % store degree information in 'DEGREES' cell array
                            CLUSTERS{ic,irep}=[CLUSTERS{ic,irep};gq]; % store cluster information in 'CLUSTERS' variable  
                        end
                    end
                    % ---
                end
            end
            % ----

            % --- PBC POSITION WRAPPING
            % Wrap positions by MIC to canonical cell (if PBC)
            if S.bc==2 || S.bc==3
                p(:,8:10)=mic_wrap_positions(p(:,8:10), S);
            end

            % --- EQUIPARTIION ANALYSIS -----------------------------------
            if P.equipartition==1
                EQUIP = equipartitionEngine('accumulate', EQUIP, S, P, p);
            end
            % -------------------------------------------------------------

            % --- PROMOTION/DEMOTION IN SBC - reset in BB/PBC -------------
            if S.bc==1
                pgp0.p(pgp(:,4),:)=pgp(:,8:10);
                p=p(1:S.N,8:10);
                pgp=pgp(:,8:10);
                idxrgswap=vecnorm(p,2,2)>S.br;
                % initialize start position array for the next timestep 
                realtoswap=p(idxrgswap,:);
                ghoststoswap=pgp0.p(idxrgswap,:);                
                pgp0.p(idxrgswap,:)=realtoswap;
                p(idxrgswap,:)=ghoststoswap;
            else
                p=p(1:S.N,8:10);
            end
            % -----------------------------------------------------------

            % --- DENSITY ANALYSIS ----------------------------------------
            if P.dens==1
                if S.bc==1
                    DCOMP = densityEngine('accumulate', DCOMP, S, P, [p;pgp]);
                else
                    DCOMP = densityEngine('accumulate', DCOMP, S, P, p);
                end
            end
            % ------------------------------------------------------------

            % ---- PDF ANALYSIS -------
            Nreal=size(p,1);
            if P.pdf==1
                if mod(qs,P.pdfthermints)==0 | qs==1
                    if S.bc==1
                        ppdf=p;
                    else
                        ppdf=p;
                    end
                    clear PDFD
                    % GET DISTANCE DISTRIBUTION
                    % calculate all distance vectors
                    if S.bc==1 || S.bc==4
                        PDFD(:,:,1)=(ppdf(:,1)-ppdf(:,1)'); 
                        PDFD(:,:,2)=(ppdf(:,2)-ppdf(:,2)');
                        PDFD(:,:,3)=(ppdf(:,3)-ppdf(:,3)');
                        PDFD=reshape(PDFD,[],3);
                        % eliminate distances between identical particles
                        PDFD(PDFD(:,1)==0,:)=[];
                    else
                        PDFD = mic_all_pair_displacements(p(:,1:3), S); % returns Mx3 array of minimum-image displacement vectors (M = N*(N-1))
                    end
                    
                    if S.bc==3
                        fccrotatev=[1 1 1] / norm([1 1 1]);
                        PDFD=FCCrotate(PDFD,fccrotatev);
                    end
                    % convert to spherical coordinates
                    [pdfaz,pdfel,pdfrho]=cart2sph(PDFD(:,1),PDFD(:,2),PDFD(:,3));
                    PDFD=[pdfaz,pdfel,pdfrho];
                    % extract azimuth and elevation from all distances smaller
                    % than sigma treshold
                    PDFDsigma=PDFD(PDFD(:,3)<(S.rp*P.pdfsigmaplusthresh),1:2);
                    % cumulative add to histograms of az, el for PDFsigma and
                    % az,el,rho for PDF
                    [temppdf,PDF.pdfedges{1}]=histcounts(PDFD(:,1),PDF.pdfedges{1});                   
                    PDFT{floor(qs/P.pdfthermints)+1,1}=uint32(temppdf'./2);
                    [temppdf,PDF.pdfedges{2}]=histcounts(PDFD(:,2),PDF.pdfedges{2});                    
                    PDFT{floor(qs/P.pdfthermints)+1,2}=uint32(temppdf'./2);
                    [temppdf,PDF.pdfedges{3}]=histcounts(PDFD(:,3),PDF.pdfedges{3});                    
                    PDFT{floor(qs/P.pdfthermints)+1,3}=uint32(temppdf'./2);                    
                    [temppdf,PDF.pdfedges{1}]=histcounts(PDFDsigma(:,1),PDF.pdfedges{1});                    
                    PDFT{floor(qs/P.pdfthermints)+1,4}=uint32(temppdf'./2);                    
                    [temppdf,PDF.pdfedges{2}]=histcounts(PDFDsigma(:,2),PDF.pdfedges{2});                    
                    PDFT{floor(qs/P.pdfthermints)+1,5}=uint32(temppdf'./2);                    
                end
            end
            % --------------------------

            % ---- STATIC STRUCTURE FACTOR ANALYSIS -------
            if P.ssf==1
                if mod(qs,P.pdfthermints)==0 | qs==1
                    ssf_ind=floor(qs/P.pdfthermints)+1;
                    tempp=p(1:S.N,1:3);
                    SSF=ssf_accumulation(tempp',SSF,ssf_ind,S.N);
                end
            end
            % --------------------------

            % --- EXCLUDED VOLUME CALCULATION BY MONTECARLO ---           
            if P.exvol==1
                if qc>P.evinterval*qev
                    axv=excludedVolume(S,P,p);
                    AV{ic,irep}(qev,:)=[qs,qc,axv];
                    qev=qev+1;
                end
            end
            % ---

            % ----- RESHUFFLE DISPLACEMENT LIBRARIES if they run out----
            if qd+2*S.N>P.nodisp
                qd=1;
                DISP=DISP(randperm(P.nodisp),:);
            end
            % -----

            % --- CONVERGENCE CHECK AFTER MINIMUM # OF COLLISIONS IS ACHIEVED ---
            if P.convergencemode==3 && qc>P.maxcoll && mod(qs,10^4)==0
                temp=DATA{ic,1}; % first two columns are step number and number of valid 2-body colls
                temp(:,3)=temp(:,1).*S.kt; % time corresponding to the step
                temp(:,4)=cumsum(temp(:,2)); % cumulative sum of collisions
                temp(:,5)=temp(:,4)./(temp(:,3).*S.bv); % collision rate: total number of collisions divided by time and volume
                tempconc=S.N/S.bv; % particle concentration;
                temp(:,6)=temp(:,5)./(tempconc^2); % rate constant (assuming second order): collision rate/concentration^2
                % tempint is an interpolation: time and rate constant are
                % the sample points; the queried points are 10000 equispaced 
                % points between the minimum time and the maximum time.
                % This is done to have an equispaced function of
                % rateconstant(time)
                tempint=interp1(temp(:,3),temp(:,6),min(temp(:,3)):(max(temp(:,3))-min(temp(:,3)))/(P.convergencewindows(1)-1):max(temp(:,3)))';
                tempint=tempint(end-P.convergencewindows(2)+1:end,:); % extract only the last P.convwindow points (usually 1000)                
                % the convergence check is by comparing the beginning
                % values of the window with the end values and see whether
                % they have changed by more than 1%. If not the series is
                % considered converged.
                beginning=mean(tempint(1:100)); 
                ending=mean(tempint(end-100:end));
                disp((abs(ending-beginning)/beginning)*100)
                if (abs(ending-beginning)/beginning)<0.01
                    kconverged=1;
                end
            elseif P.convergencemode==1
                if qs>=P.maxsteps
                    kconverged=1;
                end
            elseif P.convergencemode==2
                if qc>=P.maxcoll
                    kconverged=1;
                end
            end
            % ---

            % --- COUNTERS ---            
            qs=qs+1;
            if counterflag==0 && mod(qs,1e2)==0
                disp({filename,ic,irep,log10(qs),log10(qc)})
                counterflag=1; 
            end
            % ---

            % --- SAVING TEMP FILE ---
            if mod(log2(qs),1)==0
                if ~isfile(tempfilename)
                    % First time: Create the file with the core struct
                    save([output_folder,'\',tempfilename], 'C', 'EDGES', 'P', 'V', 'S', '-v7.3');
                else
                    % File exists: Update core variables (Overwrite them in place)
                    save([output_folder,'\',tempfilename], 'C', 'EDGES', 'P', 'V', 'S', '-append');
                end
                if P.pdf==1
                    save([output_folder,'\',tempfilename],'PDFT','PDF','-append')
                end
                if P.ssf==1
                    save([output_folder,'\',tempfilename],'SSF','-append')
                end
                if P.cluster==1
                    save([output_folder,'\',tempfilename],'CLUSTERS','DEGREES','-append')
                end
                if P.dens==1
                    save([output_folder,'\',tempfilename],'DCOMP','-append')
                end
                if P.equipartition==1
                    save([output_folder,'\',tempfilename],'EQUIP','-append')
                end
                if P.exvol==1
                    save([output_folder,'\',tempfilename],'AV','-append')
                end
            end
            % ---
            

        end
    end
    save([output_folder,'\',filename],'C','EDGES','P','V','S','-v7.3')
    if P.pdf==1
        save([output_folder,'\',filename],'PDFT','PDF','-append')
    end
    if P.ssf==1
        save([output_folder,'\',filename],'SSF','-append')
    end
    if P.cluster==1
        save([output_folder,'\',filename],'CLUSTERS','DEGREES','-append')
    end
    if P.dens==1
        save([output_folder,'\',filename],'DCOMP','-append')
    end
    if P.equipartition==1
        save([output_folder,'\',filename],'EQUIP','-append')
    end
    if P.exvol==1
        save([output_folder,'\',filename],'AV','-append')
    end
    if isfile([output_folder,'\',tempfilename])
        delete([output_folder,'\',tempfilename]);
    end
    clear DATA DEGREES AV CLUSTERS PDF* SSF* DS
end
clearvars -except S V P C