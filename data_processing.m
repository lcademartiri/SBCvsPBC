%% CLEANUP

clearvars -except ISO
close all
rng('shuffle')
warning('off', 'all');

%% SERIES NAME

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
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)

filenamecollisionseries='SBCvsPBC_%d.mat';
plottingenabled=1;

%% ANALYSIS

conditionstoprocess=[];

for ic=27
    filename=sprintf(filenamecollisionseries,ic);
    if exist(filename,'file')>0
        conditionstoprocess=[conditionstoprocess,ic];
    end
end

for ic=conditionstoprocess
    clearvars -except ISO conditionstoprocess ic filenamecollisionseries plottingenabled
    filename=sprintf(filenamecollisionseries,ic);
    load(filename,'C','P','V','S','EDGES','PDFT','SSF','DCOMP','DS','H','PDF','ASYMCORR','EQUIP')
    
    % --- ORGANIZE DATA ----------------------------------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - organizing data\n', ic, S.bc, S.phi);
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
    MEAN_SSF100=mean(SSF.SSF100,1);
    MEAN_SSF110=mean(SSF.SSF110,1);
    MEAN_SSF111=mean(SSF.SSF111,1);
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
    MMEAN_SSF100=movmean(SSF.SSF100,100,1);
    MMEAN_SSF100=MMEAN_SSF100./sum(MMEAN_SSF100,2);
    MMEAN_SSF110=movmean(SSF.SSF110,100,1);
    MMEAN_SSF110=MMEAN_SSF110./sum(MMEAN_SSF110,2);
    MMEAN_SSF111=movmean(SSF.SSF111,100,1);
    MMEAN_SSF111=MMEAN_SSF111./sum(MMEAN_SSF111,2);
    % ---
    % --- variance of moving averages
    STD_MMEAN_AZ=std(MMEAN_AZ)';
    STD_MMEAN_EL=std(MMEAN_EL)';
    STD_MMEAN_RHO=std(MMEAN_RHO)';
    STD_MMEAN_AZS=std(MMEAN_AZS)';
    STD_MMEAN_ELS=std(MMEAN_ELS)';
    STD_MMEAN_SSF100=std(MMEAN_SSF100,0,2)';
    STD_MMEAN_SSF110=std(MMEAN_SSF110,0,2)';
    STD_MMEAN_SSF111=std(MMEAN_SSF111,0,2)';
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
    % -------------------------------------------------------------------

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
    MEAN_SSF100=mean(SSF.SSF100(idtherm:end,:),1);
    MEAN_SSF100=MEAN_SSF100./sum(MEAN_SSF100);
    MEAN_SSF110=mean(SSF.SSF110(idtherm:end,:),1);
    MEAN_SSF110=MEAN_SSF110./sum(MEAN_SSF110);
    MEAN_SSF111=mean(SSF.SSF111(idtherm:end,:),1);
    MEAN_SSF111=MEAN_SSF111./sum(MEAN_SSF111);
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

    
    % --- PDF AND SSF ANALYSIS --------------------------------------------------
       
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
    % -------------------------------------------------------------------

    % --- WAITING TIME ANALYSIS --------------------------------------------------
       
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
    % --- consecutive ----
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

    % ---------------------------------------------------------------------

    % --- COLLATING DATA ------------------------------------------------
    fprintf('condition: %d - boundary condition: %d -  phi: %.3f - collating data\n', ic, S.bc, S.phi);
    ISO(ic).g=g;
    ISO(ic).pdfedges=PDF.pdfedges;
    ISO(ic).thermalized_mAZ=[PDF.pdfedges{1}(1:end-1),MEAN_AZ];
    ISO(ic).thermalized_mAZS=[PDF.pdfedges{1}(1:end-1),MEAN_AZS];
    ISO(ic).thermalized_mEL=[PDF.pdfedges{2}(1:end-1),MEAN_EL];
    ISO(ic).thermalized_mELS=[PDF.pdfedges{2}(1:end-1),MEAN_ELS];
    ISO(ic).thermalized_mRHO=[PDF.pdfedges{3}(1:end-1),MEAN_RHO];
    ISO(ic).thermalized_mSSF100=MEAN_SSF100';
    ISO(ic).thermalized_mSSF110=MEAN_SSF110';
    ISO(ic).thermalized_mSSF111=MEAN_SSF111';
    ISO(ic).stdmmAZS=[[1:numel(STD_MMEAN_AZS)]'.*S.kt,STD_MMEAN_AZS];
    ISO(ic).stdmmELS=[[1:numel(STD_MMEAN_ELS)]'.*S.kt,STD_MMEAN_ELS];
    ISO(ic).thermalizationtime=idtherm*S.kt;
    ISO(ic).thermalizationfit=bestfit;
    ISO(ic).S=S;
    ISO(ic).P=P;
    ISO(ic).V=V;
    ISO(ic).C=C;
    ISO(ic).WT.ppart=wtppart;
    ISO(ic).WT.ppair=wtppair;
    ISO(ic).CC=cc;
    % -------------------------------------------------------------------
    
end