clear all
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


load('SBCvsPBC_25.mat');
clear POS PDF PDFT SSF POS
colls=EDGES{25,1};
colls(colls(:,1)==0,:)=[];
clear EDGES
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
%%
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
    for ic=1:ncells
        dstemp=diff(PPAIR{cells(ic)});
        dspp(qdspp:qdspp-2+numel(PPAIR{cells(ic)}),1)=dstemp;
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