function [Vs,Ps,Es]=startintcomm(randkey,spnum,syssz,landprm,spsprm,res)
% setup definition and parameters for a community with spatial interactions
if(nargin<4)  landprm = [1 0.2]; end; % default parameters for landscape
if(nargin<5)  spsprm  = []; end;
if(nargin<6)  res  = 1; end;
if(length(syssz)>1)
    szs = syssz;
else
    szs = [syssz syssz]; % assume both axes have the same length, if only one number is given
end;

% defspsprm includes 9 differnet species parameters:
% [1-4: niche-cent-span niche-cent-offset niche-wdth-span niche-wdth-offset ...]
% [5-9: self-reg-strength dispersal interaction-size interaction-cutoff interaction-normalized-share]

% default values for species parameters
defspsprm=[0.6 0.2 0.02 0.1 0.5 1e-4 1 0.01 0.5];
% use default values if not given
spsprm(length(spsprm)+1:length(defspsprm))=defspsprm(length(spsprm)+1:length(defspsprm));

% random key
rng(randkey);

% define landscape (see more details in the randlandscape function)
land = randlandscape(szs,landprm,randkey);
% set minimum value of landscape to 0
land = land-min(land(:));
% normalize the landscape so max equals 1, and reshape it to a vector format
land = reshape(land/max(land(:)),prod(szs),1);

% define center of niche per species
cent = rand(1,spnum)*spsprm(1)+spsprm(2);
% define width of niche per species
wdth = randn(1,spnum)*spsprm(3)+spsprm(4);
% define the local conditions each species experiences, according to the landscape
locr = exp(-0.5*((repmat(land,1,spnum)-repmat(cent,prod(szs),1))./repmat(wdth,prod(szs),1)).^2);

% define the interaction values using a uniform random distribution
intmat = -rand(spnum)*spsprm(5);
% set self-interactions to -1 
intmat(logical(eye(spnum)))=-1;

% define the model parameters and functions 
Ps=struct('LocFunc',@L_GLVwSI,'SpaFunc',@S_RD,'IntegFunc',@I_FDE,'r',locr,'A',intmat,'Ds',spsprm(6)*ones(1,spnum),'id',spsprm(7:9),'VarNum',spnum,'Nx',szs(1)*res,'Lx',szs(1),'Ny',szs(2)*res,'Ly',szs(2));
% define auxiliary parameters
Es=struct('TsSize',0.01,'StSmall',1e-3,'NonNeg',1,'RandSeed',randkey,'Landscape',land);

% random initial distrubution of species biomass
Vs = rand(size(locr));

end