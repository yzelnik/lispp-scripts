function res=T_SAR(Vs,Ps,Es,varargin)
% Calculate a SAR (Species-Area Relationship) curve
% res=T_SAR(Vs,Ps,Es)
% Es.SarPrm = [maxrep pntnum] (default is [100 20])
% where maxrep gives the maximum number of tests per area checked
% and pntnum gives the number of different areas (curve resolution)
% res returns a pntnum x 2 matrix of areas and their species-count

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Put in some default values of Es
Es=InsertDefaultValues(Es,'SarDiversityInput',[],'SarSurviveThresh',[],'SarPrm',[100 20]);

% Setup a nearst-neighbors spatial matrix for later use
Ps.NnSm=NeighborSM(1,Ps,Es);

% As default of diversity vector, just use all variables
if(isempty(Es.SarDiversityInput))
    Es.SarDiversityInput=1:Ps.VarNum;
end;
if(isempty(Es.SarSurviveThresh))
	Es.SarSurviveThresh=Es.StSmall/10;
end;

maxrep=Es.SarPrm(1);
pntnum=Es.SarPrm(2);

len = Ps.Nx*Ps.Ny;

% Define size-axis
xax=[1 unique(round(exp((log(2):(log(len)-log(2))/(pntnum-1):log(len)))))];
repnum =  min(maxrep,ceil(len./xax)*2);

% Calculate the first two resolutions (1 and 2 pixels)
yax = zeros(size(xax));
yax(1)=calcavgdiversity(Vs(:,Es.SarDiversityInput),Es.SarSurviveThresh);
yax(2)=calcavgdiversity(cat(3,Vs(1:2:end,Es.SarDiversityInput),Vs(2:2:end,Es.SarDiversityInput)),Es.SarSurviveThresh);
% Calculate all other reg-sizes besides the first two
for ind=3:length(xax)
    rndlocs = randi(len,1,repnum(ind));
	spvects = zeros(repnum(ind),Ps.VarNum,xax(ind));
    % Run over different choices of reg for a given size
	for kk=1:repnum(ind)
        reg=FindLocalRegion([xax(ind) rndlocs(kk)],Ps,Es);
        spvects(kk,:,:)=Vs(reg,Es.SarDiversityInput)';
        %divlist(kk)=sum(sum(Vs(reg,Es.SarDiversityInput)>Es.SarSurviveThresh,1)>0);
	end;
	yax(ind) = calcavgdiversity(spvects,Es.SarSurviveThresh);
end;

res = [xax(:) yax(:)];

end

function divval=calcavgdiversity(spbiomass,thresh)

if(thresh>0)
    % Default - species count over threshold
    divval=mean(sum(max(spbiomass,[],3)>thresh,2));
elseif(thresh<=0) 
    % Shannon diversity, Es.BefDiversityInput gives location of variables to measure diversity on
    spbiomass  = sum(spbiomass,3);
    sumbiomass = sum(spbiomass,2);
    goodsites  = sumbiomass>abs(thresh);
    
    relbiomass = spbiomass(goodsites)./repmat(sumbiomass(goodsites),1,size(spbiomass,2));
    divval = mean(goodsites)*mean(-sum(relbiomass.*log(relbiomass),2));     
else
	error('Es.SarSurviveThresh not properly defined');
end;

end