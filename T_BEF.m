function [res,morestats,allpoints]=T_BEF(Vs,Ps,Es,varargin)
% Calculate a BEF (Biodiversity-Ecosystem Function) curve
% [res,morestats,allpoints]=T_BEF(Vs,Ps,Es)
% - res returns avg (binned) points for biodiversity vs. function 
%   (one column for biodiversity, one or more columns for function, if more than one scale/input)
% - morestats returns also the std and relative ratio of these points (in sets of two columns)
% - allpoints returns unbinned data (two columns)
% Es.BefDiversityInput can define range of variables to calculate diversity on, or a function that does this
% Es.BefFunctionInput can define range of variables to calculate function on (via biomass), or a function that does this
% It is possible to gbinsive both of these as cell arrays, to define more than one option
% Es.BefSurviveThresh>0 gives the threshold for diversity measurements
% (number of species), or if <0 a threshold used for Shannon diversity
% Es.BefScale gives one (or more) scale(s) at which to measure BEF 
% (0/1 means every point seperately, otherwise uses relative size or pixel number)
% Es.BefBins can specify the binning (for diversity axis), 
% which is required for more than one scale/input ( Es.BefScale>0 or Es.BefDiversityInput>0)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Put in some default values of Es
Es=InsertDefaultValues(Es,'BefDiversityInput',[],'BefFunctionInput',[],'BefSurviveThresh',[],'BefScale',0,'BefBins',[]);

% As default of both diversity and function, just use all variables
if(isempty(Es.BefDiversityInput))
    Es.BefDiversityInput=1:Ps.VarNum;
end;
if(isempty(Es.BefFunctionInput))
    Es.BefFunctionInput=1:Ps.VarNum;
end;

% At which scale(s) to measure BEF?
len=Ps.Nx*Ps.Ny;
for ii=1:length(Es.BefScale)
    if(Es.BefScale(ii)<1) % if given at values of 0 to 1 (ratios)
        Es.BefScale(ii)=ceil(Es.BefScale(ii)*len+1e-20);
    end;
    
    if(Es.BefScale(ii)>1) % If not in local scale
        % Setup a nearest-neighbors spatial matrix for later use
        Ps.NnSm=NeighborSM(1,Ps,Es);
        sampnum = ceil(len/Es.BefScale(ii)); % cover approx. the size of the system
        rndlocs = randi(len,1,sampnum);
        neigregs{ii}=sparse(sampnum,len);
        for jj=1:sampnum
            neigregs{ii}(jj,:)=(1/Es.BefScale(ii))*FindLocalRegion([Es.BefScale(ii) rndlocs(jj)],Ps,Es);
        end;
    else
        neigregs{ii}=[];
    end;
    
end;

% Wrap in cell-arrays 
if(~iscell(Es.BefDiversityInput)||~iscell(Es.BefFunctionInput))
    if(~iscell(Es.BefDiversityInput)&&~iscell(Es.BefFunctionInput))
        Es.BefDiversityInput={Es.BefDiversityInput};
        Es.BefFunctionInput={Es.BefFunctionInput};
    else
        error('Bef Inputs sholuld be both cell arrays (or none at all)');
    end;
end;

% Try to avoid some strage things
if(((length(Es.BefDiversityInput)>1) || length(Es.BefScale)>1) && isempty(Es.BefBins))
    error('If Es.BefBins is not defined, only one input and one scale can be used.');
end;

% Go over different inputs and different scales for BEF measurements
for inputind=1:length(Es.BefDiversityInput)
  for scaleind=1:length(Es.BefScale)
    divinput = Es.BefDiversityInput{inputind};
    funinput = Es.BefFunctionInput{inputind};
    Ps2=Ps;
    if(isempty(neigregs{scaleind}))
        state=Vs;  
    else
        state=neigregs{scaleind}*Vs;
        Ps2.Nx=size(state,1);
        Ps2.Ny=1;
    end;
    
    % Get divertisy vector (x-axis)
    if(isa(divinput,'function_handle'))
      xax = divinput(state,Ps2,Es);
    else
      if(min(divinput)<0)
        % We allow negative values to mean that the diversity is preselected per site (input diversity)
        if(~isempty(neigregs{scaleind}))
            xax = abs(divinput);
        else
            xax = neigregs{ii}*abs(divinput);
        end;
      else
        if(isempty(Es.BefSurviveThresh))
            Es.BefSurviveThresh=Es.StSmall/10;
        end;
        if(Es.BefSurviveThresh>0)
            % Normal case of species richness, Es.BefDiversityInput gives location of variables to measure diversity on
            xax = sum(state(:,divinput)>Es.BefSurviveThresh,2);
        elseif(Es.BefSurviveThresh<0) 
            % Shannon diversity, Es.BefDiversityInput gives location of variables to measure diversity on
            normconst=sum(state(:,divinput),2);
            normconst(normconst<abs(Es.BefSurviveThresh))=inf;
            relbiomass = state(:,divinput)./repmat(normconst,1,length(divinput));
            xax = -sum(relbiomass.*log(relbiomass),2);     
        else
            error('Es.BefSurviveThresh Not properly defined');
        end;
      end;
    end;

    % Get function vector (y-axis)
    if(isa(funinput,'function_handle'))
      yax = funinput(state,Ps2,Es);
    else
      % Default is just the biomass of a range of variables
      yax = sum(state(:,funinput),2);
    end;

    % Calcualte bins
    stats=makebins(xax,yax,Es.BefBins);

    if(inputind==1 && scaleind==1)
        res(:,1)=stats(:,1); % x-axis (biodiversity axis)
        allpoints = [xax(:) yax(:)]; % Save all data points, but just for the first instance of scale&input
    end;
    % For each instance of scale and input, save avg data, and also std&relative ratio
    res(:,1+scaleind+(inputind-1)*length(Es.BefScale)) = stats(:,2);
    morestats(:,(scaleind+(inputind-1)*length(Es.BefScale))*2+[-1 0]) = stats(:,[3 6]);
    
  end;
end;

end