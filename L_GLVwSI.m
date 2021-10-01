function VsOut=L_GLVwSI(Vs,Ps,Es)
% Generalized Lotka Volterra equations with Spatial Interactions
% VsOut=L_GLVwSI(Vs,Ps,Es)
% equation per variable N_i is: dN_i/dt = r_i*N_i + SUM{A_ij*N_j}N_i + d_i* D^2(N_i)
% Parameters are the growth,interaction-distance, and diffusion vectors (r,id,d) and interaction matrix A.
% For example, for a two-species symmetric competition:
% r=[0.1,0.1], A=-[1 2;2 1]; d=[1,2];
% Where each row in A corresponds to all the effects enacted on a given species
% Note that r can be a scalar, vector or matrix (to allow different heterogenities)

len=size(Vs,1);

if(isfield(Es,'SetupMode') && Es.SetupMode) % setup mode
    kerwdth=Ps.id(1);
    if(length(Ps.id)>1)
        cutoff=Ps.id(2);
    else
        cutoff=Es.StSmall;
    end;
    % create a grid to calculate the distance from the center
    [xx,yy] = meshgrid(-Ps.Nx/2:Ps.Nx/2-1,-Ps.Ny/2:Ps.Ny/2-1);
    kernel  = exp(-0.5*(xx.^2+yy.^2)/kerwdth);
    kernel(kernel<cutoff)=0; % set tail of 2d kernel to zero (below threshold)
    submat  = (kernel(find(kernel(Ps.Nx/2+1,:)>0),find(kernel(:,Ps.Ny/2+1)>0)));
    if(length(Ps.id)>2) && (Ps.id(3)>0) % do we have a mixture of local and spatial interactions? 
        totint = sum(submat(:)); % sum up the kernel for normalization
        submat = submat*Ps.id(3)/totint; % normalize
        cent   = round((length(submat)-1)/2)+1; % find the center
        submat(cent,cent)=submat(cent,cent)+(1-Ps.id(3)); % add the local-only interactions
    end;
    Ps.IntSM= StencilToSM(submat,Ps.Nx,Ps.Ny,0); % calculate the spatial matrix
    VsOut   = Ps;

else
  if(Es.JacMode==0)      % Model equations
    baseint=Ps.IntSM*Vs;
    if(size(Ps.r,1)==1) && (size(Ps.r,2)>1)
        VsOut = Vs.*(repmat(Ps.r,Ps.Nx,1) + baseint*(Ps.A'));
    else % Ps.r is either scalar or a matrix
        VsOut = Vs.*(Ps.r + baseint*(Ps.A'));
    end;
  else             % Jacobian of equations
    jac = zeros(len*Ps.VarNum,Ps.VarNum);
    if(length(Ps.r)==1)
        Ps.r = repmat(Ps.r,len,Ps.VarNum);
    elseif(min(size(Ps.r))==1)
        Ps.r = repmat(Ps.r(:)',len,1);
    end;
    for ii=1:Ps.VarNum
        jac((1:len)+(ii-1)*len,ii)=Ps.r(:,ii)+Vs*Ps.A(ii,:)';
        for jj=1:Ps.VarNum
            jac((1:len)+(ii-1)*len,jj)=jac((1:len)+(ii-1)*len,jj)+Ps.A(ii,jj)*Vs(:,ii);
        end;
    end;
    % write it in a large sparse matrix format 
    VsOut = ArrangeJacobian(jac,Ps,Es);
  end;
end;

end
