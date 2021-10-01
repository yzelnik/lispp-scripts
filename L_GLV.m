function VsOut=L_GLV(Vs,Ps,Es)
% Generalized Lotka Volterra equations - Local terms
% VsOut=L_GLV(Vs,Ps,Es)
% equation per variable N_i is: dN_i/dt = r_i*N_i + SUM{A_ij*N_j}N_i + d_i* D^2(N_i)
% Parameters are the growth and diffusion vectors (r,d) and interaction matrix A.
% For example, for a two-species symmetric competition:
% r=[0.1,0.1], A=-[1 2;2 1]; d=[1,2];
% Where each row in A corresponds to all the effects enacted on a given species
% Note that r can be a scalar, vector or matrix (to allow different heterogenities)

len=size(Vs,1);

if(Es.JacMode==0)      % Model equations
    if(size(Ps.r,1)==1) && (size(Ps.r,2)>1)
        VsOut = Vs.*(repmat(Ps.r,Ps.Nx,1) + Vs*(Ps.A'));
    else % Ps.r is either scalar or a matrix
        VsOut = Vs.*(Ps.r + Vs*(Ps.A'));
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

end
