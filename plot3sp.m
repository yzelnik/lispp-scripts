function plot3sp(Vs,Ps,Es)
for ii=1:3
tmp(:,:,ii) = reshape(Vs(:,ii),Ps.Ny,Ps.Nx)';
end;
imagesc(tmp*10);
axis xy;
%permute(reshape(Vs(:,1:3),[Ps.Ny Ps.Nx 3])*5,[2 1 3]))
end