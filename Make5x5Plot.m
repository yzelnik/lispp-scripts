% define basic parameters

alllands = [0.95 0.04; 0.5 0.1]; % E values of: 32, 10, respectively.
spnum = 20;  % number of species in the system
syssz = 320; % number of pixel along the x and y axes
landind=1;   % either 1 or 2. choose 1 for E=32, or 2 for E=10

landprm=alllands(landind,:);  % choose value for E
% 5 values for I and D
interdist=[1.0 3.2 10.0 32.0 100.0]; 
diffusion=[0.01 0.1 1.0 10.0 100.0];
% time-step size values corresponding to the diffusion coefficients (i.e. that do not lead to divergence)
timesteps=[0.1;0.1;0.05;0.02;0.002];

% define 9 parameters related to species properties and interactions
% [1-4: niche-cent-span niche-cent-offset niche-wdth-span niche-wdth-offset ...]
% [5-9: self-reg-strength dispersal-strength interaction-distance interaction-cutoff interaction-normalized-share]
spsprm=[0.6 0.2 0.02 0.2  1 0.1 1 0.02 0.9];

maxtime=1000;  % maximum simulation time
ssthresh=1e-5; % threshold for reaching steady-state
rndind=1;      % randomization key

% NOTE: With the above parameters (used in the main results), the
% simulations take an EXTREMELY long time to run, possibly weeks.
% un-commenting the defintions below will lead to more a reasonable run-time
%syssz   = 120;
%maxtime = 100;
%ssthresh= 1e-3;
%interdist=[1.0 3.2 10.0];

%% run simulations for different values of I and D
for ind1=1:length(diffusion)
  for ind2=1:length(interdist)
    disp([ind1 ind2])
    
    % change values of D and I according to the two indices
    spsprm(6:7)=[diffusion(ind1) interdist(ind2)]; 
    % create community (see more explanations inside the startintcomm script file)
    [Vs,Ps,Es]=startintcomm(rndind,spnum,syssz,landprm,spsprm); 
    
    % run simulation (change Es.OlDraw to 1 to see the simulation white it is running)
    outs{ind1,ind2}  = run2ss(Vs,Ps,Es,'Es.OlDraw',0,'Es.PlotFunc',@plot3sp,'Es.SsThresh',ssthresh,'Es.TimeMax',maxtime,'Es.TsSize',timesteps(ind1));
  end;
end;

%% plot out spatial profiles of states
plotmode=0; % 1 for biomass, 0 for species number
thresh=1e-3; % how much biomass to be considered an extant species

% define values of I and D
dvals = 10*sqrt(diffusion);
ivals = interdist;

% setup the plot structure
clf;
ha = tight_subplot(length(interdist),length(diffusion),[0.01 0.01],[0.01 0.12],[0.12 0.12]);

% go through all simulation results
for ii=1:length(diffusion)
  for jj=1:length(interdist)
    axes(ha(ii+(jj-1)*length(diffusion)))
    % plot out the state of the system
	if(plotmode) % biomass
      plotst(sum(outs{ii,jj},2),Ps,Es,'Es.PlotBare',1,'Es.St2Colorbar',0);
    else  % biodiversity
      plotst(sum(outs{ii,jj}>thresh,2),Ps,Es,'Es.PlotBare',1,'Es.St2Colorbar',0);
    end;
    caxis([0 20-17*plotmode]); % set the color spectrum to max of 3 or 20, depending on what is plotted
    axis square;
    % write out the I and D values 
    if(jj==1) xlabel(sprintf('D=%.1f',dvals(ii)),'fontSize',20); set(gca,'XAxisLocation','top'); end;
    if(ii==1) ylabel(sprintf('I=%.1f',ivals(jj)),'fontSize',20); end;
  end;
end;
colormap jet;

axes(ha(1))
temp = annotation('arrow', [0.2 0.88], [0.94 0.94],'color',[0 0 0],'lineWidth',10,'headWidth',35,'headLength',75);
text(450,520,'stronger dispersal','fontSize',28);
temp = annotation('arrow', [0.06 0.06], [0.8 0.12],'color',[0 0 0],'lineWidth',10,'headWidth',35,'headLength',75);
text(-210,-1000,'longer-distance interaction','fontSize',28,'rotation',90);

colorbar('location','manual','position',[0.9 0.01 0.035 0.87],'xTick',0:3:20,'fontSize',12)

%% compare to baseline (no spatial interactions or dispersal)

baseout = run2ss(Vs,Ps,Es,'Es.SsThresh',ssthresh,'Es.TimeMax',maxtime,'Es.TsSize',timesteps(1),'Ps.LocFunc',@L_GLV,'Ps.Ds',Ps.Ds*0);

for ii=1:length(diffusion)
  for jj=1:length(interdist)
    if(~isempty(outs{ii,jj}))
        difmat(ii,jj)=mean(sum((outs{ii,jj}-baseout).^2,2));
    end;
  end;
end;

% show overall difference of different states (under different values of I and D) compared to the baseline case
clf;
imagesc(difmat')
colorbar

%% other summary plots
clf;
thresh=1e-3;

% make calculations of averages and sums
for ii=1:length(diffusion)
  for jj=1:length(interdist)
    if(~isempty(outs{ii,jj}))
        sumbio(ii,jj)=mean(sum(outs{ii,jj},2));
        avgdiv(ii,jj)=mean(sum((outs{ii,jj}>thresh),2));
        totdiv(ii,jj)=max(sum((outs{ii,jj}>thresh),2));
    end;
  end;
end;


ha = tight_subplot(1,3,0.04,[0.06 0.05],[0.05 0.02]);

% plot average biomass per parameter set
axes(ha(1))
imagesc(sumbio');
set(gca,'xTick',[],'yTick',[])
caxis([1.4 1.8]); colorbar;
title('average community biomass','fontSize',14)

% plot average biodiversity per parameter set
axes(ha(2))
imagesc(avgdiv');
set(gca,'xTick',[],'yTick',[])
caxis([0 13]); colorbar;
title('average species richness','fontSize',14)

% plot total biodiversity per parameter set
axes(ha(3))
imagesc(totdiv');
set(gca,'xTick',[],'yTick',[])
caxis([0 18]); colorbar;
title('total species richness','fontSize',14)

