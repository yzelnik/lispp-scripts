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

% which points in the I-D parameter space to choose to look at?
spcpnts = [1 3 5 1 1; 1 1 1 3 5];

% define 9 parameters related to species properties and interactions
% [1-4: niche-cent-span niche-cent-offset niche-wdth-span niche-wdth-offset ...]
% [5-9: self-reg-strength dispersal-strength interaction-distance interaction-cutoff interaction-normalized-share]
spsprm=[0.6 0.2 0.02 0.2  1 0.1 1 0.02 0.9];

maxtime=1000;  % maximum simulation time
ssthresh=1e-5; % threshold for reaching steady-state
randrep=20;    % number of different randomizations to simulate

% NOTE: With the above parameters (used in the main results), the
% simulations may take an EXTREMELY long time to run, possibly weeks.
% un-commenting the defintions below will lead to more a reasonable run-time
%randrep=3;  % only 3 reps instead of 20
%spcpnts=[1 2 3 1 1; 1 1 1 2 3]; % not using high D or high I
%syssz = 200; % number of pixel along the x and y axes

%% run simulations

for rndind=1:randrep % go over randomization keys
  for pntind=1:size(spcpnts,2) % go over different parameter sets
    disp([rndind pntind])
    spsprm(6:7)=[diffusion(spcpnts(1,pntind)) interdist(spcpnts(2,pntind))];
    % create community (see more explanations inside the startintcomm script file)
    [Vs,Ps,Es]=startintcomm(rndind,spnum,syssz,landprm,spsprm);
    
    % run simulation (change Es.OlDraw to 1 to see the simulation white it is running)
    outs{rndind,pntind}  = run2ss(Vs,Ps,Es,'Es.OlDraw',0,'Es.PlotFunc',@plot3sp,'Es.SsThresh',ssthresh,'Es.TimeMax',maxtime,'Es.TsSize',timesteps(spcpnts(1,pntind)));
  end;
end;


%% calculate SAR and BEF curves
thresh=1e-3;
befprm=[1 100]; % numer of pixels for BEF calculations (1 for local, 100 for regional)
befbins=1:Ps.VarNum; % for BEF calculations, consider the different number of species
sarprms=[100 40]; % for creating the SAR curves: number of repetitions per size, number of sizes to run through

for ii=1:size(outs,1)
  for jj=1:size(outs,2)
    if(~isempty(outs{ii,jj}))
      % calculate the SAR curve
      sarres=T_SAR(outs{ii,jj},Ps,Es,'Es.SarPrm',sarprms,'Es.SarSurviveThresh',thresh);
      % calculate the BEF curve
      [befres,extra]=T_BEF(outs{ii,jj},Ps,Es,'Es.BefBins',befbins,'Es.BefSurviveThresh',thresh,'Es.BefScale',befprm);
      allsar(:,jj,ii)=sarres(:,2);
      allbef(:,jj,ii,:)=[befres(:,2:3) extra];
    end;
  end;
end;

%% plot curves out
showchs = 1:5; % which parameter sets to show
ylms=[0.6 1.3; 0.001 2.25; 0.001 2.25];
xlms=[0 4.75; 0 18; 0 18];
ttls={'SAR curve','local BEF','regional BEF'};
xlbs={'log10 area','species richness','species richness'};
ylbs={'log10 species richness','community biomass','community biomass'};
xjmp=[1 5 5];
yjmp=[0.2 1 1];

clrs='kbcmrg';
sharethresh=0.005; % what fraction of the landscape should have this value of biodiversity to be shown in the BEF curve?
mintimes = 2; % how many randomization sets need to have this value of biodiversity to be used for the average of BEF?
tmp=(allbef(:,:,:,1:2).*(allbef(:,:,:,5:6)>sharethresh)); % consider only values where there is over a fraction of sharethresh
tmp=tmp.*repmat(sum(tmp>0,3)>=mintimes,[1 1 size(allbef,3)]);
tmp(isnan(tmp))=0; 

clf; ha = tight_subplot(1,3,0.07,[0.13 0.08],[0.06 0.01]);
for plotind=1:3
    axes(ha(plotind))
  if(plotind==1)
    xax=log10(sarres(:,1));
  else
    xax=befres(:,1);
  end;

  hold on;
  for ind=1:length(showchs)
    % setup the data for the plot
    if(plotind>1)
	  data=squeeze(tmp(:,showchs(ind),:,plotind-1));
    else
	  data=log10(squeeze(allsar(:,showchs(ind),:)));
    end;
    % plot the average values
    plot(xax,sum(data,2)./sum(data>0,2),clrs(ind),'lineWidth',2); % deals with empty values
  end;
  
  for ind=1:length(showchs)
    % and now plot all the different randomizations
    if(plotind>1)
        data=squeeze(allbef(:,showchs(ind),:,plotind-1));
    else
        data=log10(squeeze(allsar(:,showchs(ind),:)));
    end;
    plot(xax,data,[':' clrs(ind)])
    
    xlim(xlms(plotind,:))
    ylim(ylms(plotind,:))
    set(gca,'xTick',0:xjmp(plotind):100,'yTick',0:yjmp(plotind):100,'xTickLabel',0:xjmp(plotind):100,'yTickLabel',0:yjmp(plotind):100)
  end;

  hold off;
  title(ttls{plotind},'fontSize',20)
  xlabel(xlbs{plotind},'fontSize',20)
  ylabel(ylbs{plotind},'fontSize',20)
end;
axes(ha(2))

%% plot biomass, diversity and three-species
thresh=1e-3; % biomass threshold for extinction
showchs=1:5; % which parameter sets to show
rndchs=1;    % which randomization to show

ttls={}; 
for ii=1:size(spcpnts,2) ttls{ii}=sprintf('D=%.1f, I=%.1f',0.1*round(10*sqrt(diffusion(spcpnts(1,ii)))*10),0.1*round(10*(interdist(spcpnts(2,ii))))); end;
ylbs={'total biomass','species richness','three species'};
panelen=length(showchs);
clf; ha = tight_subplot(3,panelen,[0.02 0.01],[0.01 0.05],[0.03 0.08]);

for ii=1:panelen
    % plot biomass
    axes(ha(ii))
    plotst(sum(outs{rndchs,showchs(ii)},2),Ps,Es,'Es.PlotBare',1,'Es.St2Colorbar',0);
    title(ttls{ii},'fontSize',20)
    caxis([0 3])
    axis square;
    
    % plot biodiversity
    axes(ha(ii+panelen))
    plotst(sum(outs{rndchs,showchs(ii)}>thresh,2),Ps,Es,'Es.PlotBare',1,'Es.St2Colorbar',0);
    caxis([0 18])
    axis square;
    
    % plot biomass distribution of 3 species
    whichspcs = [9 6 2]; % which 3 species to plot out?
    axes(ha(ii+panelen*2))
    plot3sp(outs{rndchs,showchs(ii)}(:,whichspcs),Ps,Es);
    set(gca,'xTick',[],'yTick',[]);
    axis square;
end;
colormap jet;
for ii=1:3
    axes(ha((ii-1)*panelen+1));
    ylabel(ylbs{ii},'fontSize',20)
end;

axes(ha(length(showchs)))
colorbar('location','manual','position',[0.93 0.665 0.035 0.27],'xTick',0:1:20,'fontSize',12)

axes(ha(length(showchs)*2))
colorbar('location','manual','position',[0.93 0.345 0.035 0.27],'xTick',0:3:20,'fontSize',12)

basebar=(64:-1:1)'/64;
tmpcb=axes('Units','normalized', 'Position',[0.93 0.025 0.035 0.27]);
imagesc(cat(3,[basebar;basebar*0;basebar*0],[basebar*0;basebar;basebar*0],[basebar*0;basebar*0;basebar]));
set(tmpcb,'xTick',[],'yTick',[]);
text(1.8,35,'S_1','fontSize',14)
text(1.8,100,'S_2','fontSize',14)
text(1.8,166,'S_3','fontSize',14)

%% spatial correlations
sptest=[1 3 5]; % which parameter sets to show
rndchs=1; % which randomization to show?

% go through randomizations and parameter setrs
for ii=1:randrep
    for jj=1:size(spcpnts,2)
      if(~isempty(outs{ii,jj}))
        for spind=1:spnum  % go through each species
            % calcualte spatial correlation
            mcor(:,spind,ii,jj)=spacorr2d(reshape(outs{ii,jj}(:,spind),syssz,syssz)-mean(outs{ii,jj}(:,spind))); 
            tmval=(max(mcor(:,spind,ii,jj))+min(mcor(:,spind,ii,jj)))/2; 
            [~,tbind]=min(abs(mcor(:,spind,ii,jj)-tmval)); 
            matmids(spind,ii,jj)=tbind; 
        end; 
      end;
    end; 
end;

ttls={}; 
for ii=1:length(sptest) ttls{ii}=sprintf('D=%d, I=%d',round(sqrt(diffusion(spcpnts(1,sptest(ii))))*10),round(interdist(spcpnts(2,sptest(ii))))); end;

clf; ha = tight_subplot(1,3,0.07,[0.15 0.08],[0.07 0.01]);
for ii=1:3
    axes(ha(ii))
    plot(mcor(:,:,rndchs,sptest(ii)))
    hold on;
    plot(mean(mcor(:,:,rndchs,sptest(ii)),2),'k','lineWidth',2);
	hold off;
    title(ttls{ii},'fontSize',20)
    ylim([-0.4 1])
    mm=mean(matmids(:,rndchs,sptest(ii)));
    hold on;
    plot([mm mm],[-0.4 1],'--','color',[0.5 0.5 0.5],'lineWidth',2);
    hold off;
    text(mm+15,0.7,sprintf('$X = %.1f$',mm),'fontSize',20,'color',[0.5 0.5 0.5],'interpreter','latex');
    text(mm+15,0.85,sprintf('$\\bar{X} = %.1f$',mean(mean(matmids(:,:,sptest(ii))))),'fontSize',20,'color',[0.5 0.5 0.5],'interpreter','latex');
end;

axes(ha(2))
xlabel('distance from origin','fontSize',20);
axes(ha(1))
ylabel('spatial correlation','fontSize',20);

