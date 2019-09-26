close all; 
warning off;
%_____________________load mooring data______________________
if 1; %{{{
clear all;

   % path to m_map here
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/seawater/'));
   load cmap;

%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');

% load bg state
   load('../data/Tsort.mat');

% load bottom GusTs
  load('../data/G_bbl.mat')
end %}}}

if 1 


% define grid 
xl = M(1).T.time([1 end]);
TG.time = xl(1):1:xl(2);
%TG.T    = 6:.1:19;
TG.T    = 9:.1:19;


%__________________ sumup INTERIOR chipods
TG.N_int   = zeros(length(TG.T), length(TG.time));
TG.chi_int = zeros(length(TG.T), length(TG.time));
TG.N_abbl  = zeros(length(TG.T), length(TG.time));
TG.chi_abbl = zeros(length(TG.T), length(TG.time));
   TG.Kt_abbl = zeros(length(TG.T), length(TG.time));
   TG.eps_abbl = zeros(length(TG.T), length(TG.time));
TG.N_surf    = zeros(length(TG.T), length(TG.time));
TG.chi_surf = zeros(length(TG.T), length(TG.time));
TG.N_100  = zeros(length(TG.T), length(TG.time));
TG.chi_100 = zeros(length(TG.T), length(TG.time));
TG.N_50    = zeros(length(TG.T), length(TG.time));
TG.chi_50 = zeros(length(TG.T), length(TG.time));
TG.N_30    = zeros(length(TG.T), length(TG.time));
TG.chi_30 = zeros(length(TG.T), length(TG.time));
for m  = 1:length(M)
   for g = 1:sum(M(m).G.mab>1)  
      M(m).G.ii_timegrid(g,:) = round(interp1( TG.time, 1:length(TG.time),  M(m).G.time, 'nearest'));
      M(m).G.ii_Tgrid(g,:)    = round(interp1( TG.T, 1:length(TG.T),  M(m).G.T(g,:), 'nearest'))  ;

      %%count for grid
      for t = find(~isnan(M(m).G.ii_Tgrid(g,:)))
         if ~isnan(M(m).G.ii_timegrid(g,t)) & ~isnan(M(m).G.chi(g,t))
           TG.N_int( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_int( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;
           TG.chi_int( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_int( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);

           % seperately also count whats just above the bbl
           if M(m).G.mab(g)<8
              TG.N_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;
              TG.chi_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);
              TG.Kt_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.Kt_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.Kt(g,t);
              TG.eps_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.eps_abbl( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.eps(g,t);
           end
           % whats close to the surface >17
           if (M(m).G.mab(g)-M(m).waterdepth) >-17
              TG.N_surf( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_surf( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;
              TG.chi_surf( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_surf( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);
           end

           % separte differnt water depths
           if M(m).waterdepth>60
              TG.N_100( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_100( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;
              TG.chi_100( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_100( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);
           elseif M(m).waterdepth>35
              TG.N_50( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_50( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;;
              TG.chi_50( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_50( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);
           else
              TG.N_30( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.N_30( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+1;;
              TG.chi_30( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t)) = TG.chi_30( M(m).G.ii_Tgrid(g,t), M(m).G.ii_timegrid(g,t))+M(m).G.chi(g,t);
           end
         end
      end

      %plot(ax(a), M(m).G.time, M(m).G.T(g,:), '.', 'Linewidth', 1);
   end
end

% __________________sumup BBL chipods
TG.N_bbl   = zeros(length(TG.T), length(TG.time));
TG.chi_bbl = zeros(length(TG.T), length(TG.time));
TG.eps_bbl = zeros(length(TG.T), length(TG.time));
TG.Kt_bbl = zeros(length(TG.T), length(TG.time));
   for g = 1:length(G_bbl)  
      ii_timegrid = round(interp1( TG.time, 1:length(TG.time),  G_bbl(g).time, 'nearest'));
      ii_Tgrid    = round(interp1( TG.T, 1:length(TG.T),  G_bbl(g).T, 'nearest'))  ;

      %%count for grid
      for t = find(~isnan(ii_Tgrid))
         if ~isnan(ii_timegrid(t)) & ~isnan(G_bbl(g).chi(t))
           TG.N_bbl( ii_Tgrid(t), ii_timegrid(t)) = TG.N_bbl( ii_Tgrid(t), ii_timegrid(t))+1;
           TG.chi_bbl( ii_Tgrid(t), ii_timegrid(t)) = TG.chi_bbl( ii_Tgrid(t), ii_timegrid(t))+G_bbl(g).chi(t);
           TG.eps_bbl( ii_Tgrid(t), ii_timegrid(t)) = TG.eps_bbl( ii_Tgrid(t), ii_timegrid(t))+G_bbl(g).eps(t);
           TG.Kt_bbl( ii_Tgrid(t), ii_timegrid(t)) = TG.Kt_bbl( ii_Tgrid(t), ii_timegrid(t)) ...
                     + sqrt(2.3e-3)*abs(G_bbl(g).spd(t))*.4*1.25;
         end
      end

   end
% NORMALIZATION OF SUM to get Average
TG.chi = TG.chi_int + TG.chi_bbl;

TG.chi_int( TG.chi_int==0 ) = nan;
TG.chi_int = TG.chi_int./TG.N_int;

TG.chi_bbl( TG.chi_bbl==0 ) = nan;
TG.chi_bbl = TG.chi_bbl./TG.N_bbl;
   TG.Kt_bbl( TG.chi_bbl==0 ) = nan;
   TG.Kt_bbl = TG.Kt_bbl./TG.N_bbl;
   TG.eps_bbl( TG.chi_bbl==0 ) = nan;
   TG.eps_bbl = TG.eps_bbl./TG.N_bbl;

TG.chi( TG.chi==0 ) = nan;
TG.N   = TG.N_int+TG.N_bbl;
TG.chi = (TG.chi)./TG.N;


TG.chi_abbl( TG.chi_abbl==0 ) = nan;
TG.chi_abbl = TG.chi_abbl./TG.N_abbl;
   TG.Kt_abbl( TG.chi_abbl==0 ) = nan;
   TG.Kt_abbl = TG.Kt_abbl./TG.N_abbl;
   TG.eps_abbl( TG.chi_abbl==0 ) = nan;
   TG.eps_abbl = TG.eps_abbl./TG.N_abbl;
TG.chi_surf( TG.chi_surf==0 ) = nan;
TG.chi_surf = TG.chi_surf./TG.N_surf;
TG.chi_100( TG.chi_100==0 ) = nan;
TG.chi_100 = TG.chi_100./TG.N_100;
TG.chi_50( TG.chi_50==0 ) = nan;
TG.chi_50 = TG.chi_50./TG.N_50;
TG.chi_30( TG.chi_30==0 ) = nan;
TG.chi_30 = TG.chi_30./TG.N_30;
end


%_____________________grid background state to chi grid______________________
   % interpolate each teim step on common T-grid;
   tmpTz = nan( length(TG.T), length(TS.time));
   tmpzT = nan( length(TG.T), length(TS.time));
   tmpzs = nan( length(TG.T), length(TS.time));
   for t = 1:length(TS.time)
      tmpTz(:,t) = interp1( TS.T_dzdT, (TS.dzdT(:,t)).^(-1), TG.T);
      tmpzT(:,t) = interp1( TS.T_dzdT, (TS.dzdT(:,t)), TG.T);
      tmpzs(:,t) = interp1( TS.T_vec, (TS.zs(:,t)), TG.T);
   end
   % now interpolate on common time grid
   TG.dTdz = nan(size(TG.chi));
   TG.dzdT = nan(size(TG.chi));
   TG.zs   = nan(size(TG.chi));
   for l = 1:length(TG.T)
      TG.dTdz(l,:)= clever_interp( TS.time, tmpTz(l,:), TG.time);
      TG.dzdT(l,:)= clever_interp( TS.time, tmpzT(l,:), TG.time);
      TG.zs(l,:)= clever_interp( TS.time, tmpzs(l,:), TG.time);
   end
   
   
   % maybe this should be changed 
   TG.chi(TG.N<200) = nan;

   TG.Phi_d = .5*TG.chi.*TG.dzdT;
   TG.Keff  = .5*TG.chi.*TG.dzdT.^2;
   TG.Jq    = 1025*4200*TG.Phi_d;

   TG.Phi_d_int = .5*TG.chi_int.*TG.dzdT;
   TG.Keff_int  = .5*TG.chi_int.*TG.dzdT.^2;
   TG.Jq_int    = 1025*4200*TG.Phi_d_int;
 
   TG.Phi_d_bbl = .5*TG.chi_bbl.*TG.dzdT;
   TG.Keff_bbl  = .5*TG.chi_bbl.*TG.dzdT.^2;
   TG.Jq_bbl    = 1025*4200*TG.Phi_d_bbl;

   TG.Phi_d_abbl = .5*TG.chi_abbl.*TG.dzdT;
   TG.Keff_abbl  = .5*TG.chi_abbl.*TG.dzdT.^2;
   TG.Jq_abbl    = 1025*4200*TG.Phi_d_abbl;

   TG.Phi_d_surf = .5*TG.chi_surf.*TG.dzdT;
   TG.Keff_surf  = .5*TG.chi_surf.*TG.dzdT.^2;
   TG.Jq_surf    = 1025*4200*TG.Phi_d_surf;

   TG.Phi_d_100 = .5*TG.chi_100.*TG.dzdT;
   TG.Keff_100  = .5*TG.chi_100.*TG.dzdT.^2;
   TG.Jq_100    = 1025*4200*TG.Phi_d_100;

   TG.Phi_d_50 = .5*TG.chi_50.*TG.dzdT;
   TG.Keff_50  = .5*TG.chi_50.*TG.dzdT.^2;
   TG.Jq_50    = 1025*4200*TG.Phi_d_50;

   TG.Phi_d_30 = .5*TG.chi_30.*TG.dzdT;
   TG.Keff_30  = .5*TG.chi_30.*TG.dzdT.^2;
   TG.Jq_30    = 1025*4200*TG.Phi_d_30;

   save('../data/Turb_Tspace.mat', 'TG');


%_____________________convert to Zgrid______________________
ZG.time    = TG.time;
ZG.zs      = [-100:1:0];
ZG.chi     = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_d   = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_dz  = nan( length(ZG.zs), length(ZG.time));
ZG.T       = nan( length(ZG.zs), length(ZG.time));
ZG.N       = nan( length(ZG.zs), length(ZG.time));
ZG.N_int   = nan( length(ZG.zs), length(ZG.time));
ZG.N_bll   = nan( length(ZG.zs), length(ZG.time));
ZG.chi_int     = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_d_int   = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_dz_int  = nan( length(ZG.zs), length(ZG.time));
ZG.chi_bbl     = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_d_bbl   = nan( length(ZG.zs), length(ZG.time));
ZG.Phi_dz_bbl  = nan( length(ZG.zs), length(ZG.time));

for t =1:length(ZG.time)
   tmpz = TG.zs(:,t);
   ii_nnan_tmpz = find(~isnan(tmpz));
   if ~isempty(ii_nnan_tmpz)
      ii = ii_nnan_tmpz;
      ZG.chi(:,t) = interp1( tmpz(ii), TG.chi(ii,t), ZG.zs);
      ZG.Phi_d(:,t) = interp1( tmpz(ii), TG.Phi_d(ii,t), ZG.zs);
      ZG.Phi_d_int(:,t) = interp1( tmpz(ii), TG.Phi_d_int(ii,t), ZG.zs);
      ZG.Phi_d_bbl(:,t) = interp1( tmpz(ii), TG.Phi_d_bbl(ii,t), ZG.zs);

      tmpz2 = tmpz(ii(2:end))-.5*diff(tmpz(ii));
      tmpPz =  diff(smooth(TG.Phi_d(ii,t), 11))./diff(tmpz(ii));
      ZG.Phi_dz(:,t) = interp1( tmpz2, tmpPz, ZG.zs);
      tmpPz =  diff(smooth(TG.Phi_d_int(ii,t), 11))./diff(tmpz(ii));
      ZG.Phi_dz_int(:,t) = interp1( tmpz2, tmpPz, ZG.zs);
      tmpPz =  diff(smooth(TG.Phi_d_bbl(ii,t), 11))./diff(tmpz(ii));
      ZG.Phi_dz_bbl(:,t) = interp1( tmpz2, tmpPz, ZG.zs);

      ZG.T(:,t)   =  interp1( tmpz(ii), TG.T(ii), ZG.zs);
      ZG.N(:,t)   = interp1( tmpz(ii), TG.N(ii), ZG.zs);
      ZG.N_int(:,t)   = interp1( tmpz(ii), TG.N_int(ii), ZG.zs);
      ZG.N_bbl(:,t)   = interp1( tmpz(ii), TG.N_bbl(ii), ZG.zs);
   end


end

   ZG.Tt = diff(ZG.T,1,2)./(diff(ZG.time)*3600*24);


save('../data/Turb_Zsspace.mat', 'ZG');





    
