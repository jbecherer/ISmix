close all;

if 0; %{{{
clear all;
   % path to m_map here
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/seawater/'));
   load cmap;

%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');
end %}}}


%% define T_grid

T_grid = 10:.02:19;

tl     = [datenum( 2017,9,8) datenum(2017,11,02)];
AW     = 72/24;   % average window 36 hours;
time   =  tl(1):AW/10:tl(2);   % time vector with 5 window width AW overlap    


%% loop through moorings
for m = 1:length(M)   % mooring loop

   M(m).Zs.time = time;
   M(m).Zs.T    = T_grid;
   M(m).Zs.zs      = nan( length(T_grid), length(time));

   for t = 1:length(time)  % Average time loop
      ii_t = find(  M(m).T.time >= time(t)-.5*AW & M(m).T.time <= time(t)+.5*AW  );
      
     %zs_tmp = nan( length(T_grid), length(ii_t));
     %for i = 1:length(ii_t)  % high freq time loop
     %   Ttmp = sort(M(m).T.T( :, ii_t(i)));
     %   if sum(~isnan(Ttmp))
     %      zs_tmp(:,i) = interp1( Ttmp, M(m).T.mab-M(m).waterdepth, T_grid  );
     %   end
     %end

      Ttmp = nanmedian(M(m).T.T( :, ii_t),2);
         % make sure stuff is good for interp1
         [Ttmp,ii_unique,~] = unique(Ttmp);
         ii_unique = ii_unique(~isnan(Ttmp));
         Ttmp      = Ttmp(~isnan(Ttmp));
      if sum(~isnan(Ttmp))
          M(m).Zs.zs(:, t) = interp1( Ttmp, M(m).T.mab(ii_unique)-M(m).waterdepth, T_grid  );
          %extrapolate to bottom and surface by constant strat
          [~,ii_minZs] = min(M(m).Zs.zs(:, t));
            M(m).Zs.zs(ii_minZs, t) = -M(m).waterdepth;
          [~,ii_maxZs] = max(M(m).Zs.zs(:, t));
            M(m).Zs.zs(ii_maxZs, t) = 0;
      end

   end

end

%% group moorings according to depth

Zs(1).ii_m = 1;
Zs(2).ii_m =  find_names_inM( M ,'50');
Zs(3).ii_m =  find_names_inM( M ,'40');
Zs(4).ii_m =  find_names_inM( M ,'3');
Zs(5).ii_m =  find_names_inM( M ,'25');

for iz = 1:length(Zs)
 if ~isempty(Zs(iz).ii_m)
   Zs(iz).time = time;
   Zs(iz).T    = T_grid;
   Zs(iz).zs   = nan(size(M(Zs(iz).ii_m(1)).Zs.zs));
   zstmp = Zs(iz).zs;
   for m = Zs(iz).ii_m
      zstmp = cat( 3,  zstmp, M(m).Zs.zs);
     %Zs(iz).zs = nansum(cat( 3,  zstmp, M(m).Zs.zs), 3);
     %keyboard
     %Zs(iz).zs( isnan(zstmp)& isnan(M(m).Zs.zs)) =nan;
   end
   Zs(iz).zs  = nanmean(zstmp ,3);
  end
end

      zstmp = nan(size(M(Zs(iz).ii_m(1)).Zs.zs));
      for iz = 1:length(Zs)
           zstmp = cat( 3,  zstmp, Zs(iz).zs);
      end
      zs_tot  = nanmean(zstmp ,3);
	   for t =1:2:length(Zs(1).time);
         ii_nnan = find(zs_tot(:,t));
         zs_tot(ii_nnan,t) = sort(zs_tot(ii_nnan,t));
	   end


fig = figure
   [ax, ~] = create_axes(fig,1,1, 0);

		col = get(groot,'DefaultAxesColorOrder');
   a=1;
	for t =1:2:length(Zs(1).time);
      for iz = 1:length(Zs)
          plot(ax(a), Zs(iz).T, nanmean(Zs(iz).zs(:,t),2), 'color', col(iz,:), 'Linewidth', 2);
      end
      plot(ax(a), Zs(iz).T, zs_tot(:,t), 'color', [0 0 0], 'Linewidth', 2);
		pause
      for iz = 1:length(Zs)
          plot(ax(a), Zs(iz).T, nanmean(Zs(iz).zs(:,t),2), 'color', col(iz,:)*.3+.7, 'Linewidth', 2);
      end
      plot(ax(a), Zs(iz).T, zs_tot(:,t), 'color', [0 0 0]+.7, 'Linewidth', 2);
	end
   
