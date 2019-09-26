
close all;

if 1 %{{{

clear all;
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/mixingsoftware/seawater/'));
   

   load('../data/IS_grid.mat');


 % load('~/gdrive/IS17/tau_map/data/bat_all_arrays.mat')
 % bat_f = movmean(movmean( B.bat, 10, 1), 10, 2);
 % B.bat_f = bat_f(1:10:end, 1:10:end)';
 % B.lat_f = B.lat(1:10:end);
 % B.lon_f = B.lon(1:10:end);
 % B.bat_ff = movmean(movmean(B.bat_f,30,1),30,2);

 % %_____________________all epsilon data______________________
 % load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');
 % load cmap;

end %}}}

sort_T       = 1;
save_profile = 1;


if sort_T 


DT         = 36/24; 
TS.DT   = DT; 
TS.time = (Grid.time(1)+.5*DT):.5*DT:(Grid.time(end)-.5*DT);
TS.T_vec  = 9.5:.1:19;
TS.zs_vec = -100:1:0;

TS.Ts   = nan( length(TS.zs_vec), length(TS.time));
TS.zs   =  nan( length(TS.T_vec), length(TS.time));


Nk = zeros( length(Grid.z), length(TS.time));
tlsub = [-.5 .5]*TS.DT + TS.time(1);
	ii_time =  find(Grid.time >= tlsub(1) & Grid.time <= tlsub(end) );
	Ntime =  length(ii_time);


dz    = Grid.dz;

tic
for t = 1:length(TS.time)

   if t >1
    ii_time = ii_time + floor(Ntime*.5);
   end
   TS.time(t) = mean(Grid.time(ii_time));

   % sort temperature
   T_tmp = Grid.Tcube(:,:,:, ii_time);
  T_sort = sort(T_tmp(~isnan(T_tmp)));
  % T_sort = sort(T_tmp(:));


   %% count number of boxes for each z-layer
   l_k = [];
   for k = 1:size(Grid.Tcube,1)
       Nk(k,t) = sum(sum(sum(~isnan(Grid.Tcube(k,:,:, ii_time))))); 
       l_k  = cat(2, l_k, ones( 1, Nk(k,t))*k );
   end
   
   %zs    = nan(size(T_sort));
   zs    = nan(size(l_k));
   zs(1) = Grid.z(1)+dz/Nk(1,t)-dz*.5;
   for l = 2:length(l_k)
      zs(l) = zs(l-1) + dz/Nk(l_k(l),t);
   end

   TS.Ts(:,t) = interp1( zs, T_sort, TS.zs_vec );

   [~,ii_uni,~] = unique(T_sort);
   TS.zs(:,t) = interp1( T_sort(ii_uni), zs(ii_uni), TS.T_vec );

end
toc

%_____________________cal vertical gradients______________________
TS.dTdz = nan( length(TS.zs_vec)-1, length(TS.time));
   TS.z_dTdz = TS.zs_vec(2:end)-.5*diff(TS.zs_vec);
   TS.N2 = nan( length(TS.zs_vec)-1, length(TS.time));

TS.dzdT = nan( length(TS.T_vec)-1, length(TS.time));
   TS.T_dzdT = TS.T_vec(2:end)-.5*diff(TS.T_vec);
for t = 1:length(TS.time)
   % cal vertical gradient
   
   TS.dTdz(:,t) = diff(smooth(TS.Ts(:,t),7))./diff(TS.zs_vec)';
   TS.N2(:,t) = 9.81*TS.dTdz(:,t).*sw_alpha( 35*ones(size(TS.dTdz(:,t))), TS.dTdz(:,t), 30*ones(size(TS.dTdz(:,t))), 'temp' );
      
   TS.dzdT(:,t) = diff(smooth(TS.zs(:,t),7))./diff(TS.T_vec)';

end


end



%_____________________save sorted profile______________________
if save_profile
save('../data/Tsort.mat', 'TS');
end



%_____________________plot T-sorted______________________
fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
        'Papersize',[20 10],'PaperPosition',[0 0 20 10]);
        [ax, ~] = create_axes(fig, 1, 3, 0);
         a=1;
         ii = [1 4 8]+10;
         plot(ax(a), TS.Ts(:,ii), TS.zs_vec, 'Linewidth', 1);
         
         ylabel(ax(a),'z m');
         xlabel(ax(a),'T {^\circ}C');

   
         a=2;
         plot(ax(a), TS.dTdz(:,ii), TS.z_dTdz, 'Linewidth', 1);
         xlabel(ax(a),'dTdz {^\circ}C/m');

         a=3;
         plot(ax(a), TS.N2(:,ii), TS.z_dTdz, 'Linewidth', 1);
         xlabel(ax(a),'N2 s^-2');

         legend(ax(a), datestr(TS.time(ii(1)), 'mmm-dd'),datestr(TS.time(ii(2)), 'mmm-dd'), datestr(TS.time(ii(3)), 'mmm-dd') )

         abc='abcdefghijklmnopqrst';
         yl = [-100 0];
         for a = 1:(size(ax,1)*size(ax,2))
            tabc = text_corner(ax(a), abc(a), 7);
            tabc.BackgroundColor = [1 1 1 .5];
            set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
            ylim(ax(a), yl);
            
         end
         
         print(gcf,'../pics/Sorted_bg_State.png','-dpng','-r200','-painters')
         

%__   Potential energy in sorted state
   
   % Hypsometry of sorted state
   TS.A_z      = interp1( Grid.z,Grid.A_z(:,1), TS.zs_vec, 'nearest', 'extrap');
   TS.E_pot_bg = nan(size(TS.time));
   TS.dz = nanmedian(diff(TS.zs_vec));
   ref_dens = sw_pden(35, 9, 100, 1 );
   for t = 1:length(TS.time)
      Ttmp = reshape(TS.Ts(:,t),[],1);
      Dtmp = sw_pden( ones(size(Ttmp))*35, Ttmp, abs(TS.zs_vec'), 1 )-ref_dens;
      TS.E_pot_bg(t) =  nansum(9.81*(Dtmp)'.*(TS.zs_vec+100).*TS.A_z.*TS.dz); 
   end

   

%__   Potential energy in Tcube?

Epot     = zeros(1,size(Grid.Tcube,4));
Vol      = zeros(1,size(Grid.Tcube,4));
dV = Grid.dz*Grid.dx*Grid.dy;
tic
for t = 1:size(Grid.Tcube,4)
   for z=1:length(Grid.z)
      Ttmp = reshape(Grid.Tcube(z,:,:,t),[],1);
      Dtmp = sw_pden( ones(size(Ttmp))*35, Ttmp, abs(Grid.z(z)), 1 ) - ref_dens;
      %Nnan = sum(sum(sum(~isnan(Grid.Tcube(z,:,:,t)))));
      Nnan = Grid.NNan(z,t);
      Epot(t) = Epot(t) +  9.81*(nanmean(Dtmp))*(Grid.z(z)+100)*Nnan*dV ; 
      Vol(t) = Vol(t) + Nnan*dV ; 
   end
end
toc


figure
for a=1:2
   ax(a) = subplot( 2, 1, a);
   hold(ax(a), 'on');
end
a=1;
   plot(ax(a), Grid.time, Epot);
   hold(ax(a), 'all')
   plot(ax(a), TS.time, TS.E_pot_bg);

a=2;
TS.E_pot_tot = clever_interp(Grid.time, qbutter(Epot,1/200), TS.time);
TS.E_pot_ava = TS.E_pot_tot-TS.E_pot_bg;
plot(TS.time, TS.E_pot_ava);


   
