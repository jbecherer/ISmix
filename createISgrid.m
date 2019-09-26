
close all;

if 1 %{{{

clear all;
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
load('~/gdrive/IS17/tau_map/data/bat_all_arrays.mat')
   bat_f = movmean(movmean( B.bat, 10, 1), 10, 2);
   B.bat_f = bat_f(1:10:end, 1:10:end)';
   B.lat_f = B.lat(1:10:end);
   B.lon_f = B.lon(1:10:end);
   B.bat_ff = movmean(movmean(B.bat_f,30,1),30,2);

%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');
   load cmap;

end %}}}


plot_grid = 0;
save_grid = 1;

% time vecor for Tcube
tl = [datenum(2017,9,11) datenum(2017,10,28,23,00,00)] ;
dt = 1/24/6; % 10 min time step
time = tl(1):dt:tl(2);

% only good moorings with full record
%[mss, ~] =  getFullMList(M);
% all moorings with temperature record
mss = 1:28; %length(M);
%_____________________get activity of moorings______________________
%  NNanTimes = nan( length(mss), length(time)); 
%  for m =1:length(mss)
%     if ~isempty(M(m).T.T)
%        ii_time =  find(time >= M(m).T.time(1) & time <= M(m).T.time(end));
%        NNanTimes(m,ii_time) = interp1( M(m).T.time, round(~isnan( nanmean(M(m).T.T) )), time(ii_time), 'nearest');
%     end
%  end




for i = 1:length(mss)
   m = mss(i);
   H(i) = M(m).waterdepth;
   Lat(i) = M(m).lat;
   Lon(i) = M(m).lon;
end


tic
%________ formulate lat lon grid _____ %{{{
Lat_lim = [min(Lat) max(Lat)];
Lat_lim = [min(Lat) 35.05];
Lon_lim = [min(Lon) max(Lon)];
Lon_lim = [min(Lon) -120.62];
% grid 500x500;
dx = 500;
dy = 500;
dlat = diff(Lat_lim)*(dy/1000)/m_lldist(Lon_lim([1 1]), Lat_lim([1 2]));
dlon = diff(Lon_lim)*(dx/1000)/m_lldist(Lon_lim([1 2]), Lat_lim([1 1]));
Lat_vec = Lat_lim(1):dlat:Lat_lim(2);
Lon_vec = Lon_lim(1):dlon:Lon_lim(2);

[Lon_grid, Lat_grid] = meshgrid(Lon_vec, Lat_vec);


[B.Lon_grid, B.Lat_grid] = meshgrid(B.lon_f, B.lat_f);
B_grid = griddata(B.Lon_grid(:), B.Lat_grid(:), B.bat_ff(:),Lon_grid, Lat_grid, 'nearest');

dz = 2;
z_vec = [(-99+.5*dz):dz:(0-.5*dz)];
M_grid = nan(length(z_vec), length(Lat_vec), length(Lon_vec));
for z = 1:length(z_vec)
   b_mask = ones(size(B_grid));
   b_mask(B_grid>z_vec(z)) = nan;
   if z_vec(z)>-10 % cut-off at 10m isobath 
       b_mask(B_grid>-10) = nan;
   end
   if z_vec(z)<-49 % deep water
      M_grid(z,:,:) = 1*b_mask;  % ms100
   else
      ii_mss = find(z_vec(z)>=-H);
      M_grid(z,:,:) = b_mask.*griddata(Lon(ii_mss),Lat(ii_mss)/3,mss(ii_mss),Lon_grid, Lat_grid/3, 'nearest');
   end
end

% }}}
toc

tic
%____________plot grid %{{{
if plot_grid;
colmap = jet(max(mss)+10);
colmap = cmap.hblvil(randsample(64,64),:);
colmap = rand(max(mss)+2,3);
colmap(1,:) = [.7 .6 .3];

   xl = [-120.85 -120.6];
   yl = [34.74 35.08];
 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20]);

         %[ax, ~] = create_axes(fig, 1, 1, 0);
            ax(1) = subplot( 3, 3, [1 4]);
            ax(2) = subplot( 3, 3, [2 5]);
            ax(3) = subplot( 3, 3, [3 6]);
            ax(4) = subplot( 3, 3, [7:9]);
         for a=1:4
            hold(ax(a), 'on');
         end


         col = prism(length(mss));
         col = lines(length(mss));
         j_trans = round(size(Lon_grid,1)*.8);
         
         z= [length(z_vec) round(length(z_vec)*[.7]) round(length(z_vec)*[.4])];
        for a = 1:3
            contour(ax(a), B.lon_f, B.lat_f, B.bat_f,[-100 -50 -30 0], 'k');
            contour(ax(a), B.lon_f, B.lat_f, B.bat_f,[0], 'k', 'Linewidth',3);
            contourf(ax(a), B.lon_f, B.lat_f, B.bat_f,[0 0]);
            %colormap(ax(a), cmap.bat_low);

            
            %contour(ax(a), Lon_grid, Lat_grid, B_grid,[-100 -50 -30 0], 'm');
           %for i = 1:length(mss)
           %   ii_col = find(M_grid(z(a),:,:) == mss(i));
           %   scatter(ax(a), Lon_grid(ii_col), Lat_grid(ii_col), 50 ,col(i,:), 'filled');
           %   %contour(ax(a), Lon_grid, Lat_grid, M_grid,[1 1]*mss(i), 'fill', 'on', 'FaceColor', col(i,:) );
           %end
            pcolor(ax(a), Lon_grid, Lat_grid, squeeze(M_grid(z(a), : ,:)))
            colormap(ax(a), colmap);
               caxis(ax(a),  [0 max(mss)]);

            for i = 1:length(mss)
               %t = text( Lon(i), Lat(i), [num2str(mss(i))], ...
              %t = text( Lon(i), Lat(i), ['x'], ...
              %          'verticalalignment', 'middle', 'horizontalalignment', 'center','Parent', ax(a));
              %t.BackgroundColor = [1 1 1 .8];
              %t.FontWeight      = 'bold';
              %t.FontSize         = 12;
              %t.Color            = col(i,:);
               
              plot(ax(a),  Lon(i), Lat(i), 'x', 'Color', col(i,:)*.5, 'Linewidth',3);
              plot(ax(a),  Lon(i), Lat(i), 'x', 'Color', col(i,:)*.3+.7, 'Linewidth',1);
            end

            plot(ax(a), Lon_grid, Lat_grid, 'color', [.5 .5 .5 .5],  'Linewidth', 1);
            plot(ax(a), Lon_grid', Lat_grid', 'color',[.5 .5 .5 .5],  'Linewidth', 1);

            % transect
            plot(ax(a), Lon_grid(j_trans,:), Lat_grid(j_trans,:), 'color',[0 0 0],  'Linewidth', 2);
            ylim(ax(a), yl);
            xlim(ax(a), xl);
            t = text_corner(ax(a), ['z = ' num2str(z_vec(z(a))) 'm'], 2);
              t.FontWeight      = 'bold';
            
        end
   

        a=4;

        xlim(ax(a), Lon_grid(j_trans,[1 end]));
        pcolor(ax(a), Lon_vec, z_vec, squeeze(M_grid(:,j_trans,:)))
        colormap(ax(a), colmap);
				caxis(ax(a),  [0 max(mss)]);
        yl = [-100 0];
        ylim(ax(a), yl);
        
        [lon_g, z_g] = meshgrid(Lon_vec, z_vec);
        %plot(ax(a), lon_g, z_g-.5*dz, 'color', [.5 .5 .5 .5],  'Linewidth', 1);
        %plot(ax(a), lon_g', z_g'-.5*dz, 'color', [.5 .5 .5 .5],  'Linewidth', 1);
        plot(ax(a), Lon_grid(j_trans,:), B_grid(j_trans,:), 'k', 'Linewidth', 3);
          %px = [Lon_grid(j_trans,:) fliplr(Lon_grid(j_trans,:))]; 
          %py = [B_grid(j_trans,:) fliplr(B_grid(j_trans,:))];
          %patch([px([1 1:end end])], [py([1 1:end end])], [.3 .3 .3], ...
          %      'facealpha', .7, 'edgecolor', [.5 0 0], 'Linewidth', 1, 'parent', ax(a));
         plot(ax(a), Lon_grid(j_trans,[1 end]), [0 0], 'color', [.3 .3 .7 ], 'Linewidth', 3);

			for iz = 1:length(z)
            plot(ax(a), Lon_grid(j_trans,[1 end]), [1 1]*z_vec(z(iz)), 'color', [.2 .2 .2 ], 'Linewidth', 3);
            plot(ax(a), Lon_grid(j_trans,[1 end]), [1 1]*z_vec(z(iz)), 'color', [1 1 1 ], 'Linewidth', 1);
			end

         abc='abcdefghijklmnopqrst';
         for a = 1:(size(ax,1)*size(ax,2))
            tabc = text_corner(ax(a), abc(a), 9);
            tabc.FontWeight      = 'bold';
            tabc.BackgroundColor = [1 1 1 .5];
            set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
         end
         
        



        print(gcf,'../pics/mesh_map.png','-dpng','-r200','-painters')
        
   

end %}}}
toc

tic
%_________CAL___________temperature game____{{{

Tcube    =  nan(length(z_vec), length(Lat_vec), length(Lon_vec), length(time));
% index matrix
M_grid_index =  nan(length(z_vec), length(Lat_vec), length(Lon_vec));
   for k = 1:size(M_grid,1)
    for j = 1:size(M_grid,2)
     for i = 1:size(M_grid,3)
         M_grid_index(k,j,i) = 1e6*k +1e3*j + i;
     end
    end
   end

Z_grid = repmat( z_vec',1, size(M_grid,2), size(M_grid,3) );

for m = mss % loop through good moorings
   moor_z = M(m).T.mab-M(m).waterdepth;

   for k = 1:length(moor_z)

      z_index = interp1( z_vec, 1:length(z_vec), moor_z, 'nearest');
      % depth range for the current mooring sensor
      if k ==1 % bottom
         sensor_z_range = [ -M(m).waterdepth   ,(moor_z(k) + .5*diff(moor_z([k k+1])) )];
      elseif k == length(moor_z); % surface
         sensor_z_range = [(moor_z(k) - .5*diff(moor_z([k-1 k])) ), 0];
      else
         sensor_z_range = [-.5*diff(moor_z([k-1 k])) , .5*diff(moor_z([k k+1]))] + moor_z(k)  ;
      end

      %inds_grid = find( M_grid==m & floor(M_grid_index./1e6)==z_index(k) );
      inds_grid = find( M_grid==m &  Z_grid>sensor_z_range(1) & Z_grid<=sensor_z_range(2) );
      if ~isempty(inds_grid)
         tmpT = clever_interp( M(m).T.time, M(m).T.T(k,:), time);
         % reintroduce nans in interpolation
         %tmpT( find(interp1( M(m).T.time, double(isnan(M(m).T.T(k,:))), time, 'nearest'))) = nan;
         for ii_ind = inds_grid'
            k_ind = floor(M_grid_index(ii_ind)/1e6);
            j_ind = floor((M_grid_index(ii_ind)-1e6*floor(M_grid_index(ii_ind)/1e6))/1e3);
            i_ind = M_grid_index(ii_ind)-1e3*floor(M_grid_index(ii_ind)/1e3);
            Tcube( k_ind, j_ind, i_ind,: ) = tmpT;    
         end
      end

   end

end

toc

%_____________________create hypsometry______________________
NNan = zeros(size(Tcube,1), size(Tcube,4));
for t = 1:size(Tcube,4)
   for z = 1:size(Tcube,1)
      NNan(z,t)  = sum(sum(~isnan(Tcube(z,:,:,t))));
   end
end
A_z  = NNan*dx*dy;



% fluxes on boundaries --- U-slices
% to be implemented ......


if save_grid
   Grid.lon       = Lon_vec;
   Grid.lat       = Lat_vec;
   Grid.z         = z_vec;
   Grid.NNan      = NNan;
   Grid.A_z       = A_z;
   Grid.time      = time;
   Grid.mooring   = M_grid;
   Grid.Tcube     = Tcube;
   Grid.dz        = dz;
   Grid.dx        = dx;
   Grid.dy        = dy;
   Grid.dlat      = dlat;
   Grid.dlon      = dlon;

   %save('~/gdrive/IS17/mix_data/IS_grid.mat', 'Grid', '-v7.3');
   save('../data/IS_grid.mat', 'Grid', '-v7.3');
end %}}}

