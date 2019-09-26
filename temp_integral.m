close all; 
warning off;
%_____________________load mooring data______________________
if 0; %{{{
clear all;

   % path to m_map here
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/seawater/'));
   load cmap;


load('../data/Turb_Tspace.mat');
load('../data/Turb_Zsspace.mat');
% load bg state
   load('../data/Tsort.mat');


% load the monster grid
load ../data/IS_grid.mat

 % mean temperature in cube
 Tint  = nanmean(reshape(Grid.Tcube, [], size(Grid.Tcube,4)));


 
%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');


%_______________________________ cal boundary______________________

% find moorings on the southern boundary (SB), northern (NB), and Western (WB)
   Bdy.SBmoors = unique(Grid.mooring(:,1,:));
   Bdy.SBmoors = Bdy.SBmoors(~isnan(Bdy.SBmoors));
   Bdy.NBmoors = unique(Grid.mooring(:,end,:));
   Bdy.NBmoors = Bdy.NBmoors(~isnan(Bdy.NBmoors));
   Bdy.WBmoors = 1;  % only ms100 on western boundary
   Bdy.ALLmoors = unique([Bdy.SBmoors; Bdy.NBmoors; Bdy.WBmoors]);

% cal T-flauxes for all boundary moorings
      lp_cutoff= 3/24;  % pre filter
      hp_cutoff= 36/24; % secondary filter
 for i = 1:length(Bdy.ALLmoors)
    m = Bdy.ALLmoors(i);
    %disp(M(m).mooringName)
    %tic
    if ~isempty(M(m).U.time)
       [ Bdy.M(i).time, Bdy.M(i).z,  Bdy.M(i).uT, Bdy.M(i).T, Bdy.M(i).U] ...
            = cal_Tflux( M(m).U.time, M(m).U.mab, M(m).U.U, ...
                        M(m).T.time, M(m).T.mab, M(m).T.T, lp_cutoff, hp_cutoff, M(m).waterdepth);  
        Bdy.M(i).waterdepth    = M(m).waterdepth;
      
        Bdy.M(i).uT( imag(Bdy.M(i).uT)==0 )= nan + 1i*nan;
     else
        Bdy.M(i).time = [];
        Bdy.M(i).z    = [];
        Bdy.M(i).uT   = [];
        Bdy.M(i).T    = [];
        Bdy.M(i).U    = [];
        Bdy.M(i).waterdepth    = M(m).waterdepth;
     end
   %toc
 end

Bdy.SBarea = sum(sum(~isnan(Grid.mooring(:,1,:))))*Grid.dz*Grid.dx;
Bdy.NBarea = sum(sum(~isnan(Grid.mooring(:,end,:))))*Grid.dz*Grid.dx;
Bdy.WBarea = sum(sum(~isnan(Grid.mooring(:,:,1))))*Grid.dz*Grid.dy;
Bdy.Cubevol= sum(sum(sum(~isnan(Grid.mooring(:,:,:)))))*Grid.dz*Grid.dy*Grid.dx;

daysec = 3600*24;

% average data over particular boundary
Bdy.time = TS.time;

cnt1= 1;
cnt2= 1;
cnt3= 1;
for i = 1:length(Bdy.ALLmoors)
    m = Bdy.ALLmoors(i);
    % choose a plot
    if ~isempty( Bdy.M(i).time)
       if ismember(m, Bdy.WBmoors)
          bdyWB(cnt1,:) = clever_interp( Bdy.M(i).time, real(nanmean(Bdy.M(i).uT)), Bdy.time);
          cnt1 = cnt1+1;
       elseif ismember(m, Bdy.NBmoors)
          bdyNB(cnt2,:) = clever_interp( Bdy.M(i).time, -imag(nanmean(Bdy.M(i).uT)), Bdy.time);
          cnt2 = cnt2+1;
       else
          bdySB(cnt3,:) = clever_interp( Bdy.M(i).time, imag(nanmean(Bdy.M(i).uT)), Bdy.time);
          cnt3 = cnt3+1;
       end
    end
 end
 Bdy.uT_w = nanmean(bdyWB,1);
 Bdy.uT_n = nanmean(bdyNB,1);
 Bdy.uT_s = nanmean(bdySB,1);




end %}}}





 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
 
         [ax, ~] = create_axes(fig, 5, 1, 0);
         
         col = get(groot,'DefaultAxesColorOrder');
         

cnt1= 1;
cnt2= 1;
cnt3= 1;
yl = [-8 8];
tl = [datenum(2017,9,5) datenum(2017,11,5)];
         for i = 1:length(Bdy.ALLmoors)

            m = Bdy.ALLmoors(i);

            % choose a plot
            if ~isempty( Bdy.M(i).time)
               if ismember(m, Bdy.WBmoors)
                  a=1;
                  plot(ax(a), Bdy.M(i).time, real(nanmean(Bdy.M(i).uT))*Bdy.WBarea/Bdy.Cubevol*daysec, 'color', col(cnt1,:), 'Linewidth', 1);
                  txt = text( tl(1) + diff(tl)/7*cnt1,  yl(1), [M(m).mooringName], ...
                            'verticalalignment', 'bottom', 'horizontalalignment', 'center','Parent', ax(a));
                     txt.Color            = col(cnt1,:);
                  cnt1 = cnt1+1;
                  
                  
               elseif ismember(m, Bdy.NBmoors)
                  a=2;
                  plot(ax(a), Bdy.M(i).time, -imag(nanmean(Bdy.M(i).uT))*Bdy.NBarea/Bdy.Cubevol*daysec, 'color', col(cnt2,:), 'Linewidth', 1);
                  txt = text( tl(1) + diff(tl)/7*cnt2, yl(1), [M(m).mooringName], ...
                            'verticalalignment', 'bottom', 'horizontalalignment', 'center','Parent', ax(a));
                     txt.Color            = col(cnt2,:);
                  cnt2 = cnt2+1;
               else
                  a=3;
                  plot(ax(a), Bdy.M(i).time, imag(nanmean(Bdy.M(i).uT))*Bdy.SBarea/Bdy.Cubevol*daysec, 'color', col(cnt3,:), 'Linewidth', 1);
                  txt = text( tl(1) + diff(tl)/7*cnt3,  yl(1), [M(m).mooringName], ...
                            'verticalalignment', 'bottom', 'horizontalalignment', 'center','Parent', ax(a));
                     txt.Color            = col(cnt3,:);
                  cnt3 = cnt3+1;
               end
            end
            
         end


         a=5;
         plot(ax(a), Grid.time, Tint, 'Linewidth', 1);
            Grid.dt = nanmean(diff(Grid.time));
            Tint_lp = qbutter(Tint, Grid.dt/hp_cutoff);
         plot(ax(a), Grid.time, Tint_lp, 'Linewidth', 2);
            ylabel(ax(a), '^\circ C')
             

         a=4;
        pl(1) =  plot(ax(a), Grid.time([2:end]), 5* diff(Tint_lp)./(Grid.dt), 'color', [.5 0 0], 'Linewidth', 2);
        pl(2) =  plot(ax(a), Bdy.time, (Bdy.uT_s+Bdy.uT_s+Bdy.uT_n)*Bdy.SBarea/Bdy.Cubevol*daysec, 'color', [0 0 0], 'Linewidth', 2);
        pl(3) =  plot(ax(a), Bdy.time, Bdy.uT_w*Bdy.SBarea/Bdy.Cubevol*daysec, 'Linewidth', 1);
        pl(4) =  plot(ax(a), Bdy.time, Bdy.uT_n*Bdy.SBarea/Bdy.Cubevol*daysec, 'Linewidth', 1);
        pl(5) =  plot(ax(a), Bdy.time, Bdy.uT_s*Bdy.SBarea/Bdy.Cubevol*daysec, 'Linewidth', 1);
          
        
         


         abc='abcdefghijklmnopqrst';
         for a = 1:(size(ax,1)*size(ax,2))
            tabc = text_corner(ax(a), abc(a), 7);
            tabc.BackgroundColor = [1 1 1 .5];
            set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
            
            if a<5
               ylim(ax(a), yl);
               ylabel(ax(a), ['[K/day]'])
               plot(ax(a), tl, [0 0], 'k--','Linewidth', 1);
            end
            
            xlim(ax(a), tl);
            
         end
         
        t = text_corner(ax(1), 'Advective FLuxes into the BOX', -2);
         t.FontWeight      = 'bold';
         t.FontSize         = 12;
         t.Color            = [0 0 0];
         
        t = text_corner(ax(1), 'West Bdy', 1);
        t = text_corner(ax(2), 'North Bdy', 1);
        t = text_corner(ax(3), 'South Bdy', 1);
        t = text_corner(ax(4), 'Total flux', 1);
        t = text_corner(ax(5), '\langle T \rangle_{cube}', 1);

        lg = legend(pl, '5xdTdt', 'uT_{all}*A_{bdy}', 'uT_{w}*A_{bdyW}', 'uT_{n}*A_{bdyN}', 'uT_{s}*A_{bdyS}');
           axpos = get(ax(4), 'Position');
           lg.Position =  [axpos(1)+axpos(3)*.95 axpos(2) .02 axpos(4) ];
        
        
         linkaxes(ax, 'x');
         datetick(ax(a), 'x', 'mmm-dd ',  'keeplimits');
         
         print(gcf,'../pics/T_fluxes.png','-dpng','-r200','-painters')
         
         
         
            



