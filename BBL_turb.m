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


% load bottom GusTs
  load('../data/G_bbl.mat')

load('../data/Turb_Tspace.mat');
load('../data/Turb_Zsspace.mat');
% load bg state
   load('../data/Tsort.mat');


end %}}}

%_____________________depth integral______________________
 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[20 12]*1,'PaperPosition',[0 0 20 12]*1);
 
         [ax, ~] = create_axes(fig, 2, 1, 0);
         squeeze_axes(ax,.8,1);
         a=0;

         x = TG.time;

         collection{1} = '_bbl';
         collection{2} = '_abbl';

         col = get(groot,'DefaultAxesColorOrder');
         
         %col(3,:) = col(1,:);
         col(3,:) = [0 0 0];

         for c = 1:2
            a=1;
            y = eval(['TG.eps' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [1e-8 1e-5];
            ylim(ax(a), yl);
            
            
            a=2;
            y = eval(['TG.Kt' collection{c}]);
            %y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            y1 = nanmean(y );
            y = eval(['TG.Keff' collection{c}]);
            y2 = nanmean(y );
            pl1(1+(c-1)*2) = plot(ax(a), x, y1,'--', 'color', col(c,:), 'Linewidth', 1);
            pl1(2+(c-1)*2) =plot(ax(a), x, y2,  'color', col(c,:), 'Linewidth', 1);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-6 1e-2];
            ylim(ax(a), yl);

          end

          pl1(5) = plot(ax(a), x, nanmean(TG.Keff), 'color', col(3,:), 'Linewidth', 1);
          lg = legend(pl, '\epsilon (bbl)', '\epsilon (3-8mab)');
          axpos = get(ax(1), 'Position');
          lg.Position = [.82 axpos(2) .1 axpos(4)];
          
          lg = legend(pl1, 'Pr_tK_T^{local} (bbl)', 'K_{eff} (bbl)' ,'K_{T}^{local} (3-8mab)','K_{eff} (3-8mab)', 'K_{eff} (all)');
          axpos = get(ax(2), 'Position');
          lg.Position = [.82 axpos(2) .1 axpos(4)];

          txt(1) = text_corner(ax(1), ['\epsilon'], 1);
            ylabel(ax(1),'m^2/s^3')
          txt(2) = text_corner(ax(2), ['K_T'], 1);
            ylabel(ax(2),'m^2/s')
            
          ttxt =  text_corner(ax(1), ['mixing in the BBL'], -2);
            ttxt.BackgroundColor = [1 1 1];
            ttxt.FontWeight      = 'bold';
            ttxt.FontSize         = 12;
            ttxt.Color            = [0 0 0];
          

          abc='abcdefghijklmnopqrst';
          for a = 1:(size(ax,1)*size(ax,2))
             tabc = text_corner(ax(a), abc(a), 7);
             tabc.BackgroundColor = [1 1 1 .5];
             txt(a).BackgroundColor = [1 1 1];
             txt(a).FontWeight      = 'bold';
             txt(a).FontSize         = 12;
             txt(a).Color            = [0 0 0];
             xl = TS.time([1 end]);
             xlim(ax(a), xl);
             
             
             set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
          end
          
          
         datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
         linkaxes(ax, 'x');
         
         print(gcf,'../pics/mix_bbl.png','-dpng','-r200','-painters')
         
