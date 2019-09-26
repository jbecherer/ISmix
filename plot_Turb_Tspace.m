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


end %}}}

%_____________________depth integral______________________
 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20]*.7,'PaperPosition',[0 0 30 20]*.7);
 
         [ax, ~] = create_axes(fig, 3, 1, 0);
         a=0;

         x = TG.time;

         collection{1} = '_int';
         collection{2} = '_bbl';
         collection{3} = '_abbl';
         collection{4} = '_surf';

         label{1} = 'all interior';
         label{2} = 'lower bbl';
         label{3} = 'upper bbl';
         label{4} = 'top ';

         col = get(groot,'DefaultAxesColorOrder');
         
         col(3,:) = col(1,:);
         col(1,:) = [0 0 0];

         for c = 1:4
            a=1;
            y = eval(['TG.chi' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-8 1e-5];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
               
               
            
            
            a=2;
            y = eval(['TG.Phi_d' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-7 1e-4];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);

            a=3;
            y = eval(['TG.Keff' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-6 1e-3];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
          end

         % legend(pl, 'all', 'bbl', 'top of BBL ', 'surface');


          txt(1) = text_corner(ax(1), ['\chi'], 1);
            ylabel(ax(1),'K^2/s')
          txt(2) = text_corner(ax(2), ['\Phi_d'], 1);
            ylabel(ax(2),'Km/s')
          txt(3) = text_corner(ax(3), ['K_{eff}'], 1);
            ylabel(ax(3),'m^2/s')
            
          ttxt =  text_corner(ax(1), ['Depth dependence'], -2);
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
         
         print(gcf,'../pics/mix_depthdep.png','-dpng','-r200','-painters')
         
            
         

%_____________________horizontal______________________
 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20]*.7,'PaperPosition',[0 0 30 20]*.7);
 
         [ax, ~] = create_axes(fig, 3, 1, 0);
         a=0;

         x = TG.time;

         collection{1} = '_int';
         collection{2} = '_100';
         collection{3} = '_50';
         collection{4} = '_30';

         label{1} = 'all interior';
         label{2} = 'ms100';
         label{3} = '50-40m';
         label{4} = '<35m';

         col([2:4],:) = col([5:7],:);

         for c = 1:4
            a=1;
            y  = eval(['TG.chi' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
           % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-8 1e-5];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
            
            
            a=2;
            y = eval(['TG.Phi_d' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
           % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-7 1e-4];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);

            a=3;
            y = eval(['TG.Keff' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            plot(ax(a), x, y1, 'color', col(c,:), 'Linewidth', 1);
           % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            yl = [5e-6 1e-3];
            ylim(ax(a), yl);
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
          end

          legend(pl, 'all', '100', '50-40m ', '<35m');

          txt(1) = text_corner(ax(1), ['\chi'], 1);
            ylabel(ax(1),'K^2/s')
          txt(2) = text_corner(ax(2), ['\Phi_d'], 1);
            ylabel(ax(2),'Km/s')
          txt(3) = text_corner(ax(3), ['K_{eff}'], 1);
            ylabel(ax(3),'m^2/s')
            
          ttxt =  text_corner(ax(1), ['isobath dependence'], -2);
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
         
         print(gcf,'../pics/mix_horzdep.png','-dpng','-r200','-painters')
         



%_____________________colored plot______________________
 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20]);
 
         [ax, ~] = create_axes(fig, 4, 1, 0);
         [ax(5:8), ~] = create_axes(fig, 4, 1, 0);
          squeeze_axes(ax, .47, 1)
          shift_axes(ax(1:4), -.05 , 0)
          shift_axes(ax(5:8), .43 , 0)
         
         a=0;

      % T space
            
         a=a+1;
         pcolor(ax(a), TG.time, TG.T, log10(TG.chi))
            txt(a) = text_corner(ax(a), ['log_{10}\langle \chi \rangle'], 1);

         a=a+1;
         pcolor(ax(a), TG.time, TG.T, log10(TG.dTdz))
            caxis(ax(a), [-2.2 -.5]);
            txt(a) = text_corner(ax(a), ['log_10 T_z  K/m'], 1);

         a=a+1;
         pcolor(ax(a), TG.time, TG.T, log10(TG.Jq))
            caxis(ax(a), [0 2.4]);
            txt(a) = text_corner(ax(a), ['log_{10} J_q  [W/m^2]'], 1);

         a=a+1;
         pcolor(ax(a), TG.time, TG.T, log10(TG.Keff))
            caxis(ax(a), [-6 -3]);
            txt(a) = text_corner(ax(a), ['log_{10} K_{eff} [m^2/s]'], 1);

            datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');

  
      % zs -space
            [TG.timeG,~] = meshgrid(TG.time, TG.T); 
            
         a=a+1;
         pcolor(ax(a), TG.timeG, TG.zs, log10(TG.chi))
            txt(a) = text_corner(ax(a), ['log_{10}\langle \chi \rangle'], 1);

         a=a+1;
         pcolor(ax(a), TG.timeG, TG.zs, log10(TG.dTdz))
            caxis(ax(a), [-2.2 -.5]);
            txt(a) = text_corner(ax(a), ['log_{10} T_z  K/m'], 1);

         a=a+1;
         pcolor(ax(a), TG.timeG, TG.zs, log10(TG.Jq))
            caxis(ax(a), [0 2.4]);
            txt(a) = text_corner(ax(a), ['log_{10} J_q  [W/m^2]'], 1);

         a=a+1;
         pcolor(ax(a), TG.timeG, TG.zs, log10(TG.Keff))
            caxis(ax(a), [-6 -3]);
            txt(a) = text_corner(ax(a), ['log_{10} K_{eff} [m^2/s]'], 1);

            datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
            
            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               tabc = text_corner(ax(a), abc(a), 7);
               tabc.BackgroundColor = [1 1 1 .5];
               set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
               txt(a).BackgroundColor = [1 1 1 .3];
               txt(a).FontWeight      = 'bold';
               txt(a).FontSize         = 10;
               txt(a).Color            = [0 0 0 ];
               shading(ax(a), 'flat');
               colormap(ax(a),  cmap.hblvil)
               xl = TG.time([1 end]);
               xlim(ax(a), xl);
               
               if a<5
                  ylabel(ax(a), 'T [^{\circ}C]')
                  yl = TG.T([1 end ]);
                  ylim(ax(a), yl);
                  
                 % plot(ax(a), M(1).T.time, movmean(M(1).T.T([1 end],:),3600*24/30 ,2), 'k', 'Linewidth', 1);
               else
                  ylabel(ax(a), 'z_* [m]')
                  cb = colorbar('peer', ax(a));
                     axpos = get(ax(a), 'Position');
                     cb.Position =  [axpos(1)+axpos(3)+.01 axpos(2) .01 axpos(4) ];
                  
               end
            end
            
            


         print(gcf,'../pics/Turb_Tspace.png','-dpng','-r200','-painters')
         

%_____________________flux divergence______________________

   ZG.Tt = ZG.Tt - repmat(nanmean(ZG.Tt,1), length(ZG.zs),1) ;   

 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
 
[ax, ~] = create_axes(fig, 3, 1, 0);
shift_axes(ax, -.04,0);

   a=0;

   a=a+1;
   pcolor(ax(a), ZG.time, ZG.zs, ZG.T);
   txt(a) =  text_corner(ax(a), ['T [^\circ C]'], 1);


   a=a+1;
   %pcolor(ax(a), ZG.time, ZG.zs, log10(abs(ZG.Phi_dz)));
   pcolor(ax(a), ZG.time, ZG.zs, ZG.Phi_dz*3600*24);
   caxis(ax(a), [-1 1]*.2)
   txt(a) =  text_corner(ax(a), ['d\Phi_d/dz_* [K/day]'], 1);

   a=a+1;
   %pcolor(ax(a), ZG.time(2:end), ZG.zs, log10(abs(diff(ZG.T,1,2)./(diff(ZG.time)*3600*24))));
   pcolor(ax(a), ZG.time(2:end), ZG.zs, ZG.Tt*3600*24);
   caxis(ax(a), [-1 1]*.5)
   txt(a) =  text_corner(ax(a), ['dT/dt [K/day]'], 1);

   datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
   abc='abcdefghijklmnopqrst';
   for a = 1:(size(ax,1)*size(ax,2))
      tabc = text_corner(ax(a), abc(a), 7);
      tabc.BackgroundColor = [1 1 1 .5];
      set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
      shading(ax(a), 'flat');
      colormap(ax(a), cmap.vel2)
      cb = colorbar('peer', ax(a));
         axpos = get(ax(a), 'Position');
         cb.Position =  [axpos(1)+axpos(3)+.01 axpos(2) .02 axpos(4) ];
         ylabel(ax(a), 'z_* [m]')
         
       txt(a).BackgroundColor = [1 1 1 .5];
       txt(a).FontWeight      = 'bold';
       txt(a).FontSize         = 10;
       txt(a).Color            = [0 0 0];
    
   end
   
   print(gcf,'../pics/Tt_Phidz.png','-dpng','-r200','-painters')
   

   
   
