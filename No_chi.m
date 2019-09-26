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
   load('../data/Turb_Tspace.mat');
   load('../data/Turb_Zsspace.mat');
end %}}}

tl =  [datenum(2017,9,8) datenum(2017,11,10)];


TG.SumN = nansum(TG.N,1);  
ZG.SumN = nansum(TG.N,1);  
   ZG.N = ZG.N.*repmat(TG.SumN,length(ZG.zs),1)./repmat(ZG.SumN,length(ZG.zs),1);

 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20]*.8,'PaperPosition',[0 0 30 20]*.8);
 
         [ax, ~] = create_axes(fig, 3, 1, 0);
         shift_axes(ax, -.03,-.01);
         a=0;

         tittxt = text_corner(ax(1), ['N_\chi = ' num2str(nansum(nansum(TG.N)))], -2);
               tittxt.FontWeight      = 'bold';
               tittxt.FontSize         = 10;
               tittxt.Color            = [0 0 0 ];
         
      % T space
         a=a+1;
         pcolor(ax(a), TG.time, TG.T, log10(TG.N))
            txt(a).t = text_corner(ax(a), ['log_{10}N_\chi [1/bin]'], 1);
               ylabel(ax(a), 'T [^{\circ}C]')
                     axpos = get(ax(a), 'Position');
                     ax(end+1) =  axes('Position', [axpos(1)+axpos(3)-.1 axpos(2) .1 axpos(4) ]);
                     plot(ax(end), nansum(TG.N,2) , TG.T, 'Linewidth', 1);
                     hold(ax(end),'on');
                     plot(ax(end), nansum(TG.N_bbl,2) , TG.T, 'Linewidth', 1);
                     xl = [1e1 8e4];
                     xlim(ax(end), xl);
                     set(ax(end),'XaxisLocation', 'top', 'Xscale', 'log', 'box', 'on', 'TickDir', 'out', 'Layer', 'top', 'YTickLabel',{});
                     
                     
                     

      % zs -space
         a=a+1;
         pcolor(ax(a), ZG.time, ZG.zs, log10(ZG.N))
            txt(a).t = text_corner(ax(a), ['log_{10}N_\chi [1/bin]'], 1);
            ylabel(ax(a), 'z* [m]')
                     axpos = get(ax(a), 'Position');
                     ax(end+1) =  axes('Position', [axpos(1)+axpos(3)-.1 axpos(2) .1 axpos(4) ]);
                     plot(ax(end), nansum(ZG.N,2) , ZG.zs, 'Linewidth', 1);
                     hold(ax(end),'on');
                     plot(ax(end), nansum(ZG.N_bbl,2) , ZG.zs, 'Linewidth', 1);
                     xlim(ax(end), xl);
                     set(ax(end), 'Xscale', 'log', 'box', 'on', 'TickDir', 'out', 'Layer', 'top','XTickLabel',{}, 'YTickLabel',{});

      % vertical integral
         a=a+1;
         plot(ax(a), TG.time, TG.SumN, 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_bbl,1), 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_abbl,1)+nansum(TG.N_bbl,1), 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_surf,1), 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_100,1), 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_50,1), 'Linewidth', 1);
         plot(ax(a), TG.time, nansum(TG.N_30,1), 'Linewidth', 1);
         legend(ax(a), 'all', 'bbl (1.5 mab)', '<6mab', 'surface (z>-15m)', 'ms100', 'wd =50-40m', 'wd<35m')
            txt(a).t = text_corner(ax(a), ['N_\chi [1/day]'], 1);
            set(ax(a), 'Yscale', 'log');
         

            datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
            
                     a=a+1;
                     txt(a).t = text_corner(ax(a), ['N_\chi [10/K]'], 1);
                     a=a+1;
                     txt(a).t = text_corner(ax(a), ['N_\chi [1/m]'], 1);

            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               tabc = text_corner(ax(a), abc(a), 7);
               tabc.BackgroundColor = [1 1 1 .5];
               set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
               txt(a).t.BackgroundColor = [1 1 1 .3];
               txt(a).t.FontWeight      = 'bold';
               txt(a).t.FontSize         = 10;
               txt(a).t.Color            = [0 0 0 ];
               if a<3
                  shading(ax(a), 'flat');
                  cb = colorbar('peer', ax(a));
                     axpos = get(ax(a), 'Position');
                     cb.Position =  [axpos(1)+axpos(3)+.01 axpos(2) .02 axpos(4) ];
                  colormap(ax(a),  cmap.hblvil)
               end
               if a<4
                xlim(ax(a), tl);
               end
               
            end
            
            


         print(gcf,'../pics/No_chi.png','-dpng','-r200','-painters')
         
            
                  
