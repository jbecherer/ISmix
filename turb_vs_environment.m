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

%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');


% wind
load /home/johannes/gdrive/IS17/collection_for_Jim/data/Mini_Met_Innershelf.mat;

Met = MiniMet;
Met.Van = load('/home/johannes/gdrive/IS17/collection_for_Jim/data/meteo_vandenberg.mat');

Met.sw  = clever_interp( Met.Van.M.time, Met.Van.M.radiation, Met.TimeOR2);
Met.speedOR2  = clever_interp( Met.TimeOR1, Met.WindSpeedOR1, Met.TimeOR2);
Met.SSTOR2  = clever_interp( Met.TimeSST, Met.SST, Met.TimeOR2);

%    [Met.Jb,Met.Jq]=surfaceflux( Met.AirTempOR2, zeros(size(Met.AirTempOR2)), Met.AirPressureOR2 ...
%                     ,Met.SSTOR2 , 35*ones(size(Met.TimeOR2)), Met.speedOR2 ...
%                     , Met.RelativeHumidityOR2 ,   Met.sw, 0, 0);

         Met.U = Met.Van.M.u10.*exp(1i*(-(Met.Van.M.udir-270))/180*pi);;

% pressure record
    load ~/gdrive/IS17/tau_map/data/pressure_record_PS50.mat;
    P.press = movmean(P.press,10);


   [mss, mshalf] = getFullMList(M);
   %mss = [mss mshalf];
   mss100 = mss(1);
   mss50 = []; for m=mss; if( M(m).waterdepth>45 & M(m).waterdepth<55); mss50 = [mss50 m]; end;end
   mss40 = []; for m=mss; if( M(m).waterdepth>35 & M(m).waterdepth<45); mss40 = [mss40 m]; end;end
   mss30 = []; for m=mss; if( M(m).waterdepth<35); mss30 = [mss30 m]; end;end

   IMs{1} = mss100;
   IMs{2} = mss50;
   IMs{3} = mss40;
   IMs{4} = mss30;

      anything = 'Ef.int_uppp';
      prefix   =  'log10(real(';
      suffix   =  '))';

       for i = 1:length(IMs);
          [any_time, any_avg, any_std] = mooring_average_anything(M, IMs{i}, anything, prefix, suffix);
          if(i>1)
           ii_notnan = find(~isnan(any_avg));
           any_time   = any_time(ii_notnan);
           any_avg = any_avg(ii_notnan);
           any_std    = any_std(ii_notnan);
          end
          Fe.time{i} = any_time;
          Fe.Fe{i}   = 10.^any_avg;
          Fe.Festd{i}   = any_std ;
       end


 %[N2_time, N2_maxavg, N2_std] = mooring_average_N2(M, mss);
 %      ii_notnan = find(~isnan(N2_maxavg));
 %      N2_time   = N2_time(ii_notnan);
 %      N2_maxavg = N2_maxavg(ii_notnan);
 %      N2_std    = N2_std(ii_notnan);

load('../data/Turb_Tspace.mat');
load('../data/Turb_Zsspace.mat');
% load bg state
   load('../data/Tsort.mat');


end  %}}}

%___________________set the clock___________________


 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 30]*.8,'PaperPosition',[0 0 30 30]*.8);

         [ax, ~] = create_axes(fig, 7, 1, 0);
         shift_axes(ax, -0.02, 0)
        
         col = get(groot,'DefaultAxesColorOrder');
         col(1,:) = [0 0 .4];
         col(3,:) = [.3 .6 0];


         xl = TS.time([1 end]);

         
         a=0;


 %_____________________surface______________________{{{
 if 1
         a=a+1;
         plot(ax(a), P.time, sqrt(2)*qbutter(movstd(P.press-nanmedian(P.press),60*36), 1/(60*36)),'k', 'Linewidth', 2);

            t = text_corner(ax(a), ['2^{1/2}  std'], 14);
               t.BackgroundColor = [1 1 1 .5];
               t.FontWeight      = 'bold';
       xlim(ax(a), xl);
       ylabel(ax(a), 'm');
       ylim(ax(a), [.3 1.1]);
       

         [axp, axmoon] = copy_axes(ax(a));
         plot(axp, P.time, P.press-nanmedian(P.press), 'color', [.5 .5 .5 .5], 'Linewidth', 1);
         plot(axp, P.time, qbutter(movmin(P.press-nanmedian(P.press),60*36), 1/(60*36)),'color',[.5 .5 .5], 'Linewidth', 2);
         plot(axp, P.time, qbutter(movmax(P.press-nanmedian(P.press),60*36), 1/(60*36)),'color',[.5 .5 .5], 'Linewidth', 2);

         plot(ax(a), P.time, .5*qbutter(movmax(P.press-nanmedian(P.press),60*36) - movmin(P.press-nanmedian(P.press),60*36), 1/(60*36)),'color',[.7 .2 .2], 'Linewidth', 2);
            xlim(axp, xl);
            ylim(axp, [-1.5 1.5]);
            set(axp, 'Yaxislocation', 'right', 'Ycolor', [.5 .5 .5], 'Ytick', [-1 0 1])
            ylabel(axp, 'm')
            t = text_corner(ax(a), ['surface elevation'], 3);
               t.BackgroundColor = [1 1 1 .5];
               t.Color = [1 1 1]*.5;
               t.FontWeight      = 'bold';
            t = text_corner(ax(a), ['1/2(max-min)'], 1);
               t.BackgroundColor = [1 1 1 .5];
               t.Color = [1 0 0 ]*.5;
               t.FontWeight      = 'bold';

         %__________moon
         plot(axmoon, [datenum(2017,9,20,6,0,0) datenum(2017,10,19, 20,0,0)], [1 1], 'o', 'color', [.7 .7 .7], 'MarkerSize', 12,'Linewidth', 2, 'MarkerFaceColor', [0 0 0]);
         plot(axmoon, [datenum(2017,9,6,8,0,0) datenum(2017,10,5, 20,0,0) datenum(2017,11,4, 6,0,0)], [1 1 1], 'o', 'color', [.3 .3 .3], 'MarkerSize', 12,'Linewidth', 2, 'MarkerFaceColor', [1 1 .8]);
         set(axmoon, 'visible', 'off', 'clipping', 'off');
         xlim(axmoon, xl);
         ylim(axmoon, [0 .9]);
            txt = text_corner(axmoon, ['moon phase'], -1);
               txt.BackgroundColor = [1 1 1 .5];
               txt.FontWeight      = 'bold';
         
 end
%}}}

 %_____________________wind______________________{{{
  if 1
         a=a+1;
      txt = text_corner(ax(a), ['wind'], 3);
      txt.BackgroundColor = [1 1 1];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
            yl = [-6.5 2];
            dt1 = (3600*2)/(24*3600);
            time1 = xl(1):dt1:xl(2);
            plot_wind( time1,  qbutter(clever_interp(Met.Van.M.time, Met.U, time1),.1), ax(a), fig, xl, yl, [.5 .5 .5 .3], 1 );
            plot_wind( time1(1:20:end),  clever_interp(Met.Van.M.time, Met.U, time1(1:20:end)), ax(a), fig, xl, yl, [0 0 0 1], 1.5 );
               xlim(ax(a), xl);
               plot(ax(a), xl, [0 0],'k--', 'Linewidth', 1);

            ylabel(ax(a), 'm/s')
            ylim(ax(a), yl);

        [axw, ~]  = copy_axes(ax(a));
         plot(axw,  time1,  abs(qbutter(clever_interp(Met.Van.M.time, Met.U, time1), nanmean(diff(time1))*24/36 )), 'color', [0 0 .5], 'Linewidth', 2);
            ylabel(axw, 'speed m/s')
            ylim(axw, [0 10]);
            set(axw, 'YAxisLocation','right', 'Ycolor', [0 0 .5], 'Ytick', [0 5 10], 'Tickdir', 'out');
            xlim(axw, xl);
        
   end
%}}}

%_____________________N2_____________________{{{_
if 0
         a=a+1;
         px = [N2_time fliplr(N2_time)];
         py = [ movmean(N2_maxavg-N2_std,3) , movmean(fliplr(N2_maxavg+N2_std),3)]*1e4;
         patch(px, py, [.5 .5 .5], 'facealpha', .5, 'edgecolor', [.5 .5 .5], 'Linewidth', 1, 'parent', ax(a));
         plot(ax(a), N2_time, movmean(N2_maxavg, 3)*1e4, 'k',  'Linewidth', 2);
        xlim(ax(a), xl);
        ylabel(ax(a), '10^{-4}s^{-2}');
         t = text_corner(ax(a), ['N^2_{max}'], 3);
               t.BackgroundColor = [1 1 1 .5];
               t.FontWeight      = 'bold';
end
%}}}
         
%_____________________Fe______________________{{{
if 1
         a=a+1;
      txt = text_corner(ax(a), ['F_E'], 3);
      txt.BackgroundColor = [1 1 1];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
        

      i_fe=2;
      px = [Fe.time{i_fe} fliplr(Fe.time{i_fe})];
      py = 10.^[ movmean(log10(Fe.Fe{i_fe})-Fe.Festd{i_fe},3) , fliplr(movmean(log10(Fe.Fe{i_fe})+Fe.Festd{i_fe},3))];
      patch(px, py, col(i_fe,:), 'facealpha', .5, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
       plot(ax(a), Fe.time{i_fe} , Fe.Fe{i_fe} , 'color', col(i_fe,:)*.6,  'Linewidth', 2);
       
       ylabel(ax(a), 'W/m');
       set(ax(a), 'Yscale', 'log');
       yl = [5e0 1e2];
       ylim(ax(a), yl);
end
%}}}


%_____________________DT______________________{{{
if 1
         a=a+1;
      txt = text_corner(ax(a), ['{\Delta}T = T(z_*=-5m) - T(z_*=-40m)'], 3);
      txt.BackgroundColor = [1 1 1 .5];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];

      ii = [find(TS.zs>-40,1) find(TS.zs>-5,1)  ];
        
      plot(ax(a), TS.time, smooth(diff(TS.Ts(ii,:)),3), 'k', 'Linewidth', 2);
      
       ylabel(ax(a), 'K');
       
      yl = [1 5];
      ylim(ax(a), yl);
end
%}}}


%_____________________chi______________________{{{
if 1
      a=a+1;
      txt = text_corner(ax(a), ['\chi'], 3);
      txt.BackgroundColor = [1 1 1];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
         x = TG.time;

         collection{1} = '';
         collection{2} = '_bbl';
         collection{3} = '_abbl';
         collection{4} = '_surf';

         label{1} = 'all';
         label{2} = 'lower bbl';
         label{3} = 'upper bbl';
         label{4} = 'top ';

         col = get(groot,'DefaultAxesColorOrder');
         
         col(3,:) = col(1,:);
         col(1,:) = [0 0 0];

         for c = 1%:4
            y = eval(['TG.chi' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x,smooth( y1,5), 'color', col(c,:), 'Linewidth', 2);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
         end
               
        %yl = [5e-8 1e-5];
        yl = [2e-7 3e-6];
        ylim(ax(a), yl);
        ylabel(ax(a), '\chi K^2/s')
       %       
       % a=a+1;
       % collection{1} = '_int';
       % collection{2} = '_100';
       % collection{3} = '_50';
       % collection{4} = '_30';

       % label{1} = 'all interior';
       % label{2} = 'ms100';
       % label{3} = '50-40m';
       % label{4} = '<35m';

       % col([2:4],:) = col([5:7],:);

       % for c = 1:4
       %    y  = eval(['TG.chi' collection{c}]);
       %    y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
       %    pl(c) = plot(ax(a), x, smooth(y1,5), 'color', col(c,:), 'Linewidth', 1);
       %   % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
       %    set(ax(a), 'Yscale', 'log');
       %    yl = [5e-8 1e-5];
       %    ylim(ax(a), yl);
       %    txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
       %                'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
       %       txt.Color            = col(c,:);
       % end
end            
%}}}

%_____________________Phid______________________{{{
if 1
      a=a+1;
      txt = text_corner(ax(a), ['\Phi_d'], 3);
      txt.BackgroundColor = [1 1 1];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
      
      
         x = TG.time;

         collection{1} = '';
         collection{2} = '_bbl';
         collection{3} = '_abbl';
         collection{4} = '_surf';

         label{1} = 'all';
         label{2} = 'lower bbl';
         label{3} = 'upper bbl';
         label{4} = 'top ';

         col = get(groot,'DefaultAxesColorOrder');
         
         col(3,:) = col(1,:);
         col(1,:) = [0 0 0];

         for c = 1%:4
            y = eval(['TG.Phi_d' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x,smooth( y1,5), 'color', col(c,:), 'Linewidth', 2);
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
         end
               
        yl = [1e-6 2e-5];
        ylim(ax(a), yl);
        ylabel(ax(a), 'K m/s')
       %       
       % a=a+1;
       % collection{1} = '_int';
       % collection{2} = '_100';
       % collection{3} = '_50';
       % collection{4} = '_30';

       % label{1} = 'all interior';
       % label{2} = 'ms100';
       % label{3} = '50-40m';
       % label{4} = '<35m';

       % col([2:4],:) = col([5:7],:);

       % for c = 1:4
       %    y  = eval(['TG.chi' collection{c}]);
       %    y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
       %    pl(c) = plot(ax(a), x, smooth(y1,5), 'color', col(c,:), 'Linewidth', 1);
       %   % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
       %    set(ax(a), 'Yscale', 'log');
       %    yl = [5e-8 1e-5];
       %    ylim(ax(a), yl);
       %    txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
       %                'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
       %       txt.Color            = col(c,:);
       % end
end            
%}}}

%_____________________Keff______________________{{{
if 1
      a=a+1;
      txt = text_corner(ax(a), ['K_{eff}'], 3);
      txt.BackgroundColor = [1 1 1];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
         x = TG.time;

         collection{1} = '';
         collection{2} = '_bbl';
         collection{3} = '_abbl';
         collection{4} = '_surf';

         label{1} = 'all';
         label{2} = 'lower bbl';
         label{3} = 'upper bbl';
         label{4} = 'top ';

         col = get(groot,'DefaultAxesColorOrder');
         
         col(3,:) = col(1,:);
         col(1,:) = [0 0 0];

         for c = 1%:4
            y = eval(['TG.Keff' collection{c}]);
            y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
            pl(c) = plot(ax(a), x,smooth(y1,5), 'color', col(c,:), 'Linewidth', 2);;
            %plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');
            txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
                        'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
               txt.Color            = col(c,:);
         end
               
         yl = [5e-6 3e-4];
         ylim(ax(a), yl);
         ylabel(ax(a), 'm^2/s')
               
       % a=a+1;
       % collection{1} = '_int';
       % collection{2} = '_100';
       % collection{3} = '_50';
       % collection{4} = '_30';

       % label{1} = 'all interior';
       % label{2} = 'ms100';
       % label{3} = '50-40m';
       % label{4} = '<35m';

       % col([2:4],:) = col([5:7],:);

       % for c = 1:4
       %    y  = eval(['TG.Keff' collection{c}]);
       %    y1 = eval(['weighted_mean( TG.N' collection{c} ', y)']);
       %    pl(c) = plot(ax(a), x, smooth(y1,5), 'color', col(c,:), 'Linewidth', 1);
       %   % plot(ax(a), x, nanmean(y), '--', 'color', col(c,:), 'Linewidth', 1);
       %    set(ax(a), 'Yscale', 'log');
       %    yl = [5e-8 1e-5];
       %    ylim(ax(a), yl);
       %    txt = text( 0 + .2*c, 1, [label{c} ': ' num2str(nanmean(y1), '%3.1e')], ...
       %                'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
       %       txt.Color            = col(c,:);
       % end
       % yl = [5e-6 1e-3];
       % ylim(ax(a), yl);
       % ylabel(ax(a), 'K_{eff} m^2/s')
end
            
%}}}


%_____________________Fe relative to 50______________________{{{
if 0
         a=a+1;
        
      %for i = 1:length(i_ms)
      %    plot(ax(a), M(i_ms(i)).Ef.time, M(i_ms(i)).Ef.int_uppp, 'color', col(i,:),  'Linewidth', 1);
      %end
      anything = 'Ef.int_uppp';
      prefix   =  'log10(real(';
      suffix   =  '))';
       [any_time, any_avg_50, ~] = mooring_average_anything(M, IMs{2}, anything, prefix, suffix);
         any_avg_50 = log10(qbutter(10.^any_avg_50,1/3));

       for i = 1:length(IMs);
          [any_time, any_avg, any_std] = mooring_average_anything(M, IMs{i}, anything, prefix, suffix);
         %if(i>1)
         % ii_notnan = find(~isnan(any_avg));
         % any_time   = any_time(ii_notnan);
         % any_avg = any_avg(ii_notnan);
         % any_std    = any_std(ii_notnan);
         % px = [any_time fliplr(any_time)];
         % py = 10.^[ movmean(any_avg-any_std,3) , movmean(fliplr(any_avg+any_std),3)];
         % patch(px, py, col(i,:), 'facealpha', .5, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
         % end
         pl(i) =  plot(ax(a), any_time, 10.^any_avg./10.^any_avg_50, 'color', col(i,:),  'Linewidth', 2);
       end
       xlim(ax(a), xl);
       ylabel(ax(a), 'F_E/F_E^{OC50}');
       set(ax(a), 'Yscale', 'log');
       yl = [1e-2 1e2];
       ylim(ax(a), yl);
end
%}}}
       
%_____________________eps BBL______________________{{{
if 0
         a=a+1;
        
      %for i = 1:length(i_ms)
      %    plot(ax(a), M(i_ms(i)).Ef.time, M(i_ms(i)).Ef.int_uppp, 'color', col(i,:),  'Linewidth', 1);
      %end
      anything = 'Taub.epsb';
      prefix   =  'log10(real(';
      suffix   =  '))';

       for i = 1:length(IMs);
          [any_time, any_avg, any_std] = mooring_average_anything(M, IMs{i}, anything, prefix, suffix);
          if(i>1)
           ii_notnan = find(~isnan(any_avg));
           any_time   = any_time(ii_notnan);
           any_avg = any_avg(ii_notnan);
           any_std    = any_std(ii_notnan);
           px = [any_time fliplr(any_time)];
           py = 10.^[ movmean(any_avg-any_std,3) , movmean(fliplr(any_avg+any_std),3)];
           patch(px, py, col(i,:), 'facealpha', .2, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
           end
         pl(i) =  plot(ax(a), any_time, qbutter(10.^any_avg,1/3), 'color', col(i,:),  'Linewidth', 2);
       end
       xlim(ax(a), xl);
       ylabel(ax(a), '\epsilon_{bbl} [m^2/s^3]');
       set(ax(a), 'Yscale', 'log');
       yl = [6e-8 3e-6];
       ylim(ax(a), yl);
end
%}}}

%_____________________eps BBL relative to 50______________________{{{
if 0
         a=a+1;
        
      %for i = 1:length(i_ms)
      %    plot(ax(a), M(i_ms(i)).Ef.time, M(i_ms(i)).Ef.int_uppp, 'color', col(i,:),  'Linewidth', 1);
      %end
      anything = 'Taub.epsb';
      prefix   =  'log10(real(';
      suffix   =  '))';

       [any_time, any_avg_50, ~] = mooring_average_anything(M, IMs{2}, anything, prefix, suffix);
         any_avg_50 = log10(qbutter(10.^any_avg_50,1/3));
       for i = 1:length(IMs);
          [any_time, any_avg, any_std] = mooring_average_anything(M, IMs{i}, anything, prefix, suffix);
         %if(i>1)
         % ii_notnan = find(~isnan(any_avg));
         % any_time   = any_time(ii_notnan);
         % any_avg = any_avg(ii_notnan);
         % any_std    = any_std(ii_notnan);
         % px = [any_time fliplr(any_time)];
         % py = 10.^[ movmean(any_avg-any_std,3) , movmean(fliplr(any_avg+any_std),3)];
         % patch(px, py, col(i,:), 'facealpha', .2, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
         % end
         pl(i) =  plot(ax(a), any_time, qbutter(10.^any_avg,1/3)./10.^any_avg_50, 'color', col(i,:),  'Linewidth', 2);
       end
       xlim(ax(a), xl);
       ylabel(ax(a), '\epsilon_{bbl}/\epsilon_{bbl}^{OC50}');
       set(ax(a), 'Yscale', 'log');
       yl = [1e-1 1e1];
       ylim(ax(a), yl);
end
%}}}
      
%_____________________eps int______________________{{{
if 0
         a=a+1;
        
      %for i = 1:length(i_ms)
      %    plot(ax(a), M(i_ms(i)).Ef.time, M(i_ms(i)).Ef.int_uppp, 'color', col(i,:),  'Linewidth', 1);
      %end
      anything = 'G.eps(find(M(m).G.mab>1),:)';
      prefix   =  'log10(nanmean(';
      suffix   =  ',1))';
      

       for i = 1:length(IMs);

          i_ms_gust = [];
          for j = 1:length(IMs{i})
             if sum(M(IMs{i}(j)).G.mab>1) >1
                i_ms_gust = [i_ms_gust  IMs{i}(j)];
             end

          end

          
          [any_time, any_avg, any_std] = mooring_average_anything(M, i_ms_gust, anything, prefix, suffix);
          if(i>1)
           ii_notnan = find(~isnan(any_avg));
           any_time   = any_time(ii_notnan);
           any_avg = any_avg(ii_notnan);
           any_std    = any_std(ii_notnan);
           px = [any_time fliplr(any_time)];
           py = 10.^[ movmean(any_avg-any_std,3) , movmean(fliplr(any_avg+any_std),3)];
           patch(px, py, col(i,:), 'facealpha', .2, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
           end
         pl(i) =  plot(ax(a), any_time, movmean(10.^any_avg,3), 'color', col(i,:),  'Linewidth', 2);
       end
       xlim(ax(a), xl);
       ylabel(ax(a), '\epsilon_{int} [m^2/s^3]');
       set(ax(a), 'Yscale', 'log');
       yl = [2e-9 3e-7];
       ylim(ax(a), yl);
end
%}}}

%_____________________eta______________________{{{
if 0
         a=a+1;
      anything = 'Ef.eta';
     % anything = 'T.eta_IT';
      prefix   =  '';
      suffix   =  '';
       for i = 1:length(IMs);
          [any_time, any_avg, any_std] = mooring_average_anything(M, IMs{i}, anything, prefix, suffix);
           ii_notnan = find(~isnan(any_avg));
           any_time   = any_time(ii_notnan);
           any_avg = any_avg(ii_notnan);
           any_std    = any_std(ii_notnan);
           px = [any_time fliplr(any_time)];
           py = [ movmean(any_avg-any_std,3) , movmean(fliplr(any_avg+any_std),3)];
           patch(px, py, col(i,:), 'facealpha', .2, 'edgecolor', 'none', 'Linewidth', 1, 'parent', ax(a));
          plot(ax(a), any_time, any_avg, 'color', col(i,:),  'Linewidth', 2);
       end
       xlim(ax(a), xl);
       ylabel(ax(a), '\eta_0 [m]');
       yl = [0 15];
       ylim(ax(a), yl);
end
%}}}

       datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');

      % legend(pl, 'ms100', '50m', '40m', '<=30m')


       abc='abcdefghijklmnopqrst';
       for a = 1:(size(ax,1)*size(ax,2))
          tabc = text_corner(ax(a), abc(a), 7);
          tabc.BackgroundColor = [1 1 1 .5];
          tabc.FontWeight = 'bold';
          if a~=2
             set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
          else
            set(ax(a), 'TickDir', 'out', 'Layer', 'top');
          end

          xlim(ax(a), xl);
       end
       


       print(gcf,'../pics/turb_vs_env.png','-dpng','-r200','-painters')
       



%_____________________Correlation plot______________________


 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20]*.7,'PaperPosition',[0 0 30 20]*.7);
 
         [ax, ~] = create_axes(fig, 3, 3, 0);
         

         a=0
         dt = .7;

          xt_Fe  = Fe.time{2};
          x_Fe   = log10(real(Fe.Fe{2}));
          xl_Fe  = ([5 100]);;

         xt_wind = TG.time(1):.5*nanmean(diff(TG.time)):TG.time(end);
         %x_wind  = log10(abs(clever_interp(Met.Van.M.time, Met.U, xt_wind)).^2);
         x_wind  = abs(clever_interp(Met.Van.M.time, Met.U, xt_wind));
         xl_wind = ([0 15]);


         xt_DT = TS.time;
         x_DT  = diff(TS.Ts([find(TS.zs>-40,1) find(TS.zs>-5,1)],:));
         xl_DT = ([1 4]);
         

      %_____________________Fe correlation______________________
            % Fe
            xt =  xt_Fe;
            x  =  x_Fe;
            xl = xl_Fe;

         a=a+1;

         % chi
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.chi));
             yl = ([8e-8 3e-6]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;

               marker_DT     = clever_interp(xt_DT, x_DT, time);
               marker_wind   = clever_interp(xt_wind, x_wind, time);

            scatter(ax(a), 10.^(Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
            ii_what = find(marker_wind>1.3);;
            scatter(ax(a), 10.^(Xg(ii_what)), 10.^Yg(ii_what), 16, [.5 0 0], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
               ylabel(ax(a),'\chi  K^2/s')
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f   %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);
         
             xfit = log10(xl(1)):(diff(log10(xl))*.01):log10(xl(2));
             [pfit,S,mu] = polyfit((Xg), Yg,1);
             [yfit] = polyval(pfit,xfit);
             plot(ax(a), xfit, yfit, 'Linewidth', 1);

             
             
         a=a+1;
         % Phi_d
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Phi_d));
            y(isinf(y))= nan;
             yl = ([8e-7 2e-5]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), 10.^(Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
            scatter(ax(a), 10.^(Xg(ii_what)), 10.^Yg(ii_what), 16, [.5 0 0], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
               ylabel(ax(a),'\Phi_d  K m/s')
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

         a=a+1;
         % Keff
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Keff));
            y(isinf(y))= nan;
             yl = ([8e-6 8e-4]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), 10.^(Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
            scatter(ax(a), 10.^(Xg(ii_what)), 10.^Yg(ii_what), 16, [.5 0 0], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
               xlabel(ax(a),'F_E   W/m')
               ylabel(ax(a),'K_{eff}  m^2/s')
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

      %_____________________wind cor______________________
        xt =  xt_wind;
        x = x_wind;
       xl = xl_wind; 

         a=a+1;
            % chi
            yt = TG.time;
            y = log10(weighted_mean( TG.N, TG.chi));
                yl = ([8e-8 3e-6]);
            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), (Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

         a=a+1;
         % Phi_d
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Phi_d));
            y(isinf(y))= nan;
             yl = ([8e-7 2e-5]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a),(Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

         a=a+1;
         % Keff
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Keff));
            y(isinf(y))= nan;
             yl = ([8e-6 8e-4]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), (Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

               xlabel(ax(a),'(wind speed)^2   m^2/s^2')

      %_____________________DT cor______________________
            xt = xt_DT;
            x  = x_DT;
            xl = xl_DT;

         a=a+1;
            % chi
            yt = TG.time;
            y = log10(weighted_mean( TG.N, TG.chi));
                yl = ([8e-8 3e-6]);
            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), (Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

         a=a+1;
         % Phi_d
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Phi_d));
            y(isinf(y))= nan;
             yl = ([8e-7 2e-5]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), (Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

         a=a+1;
         % Keff
         yt = TG.time;
         y = log10(weighted_mean( TG.N, TG.Keff));
            y(isinf(y))= nan;
             yl = ([8e-6 8e-4]);

            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            scatter(ax(a), (Xg), 10.^Yg, 16, [0 0 .5], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log');
             text_corner(ax(a), ['r = ' num2str(r,'%1.2f') ' ' num2str(rb, '[%0.2f - %0.2f]') ], 1);
             text_corner(ax(a), ['p = ' num2str(p,'%1.4f')], 3);

               xlabel(ax(a),'{\Delta}T  K')


               abc='abcdefghijklmnopqrst';
               for a = 1:(size(ax,1)*size(ax,2))
                  tabc = text_corner(ax(a), abc(a), 7);
                  tabc.BackgroundColor = [1 1 1 .5];
                  set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
               end
               

            print(gcf,'../pics/turb_cor.png','-dpng','-r200','-painters')
            

% Fe vs DT

 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[20 10]*.8,'PaperPosition',[0 0 20 10]*.8)

         [ax, ~] = create_axes(fig, 1, 2, 0);
            squeeze_axes(ax, .95 ,.95);
            shift_axes(ax,0 ,.04);
            shift_axes(ax(2), .05,0);
         

         %{{{ Fe vs DT
         
         a=1;
         
 
            xt = TS.time;
            x = diff(TS.Ts([find(TS.zs>-40,1) find(TS.zs>-5,1)],:));
               xl = ([1 4]);

               yl = ([1e-1 1e3]);
            col = get(groot,'DefaultAxesColorOrder');
             text( xl(1), yl(2), ['  r = '], ...
                      'verticalalignment', 'top', 'horizontalalignment', 'left','Parent', ax(a));

            Label{1} = 'ms100';
            Label{2} = '50m';
            Label{3} = '40m';
            Label{4} = '<35m';

        for i=1:4 
            yt = Fe.time{i};
            %y = log10(real(Fe.Fe{i}));
            y = (real(Fe.Fe{i}));
            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            %pl(i) = scatter(ax(a), (Xg), Yg, 16, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
            pl(i) = scatter(ax(a), (Xg), Yg, 16, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             txt = text( i*1/4 -.05, 1, [Label{i}], ...
                      'units', 'normalized','verticalalignment', 'bottom', 'horizontalalignment', 'center','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      txt.Color            = col(i,:)*.8;
             txt = text( i*1/4 -.05, 1, [num2str(r, '%1.2f')], ...
                      'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      %txt.FontWeight      = 'bold';
                      %txt.FontSize         = 12;
                      txt.Color            = col(i,:)*.8;
                      
             txt = text( i*1/4 -.05, .92, num2str(rb, '[%0.2f - %0.2f]') , ...
                      'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      %txt.FontWeight      = 'bold';
                      txt.FontSize         = 7;
                      txt.Color            = col(i,:)*.8;
          end
          plot(ax(a), xl, xl.^2,'k', 'Linewidth', 1);
          plot(ax(a), xl, xl,'k', 'Linewidth', 1);
          t = text( xl(2), xl(2), ['{\Delta}T'], ...
                  'verticalalignment', 'middle', 'horizontalalignment', 'left','Parent', ax(a));
          t = text( xl(2), xl(2)^2, ['{\Delta}T^2'], ...
                  'verticalalignment', 'middle', 'horizontalalignment', 'left','Parent', ax(a));
          
          

             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             ylabel(ax(a),'F_E   W/m')

             xlabel(ax(a),'{\Delta}T  K')

             %}}}


      %{{{ Fe vs neap spring

         a=2;
            xt_T = TS.time;
            x_T = diff(TS.Ts([find(TS.zs>-40,1) find(TS.zs>-5,1)],:));
            xt =  P.time;
            x = qbutter(movmax(P.press-nanmedian(P.press),60*36) - movmin(P.press-nanmedian(P.press),60*36), 1/(60*36));
             xl = ([1 2]);
             dt =1;
            yl = ([1e-1 1e3]);

        for i=1:4 
            yt = Fe.time{i};
            y = (real(Fe.Fe{i}));
            logit = 0;
            plot_correlation( ax(a), xt, x, yt, y, dt, Label{1}, i, 4, xl, yl, logit,  col, 16 )
        end
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             %ylabel(ax(a),'F_E   W/m')
             xlabel(ax(a),'max-min surface [m]')


          %}}}


             print(gcf,'../pics/FeVsDT.png','-dpng','-r200','-painters')


%_____________________Conditional sampling______________________%{{{
 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20]);
 

         [ax, ~] = create_axes(fig, 3, 3, 0);
            squeeze_axes(ax, .92 ,1);
            for a = 1:numel(ax)
               squeeze_axes(ax(a), .9 ,.9);
            end

            shift_axes(ax, .02, 0);
            shift_axes(ax(7:9), .03, 0);

            a=0;
         


      % X axis
            Label{1} = 'ms100';
            Label{2} = '50m';
            Label{3} = '40m';
            Label{4} = '<35m';
            xt_T = TS.time;
            x_T = diff(TS.Ts([find(TS.zs>-40,1) find(TS.zs>-5,1)],:));
            xt =  P.time;
            x = qbutter(movmax(P.press-nanmedian(P.press),60*36) - movmin(P.press-nanmedian(P.press),60*36), 1/(60*36));
             xl = ([1 2]);
             dt =1;
            yl = ([1e-1 1e3]);


      %{{{

         a=a+1;

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            logit = 2;
            plot_correlation( ax(a), xt, x, yt, y, dt, Label{i}, i, 4, xl, yl, logit,  col, 16 )
        end
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             ylabel(ax(a),'F_E   W/m')
             %xlabel(ax(a),'max-min surface [m]')

          %}}}
          txt = text_corner(ax(a), ['all'], 7);
            txt.BackgroundColor = [1 1 1];
            txt.FontWeight      = 'bold';
            txt.FontSize         = 12;
            

          thresh_wind = 2;
      %{{{ 
         a=a+1;

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            logit = 2;
            plot_correlation( ax(a), xt, x, yt, y, dt, '', i, 4, xl, yl, logit,  col, 16, ...
                  xt_wind, x_wind > thresh_wind)
        end
             ylabel(ax(a),' F_E   W/m')
             %xlabel(ax(a),'max-min surface [m]')

          %}}}
          txt = text_corner(ax(a), ['wind > ' num2str(thresh_wind, '%2.1f') 'm/s'], 7);
            txt.BackgroundColor = [1 1 1];
            txt.FontWeight      = 'bold';
            txt.FontSize         = 12;
         
      %{{{  

         a=a+1;

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            logit = 2;
            plot_correlation( ax(a), xt, x, yt, y, dt, '', i, 4, xl, yl, logit,  col, 16, ...
                  xt_wind, x_wind < thresh_wind)
        end
             ylabel(ax(a),' F_E   W/m')
             xlabel(ax(a),'max-min surface [m]')

          %}}}
          txt = text_corner(ax(a), ['wind < ' num2str(thresh_wind, '%2.1f') 'm/s'], 7);
            txt.BackgroundColor = [1 1 1];
            txt.FontWeight      = 'bold';
            txt.FontSize         = 12;


          thresh_DT = 2.5;
      %{{{ 
         a=a+1;
         a=a+1;

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            logit = 2;
            plot_correlation( ax(a), xt, x, yt, y, dt, '', i, 4, xl, yl, logit,  col, 16, ...
                  xt_DT, x_DT > thresh_DT)
        end
            % ylabel(ax(a),' F_E   W/m')
             %xlabel(ax(a),'max-min surface [m]')

          %}}}
          txt = text_corner(ax(a), ['\Delta T > ' num2str(thresh_DT, '%2.1f') 'K'], 7);
            txt.BackgroundColor = [1 1 1];
            txt.FontWeight      = 'bold';
            txt.FontSize         = 12;
         
      %{{{  

         a=a+1;

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            logit = 2;
            plot_correlation( ax(a), xt, x, yt, y, dt, '', i, 4, xl, yl, logit,  col, 16, ...
                  xt_DT, x_DT < thresh_DT)
        end
            % ylabel(ax(a),' F_E   W/m')
             xlabel(ax(a),'max-min surface [m]')

          %}}}
          txt = text_corner(ax(a), ['\Delta T < ' num2str(thresh_DT, '%2.1f') 'K'], 7);
            txt.BackgroundColor = [1 1 1];
            txt.FontWeight      = 'bold';
            txt.FontSize         = 12;

      %{{{  

         a=a+1;

        for i=1 
            logit = 0;
            plot_correlation( ax(a), xt, x, xt_wind, x_wind, dt, '', i, 4, xl, [0 10], logit,  hot(4), 16)
        end
            ylabel(ax(a),' wind m/s')
             xlabel(ax(a),'max-min surface [m]')

          %}}}

         print(gcf,'../pics/cont_sample.png','-dpng','-r200','-painters')

         
         %}}}
return
             
%_Fe vs DT2

 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[12 10]*.8,'PaperPosition',[0 0 12 10]*.8)
         
         a=1;
            ax(a) = subplot( 1, 1, a);
            hold(ax(a), 'on');
         
 
            xt = TS.time;
            x = log10((diff(TS.Ts([find(TS.zs>-40,1) find(TS.zs>-5,1)],:))).^2);
               xl = ([1e-1 1e1]);

               yl = ([1e-1 1e3]);
            col = get(groot,'DefaultAxesColorOrder');
             text( xl(1), yl(2), ['  r = '], ...
                      'verticalalignment', 'top', 'horizontalalignment', 'left','Parent', ax(a));

            Label{1} = 'ms100';
            Label{2} = 'H=50m';
            Label{3} = ' 40m';
            Label{4} = '<35m';

        for i=1:4 
            yt = Fe.time{i};
            y = log10(real(Fe.Fe{i}));
            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;
            pl(i) = scatter(ax(a), 10.^(Xg), 10.^Yg, 16, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             txt = text( xl(1)+i*diff(xl)/5, yl(2), [Label{i}], ...
                      'verticalalignment', 'bottom', 'horizontalalignment', 'left','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      txt.Color            = col(i,:)*.8;
             txt = text( xl(1)+i*diff(xl)/5, yl(2), [num2str(r, '%1.2f')], ...
                      'verticalalignment', 'top', 'horizontalalignment', 'left','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      %txt.FontWeight      = 'bold';
                      %txt.FontSize         = 12;
                      txt.Color            = col(i,:)*.8;
                      
          end

             ylim(ax(a), yl);
             xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             ylabel(ax(a),'F_E   W/m')

             xlabel(ax(a),'({\Delta}T  K)^2')
