close all; 
warning off;
%_____________________load mooring data______________________
if 0; %{{{;
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


 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[20 20],'PaperPosition',[0 0 20 20]);
 
   [ax, ~] = create_axes(fig, 3, 1, 0);
   squeeze_axes(ax, 1, .6)
   shift_axes(ax, 0, .35);
   
   
  m=2;


  a=1;
  %plot(ax(a), M(m).G.time, max(movmean(M(m).G.chi, 5,2, 'omitnan')), 'Linewidth', 1);
  %plot(ax(a), M(m).G.time, min(movmean(M(m).G.chi, 5,2, 'omitnan')), 'Linewidth', 1);
  plot(ax(a), M(m).G.time, nanmean(movmean(M(m).G.chi, 5,2, 'omitnan')), 'Linewidth', 1);
  %plot(ax(a), M(m).G.time, nanmedian(movmean(M(m).G.chi, 5,2, 'omitnan')), 'Linewidth', 1);
   plot(ax(a), M(m).G.time, nanmean(movmean(M(m).G.chi, 20,2, 'omitnan')), 'Linewidth', 3);
   set(ax(a), 'Yscale', 'log');
   txt = text_corner(ax(a), ['\langle\chi\rangle_z @ OC50'], 1);
      txt.BackgroundColor = [1 1 1 .5];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
   ylabel(ax(a), ['K^2/s']);
   
   

  a=2;
  X = nanmean(M(m).T.T);
  DX = max(M(m).T.T)-min(M(m).T.T);
   plot(ax(a), M(m).T.time, movstd(X,50 )./movmean(DX,50), 'Linewidth', 1);
   txt = text_corner(ax(a), ['std{\langle}T{\rangle_{z}}/{\Delta}T (non-dim IW amplitude)'], 1);
      txt.BackgroundColor = [1 1 1 .5];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
   ylabel(ax(a), ['']);

  a=3;
   plot(ax(a), M(m).T.time, nanmean(movmean(M(m).T.T, 5,2, 'omitnan')), 'Linewidth', 1);
   txt = text_corner(ax(a), ['T'], 1);
      txt.BackgroundColor = [1 1 1 .5];
      txt.FontWeight      = 'bold';
      txt.FontSize         = 12;
      txt.Color            = [0 0 0];
   ylabel(ax(a), ['{^\circ}C']);
   

   linkaxes(ax, 'x');

   for a=1:3
      xl = [0 4]+10 + datenum(2017,9,12);
      xlim(ax(a), xl);
   end
   
   datetick(ax(a),  'x', 'mmm-dd',  'keeplimits');
   
   
   
 

[ax(4:5), ~] = create_axes(fig, 1, 2, 0);
squeeze_axes(ax(4:5), .9, .3)
shift_axes(ax(5), 0.05, 0);

%_____________________correlation______________________
a=4;
dt=1/24/2;
xt = M(m).T.time+10/60/24;
x = movstd(X,50 )./movmean(DX,50);
%x = movstd(X,100 );
yt = M(m).G.time;
y = nanmean(movmean(M(m).G.chi, 10,2, 'omitnan'));
xl = [1e-3 2e-1];
yl = [1e-8 1e-4];


            [Xg, Yg, time, r, p, rb] = xy_correlate( xt, log10(x), yt, log10(y), dt) ;
            scatter(ax(a), 10.^Xg, 10.^Yg, 16, [0 .5 0], 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
             txt = text( 0, 1, ['r = ' num2str(r, '%1.2f') num2str(rb, ' [%1.2f   %1.2f]')], 'units', 'normalized', ...
                      'verticalalignment', 'top', 'horizontalalignment', 'left','Parent', ax(a));
                      txt.BackgroundColor = [1 1 1 0];
                      

            ylim(ax(a), yl);
            xlim(ax(a), xl);
             set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
             %set(ax(a), 'Yscale', 'log');
             ylabel(ax(a),'\langle \chi \rangle_z')

             xlabel(ax(a),'std T/ \Delta T')

   
%_____________________sorted average______________________

[Xg_sort, ii_sort] = sort(Xg, 'descend');
ii_sort = ii_sort(~isnan(Xg_sort) & ~isnan(Yg(ii_sort)));
Xg_sort = Xg(ii_sort);

plot(ax(a), 10.^Xg_sort, cumsum(10.^Yg(ii_sort))./cumsum(ones(length(ii_sort),1)), 'k', 'Linewidth', 3);

%_____________________sorted average______________________
a=5;

 
plot(ax(a), [1:length(ii_sort)]/length(ii_sort), cumsum((10.^Yg(ii_sort)))/nansum(10.^Yg(ii_sort)), 'k', 'Linewidth', 3);
plot(ax(a), [0 1], [0 1], 'Linewidth', 1);
plot(ax(a), [0 .5], [0 1], 'Linewidth', 1);
plot(ax(a), [0 .333], [0 1], 'Linewidth', 1);
             ylabel(ax(a),'\Sigma_{i=1}^{Isort} \chi_{i}/\Sigma_1^{N}\chi_{i}')
             xlabel(ax(a),'I_{sort}/N')

   abc='abcdefghijklmnopqrst';
   for a = 1:(size(ax,1)*size(ax,2))
      tabc = text_corner(ax(a), abc(a), 9);
      tabc.BackgroundColor = [1 1 1 .5];
      set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
   end
   
             
   print(gcf,'../pics/IW_vs_chi.png','-dpng','-r200','-painters')
