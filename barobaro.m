
%_____________________load mooring data______________________
if 1
clear all;
   % path to m_map here
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/ganges/work/git_software/mixingsoftware/general/'));
   addpath(genpath('~/ganges/work/git_software/mixingsoftware/seawater/'));
   addpath(genpath('~/arbeit/matlab_tbx/mixingsoftware/general/'));
   addpath(genpath('~/arbeit/matlab_tbx/mixingsoftware/seawater/'));
   
   
   load cmap;

%_____________________all epsilon data______________________
  % load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');
% wind
load /home/johannes/gdrive/IS17/collection_for_Jim/data/Mini_Met_Innershelf.mat;

Met = MiniMet;
Met.Van = load('/home/johannes/gdrive/IS17/collection_for_Jim/data/meteo_vandenberg.mat');

Met.sw  = clever_interp( Met.Van.M.time, Met.Van.M.radiation, Met.TimeOR2);
Met.speedOR2  = clever_interp( Met.TimeOR1, Met.WindSpeedOR1, Met.TimeOR2);
Met.SSTOR2  = clever_interp( Met.TimeSST, Met.SST, Met.TimeOR2);

    [Met.Jb,Met.Jq]=surfaceflux( Met.AirTempOR2, zeros(size(Met.AirTempOR2)), Met.AirPressureOR2 ...
                     ,Met.SSTOR2 , 35*ones(size(Met.TimeOR2)), Met.speedOR2 ...
                     , Met.RelativeHumidityOR2 ,   Met.sw, 0, 0);


% pressure record
    load ~/gdrive/IS17/tau_map/data/pressure_record_PS50.mat;
    P.press = movmean(P.press,10);
end


 %moorname = 'oc40s'
 %m = find_names_inM(M,moorname) ;

	[ms, mshalf] = getFullMList(M) 

  tl  = [datenum(2017,9,8) datenum(2017,11,2)];

 fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20])

for m = ms

         [ax, ~] = create_axes(fig, 5, 1, 0);
         
         flp = M(m).U.dt/(3600*36);
         fti = M(m).U.dt/(3600*6);
         yl = [-.2 .2];

         col = get(groot,'DefaultAxesColorOrder');
         
         Met.U = Met.WindSpeedOR1.*exp(1i*(-(Met.WindDirectionOR1-90))/180*pi);
         a=1;
         plot(ax(a), Met.TimeAvg, qbutter(Met.AvgWindSpeed,1), 'k', 'Linewidth', 1);
         plot(ax(a), Met.TimeAvg, qbutter(Met.AvgWindSpeed,1/36), 'k', 'Linewidth', 3);
         plot(ax(a), Met.TimeOR1, real(qbutter(Met.U,1/3000)), 'color', col(1,:), 'Linewidth', 1);
         plot(ax(a), Met.TimeOR1, imag(qbutter(Met.U,1/3000)), 'color', col(2,:), 'Linewidth', 1);
               xlim(ax(a), tl);

         a=2;
            plot(ax(a), M(m).U.time, real(qbutter( nanmean(M(m).U.U), flp )), 'Linewidth', 1);
            plot(ax(a), M(m).U.time, imag(qbutter( nanmean(M(m).U.U), flp )), 'Linewidth', 1);
               plot(ax(a), tl, [0 0] , '--k', 'Linewidth', 1);
               legend(ax(a), 'east', 'north');
               t = text_corner(ax(a), ['\langle u\rangle_z^{lp 36h} [m/s]'], 1);
               xlim(ax(a), tl);
					ylim(ax(a), yl);
					
               
         a=3;
         plot(ax(a), P.time, (P.press-nanmedian(P.press)),'k', 'Linewidth', 2);
            t = text_corner(ax(a), ['\eta_{surf} [m]'], 1);
               t.BackgroundColor = [1 1 1];
               t.FontWeight      = 'bold';
               xlim(ax(a), tl);

         a=4;
            plot(ax(a), M(m).U.time, real(qbutter( nanmean(M(m).U.U), fti, flp )), 'Linewidth', 1);
            plot(ax(a), M(m).U.time, imag(qbutter( nanmean(M(m).U.U), fti, flp )), 'Linewidth', 1);
               plot(ax(a), tl, [0 0] , '--k', 'Linewidth', 1);
               legend(ax(a), 'east', 'north');
               t = text_corner(ax(a), ['\langle u\rangle_z^{bp 6-36h} [m/s]'], 1);
               xlim(ax(a), tl);
					ylim(ax(a), yl);


         a=5;
            plot(ax(a), M(m).U.time, real(qbutter( nanmean(M(m).U.U(1:2,:)), fti)), 'Linewidth', 1);
            plot(ax(a), M(m).U.time, imag(qbutter( nanmean(M(m).U.U(1:2,:)), fti)), 'Linewidth', 1);
               plot(ax(a), tl, [0 0] , '--k', 'Linewidth', 1);
               legend(ax(a), 'east', 'north');
               t = text_corner(ax(a), ['u_{bottom}^{lp 6p} [m/s]'], 1);
               xlim(ax(a), tl);
					ylim(ax(a), yl);

      linkaxes(ax, 'x');
      datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
      t = text_corner(ax(1), [M(m).mooringName], -2);


      print(gcf,['../pics/more_all/' M(m).mooringName  '.png'],'-dpng','-r200','-painters')
      clf
      
      
end
      
            
            



