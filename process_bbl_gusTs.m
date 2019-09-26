clear all;
close all;

addpath(genpath('~/ganges/work/chipod_gust/software/'));

do_proc             = 0;
do_parallel         = 0;
do_combine_and_mask = 0;
do_plot             = 0;


Glist = readtable('../data/IS17_bottom_gusts.csv');
for g = 1:size(Glist,1)
   name{g}  = Glist{g,1};
   moor{g}  = Glist{g,2};
   dpath{g} = Glist{g,3};
   start(g) = datenum(Glist{g,4});
   stop(g)  = datenum(Glist{g,5});
end


% do the acrtual processing
if do_proc
   for g = 1:length(name)
      g

      time_range = [start(g) stop(g)];
      basedir = ['/home/johannes/ganges/' char(dpath{g})];
      load([basedir '/input/vel_m.mat']);
      mab = 2; % actual height of the adcp-bin ???

       S.time  =  vel_m.time;
       S.spd   =  vel_m.spd;
       Eps.time   =  vel_m.time;
       Eps.eps    =  (2.3e-3*abs(vel_m.spd).^2).^1.5/.4/mab;
       generate_eps_chi(basedir, Eps, S,  do_parallel,  time_range, '_ISmix');
   end
end



if do_combine_and_mask

  for g = 1:length(name) 
      fid = ['/home/johannes/ganges/' char(dpath{g}) '/proc/eps_chi_ISmix.mat'];
      load(fid);
      
      eps_min = 1e-8;
      spd_min = 0.05;
      spec_floor = .5*nanmean(chi.spec_floor)*nanmean(chi.nfft);
      mask = ~(chi.spec_area< spec_floor | chi.eps<eps_min | chi.spd<.05 ...
               | chi.time<start(g) | chi.time>stop(g));

      [G_bbl(g).time, G_bbl(g).chi] = average_60sec( chi.time, chi.chi, mask);
      [~, G_bbl(g).eps] = average_60sec( chi.time, chi.eps, mask);
      [~, G_bbl(g).spd] = average_60sec( chi.time, chi.spd, ones(size(chi.time)));
      [~, G_bbl(g).T] = average_60sec( chi.time, chi.T, ones(size(chi.time)));

      G_bbl(g).GusT_name = char(name{g});
      G_bbl(g).Moor_name = char(moor{g});
      G_bbl(g).depth = nanmedian(chi.depth);



       
       if do_plot
           fig = figure('Color',[1 1 1],'visible','off','Paperunits','centimeters',...
                   'Papersize',[30 20],'PaperPosition',[0 0 30 20])
           [ax, ~] = create_axes(fig, 5, 1, 0);
           
           a=1;
           plot(ax(a), chi.time, chi.chi, 'color', [.6 .6 .6],  'Linewidth', 1);
           plot(ax(a), chi.time(mask), chi.chi(mask), 'color', [0 0 0],  'Linewidth', 1);
           plot(ax(a), G_bbl(g).time, G_bbl(g).chi, 'color', [1 0 0],  'Linewidth', 2);
           ylabel(ax(a), '\chi');
           set(ax(a), 'Yscale', 'log');
           yl = [1e-10 1e-3];
           ylim(ax(a), yl);
           

           a=2;
           plot(ax(a), chi.time, chi.eps, 'color', [.6 .6 .6],  'Linewidth', 1);
           plot(ax(a), chi.time(mask), chi.eps(mask), 'color', [0 0 0],  'Linewidth', 1);
           plot(ax(a), G_bbl(g).time, G_bbl(g).eps, 'color', [1 0 0],  'Linewidth', 2);
           plot(ax(a), chi.time([1 end]), [1 1]*eps_min , 'Linewidth', 1);
           
           ylabel(ax(a), '\epsilon');
           set(ax(a), 'Yscale', 'log');
           yl = [1e-10 1e-3];
           ylim(ax(a), yl);

           
           a=3;
           plot(ax(a), chi.time, chi.spd, 'color', [.6 .6 .6],  'Linewidth', 1);
           plot(ax(a), chi.time(mask), chi.spd(mask), 'color', [0 0 0],  'Linewidth', 1);
           plot(ax(a), chi.time([1 end]), [1 1]*spd_min , 'Linewidth', 1);
           plot(ax(a), G_bbl(g).time, G_bbl(g).spd, 'color', [1 0 0],  'Linewidth', 2);
           ylabel(ax(a), 'speed');
           yl = [0 .4];
           ylim(ax(a), yl);

           a=4;
           plot(ax(a), chi.time, chi.T, 'color', [.6 .6 .6],  'Linewidth', 1);
           plot(ax(a), chi.time(mask), chi.T(mask), 'color', [0 0 0],  'Linewidth', 1);
           plot(ax(a), G_bbl(g).time, G_bbl(g).T, 'color', [1 0 0],  'Linewidth', 2);
           ylabel(ax(a), 'T');
           yl = [9 16];
           ylim(ax(a), yl);

           a=5;
           plot(ax(a), chi.time, chi.spec_area, 'color', [.6 .6 .6],  'Linewidth', 1);
           plot(ax(a), chi.time(mask), chi.spec_area(mask), 'color', [0 0 0],  'Linewidth', 1);
           plot(ax(a), chi.time([1 end]), [1 1]*spec_floor , 'Linewidth', 1);
           ylabel(ax(a), 'var(T)');
           set(ax(a), 'Yscale', 'log');
           %yl = [1e-10 1e-3];
           %ylim(ax(a), yl);


           datetick(ax(a), 'x', 'mmm-dd',  'keeplimits');
           
           linkaxes(ax, 'x');

           t = text_corner(ax(1), ['\epsilon_{bbl} - masking for ' char(moor{g}) ], -2);
           t.BackgroundColor = [1 1 1];
           t.FontWeight      = 'bold';
           t.FontSize         = 12;
           t.Color            = [0 0 0];
           
           

           abc='abcdefghijklmnopqrst';
           for a = 1:(size(ax,1)*size(ax,2))
              tabc = text_corner(ax(a), abc(a), 7);
              tabc.BackgroundColor = [1 1 1 .5];
              set(ax(a), 'box', 'on', 'TickDir', 'out', 'Layer', 'top');
              plot(ax(a), [1 1]*start(g), get(ax(a),'Ylim'), 'Linewidth', 1);
              plot(ax(a), [1 1]*stop(g), get(ax(a),'Ylim'), 'Linewidth', 1);
              
           end
           
           
           print(gcf,['../pics/quality_eps_bbl_mask/' char(moor{g}) '.png'],'-dpng','-r200','-painters')
           
              
       end

   end

   save('../data/G_bbl.mat', 'G_bbl');
      
end

% filter scheme
%
%  mask = (chi.spec_area<.5*nanmean(chi.spec_floor)*nanmean(chi.nfft) | chi.eps<1e-8 | chi.spd<.1);


function [time_avg, y_avg] = average_60sec( time, y, mask); 
      DT = 1/24/60; % 60 sec averages
      dt = diff(time([1 2]));
      time_avg = time(1):DT:time(end);
      ww = round(DT/dt);

      % average down to 60 sec
      tmp        = y;
      tmp(~mask) = nan;
      tmpavg = movmean( tmp, ww, 'omitnan');
      tmpnan = movsum( mask,  ww, 'omitnan');
      tmpavg(tmpnan<3) = nan; % set all averages to nan that contain less than 3 points;
      y_avg = interp1( time, tmpavg, time_avg);
end
      
