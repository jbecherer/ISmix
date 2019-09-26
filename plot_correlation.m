function [] = plot_correlation( ax, xt, x, yt, y, dt, Label, i, N, xl, yl, logit,  col, Msize, ct, c )



   [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) ;

   % condition

   if nargin > 15
    true_cond = (interp1( ct, double(c), time, 'nearest') ==1);
   else
      true_cond = ~isnan(Yg);
   end
   if size(Yg,1) ~= size(true_cond,1)
       true_cond = true_cond';
   end

   tt = find(~isnan(Xg)&~isnan(Yg) & true_cond);
   [R, P, rl, ru] = corrcoef( Xg(tt), Yg(tt));
   r  = R(1,2);
   p  = P(1,2);
   rb = [rl(1,2) ru(1,2)];


   if logit == 3
      scatter(ax, 10.^Xg(tt), 10.^Yg(tt), Msize, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
      plot(ax, xl, [1 1]*nanmean(10.^Yg(tt)), '--', 'color', col(i,:), 'Linewidth', 1);
   elseif logit == 2
      scatter(ax, Xg(tt), 10.^Yg(tt), Msize, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
      plot(ax, xl, [1 1]*nanmean(10.^Yg(tt)), '--', 'color', col(i,:), 'Linewidth', 1);
   elseif logit == 1
      scatter(ax, 10.^Xg(tt), Yg(tt), Msize, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
      plot(ax, xl, [1 1]*nanmean(Yg(tt)), '--', 'color', col(i,:), 'Linewidth', 1);
   else
      scatter(ax, Xg(tt), Yg(tt), Msize, col(i,:), 'Marker', 'o','MarkerFaceColor', [.5 .5 .5], 'Linewidth', 2);
      plot(ax, xl, [1 1]*nanmean(Yg(tt)), '--', 'color', col(i,:), 'Linewidth', 1);
   end

            if ~isempty(Label)
             text( xl(1), yl(2), ['  r = '], ...
                      'verticalalignment', 'top', 'horizontalalignment', 'left','Parent', ax);
            end
             txt = text( i*1/N -.05, 1, [Label], ...
                      'units', 'normalized','verticalalignment', 'bottom', 'horizontalalignment', 'center','Parent', ax);
                      txt.BackgroundColor = [1 1 1 0];
                      txt.Color            = col(i,:)*.8;
             txt = text( i*1/N -.05, 1, [num2str(r, '%1.2f')], ...
                      'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax);
                      txt.BackgroundColor = [1 1 1 0];
                      %txt.FontWeight      = 'bold';
                      %txt.FontSize         = 12;
                      txt.Color            = col(i,:)*.8;
                      
             txt = text( i*1/N -.05, .92, num2str(rb, '[%0.2f - %0.2f]') , ...
                      'units', 'normalized', 'verticalalignment', 'top', 'horizontalalignment', 'center','Parent', ax);
                      txt.BackgroundColor = [1 1 1 0];
                      %txt.FontWeight      = 'bold';
                      txt.FontSize         = 7;
                      txt.Color            = col(i,:)*.8;


             ylim(ax, yl);
             xlim(ax, xl);

          plot(ax, xl, xl.^2,'k', 'Linewidth', 1);
          plot(ax, xl, xl,'k', 'Linewidth', 1);
          t = text( xl(2), xl(2), ['x'], ...
                  'verticalalignment', 'middle', 'horizontalalignment', 'left','Parent', ax);
          t = text( xl(2), xl(2)^2, ['x^2'], ...
                  'verticalalignment', 'middle', 'horizontalalignment', 'left','Parent', ax);


   if logit == 3
      set(ax, 'Yscale', 'log', 'Xscale', 'log');
   elseif logit == 2
      set(ax, 'Yscale', 'log');
   elseif logit == 1
      set(ax, 'Xscale', 'log');
   end


