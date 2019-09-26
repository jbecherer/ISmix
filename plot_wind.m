function [] =  plot_wind( time,  U, ax, fig, tl, yl, col, LW )
   f_pos = get(fig, 'PaperPosition');
   ax_pos = get(ax, 'Position');

   f_aspect = f_pos(3)/f_pos(4);
   ax_aspect = ax_pos(3)/ax_pos(4);

   t_factor = diff(tl)/diff(yl)/f_aspect/ax_aspect;


   for t = 1:length(time)
      plot(ax, [0 real(U(t))*t_factor]+time(t) , [0 imag(U(t))], 'color', col ,  'Linewidth', LW);
   end
   
  %  plot(ax, time, imag(U), 'Linewidth', 1);
   ylim(ax, yl);
   


end
