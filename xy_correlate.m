function [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt) 
%% [Xg, Yg, time, r, p, rb] = xy_correlate( xt, x, yt, y, dt)
%
%  This function correlates to time series 
%     
%     INPUT:
%        xt    : time vector of x
%        x     : value vector of x
%        yt    : time vector of y
%        y     : value vector of y
%        dt    : dt of comon time vector (days)
%
%     OUTPUT
%        Xg    : gridded x-values (normally averages)
%        Yg    : gridded y-values
%        time  : common time vector
%        r     : correlation
%        p     : anti-correlation (see doc corr)
%        rb    : vector with lower and upper bound for r
%
%   created by: 
%        Johannes Becherer
%        Thu Jun 13 13:00:21 PDT 2019
   

   if nargin <5
      dt = min( nanmedian(diff(xt)), nanmedian(diff(yt)));
   end

   tl    = [max( xt(1), yt(1)) min( xt(end), yt(end))];
   time  = tl(1):dt:tl(2);


   Xg = clever_interp( xt, x, time);
   Yg = clever_interp( yt, y, time);

   % make sure the vector orientation is such that corr likes it
   if size(Xg,1) == 1
      Xg = Xg';
   end
   if size(Yg,1) == 1
      Yg = Yg';
   end

   ii_notnan = find( ~isnan(Xg)& ~isnan(Yg) );

   [R, P, rl, ru] = corrcoef( Xg(ii_notnan), Yg(ii_notnan));

   r  = R(1,2);
   p  = P(1,2);
   rb = [rl(1,2) ru(1,2)];

end
