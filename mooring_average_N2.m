function [N2_time, N2_maxavg, N2_std] = mooring_average_N2(M, i_ms, tl, dt)
%  This function averages all moorings max N2 together
tl = datenum( 2017, 9, 12) + [-10 60];
dt = 1/2;
            

%generate time vector
N2_time = tl(1):dt:tl(2);

%initialize temporay fields
tmp2d_N2 = nan(length(i_ms), length(N2_time));

for i = 1:length(i_ms)
   m = i_ms(i);
   [~, tmp1d_N2max, ~,~,~,~] = N2_avgmaxmin( M(m).T, M(m).waterdepth, 5 );
   tmp2d_N2(i,:) = clever_interp( M(m).T.time, tmp1d_N2max, N2_time);
end

N2_maxavg = nanmean(tmp2d_N2);
N2_std    = nanstd(tmp2d_N2);

