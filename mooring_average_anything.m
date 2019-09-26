function [time, any_avg, any_std] = mooring_average_anything(M, i_ms, anything, prefix, suffix,  tl, dt)
% [time, any_avg, any_std] = mooring_average_anything(M, i_ms, anything, prefix, suffix,  tl, dt)
%  This function averages anything over all moorings
%   anything = 'Ef.uppp';
%

if nargin<5
   prefix   = '';  % could be 'log10('
   suffix   = '';  % ')'
end
if nargin<7
   tl = datenum( 2017, 9, 12) + [-10 60];
   dt = 1/2;
end
            
istr_dot = strfind(anything, '.');
eval_statement = ['clever_interp( M(m).' anything(1:istr_dot) 'time , ' ...
                  prefix 'M(m).' anything suffix ...
                  ', time)'];

eval_statement_nan = ['interp1( M(m).' anything(1:istr_dot) 'time , double(~isnan(' ...
                  prefix 'M(m).' anything suffix ...
                  ')), time, ''nearest'' )'];

%generate time vector
time = tl(1):dt:tl(2);

%initialize temporay fields
tmp2d = nan(length(i_ms), length(time));

for i = 1:length(i_ms)
   m = i_ms(i);
   tmp2d(i,:) = eval(eval_statement);
   ii_nan = find(eval(eval_statement_nan)<1);
   tmp2d(i, ii_nan) = nan;

end

if numel(tmp2d)>numel(time)
   any_avg    = nanmean(tmp2d);
   any_std    = nanstd(tmp2d);
else  % if only one mooring is selected
   any_avg    = tmp2d;
   any_std    = zeros(size(tmp2d));
end

