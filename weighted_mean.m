function  [y_mean] =  weighted_mean( N, Y)
   Yw = N.*Y;
   y_mean = nansum(Yw)./nansum(N);
end
