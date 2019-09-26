function [ii]  = find_names_inM( M, namePattern )
%% [ii]  = find_names_inM( M, namePattern )
%

   ii = [];

   for m = 1:length(M)
      if contains( M(m).mooringName, namePattern, 'IgnoreCase', true )
         ii = [ii m];
      end
   end


end


