function [ms, mshalf, col] = getFullMList(M)

   coln = [.2 .2 .8]; % north
   colc = [.8 0 .8]; % center
   cols = [.6 .8 0]*.7; % south

   cnt = 1;
   ms = [];
   ms =  [ms find_names_inM( M, 'MS100')];
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'oc50')];
      col(cnt,:) = coln; cnt = cnt+1;
   ms =  [ms find_names_inM( M, 'NRLN50')];
      col(cnt,:) = coln; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'ps50')];
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'vb50n')];
      col(cnt,:) = cols; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'NRLS50')];
      col(cnt,:) = cols; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'oc40s')];
      col(cnt,:) = coln; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'ps40n')];
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'ps40m')];
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'oc25na')];
      col(cnt,:) = coln; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'oc25sa')];
      col(cnt,:) = coln; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'ps30m')];
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'ps30s')]; 
      col(cnt,:) = colc; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'vb30n')];
      col(cnt,:) = cols; cnt = cnt+1;
   ms = [ms find_names_inM( M, 'vb30s')];
      col(cnt,:) = cols; cnt = cnt+1;

   mshalf = [];
	mshalf = [mshalf find_names_inM( M, 'VB50S')];  
	mshalf = [mshalf find_names_inM( M, 'MS150')];
	mshalf = [mshalf find_names_inM( M, 'oc25nb')];
	mshalf = [mshalf find_names_inM( M, 'oc40n')];

   msextra = [ms find_names_inM( M, 'vb25n')];
   

end
