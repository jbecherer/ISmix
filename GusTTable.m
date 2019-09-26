close all; 
warning off;
%_____________________load mooring data______________________
if 0; %{{{
clear all;

   % path to m_map here
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/m_map/'));
   addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));
   addpath(genpath('~/arbeit/matlab_tbx/My_tbx/seawater/'));
   addpath(genpath('./latexTable/'));
   
   load cmap;

%_____________________all epsilon data______________________
   load('/home/johannes/gdrive/IS17/collection_for_Jim/data/all_IS_data_hwg.mat');

% load bg state
   load('../data/Tsort.mat');

% load bottom GusTs
  load('../data/G_bbl.mat')
end %}}}

moorName   = {};
GustHeight = {};
Lat        = {};
Lon        = {};
water_depth= {};

cnt = 1;
for m = 1:length(M)
      moorName{cnt}   = M(m).mooringName;
      Lat{cnt}        = num2str( M(m).lat, '%3.2f');
      Lon{cnt}        = num2str( M(m).lon, '%3.2f');
      water_depth{cnt}= num2str( round(M(m).waterdepth));
      if sum((M(m).G.mab>1))>0
         GustHeight{cnt} = num2str( round(M(m).G.mab(M(m).G.mab>1)*2)/2, '%2.1f ' );
      else
         GustHeight{cnt} = ' x  '
      end
      cnt = cnt+1;
end


T = table(   moorName', water_depth', Lon', Lat',  GustHeight');

input.data = T;
%input.tableColLabels = {'mooring','water depth','GusT/\chi-pod depth'};
latex = latexTable(input);

fid=fopen('../data/GusTTable.tex','w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
