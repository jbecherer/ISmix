function [ time, z,  uT, T_a, U_a] = cal_Tflux( u_time, u_z, U, t_time, t_z, T, lp_cutoff, hp_cutoff, H)
%%
%
%     H        :  water depth +
%  
%  OUtput:
%     time     :  common time stemp
%     z        :  common vertical axis
%     uppp     :  u'p' depth and time depedent (complex = up+ 1i*vp)
%     pp       :  p'
%     up       :  u' (complex)
%     rho_a    :  low-pass desity depedence
%     U_a      :  low-pass time depdends
%
%   created by: 
%        Johannes Becherer
%        Mon Jun 17 13:56:22 PDT 2019


%  test routine
   if 0
      [ms, mshalf] = getFullMList(M);
      m=26%ms(11);
      H        = M(m).waterdepth;
      u_time   = M(m).U.time;
      u_z      = M(m).U.mab;
      U        = M(m).U.U;
      t_time   = M(m).T.time;
      t_z      = M(m).T.mab;
      T        = M(m).T.T;
      lp_cutoff= 3/24;
      hp_cutoff= 72/24; 
   else
      if nargin < 9
         H = round(max( max(u_z), max(t_z) ));
      else
         H = abs(H);
      end
   end


% remove tmeperature out lineers
 if size(T,1)>6
   T = T - repmat( nanmean(movmean(T,1),2)-qbutter(nanmean(movmean(T,1),2),.3), 1, length(t_time));
 end

% construct time vector
   dt    = 10*max( nanmedian(diff(u_time)),nanmedian(diff(t_time)) );
   tl(1) = max( u_time(1), t_time(1) );
   tl(2) = min( u_time(end), t_time(end) );
   time  =  tl(1):dt:tl(2);

% construct common z vector;
   zl(1) = 0;
   zl(2) = H;
   dz    = 1;
   z     = zl(1):dz:zl(2);

% interpolate u and T on common grid

   U_g   =  nan(length(z), length(time));
   T_g   =  nan(length(z), length(time));
   U_a   =  nan(length(z), length(time));
   uppp  =  nan(length(z), length(time));
   up    =  nan(length(z), length(time));
   pp    =  nan(length(z), length(time));
   rho_a =  nan(length(z), length(time));
   N2    =  nan(length(z), length(time));
   eta   =  nan(length(z), length(time));


   for t = 1:length(time)
      itu   =  find(u_time>=time(t),1); 
      itt   =  find(t_time>=time(t),1); 

      if ~isempty(itu)
         U_g(:,t) =  clever_interp( u_z,  U(:, itu), z );
        % top bottom extra polation
        neinnan  = find( ~isnan(U_g(:,t)) );
        janan    = find( isnan(U_g(:,t)) );
        if length(neinnan)>3 & ~isempty(janan)
           % top
           iit = find( janan > neinnan(end) );
           if ~isempty(iit)
              U_g(janan(iit),t) =  median(U_g(neinnan(end-2),t));
           end
           % bottom % log-layer
           iib = find( janan < neinnan(1) );
           if length(iib)>0
              us = sqrt(2e-3*abs(U_g(neinnan(1),t)).^2)*exp(1i*angle(U_g(neinnan(1),t)));
              U_g(janan(iib(1:end)),t) = us/.4*log((z(janan(iib(1,:)))+.5*dz)/1e-3);
           end
        end
      end
      if ~isempty(itt)
         T_g(:,t) =  clever_interp( t_z,  T(:, itt), z );
        % top bottom extra polation
        neinnan  = find( ~isnan(T_g(:,t)) );
        janan    = find( isnan(T_g(:,t)) );
        if ~isempty(neinnan) & ~isempty(janan)
           % bottom
           iib = find( janan < neinnan(1) );
           if ~isempty(iib)
              T_g(janan(iib),t) =  T_g(neinnan(1),t);
           end
           % top
           iit = find( janan > neinnan(end) );
           if ~isempty(iit)
              T_g(janan(iit),t) =  T_g(neinnan(end),t);
           end
        end
      end

   end
   U_g = movmean(U_g, round(lp_cutoff/dt), 2);
      U_g(isnan(real(U_g))) = nan + 1i*nan;
   T_g = movmean(T_g, round(lp_cutoff/dt), 2);

   uT             = movmean( T_g.*U_g, round(hp_cutoff/dt), 2, 'omitnan');
      uT(isnan(real(U_g))) = nan + 1i*nan;
   T_a            =  movmean( T_g, round(hp_cutoff/dt), 2, 'omitnan');
   U_a            = movmean( U_g, round(hp_cutoff/dt), 2, 'omitnan');
      U_a(isnan(real(U_g))) = nan + 1i*nan;

