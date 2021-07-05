% ########### DCG_simulation_band #########################################
% version 1.01, June 2021
% P. Moore
% Matlab script (developed in Matlab 2018b)
% -----
% *numerical simulation of coupled ablation and supraglacial debris
% transport.
% *uses conservative FD scheme: D,q computed at midpoints, h,T at nodes.
% *uses either Anderson linear or Roering-type nonlinear diffusion rule
% *** This simulation has a single englacial debris source on the right
% side of the domain
%
% #########################################################################

% the next line closes plot windows and clears the workspace, so make sure
% everything you need is saved, or alternatively comment the line below
% (but beware of overwriting variable names).
clearvars; close all;

% >>>>>> uncomment below to initiate movie recording utility
%v = VideoWriter('movie_name.avi');   %>>>>>> change filename as needed
%open(v);

% -------------------------------------------------------------------------

% ########### Domain set-up #############
% >>>>>>> change L, zinit, dx as needed.
L = 50;         % total domain length
zinit = 50;     % initial ice thickness, same units as L
dx = 1;         % horizontal spatial discretization
nx = L/dx;      % number of 1D nodes
X = 0:dx:L-dx;  % height & thickness computed on nodes

% define debris source patterns from ice
deb_conc = 0.1; 
deb_w = 10*dx;
deb = zeros(1,nx);
deb(1,nx-deb_w:end) = deb_conc;

% --------------- variables and parameters -----------------------------

% ########### Model parameters and initial conditions ##########
% >>>>>>>>> debris transport constants <<<<<<<<
D0 = 1.0;  % m^2/yr [per 1/1 slope]; a maximum for large thickness
Sc = 0.9;       % critical slope for landsliding, expressed as tan(alpha)

% >>>>>>>>> ablation model <<<<<<<<<<<<<<<<<<<<
% Hyperbolic Ostrem-curve approximation following Anderson and Anderson
% 2016

Li = 334000;            % latent heat
rhoi = 900;             % kg/m^3
Kheat = 1;              % thermal conductivity of debris 
sec2yr = 60*60*24*365;  % time unit conversion factor
Ts = 2;                 % debris surface temperature <<<<<<<<<<<<<<<<<<<<<
thick0 = 0.01;          % m, reference debris thickness for hyperbolic 
                            % decay, equivalent to H*

% volume adjustments for debris melt-out
n = 0.35;               % final debris porosity after meltout
pf = 1-n;

% --------------- model construction -------------------------------------

% ######## time-stepping set-up ##########
% This model currently uses a fixed time-step 
% >>>>>>>> TO-DO: adaptive TS inside model loop.
dt = 0.1*(dx^2)/(2*D0); % years. this is currently a fixed time-step
time_tot = 500; % years. may need to lengthen for thicker deb or lower Ts
nt = ceil(time_tot/dt); % number of time-steps in complete simulation
t = 0:dt:ceil(time_tot-dt); % time vector


% ######### initial conditions ############
% This simulation begins with a uniform debris surface elevation, but
% nonuniform underlying ice surface elevation. The resulting differences in
% debris thickness give rise to differential melt and lead to relief
% developing on the ice surface
% zice = elevation of the ice surface
% zdeb = elevation of the debris surface
% maxdebthick = 2; % >>>>> change as needed
% mindebthick = 0.1; % >>>>> change as needed
% % width in nx units of ramp between thick and thin
% transzonewidth = 10; % >>>>> change as needed 

% lines below create a uniform ice surface without initial debris.
zdeb0 = zinit*ones(1,nx);
zice0 = zdeb0;
zi = zice0;

% ####### initialize variable arrays ####################
zdeb = zeros(nx,nt); 
zdeb(:,1) = zdeb0;
debthick = zeros(nx,nt); 
debthick(:,1) = zdeb(:,1)-zi(:,1);
H=zdeb0.*ones(1,nx);              % debris surface elevation <===========
thick=(zdeb0-zice0).*ones(1,nx);  % debris thickness <==========
zice=zice0.*ones(1,nx);           % ice surface elevation = H-thick <==========

% ######## pre-allocate iteration arrays ############
Htemp = zeros(1,nx);             % temporary surface elevation vector
melt=zeros(1,nx);                % melt per timestep
dthick=zeros(1,nx);              % new debris added to surface by ablation
thicktemp=zeros(1,nx);           % temp storage vector for debris thickness
thicknew=zeros(1,nx);            % temp storage vector for debris thickness
Hnew=zeros(1,nx);                % debris surface elevation
zicetemp=zeros(1,nx);
D=zeros(1,nx-1);
slope=zeros(1,nx-1);
flux=zeros(1,nx-1);

% -------------------------------------------------------------------------

% ################# BEGIN MODEL LOOPS #####################################
timei = 0;
% ------- outer loop: adjust debris thickness -------

for i = 2:nt
	for j = 1:nx
		melt(j) = dt*Kheat*Ts/(rhoi*Li*(thick(j)+thick0))*sec2yr;  	% melt!
		if zice(j)>melt(j) % this catches melt exceeding ice thickness
			melt(j)=melt(j);
        else
            melt(j) = zice(j);
		end
		dthick(j) = melt(j)*deb(j)/pf;         % debris thickness released
    	thicktemp(j) = thick(j) + dthick(j);   % new nodal debris thickness
    	Htemp(j) = H(j) - melt(j) + dthick(j); % new surface elevation
    	zicetemp(j) = Htemp(j) - thicktemp(j);
    end

% ------- inner loop: evolve surface ----------------
	
% --- step 1; get midpoint diffusivities and fluxes -------

	for jj = 1:(nx-1)    % there are nx-1 midpoints; index # is LEFT node
        % this diffusivity function is based on Anderson (2000)
        D(jj) = D0*...
            (1-(exp(-(0.5*(thicktemp(jj)+thicktemp(jj+1)))/(thick0))));
        slope(jj) = (Htemp(jj+1)-Htemp(jj))/dx;
        %---- choose nonlinear or linear flux; comment out the other
		%flux(jj) = D(jj)*slope(jj)/(1-(slope(jj)/Sc).^2);  % nonlinear
        flux(jj) = D(jj)*slope(jj);   % linear
	end

% --- step 2; compute flux divergence, height changes ------	

	% --- solve for new thickness between edges ----
	for jjj = 2:(nx-1)
		thicknew(jjj) = thicktemp(jjj) + dt/dx*(flux(jjj)-flux(jjj-1));
		Hnew(jjj) = zicetemp(jjj) + thicknew(jjj);
	end

% --- enforce BCs : no-flux boundaries ---

	% ____ left ______            ____ right ______
	Hnew(1)=Hnew(2);             Hnew(nx)=Hnew(nx-1);
    thicknew(1)=thicknew(2);     thicknew(nx)=thicknew(nx-1);

	% --- update dependent variables ---------
	H = Hnew;
	thick = thicknew;
	zice = zicetemp;
	timei = timei+dt;
	zdeb(:,i)=H;  % storage array
	debthick(:,i)=thicknew;
    
    if sum(zice)==0
        break
    end
    
	% Plot result
	figure(1), clf
    hold on;
    subplot(3,1,1:2);
    axis([0 L 0 1.1*zice0(1)]);
	plot(X,H,'k-'); %hold on
%	plot(X,(H-thick),'b--');
    ylim([0,1.1*zice0(1)]);
	ylabel('elevation, m')
	title(['elevation after ',num2str(timei),' years'])
    subplot(3,1,3)
    area(X,thicktemp);
    axis([0 L 0 1.5]);
    xlabel('x [m]');
    ylabel('thickness, m')
	drawnow
    
      % Capture the plot as an image for video; uncomment below
      %frame = getframe(gcf);   ####
      %writeVideo(v,frame);     #####
end
%close(v);                      % ##### uncomment if video used

% ############### END MODEL LOOP #########################################

% Compute results parameters
HH = mean(debthick(:,1));
XX = L/HH;
Z = 50/HH;
tt = XX^2/D0;
tm = Z*Li*rhoi/(Kheat*Ts*sec2yr);
M = tm/tt;          % if Rt >> 1, melt more efficient
ditime = i*dt;

debIQR = iqr(debthick(:,i));
IQRnorm = debIQR/HH;

data = [HH, XX, Z, D0, Kheat, Ts, ditime, IQRnorm];
format shortG; disp(data)

% #### code below produces various summary plots #####################
% figure;
% %subplot(2,1,1);
% %patch([nx/2-deb_w nx/2+deb_w nx/2+deb_w nx/2-deb_w],[0 0 zice0 zice0],[0.8 0.8 0.8]); hold on;
% plot(X,zdeb(:,1:100:end),'k-'); axis([0 L 0 50]);
% %subplot(2,1,2); %area(X,zdeb(:,end),0);  hold on;  
% %plot(X,debthick(:,1:500:end),'k-'); %axis([0 400 0 0.3]);
% %figure; surf(debthick,'LineStyle','none')
% xlabel('horizontal distance'); ylabel('vertical distance');
% format short g; str = sprintf('M = %f, relative relief = %f',M,IQRnorm);
% title(str);