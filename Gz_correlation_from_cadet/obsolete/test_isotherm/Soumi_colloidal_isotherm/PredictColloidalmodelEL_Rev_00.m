
function [solution] = PredictColloidalmodelEL_Rev_00(protein_load, TimeSections, FlowSections, pH_time, pH, par)


%IS= ionic strength , first component
%ke=ka*IS^-kb+kc.*exp(kd*IS);
%bpp=bppa*IS^bppb+bppd.*exp(IS*bppc)
% kappa=1e9/(kappac*IS^kappacp+.kappacd)
% kappacd = 2.5
% ParamType ke = b0l * pow(tempC, b1l) + b2l;
% calculation of Surface diffusion in code
%ke = b0l * IS^ b1l + b2l;
%1 1 0
	%mult = mult + b6l*(ke^b3l)  + b4l* exp((b5l*ke));
    %Ds=Ds * mult;
 % mult = 1

n=1; % one component protein

 %% Step A: Construct pre-column (to model holdup)
	% ===================================================

	% General rate model unit operation
	mGrm1 = GeneralRateModel();

	% Discretization
	mGrm1.nComponents = 2; % each component data +salt
	mGrm1.nCellsColumn = 60; 
	mGrm1.nCellsParticle =30;  
    mGrm1.nBoundStates = [0 1];
	

	% Initial conditions
	mGrm1.initialBulk = [10^(-pH(1)) 0]; % [mol / m^3], also used for the particle mobile phase for component 1 and 2
	mGrm1.initialSolid = 0; % [mol / m^3] for component 2
		
	% Transport
	mGrm1.dispersionColumn          = 3.3e-6; % [m^2 / s]
	mGrm1.filmDiffusion             = [1e-5 1e-5]; % [m/s]
	

    mGrm1.diffusionParticle         = [0 0 ]; % [m^2 / s]
    mGrm1.diffusionParticleSurface  = [0]; % [m^2 / s]
   
    
    mGrm1.diffusionParticleSurfaceB0  = [0.0 1 ]; 
    mGrm1.diffusionParticleSurfaceB1  = [0.0 1 ]; 
    mGrm1.diffusionParticleSurfaceB2  = [0.0 0 ]; 
    mGrm1.diffusionParticleSurfaceB3  = [0.0 0 ]; 
    mGrm1.diffusionParticleSurfaceB4  = [0.0 0 ]; % EDIT Changed from [0, 1]
    mGrm1.diffusionParticleSurfaceB5  = [0.0 0 ]; 
    mGrm1.diffusionParticleSurfaceB6  = [0.0 1 ]; 
 
 
    % Geometry 
    columnRadius                      = 7.5E-4;% [m]
    columnCSA                         = pi*(columnRadius^2); %[m^2]
	mGrm1.columnLength                = 0.08; % [m] % 0.12 default
    mGrm1.interstitialVelocity        = [((0.34E-6)/(60*columnCSA))/1];%, ...
                                         %((0.085E-6)/(60*columnCSA))/1, ...
                                         %((0.340E-6)/(60*columnCSA))/1, ...
                                         %((0.340E-6)/(60*columnCSA))/1]; % [m/s] % EDIT - from 0.4669
    mGrm1.particleRadius              = 2.5e-5; % [m] % EDIT - from 0.4669
	mGrm1.porosityColumn              = 1; % [-]
	mGrm1.porosityParticle            = 0.300; % [-]    % EDIT - from 0.4669
	
    mGrm1.useAnalyticJacobian = false; % leave at false always, not defined in code
    
	% Adsorption
	isotherm1 = LinearBinding();
    isotherm1.kineticBinding = true; 
    isotherm1.kA = [0  0.0001];
    isotherm1.kD = [0  1];

	mGrm1.bindingModel = isotherm1;


    %% Step 1: Construct general rate model unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = GeneralRateModel();

	% Discretization
    
    
	mGrm.nComponents = 2; % each component data +salt
	mGrm.nCellsColumn = 60; 
	mGrm.nCellsParticle = 30;  
    mGrm.nBoundStates = [0 1];
	

	% Initial conditions
	mGrm.initialBulk = [10^(-pH(1)) 0]; % [mol / m^3], also used for the particle mobile phase for component 1 and 2
	mGrm.initialSolid = 0; % [mol / m^3] for component 2
		
	% Transport
	mGrm.dispersionColumn          = 3.3e-6; % [m^2 / s] - 26APR21-Changed from 3.3e-6
	mGrm.filmDiffusion             = [1e-5 1e-5]; % [m/s]
	
%IS= ionic strength , first component

% kappa=1e9/(kappac*IS^kappacp+.kappacd)
% kappacd = 2.5
% ParamType ke = b0l * pow(tempC, b1l) + b2l;
% calculation of Surface diffusion in code
%ke = b0l * IS^ b1l + b2l;
%1 1 0
	%mult = mult + b6l*(ke^b3l)  + b4l* exp((b5l*ke));
    %Ds=Ds * mult;
 % mult = 1
 
    mGrm.diffusionParticle         = [6e-11 2.9e-12 ]; % [m^2 / s]
    mGrm.diffusionParticleSurface  = 1E-13; % [m^2 / s]
   
    
    mGrm.diffusionParticleSurfaceB0  = [0.0 1 ]; % Keep 1
    mGrm.diffusionParticleSurfaceB1  = [0.0 1 ]; % Keep 1
    mGrm.diffusionParticleSurfaceB2  = [0.0 0 ]; % Keep 0
    mGrm.diffusionParticleSurfaceB3  = [0.0 0 ]; % Keep 0
    mGrm.diffusionParticleSurfaceB4  = [0.0 0.274 ]; % Change to pre-exponential term
    mGrm.diffusionParticleSurfaceB5  = [0.0 14756]; % Change to exponential coefficient
    mGrm.diffusionParticleSurfaceB6  = [0.0 0 ]; % Keep 0
 
 
    % Geometry
    % Particle Geometry
    mGrm.particleRadius              = 2.5e-5;  % [m]  
    mGrm.porosityParticle            = 0.30;    % [-] 
	% Column Geometry
    mGrm.porosityColumn              = 0.42;    % [-] (Voidage)
    mGrm.columnLength                = 0.025;   % [m] (2.5 cm column, 0.170 mL)
    mGrm.interstitialVelocity        = [(0.34/0.170)*(4.166667e-4)/mGrm.porosityColumn]; %...
                                        %(0.085/0.170)*(4.166667e-4)/mGrm.porosityColumn, ...
                                        %(0.34/0.170)*(4.166667e-4)/mGrm.porosityColumn,  ...
                                        %(0.34/0.170)*(4.166667e-4)/mGrm.porosityColumn] ; % [m/s] % EDIT - from 0.4669
   
	
    mGrm.useAnalyticJacobian =false; % leave at false always, not defined in code
    
    %ke=ka*IS^-kb+kc.*exp(kd*IS);
%bpp=bppa*IS^bppb+bppd.*exp(IS*bppc)
	% Adsorption
	isotherm = ColloidalMultiBinding();
    isotherm.kineticBinding = true; % Quasi-stationary binding
	isotherm.kkin       = [1e10  1e10];
    isotherm.kkd       = [1e10  1e10] ;
    
    isotherm.ka         =  [ 0 0];  
	isotherm.kb         =  [ 0 0];%1.72 ] ; 
    isotherm.kc         =  [ 0 24.658]; %-3.93] ; 
    isotherm.kd         =  [ 0 -16202] ; % for exponent function otherwise 0
    
    isotherm.bppe       = [ 0 1]; %not defined leave it at 1
    isotherm.bppa       = [ 0 0];%37220] ; 
    isotherm.bppb       = [ 0 0];%-2.25]; 
    isotherm.bppc       = [ 0 0]; 
    isotherm.bppd       = [ 0 3.25];%0.26] ;  
    
    isotherm.radius     = [1.59e-9 4.5e-9 ];
      
   isotherm.phaseratio = (4.68e7); % 4.68e7 for mAb
    
    
  
    isotherm.cordNum = 6;
    isotherm.linearApprox = 1e-7;
    
    %   kappa=1e9/(kappac*IS^kappacp+.kappacd)
    %   kappacd = 2.5
    isotherm.kappac =  0;
    isotherm.kappacp = 0;
    isotherm.kappacd = 2.5;
	mGrm.bindingModel = isotherm;


	%% Step 2: Construct inlet unit operation
	% ===================================================

	% Inlet unit operation
	mIn = PiecewiseCubicPolyInlet();
	mIn.nComponents = 2;
    
    % Get index for protein injection time
    time_length = length(pH_time);
    [protein_injection_delta, protein_injection_index_start]= min(abs(pH_time-TimeSections(1))); 
    [protein_injection_delta, protein_injection_index]= min(abs(pH_time-TimeSections(2)));
    [protein_injection_delta, protein_wash1_index]= min(abs(pH_time-TimeSections(3)));
    [protein_injection_delta, protein_wash_index]= min(abs(pH_time-(TimeSections(4))));
    
    
    % Get injection profiles
    protein_injection_data=zeros(time_length,1); % protein data, zeros 
    protein_injection_data(protein_injection_index_start:protein_injection_index)=protein_load/146; % protein data, injection
    
    % Adjust pH profile for time delay
    pH_delay = 100;
    pH_end_array = pH(end)*ones(pH_delay,1);
    pH = [pH;pH_end_array];
    pH(1:pH_delay) = [];
    
    % Make inlet profiles - protein and H+
    inletProfileProtein = PiecewiseCubicPolyProfile.fromUniformData(pH_time, protein_injection_data); % Get protein profile
    inletProfile = PiecewiseCubicPolyProfile.fromUniformData(pH_time,10.^(-pH)); % Get H+ ion profile
    
    % Consolidate - column 1 is H+, column 2 is protein
    % NOTE - Ensure only constant term remains! The fitting produces some
    % artifacts which seem to significantly impact outlet profile
    inletProfile.constant(:,2)  =inletProfileProtein.constant(:,1);
    inletProfile.linear(:,2)    =0;%inletProfileProtein.linear(:,1);
    inletProfile.quadratic(:,2) =0;%inletProfileProtein.quadratic(:,1);
    inletProfile.cubic(:,2)     =0;%inletProfileProtein.cubic(:,1);

    % Put it all down in the inlet unit operation "mIn" defined at the top
    mIn.constant       = inletProfile.constant;
	mIn.linear         = inletProfile.linear;
	mIn.quadratic      = inletProfile.quadratic;
	mIn.cubic          = inletProfile.cubic;
     
    % Adjust interstitial velocity vector
    flow_array = ones(time_length-1,1); % In ml/min
    flow_array(protein_injection_index_start:protein_injection_index) = FlowSections(1);
    flow_array(protein_injection_index:protein_wash1_index) = FlowSections(2);
    flow_array(protein_wash1_index:protein_wash_index) = FlowSections(3);
    flow_array(protein_wash_index:end) = FlowSections(4);
    mGrm1.interstitialVelocity = [(flow_array*1E-6)/(60*columnCSA)];
    mGrm.interstitialVelocity  = [(flow_array/0.170)*(4.166667e-4)/mGrm.porosityColumn];
    

    %% Step A: Construct pre-column (to model holdup)
	% ===================================================

	% General rate model unit operation
	mGrm2 = GeneralRateModel();

	% Discretization
	mGrm2.nComponents = 2; % each component data +salt
	mGrm2.nCellsColumn = 10; 
	mGrm2.nCellsParticle =10;  
    mGrm2.nBoundStates = [0 1];
	

	% Initial conditions
	mGrm2.initialBulk = [10^(-pH(1)) 0]; % [mol / m^3], also used for the particle mobile phase for component 1 and 2
	mGrm2.initialSolid = 0; % [mol / m^3] for component 2
		
	% Transport
	mGrm2.dispersionColumn          = 1e-10; % [m^2 / s]
	mGrm2.filmDiffusion             = [1e-10 1e-10]; % [m/s]
	

    mGrm2.diffusionParticle         = [0 0 ]; % [m^2 / s]
    mGrm2.diffusionParticleSurface  = [0]; % [m^2 / s]
   
    
    mGrm2.diffusionParticleSurfaceB0  = [0.0 1 ]; 
    mGrm2.diffusionParticleSurfaceB1  = [0.0 1 ]; 
    mGrm2.diffusionParticleSurfaceB2  = [0.0 0 ]; 
    mGrm2.diffusionParticleSurfaceB3  = [0.0 0 ]; 
    mGrm2.diffusionParticleSurfaceB4  = [0.0 0 ]; % EDIT Changed from [0, 1]
    mGrm2.diffusionParticleSurfaceB5  = [0.0 0 ]; 
    mGrm2.diffusionParticleSurfaceB6  = [0.0 1 ]; 
 
 
    % Geometry 
    columnRadius2                      = 5E-6;% [m]
    columnCSA2                         = pi*(columnRadius2^2); %[m^2]
	mGrm2.columnLength                = 0.001; % [m] % 0.12 default
    mGrm2.interstitialVelocity        = [(flow_array*1E-6)/(60*columnCSA2)];%, ...
                                         %((0.085E-6)/(60*columnCSA))/1, ...
                                         %((0.340E-6)/(60*columnCSA))/1, ...
                                         %((0.340E-6)/(60*columnCSA))/1]; % [m/s] % EDIT - from 0.4669
    mGrm2.particleRadius              = 2.5e-5; % [m] % EDIT - from 0.4669
	mGrm2.porosityColumn              = 1; % [-]
	mGrm2.porosityParticle            = 0.300; % [-]    % EDIT - from 0.4669
	
    mGrm2.useAnalyticJacobian = false; % leave at false always, not defined in code
    
	% Adsorption
	isotherm2 = LinearBinding();
    isotherm2.kineticBinding = true; 
    isotherm2.kA = [0  0.0001];
    isotherm2.kD = [0  1];

	mGrm2.bindingModel = isotherm2;
    
    %%END OF BACKMIXING
    
    
	% Step 3: Assemble system of unit operations
	% ===================================================
    mOut = OutletModel();
	mOut.nComponents = 2;
	% Construct ModelSystem and assign unit operations (order determines IDs)
	mSys = ModelSystem();
	mSys.models = [mIn, mGrm1, mGrm, mGrm2, mOut];%, mOut];% [0 1 2]

	% Define valve configurations / unit operation connections
	% Valve configuration active on entering section 0
	mSys.connectionStartSection = [0 ];% no meaning yet
	% Connect unit 0 with unit 1 (-1 = all components)
	%mSys.connections = {[0, 1, -1, -1.0, 1.0 ]}; %% first 0 first unit inlet.. 1 .. 2 inlet GRM .. then all component {[0, 1, -1, -1], [1, 2, -1, -1] }
    mSys.connections = {[0, 1, -1, -1,  1.0;    ...
 	                     1, 2, -1, -1,  1.0;    ...
                         2, 3, -1, -1,  1.0;
                         3, 4, -1, -1,  1.0]};

                     

	% Step 4: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();
    sim.nThreads = 12;

    sim.solutionTimes = pH_time; % [s], time points at which solution is computed
    sim.sectionTimes = inletProfile.breaks;
	sim.sectionContinuity = inletProfile.continuity;
            
    sim.model = mSys;

    
	% Step 5: Run the model and plot the results
	% ===================================================

    sim.absTol=1.0000e-08; % 1.0000e-08;
    sim.relTol=1.0000e-06; % 1.0000e-06;
    sim.algTol= 1.0000e-08; % 1.0000e-08;
    sim.relTolSens=1.000e-6; % 1.0000e-06;
    sim.maxSteps=1000; % 1000
    
    result = sim.run();

	% Extract solution into a matrix with time being the first column
	% Note that we need to extract the outlet of the second unit operation, which
	% is the general rate model (order in the ModelSystem matters)
    solution = [result.solution.time, result.solution.outlet{5}, [flow_array;flow_array(end)]];

        
    %solution = solution(:,3); % This line is added to return only protein values from the model
    %plotyy(solution(:,1),solution(:,3),solution(:,1),-log10(solution(:,2)));
    

    
     

end


