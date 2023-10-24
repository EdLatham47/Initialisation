%all of the three pulse deer written by me into one script. 

function OneScript()
    % System Specifications. 
    sys.magnet = 0.3333333; %T
    % SPins to be viewed.
    sys.isotopes = {'E','E'}; %Two seperated electrons, as oppose to two in one shell (E2)
    sys.disable={'trajlevel'}; % (For less than five spins, speeds up computation.)
    %           SPIN INTERACTIONS
    % defines structure of spins. 
    inter.zeeman.eigs = cell(1,2);
    inter.zeeman.euler = cell(1,2);
    inter.coordinates = cell(2,1);
    % First Spin g-values then euler angle relative to B0.
    inter.zeeman.eigs{1}=[2.0034 2.0031 2.0027];
    inter.zeeman.euler{1}=[0 0 0];
    inter.coordinates{1}=[0.00 0.00 0.00];
    % Second Spin g-values then euler angle relative to B0.
    inter.zeeman.eigs{2}=[2.005 2.007 2.009];
    inter.zeeman.euler{2}=[0 0 0];
    inter.coordinates{2}=[30.00 0.00 0.00]; % Angstrungs
    %inter.coupling.scalar={0 30*10^6; 0 0}; % Coupling
    %   -    Providinng a specific Coupling such as this, makes the
    %   Frequency Sweep Spectrum very 'clean' whereas two distances
    %   provided provides the heavy 'interference'.

    % - Relaxation
    inter.relaxation={'lindblad'}; %lindblad
    inter.lind_r1_rates=[1.0 1.0]*10^2; % Hz
    inter.lind_r2_rates=[1.0 1.0]*10^5; % Hz
    inter.temperature=10; %Kelvin. Default: Spinach uses the high-temperature approximation. 
                        % Specifying zero makes the system start at the lowest energy eigenstate of the Hamiltonian.
    inter.rlx_keep='diagonal';%Not sure what these do..... will require further investigation.
    inter.equilibrium='zero';
    %Create Spin System using the specifications above. (Before experimental parameters.)
    spin_system=create(sys,inter);
    % Basis Set
    bas.formalism = 'sphten-liouv';
    bas.approximation = 'none';
    %rotate spin sys into basis
    spin_system=basis(spin_system,bas);


    %           COMPUTATIONAL PARAMETERS
    % spherical averaging grid for solid state experiments
    parameters.grid='rep_2ang_6400pts_sph'; % ~./spinnach_2_7_6049/kernel/grids
    % MEthod. soft puse propagation method,'expv' for Krylov propagation,'expm' for exponential propagation, 'evolution' for Spinach evolution function
    parameters.method='expm'; % was evolution
    % Fidelity of the pulse. (Keep in 2^n for ease of computation.)
    parameters.npoints = 256;
    % zerofill - Fidelity of Free Induction Decay Plots. (First Figure) 
    parameters.zerofill = 2048;
    parameters.axis_units = 'GHz-labframe';
    % if set to 1, the spectrum is differentiated before plotting
    parameters.derivative=1;
    % if set to 1, the spectrum is multiplied by -1 before plotting (Left to right right to left for the Hz)
    parameters.invert_axis=0; 
    % As oppose to nmr
    parameters.assumptions='deer-zz';
    %STFU
    parameters.verbose=0;
    
%gaussleg.m for magic stuff to do with linewidth calculating.
    % linewidth 0.1 20 gauss 0.002Tesla +/-1gauss 
    %           EXPERIMENTAL PARAMETERS
    % offset - the offset of the pulse from the center of the pulse.
    parameters.offset = 0;   % a cell array giving transmitter off-sets in Hz on each of the spins listed in parameters.spins array, (Lamor Frequency from the Transmitter Frequency. 
    % http://nmrwiki.org/wiki/index.php?title=Offset#:~:text=From%20NMR%20Wiki&text=Offset%20(%CE%A9%2C%20omega%3B%20majuscule,for%20offset%20are%20in%20hertz.
    %parameters.Bw=10^9; %Hz https://spindynamics.org/wiki/index.php?search=parameters.Bw&title=Special%3ASearch&go=Go
    % Doesn't change even at 10^9

    % Sweep width in each dimension of the spectrum (Hz). - RESOLVED
    parameters.sweep = 1e8;     % Width of Sweep for First Figure.
    % Linewidth?
    parameters.fwhm=0.01;    %line FWHM, Tesla

    parameters.spins = {'E'}; % spins of the molecule that will be activlley operated on.
    % a combined initial state of both spins, alligned in the Bz direction. 
    parameters.rho0=state(spin_system,'Lz','E');
    % detection state of the output files. (pi/2 tf L+), can add more than one detection state.
    parameters.coil = state(spin_system,'L+','E');
    %       THE INDIVIDUAL PULSES
    % DEER pulse parameters - all params can be lists for multi-pulse sequences. 
    %pulse_rnk 1 would be a perfect pulse. 2 does narruto change. 3 will aid in more complex-multi rotations. 
    parameters.pulse_rnk=[2 2 2];
    % Pulse duration
    parameters.pulse_dur=[19e-9 40e-9 38e-9];       % altered 30ns suggestion - second spin takes slightly longer.
    % Phase on initial wave
    parameters.pulse_phi=[pi pi pi];
    % Power of the pulse
    parameters.pulse_pwr=2*pi*[12e6 12e6 12e6];
    % Frequency of the pulse - independantly specified is a bit odd.
    parameters.pulse_frq=[9.34e9 9.363e9 9.345e9];
    %        GAPS BETWEEN PULSES
    parameters.p1_p3_gap=1e-6;
    %fidelity of p2 
    parameters.p2_nsteps=100;
    parameters.DeerCut=.2; %ammount of the Deer data to cut the echo stack, as a decimal from 0-1. this is dedeucted from both ends. 
    % readout gap for echo
    parameters.echo_time=200e-9;
    % pfidelity of echo
    parameters.echo_nsteps=100;
    %           SIMULATION
    % End of inputs, begin simulation.
    
    %P = sphten2zeeman(spin_system);
    %Happiness=P*parameters.rho0
    %stateinfo(spin_system, Happiness, 4)
    %works but they arent set up right.
    relaxan(spin_system)
    OneScript_Plots(spin_system, parameters)
    end