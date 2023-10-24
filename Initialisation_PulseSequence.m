function echo_stack=OneScript_PulseSequence(spin_system, parameters, H, R, K)
    L=H+1i*R+1i*K;
    % Electron Pulse Operators
    Ep=operator(spin_system,'L+', parameters.spins{1});
    Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

    % Freq Offsets. 
    parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)-...
                            parameters.pulse_frq-parameters.offset;

    % First Pulse on initialise state (parameters.rho0)
    rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0, parameters.pulse_frq(1),parameters.pulse_pwr(1),...
                                                            parameters.pulse_dur(1),parameters.pulse_phi(1), ...
                                                             parameters.pulse_rnk(1),parameters.method);
    % Evolutino after first pulse. 
    stepsize=parameters.p1_p3_gap/parameters.p2_nsteps;
    %evolution of states after first pulse. - trajectory allows for toatiaons in xy plane.  
    rho_stack=evolution(spin_system, L, [], rho, stepsize, parameters.p2_nsteps, 'trajectory');

    % Second Pulse on the evolved states.
    rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(2),parameters.pulse_pwr(2),...
                                                            parameters.pulse_dur(2),parameters.pulse_phi(2), ...
                                                            parameters.pulse_rnk(2),parameters.method);
    % Second Evolution.
    rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),stepsize,parameters.p2_nsteps,'refocus');
    
    %Third Pulse on the evolved states.
    rho_stack=shaped_pulse_af(spin_system,L,Ex,Ex,rho_stack,parameters.pulse_frq(3),parameters.pulse_pwr(3),...
                                                            parameters.pulse_dur(3),parameters.pulse_phi(3), ...
                                                            parameters.pulse_rnk(3),parameters.method);
    
    % Refocus after pulse 3.
    echo_location=parameters.p1_p3_gap+parameters.pulse_dur(1)/2+...
                                       parameters.pulse_dur(2)-parameters.echo_time/2;
    rho_stack=evolution(spin_system,L,[],rho_stack,echo_location,1,'final');

    stepsize = parameters.echo_time/parameters.echo_nsteps;
    % Final Evolution
    echo_stack=evolution(spin_system,L,parameters.coil,rho_stack,stepsize,parameters.echo_nsteps,'observable');
    end