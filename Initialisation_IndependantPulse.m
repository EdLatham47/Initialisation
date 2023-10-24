function fids = OneScript_IndependantPulse(spin_system, parameters, H, R, K)
    % This function shows how soft pulses affect the magnetisation. 
    % It is a hypothetical experiment where a soft pulse specified by the user
    % is performed, immediately followed by an ideal pi/2 pulse on all spins 
    % followed by infinite-bandwidth time-domain detection.

    % Compose Liouvillian
    L = H + 1i*R + 1i*K;

    %Make the Pulse Operators. (Electron Plus on spin 1) - ' is transpose.
    Ep=operator(spin_system,'L+',parameters.spins{1});
    Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

    % Adjust Frequency from input to after experimental alterations
    parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)-parameters.pulse_frq-parameters.offset;

    %                  Conduct Pulses. (Electron Plus on spin 1)
    %Pulse 1 -'amplitude frequency' coordinates rather than xy. 

    rho1=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1),...
        parameters.pulse_dur(1),parameters.pulse_phi(1),...
        parameters.pulse_rnk(1),parameters.method);
    rho2=shaped_pulse_af(spin_system,L,Ex,Ey,rho1,parameters.pulse_frq(2),parameters.pulse_pwr(2),...
        parameters.pulse_dur(2),parameters.pulse_phi(2),...
        parameters.pulse_rnk(2),parameters.method);
    rho3=shaped_pulse_af(spin_system,L,Ex,Ey,rho2,parameters.pulse_frq(3),parameters.pulse_pwr(3),...
        parameters.pulse_dur(3),parameters.pulse_phi(3),...
        parameters.pulse_rnk(3),parameters.method);
    %Hard Pulse to ratate into xy for detection.
    parameters.rho0=step(spin_system,Ex,[parameters.rho0,rho1,rho2,rho3],pi/2);
    %acquire does the Free Induction Decay Signal. Time dynamics of Detection State (param.coil)
    fids=acquire(spin_system,parameters,H,R,K);
    % now rho0 is a list of multiple density matrices called fids.
    end