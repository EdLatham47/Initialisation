function OneScript_Plots(spin_system,parameters)
    % FIRST             -- Free induction Decay, Each pulse from a rho0 state.
    % Checks out the individual pulses wihtout the mess of evolutino in the middle. 
    %fids = powder(spin_system, @OneScript_IndependantPulse, parameters, parameters.assumptions );
    % <https://spindynamics.org/wiki/index.php?title=Apodization.m> 
    % Add more of these for longer pulse sequences :)))
    %fid_1=apodization(fids(:,1),'crisp-1d'); spectrum_1=fftshift(fft(fid_1,parameters.zerofill));
    %fid_2=apodization(fids(:,2),'crisp-1d'); spectrum_2=fftshift(fft(fid_2,parameters.zerofill));
    %fid_3=apodization(fids(:,3),'crisp-1d'); spectrum_3=fftshift(fft(fid_3,parameters.zerofill));
    %fid_4=apodization(fids(:,4),'crisp-1d'); spectrum_4=fftshift(fft(fid_4,parameters.zerofill));

    % Plotting those spectrums.
    %figure(); ylabel('Magnetisation, arb.u.'); xlabel('frequency, MHz');
    %subplot(2,2,1); plot_1d(spin_system, real(spectrum_1), parameters,'r-'); title('Freq Swp EPR')
    %subplot(2,2,2); plot_1d(spin_system, real(spectrum_2), parameters,'r-'); title('Spec aft P1')
    %subplot(2,2,3); plot_1d(spin_system, real(spectrum_3), parameters,'r-'); title('Spec aft P2')
    %subplot(2,2,4); plot_1d(spin_system, real(spectrum_4), parameters,'r-'); title('Spec aft P3')
    %drawnow();
    
    % Second Sim. complete Successive pulses with evolution in middle.
    echo_stack = powder(spin_system, @OneScript_PulseSequence, parameters, parameters.assumptions );

    echo_axis=1e9*linspace(-parameters.echo_time/2,parameters.echo_time/2, parameters.echo_nsteps+1);
    % time axis for the Deer trace.
    %deer_axis=1e6*linspace(0,parameters.p1_p3_gap, parameters.p2_nsteps+1);
    deer_axis=1e6*linspace(0,parameters.p1_p3_gap, parameters.p2_nsteps+1);
    [deer_axis_2d,echo_axis_2d]=meshgrid(deer_axis,echo_axis);

    % Plotting the echo stack.
    DeerCut=round(parameters.p2_nsteps*parameters.DeerCut+1):round(parameters.p2_nsteps*(1-parameters.DeerCut));
    figure(); surf(deer_axis_2d(:,DeerCut),echo_axis_2d(:,DeerCut),real(echo_stack(:,DeerCut)));
    title('echo stack, unphased'); ylabel('echo window, ns');
    xlabel('2nd pulse position \mu s'); axis tight; kgrid;

    %Extract and phase the echo modulation
    [deer_echoes, deer_sigmas, deer_traces] = svd(echo_stack);
    deer_echoes = deer_echoes*deer_sigmas/deer_traces(1);
    deer_traces=deer_traces*deer_sigmas/deer_traces(1);

    %Plot echo components
    figure(); plot(echo_axis,real(deer_echoes(:,1:3)));
    xlabel('echo window, ns'), axis tight;
    ylabel('echo, real channel'); kgrid; 
    title('Principal components of the echo stack');
    legend({'u_1\cdot\sigma_1',...
        'u_2\cdot\sigma_2',...
        'u_3\cdot\sigma_3'});
    
    % Plot Deer components.
    figure(); plot(deer_axis, real(deer_traces(:,1:3)));
    xlabel('2nd pulse insertion points \mu s'), axis tight;
    ylabel('echo, real channel'); kgrid;
    title('Principal components of the DEER stack');
    legend({'v_1\cdot\sigma_1',...
        'v_2\cdot\sigma_2',...
        'v_3\cdot\sigma_3'});
    end