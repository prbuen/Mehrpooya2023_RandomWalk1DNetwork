%% Adel Mehrpooya, Vivien J Challis, Pascal R Buenzli (2023)

%% Requirements
% randsample requires the "Statistics and Machine Learning toolbox"

if exist("datadir")==7
    rmdir("datadir", 's')
end
clear all
close all

t_end = 50 % Final time

plot_fluxes = 1;
solve_continuous = 0; % slow

%% Discrete and stochastic models                                                  

%% Parameters
global K Dxi Dt q_left1 q_right1 q_left2 q_right2 N_tot_disc L

% space and time scales
L = 50; % [um]
T = 50; % [s]

% network
K = 61;
q_left1 = 0.16;
q_right1 = 0.16;
q_left2 = 0.16;
q_right2 = 0.16;
q_0 = 1-q_left1-q_right1-q_left2-q_right2;

% signalling molecules
N_tot_disc = 10000; % initial
N_tot_stoch = N_tot_disc;

N_ss = 400; % reaction steady state
lambda1 = 00/T; % 0, 100/T, etc.
lambda2 = lambda1/N_ss;

% signals
N_impulse = 150;
t_signal = 2*T;
i_signal = round((K-1)/2);

% space and time discretisation
Dt = 0.0005*T % discrete model

x_max = 2*L;
x_min = -2*L;
Dxi = (xi_vs_x(x_max) - xi_vs_x(x_min))/(K-1);

x_min = x_vs_i(1);
x_max = x_vs_i(K); % i=1,2,..,K

% flux
T_flux_stoch_history = 0.1*T;
N_flux_stoch_history = round(T_flux_stoch_history/Dt);

% Node territory
Dx = zeros(1,K);
for i = 1:K
    Dx(i) = 0.5*( x_vs_i(i+1) - x_vs_i(i-1) );
end

Dtilde = (((Dxi)^2)/(2*Dt))*(q_left1 + q_right1 + 4*q_left2 + 4*q_right2) % print value for info


%% Initialisation (discrete and stochastic)
t = 0;
iter = 0; % number of time steps taken
save_disc_every = 10; % iterations
save_iter = 0; % counter: number of states saved

% allocate space to save states
N_disc = zeros(1,K); % distribution of molecules (node occupancy)
N_stoch = zeros(1,K); % distribution of molecules (node occupancy)
flux_stoch_left_history = zeros(1,N_flux_stoch_history); % past values so averages can be computed
flux_stoch_right_history = zeros(1,N_flux_stoch_history);
flux_stoch_history_index = 1; % index in flux_stoch_history[] where latest values will be written. Cyclic in {1,...,N_flux_stoch_history}
flux_disc_left = 0;
flux_disc_right = 0;

% Seed the discrete model with molecules in the middle node
i0 = (K-1)/2 + 1; % index of middle node (+1 since indices start at 1)
N_disc(i0) = N_tot_disc;
N_stoch(i0) = N_tot_stoch;


%% Time stepping (discrete and stochastic)
mkdir('datadir')
flux_stoch_left = mean(flux_stoch_left_history);
flux_stoch_right = mean(flux_stoch_right_history);
save('datadir/state_disc_0.mat', 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right') % save data for plotting later
N_disc_old = zeros(size(N_disc)); % allocate, for time stepping
N_stoch_old = zeros(size(N_stoch));
while t < t_end - 0.1*Dt
    % keep old values
    N_disc_old = N_disc;
    N_stoch_old = N_stoch;
    N_stoch_old_left = N_stoch(1); % N_stoch_old is updated mid-way so we need these
    N_stoch_old_right = N_stoch(K);

    % update discrete
    % when molecules reach i=1 or i=K, they stay there forever (absorbing nodes)
    for i=2:K-1
        % fraction jump size -1
        dN = q_left1*N_disc_old(i);
        N_disc(i-1) = N_disc(i-1) + dN;
        N_disc(i) = N_disc(i) - dN;
        
        % fraction jump size 1
        dN = q_right1*N_disc_old(i);
        N_disc(i+1) = N_disc(i+1) + dN;
        N_disc(i) = N_disc(i) - dN;

        % fraction jump size -2
        dN = q_left2*N_disc_old(i);
        N_disc(max(1,i-2)) = N_disc(max(1,i-2)) + dN;
        N_disc(i) = N_disc(i) - dN;
        
        % fraction jump size 2
        dN = q_right2*N_disc_old(i);
        N_disc(min(K,i+2)) = N_disc(min(K,i+2)) + dN;
        N_disc(i) = N_disc(i) - dN;

        % reaction
        N_disc(i) = N_disc(i) + Dt*(lambda1 - lambda2*N_disc_old(i));
    end
    % fluxes
    flux_disc_left = (N_disc(1) - N_disc_old(1))/Dt;
    flux_disc_right = (N_disc(K) - N_disc_old(K))/Dt;

    % update stochastic
    % when molecules reach i=1 or i=N, they stay there forever (absorbing nodes)
    for i=2:K-1
        % for all random molecules at i: move them left, right, or let them stay 
        j=randsample([-2,-1,0,1,2], N_stoch_old(i), true, [q_left2, q_left1, q_0, q_right1, q_right2]);
        for m=1:N_stoch_old(i)
            target = min( max(i+j(m), 1), K); % = i+j(m) capped between 1,..,K
            N_stoch(target) = N_stoch(target) + 1; % add to neighbour j or left/right boundary
            N_stoch(i) = N_stoch(i) - 1; % remove from current node
        end
    end
    % reactions
    N_stoch_old = N_stoch;
    for i=2:K-1
        c = lambda1; % source
        e = lambda2*N_stoch_old(i); % sink
        r = c + e;
        assert(Dt*r < 1, "**** Probability of reaction Dt*r = %f is too high at t = %f", Dt*r, t);

        if rand(1,1) <= Dt*r
            if rand(1,1) <= c/r
                N_stoch(i) = N_stoch(i) + 1;
            else
                N_stoch(i) = N_stoch(i) - 1;
            end
        end
    end
    % fluxes (add to history; moving mean calculated when saving)
    flux_stoch_left_history(flux_stoch_history_index) = (N_stoch(1) - N_stoch_old_left)/Dt;
    flux_stoch_right_history(flux_stoch_history_index) = (N_stoch(K) - N_stoch_old_right)/Dt;
    flux_stoch_history_index = flux_stoch_history_index + 1;
    if flux_stoch_history_index > N_flux_stoch_history
        flux_stoch_history_index = 1;
    end
    
    % update time
    iter = iter + 1;
    t = t + Dt;

    % Save states to files state_1.mat, state_2.mat, etc.
    if rem(iter, save_disc_every) == 0
        flux_stoch_left = mean(flux_stoch_left_history);
        flux_stoch_right = mean(flux_stoch_right_history);
        
        save_iter = save_iter + 1;
        save(sprintf("datadir/state_disc_%d.mat", save_iter), 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right')
    end
end

disp("Discrete and stochastic models: done")



%% Continuum model                                                                       
if solve_continuous
    
    %% Parameters
    % Also some quantities defined in the discrete/stochastic model

    dt = 0.0001*T; % continuum model. Must have dt < 0.5*dx^2*D (CFL condition)
    dx = 0.1*L; % continuum model. Initial guess, may be adapted slightly
    M = round((x_max-x_min)/dx) % number of discretisation intervals
    if rem(M, 2)==0
        M = M + 1;
        assert(rem(M, 2)==1) % M odd so there is no discrete spatial point in middle of domain, where v is discontinuous
    end
    dx = (x_max - x_min)/M;

    % times in discrete and continuum models must match every now and then. Usually dt << Dt is required for precision
    dt_CFL = (dx^2)/(2*max(D(x_min:dx:x_max))) % CFL condition timestep
    Dt_over_dt = 10; % 1,2,3,...,10,...,100 as needed for precision
    dt = Dt/Dt_over_dt
    assert(dt < dt_CFL)

    
    %% Initialisation (continuum)
    t = 0;
    iter = 0; % number of time steps taken
    save_cont_every = save_disc_every*Dt_over_dt; % iterations
    save_iter_cont = 0; % number of states saved
    i0 = (M-1)/2 + 1; % index of middle node (+1 since indices start at 1)

    n = zeros(1,M+1); % number density of molecules
    n(i0) = N_tot_disc/dx;
    n_old = zeros(size(n)); % allocate, for time stepping

    %% Time stepping (continuum)
    mkdir('datadir')
    save('datadir/state_cont_0.mat', 'n') % save data for plotting later

    while t < t_end - 0.1*dt
        % keep old values
        n_old = n;
        
        for i=2:M % Dirichlet BC: don't update boundary values (=0)
                  % space discretisations of the continuous model (convenience)
            xx_i = x_min + (i-1)*dx; % x_i (name clash with function)
            xx_ip1 = xx_i + dx; % x_{i+1}
            xx_im1 = xx_i - dx; % x_{i-1}
            xx_ip0p5 = xx_i + dx/2.; % x_{i+1/2}
            xx_im0p5 = xx_i - dx/2.; % x_{i-1/2}
            
            % advection (1st order upwinding)
            if(v(xx_i) >  0)
                dn_adv = - (v(xx_i)*n_old(i) - v(xx_im1)*n_old(i-1))/dx;
            else
                dn_adv = - (v(xx_ip1)*n_old(i+1) - v(xx_i)*n_old(i))/dx;
            end

            % diffusion (2nd order central diff, known diffusivity function)
            dn_diffus = ( D(xx_ip0p5)*( n_old(i+1)-n_old(i) ) - D(xx_im0p5)*( n_old(i)-n_old(i-1) ) )/dx^2;
            % update n
            n(i) = n_old(i) + dt*dn_adv + dt*dn_diffus;
        end
        
        % update time
        iter = iter + 1;
        t = t + dt;

        % Save states to files state_1.mat, state_2.mat, etc.
        if rem(iter, save_cont_every) == 0
            save_iter_cont = save_iter_cont + 1
            save(sprintf("datadir/state_cont_%d.mat", save_iter_cont), 'n')
        end
    end

    disp("Continuum model: done")
end


%% Plotting                                                                   
h_network = figure;

% plot initial condition
n_max=4800/L; % for plotting
% n_max=1.5*N_ss;
t=0;
x = x_vs_i([1:K]); % node positions
load('datadir/state_disc_0.mat', 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right')
hold on
bar(x, N_disc./Dx, 'c') 
plot(x, (N_tot_disc/N_tot_stoch)*N_stoch./Dx, 'k-') % rescale if N_tot_disc != N_tot_stoch
if solve_continuous
    load('datadir/state_cont_0.mat', 'n')
    plot(x_min:dx:x_max, n, 'g-')
end
fplot(@(x) 0.0, 'r-')
title(sprintf("t=%g", t))
xlim([x_min, x_max])
ylim([0,n_max])
xlabel('x')
ylabel('n(x,t)')
if solve_continuous
    legend('determ', 'stoch', 'cont', '$n_\infty$', 'Interpreter', 'latex')
else
    legend('determ', 'stoch', '$n_\infty$', 'Interpreter', 'latex')
end
hold off
drawnow

for s = 1:save_iter
    t = s*save_disc_every*Dt;
    load(sprintf("datadir/state_disc_%d.mat", s), 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right')

    % pause(0.01) % slow down
    figure(h_network)
    clf
    hold on

    % plot discrete:
    % pl_disc = bar(x, N_disc./(m_disc*Dx), 'c');
    % pl_disc = stairs(x, N_disc./(m_disc*Dx), 'c-');
    for i=1:length(x)
        rectangle('position', [x(i)-(x(i)-x_vs_i(i-1))/2, 0, (x_vs_i(i+1)-x_vs_i(i-1))/2, N_disc(i)/Dx(i)], 'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]);
    end
    pl_disc = fill(NaN, NaN, 'w', 'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]); % for the legend
 
    % plot stochastic:
    pl_stoch = plot(x, (N_tot_disc/N_tot_stoch)*N_stoch./Dx, 'b-'); % rescale if N_tot_stoch != N_tot_disc

    % plot continuous:
    if solve_continuous
        load(sprintf("datadir/state_cont_%d.mat", s), 'n')
        pl_cont = plot(x_min:dx:x_max, n, 'g-');
    end
    
    % plot exact:
    pl_exact = fplot(@(x) p_exact(x,t), [x_min, -1e-4], 'r--','LineWidth', 1.2); % fplot doesn't return a handle
    fplot(@(x) p_exact(x,t), [1e-4, x_max], 'r--', 'LineWidth', 1.2)

    % decoration
    title(sprintf("t=%g", t))
    xlim([x_min, x_max])
    ylim([0,n_max])
    xlabel('x')
    ylabel('density n(x,t)')
    if solve_continuous
        legend([pl_disc, pl_stoch, pl_cont, pl_exact], 'determ', 'stoch', 'FD', '$n_\infty$', 'Interpreter', 'latex')
    else
        legend([pl_disc, pl_stoch, pl_exact], 'determ', 'stoch', '$n_\infty$', 'Interpreter', 'latex')
    end
    hold off

    drawnow
end

if plot_fluxes
    t = 0:save_disc_every*Dt:t_end;
    Jm_disc = zeros(size(t));
    Jm_stoch = zeros(size(t));
    Jp_disc = zeros(size(t));
    Jp_stoch = zeros(size(t));

    for s = 0:save_iter
        fluxes = load(sprintf("datadir/state_disc_%d.mat", s), 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right');
        Jm_disc(s+1) = fluxes.flux_disc_left;
        Jp_disc(s+1) = fluxes.flux_disc_right;
        Jm_stoch(s+1) = fluxes.flux_stoch_left;
        Jp_stoch(s+1) = fluxes.flux_stoch_right;        
    end

    h_fluxes2 = figure;
    xlim([0,t_end])
    ylim([0,18e3/L])
    xlabel('t')
    ylabel('flux')
    hold on
    pl_Jm_disc = plot(t, Jp_disc, '-', 'Color', [0.7,0.7,0.7], 'LineWidth', 2);
    pl_Jp_disc = plot(t, Jm_disc, '--k', 'LineWidth', 2);
    pl_Jm_stoch = plot(t-T_flux_stoch_history/2., Jm_stoch);
    pl_Jp_stoch = plot(t-T_flux_stoch_history/2., Jp_stoch);
    legend([pl_Jm_disc, pl_Jp_disc, pl_Jm_stoch, pl_Jp_stoch], '$J_{-}^{\textrm{disc}}$', '$J_{+}^{\textrm{disc}}$', '$J_{-}^{\textrm{stoch}}$', '$J_{+}^{\textrm{stoch}}$', 'Interpreter', 'latex');
    hold off
end


%% Functions
function x=x_vs_i(i)
    global K Dxi L
    i0 = (K-1)/2 + 1; % index of middle node (+1 since indices start at 1)
    xi = (i-i0)*Dxi;
    % x = xi; % regular lattice
    x = L*sign(xi).*log(1+abs(xi)); % log lattice
end

function xi=xi_vs_x(x)
    global L
    % xi = x; % regular lattice
    xi = sign(x).*(exp(abs(x)/L)-1); % log lattice
end

function y=g(x) % metric
    global L
    % y = 1; % regular lattice
    y = L*exp(-abs(x)/L); % log lattice
end

function y=gg_x(x)
    global L
    % y = 0.; % regular lattice
    y = -L*sign(x).*exp(-2*abs(x)/L); % log lattice
end

function y=D(x) % Diffusion in x space
    global Dxi Dt q_left1 q_right1 q_left2 q_right2;
    y=(((Dxi)^2)/(2*Dt)) .* (q_left1 + q_right1 + 4*q_left2 + 4*q_right2) .* (g(x).^2);
end

function y=v(x) % advective velocity
    global Dxi Dt q_left1 q_right1 q_left2 q_right2

    avg_s = q_right1 - q_left1 + 2*q_right2 - 2*q_left2;
    avg_s2 = q_left1 + q_right1 + 4*q_left2 + 4*q_right2;
    
    y = (Dxi/Dt) .* g(x) .* avg_s - (Dxi^2/Dt).* gg_x(x) .* avg_s2/2.;
end

function y=p_exact(x,t)
    global Dxi Dt q_left1 q_right1 q_left2 q_right2 N_tot_disc
    Dtilde = (((Dxi)^2)/(2*Dt))*(q_left1 + q_right1 + 4*q_left2 + 4*q_right2);
    vtilde = (Dxi/Dt)*(q_right1-q_left1 + 2*(q_right2 - q_left2));
    y = N_tot_disc*(1./sqrt(4*pi*Dtilde.*t)).*exp( - (xi_vs_x(x) - vtilde.*t).^2./(4*Dtilde.*t))./g(x);
end
