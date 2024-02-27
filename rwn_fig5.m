%% Adel Mehrpooya, Vivien J Challis, Pascal R Buenzli (2023)

%% Requirements
% randsample requires the "Statistics and Machine Learning toolbox"

if exist("datadir")==7
    rmdir("datadir", 's')
end
clear all
close all

T = 2 % Final time

plot_fluxes = 0;
solve_continuous = 1; % slow

global K Dxi Dt D0 N_tot_disc x_min x_max dx xi_tabulated

%% Continuum model parameters used in the discrete model
%The inverse function xi(x) is tabulated in xi_tabulated for half integer steps dx; Here are the results for M=150 intervals (M+1 points)
xi_tabulated = [-5, -4.98611, -4.97221, -4.9583, -4.94436, -4.9304, -4.91639, -4.90234, -4.88824, -4.87407, -4.85983, -4.84551, -4.8311, -4.81659, -4.80197, -4.78722, -4.77235, -4.75733, -4.74214, -4.72679, -4.71125, -4.69551, -4.67954, -4.66333, -4.64686, -4.6301, -4.61303, -4.59562, -4.57784, -4.55965, -4.54101, -4.52187, -4.50219, -4.48189, -4.46092, -4.43918, -4.41657, -4.39298, -4.36826, -4.34223, -4.31466, -4.28525, -4.25361, -4.21923, -4.18135, -4.13889, -4.09018, -4.03246, -3.961, -3.86809, -3.75, -3.63191, -3.539, -3.46754, -3.40982, -3.36111, -3.31865, -3.28077, -3.24639, -3.21475, -3.18534, -3.15777, -3.13174, -3.10702, -3.08343, -3.06082, -3.03908, -3.01811, -2.99781, -2.97813, -2.95899, -2.94035, -2.92216, -2.90438, -2.88697, -2.8699, -2.85314, -2.83667, -2.82046, -2.80449, -2.78875, -2.77321, -2.75786, -2.74267, -2.72765, -2.71278, -2.69803, -2.68341, -2.6689, -2.65449, -2.64017, -2.62593, -2.61176, -2.59766, -2.58361, -2.5696, -2.55564, -2.5417, -2.52779, -2.51389, -2.5, -2.48611, -2.47221, -2.4583, -2.44436, -2.4304, -2.41639, -2.40234, -2.38824, -2.37407, -2.35983, -2.34551, -2.3311, -2.31659, -2.30197, -2.28722, -2.27235, -2.25733, -2.24214, -2.22679, -2.21125, -2.19551, -2.17954, -2.16333, -2.14686, -2.1301, -2.11303, -2.09562, -2.07784, -2.05965, -2.04101, -2.02187, -2.00219, -1.98189, -1.96092, -1.93918, -1.91657, -1.89298, -1.86826, -1.84223, -1.81466, -1.78525, -1.75361, -1.71923, -1.68135, -1.63889, -1.59018, -1.53246, -1.461, -1.36809, -1.25, -1.13191, -1.039, -0.967538, -0.90982, -0.861108, -0.818652, -0.780775, -0.746389, -0.714754, -0.685344, -0.657769, -0.631736, -0.607015, -0.583426, -0.560822, -0.539082, -0.518108, -0.497814, -0.47813, -0.458995, -0.440354, -0.422163, -0.404381, -0.38697, -0.369901, -0.353143, -0.336671, -0.320462, -0.304495, -0.288749, -0.273208, -0.257855, -0.242675, -0.227653, -0.212777, -0.198034, -0.183412, -0.168901, -0.15449, -0.140169, -0.12593, -0.111762, -0.0976567, -0.0836062, -0.069602, -0.0556361, -0.0417006, -0.0277878, -0.0138901, 1.54198e-16, 0.0138901, 0.0277878, 0.0417006, 0.0556361, 0.069602, 0.0836062, 0.0976567, 0.111762, 0.12593, 0.140169, 0.15449, 0.168901, 0.183412, 0.198034, 0.212777, 0.227653, 0.242675, 0.257855, 0.273208, 0.288749, 0.304495, 0.320462, 0.336671, 0.353143, 0.369901, 0.38697, 0.404381, 0.422163, 0.440354, 0.458995, 0.47813, 0.497814, 0.518108, 0.539082, 0.560822, 0.583426, 0.607015, 0.631736, 0.657769, 0.685344, 0.714754, 0.746389, 0.780775, 0.818652, 0.861108, 0.90982, 0.967538, 1.039, 1.13191, 1.25, 1.36809, 1.461, 1.53246, 1.59018, 1.63889, 1.68135, 1.71923, 1.75361, 1.78525, 1.81466, 1.84223, 1.86826, 1.89298, 1.91657, 1.93918, 1.96092, 1.98189, 2.00219, 2.02187, 2.04101, 2.05965, 2.07784, 2.09562, 2.11303, 2.1301, 2.14686, 2.16333, 2.17954, 2.19551, 2.21125, 2.22679, 2.24214, 2.25733, 2.27235, 2.28722, 2.30197, 2.31659, 2.3311, 2.34551, 2.35983, 2.37407, 2.38824, 2.40234, 2.41639, 2.4304, 2.44436, 2.4583, 2.47221, 2.48611, 2.5, 2.51389, 2.52779, 2.5417, 2.55564, 2.5696, 2.58361, 2.59766, 2.61176, 2.62593, 2.64017, 2.65449, 2.6689, 2.68341, 2.69803, 2.71278, 2.72765, 2.74267, 2.75786, 2.77321, 2.78875, 2.80449, 2.82046, 2.83667, 2.85314, 2.8699, 2.88697, 2.90438, 2.92216, 2.94035, 2.95899, 2.97813, 2.99781, 3.01811, 3.03908, 3.06082, 3.08343, 3.10702, 3.13174, 3.15777, 3.18534, 3.21475, 3.24639, 3.28077, 3.31865, 3.36111, 3.40982, 3.46754, 3.539, 3.63191, 3.75, 3.86809, 3.961, 4.03246, 4.09018, 4.13889, 4.18135, 4.21923, 4.25361, 4.28525, 4.31466, 4.34223, 4.36826, 4.39298, 4.41657, 4.43918, 4.46092, 4.48189, 4.50219, 4.52187, 4.54101, 4.55965, 4.57784, 4.59562, 4.61303, 4.6301, 4.64686, 4.66333, 4.67954, 4.69551, 4.71125, 4.72679, 4.74214, 4.75733, 4.77235, 4.78722, 4.80197, 4.81659, 4.8311, 4.84551, 4.85983, 4.87407, 4.88824, 4.90234, 4.91639, 4.9304, 4.94436, 4.9583, 4.97221, 4.98611, 5];

M = 150
x_max = 5;
x_min = -x_max;
dx = (x_max - x_min)/M; % continuum model


%% Discrete and stochastic models                                                  

%% Parameters

% network
% K = 51;
K = 101;

q_left2 = 0.;
q_right2 = 0.;
D0 = 1;

% signalling molecules
N_tot_disc = 6000; % initial
N_tot_stoch = N_tot_disc;

N_ss = 0; % reaction steady state
lambda1 = 0;
lambda2 = 0;

% signals
N_impulse = 0;
t_signal = -1;
i_signal = round((K-1)/2);

% space and time discretisation
Dt = 0.0001 % discrete model

x_max = 5;
x_min = -x_max;
xi_min = x_min;
xi_max = x_max;
Dxi = (xi_max - xi_min)/(K-1)

x_min = x_vs_i(1)
x_max = x_vs_i(K) % i=1,2,..,K

xi_min = xi_vs_x(x_min); % check same as x_min, x_max
xi_max = xi_vs_x(x_max);

% flux
T_flux_stoch_history = 0.1;
N_flux_stoch_history = round(T_flux_stoch_history/Dt);

% Node territory
Dx = zeros(1,K);
for i = 1:K
    Dx(i) = 0.5*( x_vs_i(i+1) - x_vs_i(i-1) );
end


%% Initialisation (discrete and stochastic)
t = 0;
iter = 0; % number of time steps taken
save_disc_every = 1000; % iterations
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
N_disc = Dx.*ones(1,K)*(N_tot_disc*pi/(4*x_max)).*cos(pi*x_vs_i(1:K)/(2*x_max)); % distribution of molecules (node occupancy)
N_stoch = round(Dx.*ones(1,K)*(N_tot_disc*pi/(4*x_max)).*cos(pi*x_vs_i(1:K)/(2*x_max))); % distribution of molecules (node occupancy)

mkdir('datadir')


%% Time stepping (discrete and stochastic)
flux_stoch_left = mean(flux_stoch_left_history);
flux_stoch_right = mean(flux_stoch_right_history);
save('datadir/state_disc_0.mat', 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right') % save data for plotting later
N_disc_old = zeros(size(N_disc)); % allocate, for time stepping
N_stoch_old = zeros(size(N_stoch));
while t < T - 0.1*Dt
    % keep old values
    N_disc_old = N_disc;
    N_stoch_old = N_stoch;
    N_stoch_old_left = N_stoch(1); % N_stoch_old is updated mid-way so we need these
    N_stoch_old_right = N_stoch(K);

    % update discrete
    % when molecules reach i=1 or i=K, they stay there forever (absorbing nodes)
    for i=2:K-1
        xi = xi_min + (i-1)*Dxi;
        
        % fraction jump size -1
        dN = q_left1(xi).*N_disc_old(i);
        N_disc(i-1) = N_disc(i-1) + dN;
        N_disc(i) = N_disc(i) - dN;
        
        % fraction jump size 1
        dN = q_right1(xi).*N_disc_old(i);
        N_disc(i+1) = N_disc(i+1) + dN;
        N_disc(i) = N_disc(i) - dN;

        % % fraction jump size -2
        % dN = q_left2*N_disc_old(i);
        % N_disc(max(1,i-2)) = N_disc(max(1,i-2)) + dN;
        % N_disc(i) = N_disc(i) - dN;
        
        % % fraction jump size 2
        % dN = q_right2*N_disc_old(i);
        % N_disc(min(K,i+2)) = N_disc(min(K,i+2)) + dN;
        % N_disc(i) = N_disc(i) - dN;

        % reaction
        N_disc(i) = N_disc(i) + Dt*(lambda1 - lambda2*N_disc_old(i));
    end
    % fluxes
    flux_disc_left = (N_disc(1) - N_disc_old(1))/Dt;
    flux_disc_right = (N_disc(K) - N_disc_old(K))/Dt;

    % update stochastic
    % when molecules reach i=1 or i=N, they stay there forever (absorbing nodes)
    for i=2:K-1
        xi = xi_min + (i-1)*Dxi;
        % for all random molecules at i: move them left, right, or let them stay
        j=randsample([-2,-1,0,1,2], N_stoch_old(i), true, [q_left2, q_left1(xi), q_0(xi), q_right1(xi), q_right2]);
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

    dt = 0.000002; % continuum model. Must have dt < 0.5*dx^2*D (CFL condition)

    % times in discrete and continuum models must match every now and then. Usually dt << Dt is required for precision
    % dt_CFL = (dx^2)/(2*max(D(x_min:dx:x_max))) % CFL condition timestep. Can only check after xi_tabulated is defined
    Dt_over_dt = 20; % 1,2,3,...,10,...,100 as needed for precision
    dt = Dt/Dt_over_dt
    % assert(dt < dt_CFL)

    
    %% Initialisation (continuum)
    t = 0;
    iter = 0; % number of time steps taken
    save_cont_every = save_disc_every*Dt_over_dt; % iterations
    save_iter_cont = 0; % number of states saved
    % i0 = M/2+1; % index of middle node 

    n = cos(pi*(x_min+[0:M]*dx)/(2*x_max))*N_tot_disc*pi/(4*x_max); % number density of molecules; M+1 nodes incl. boundaries
    n_old = zeros(size(n)); % allocate, for time stepping

    %% Time stepping (continuum)
    mkdir('datadir')
    save('datadir/state_cont_0.mat', 'n') % save data for plotting later

    while t < T - 0.1*dt
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
n_max=N_tot_disc/5.; % for plotting
% n_max=1.5*N_ss;
t=0;
x = x_vs_i([1:K]); % node positions
load('datadir/state_disc_0.mat', 'N_disc', 'N_stoch', 'flux_disc_left', 'flux_disc_right', 'flux_stoch_left', 'flux_stoch_right')
hold on
% bar(x, N_disc./Dx, 'c') 
for i=1:length(x)
    rectangle('position', [x(i)-(x(i)-x_vs_i(i-1))/2, 0, (x_vs_i(i+1)-x_vs_i(i-1))/2, N_disc(i)/Dx(i)], 'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]);
end
pl_disc = fill(NaN, NaN, 'w', 'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]); % for the legend
plot(x, (N_tot_disc/N_tot_stoch)*N_stoch./Dx, 'k-') % rescale if N_tot_disc != N_tot_stoch
if solve_continuous
    load('datadir/state_cont_0.mat', 'n')
    plot(x_min:dx:x_max, n, 'g-')
end
% fplot(@(x) 0.0, 'r-')
fplot(@(x) p_exact(x,0), 'r-')
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
    t = 0:save_disc_every*Dt:T;
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
    xlim([0,T])
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
    global K Dxi
    i0 = (K-1)/2 + 1; % index of middle node (+1 since indices start at 1)
    xi = (i-i0)*Dxi;
    % x = xi; % regular lattice
    % x = tanh(xi); % tanh lattice
    x = xi + (1/pi)*sin(4*pi.*xi/5); % sin lattice [=xi+A*sin(k*xi), k=4pi/5, A=1/pi]
    % x = sign(xi).*log(1+abs(xi)); % log lattice
end

function xi=xi_vs_x(x)
    global x_min dx xi_tabulated
    
    % xi = x; % regular lattice
    % xi = 0.5*log((1+x)/(1-x)); % tanh lattice

    % sin lattice:
    % Use xi_tabulated in which xi(x) is tabulated for half integer steps dx; i2 is an integer that corresponds to 2*i (i is a half integer)
    i2 = round(2*(x-x_min)/dx)+1; % +1 since indices start at 1
    xi = xi_tabulated(i2);

    % xi = sign(x).*(exp(abs(x))-1); % log lattice
end

function y=g(x) % metric
    % y = 1; % regular lattice
    % y = 1-x.^2; % tanh lattice
    y = 1 + 0.8*cos(4*pi*xi_vs_x(x)/5); % sin lattice [= 1+k*A*cos(k*xi(x)), k=4pi/5, kA=8/10]
    % y = exp(-abs(x)); % log lattice
end

function y=gg_x(x)
    % y = 0.; % regular lattice
    % y = -2*x.*(1-x.^2); % tanh lattice
    y = -0.8*(4*pi/5)*sin(4*pi*xi_vs_x(x)/5); % sin lattice [= -k^2*A*sin(k*xi(x))]
    % y = -sign(x).*exp(-2*abs(x)); % log lattice
end

function y=D(x) % Diffusion in x space
    global D0;
    y=D0;
end

function y=v(x) % advective velocity
    y=0;
end

function g = g_vs_xi(xi)
    kA = 0.8;
    k = 4*pi/5;
    g = 1 + kA*cos(k*xi);
end

function ggx = gg_x_vs_xi(xi)
    kA = 0.8;
    k = 4*pi/5;
    ggx = -k*kA*sin(k*xi);
end
    
function q0 = q_0(xi)
    global D0 Dxi Dt
    q0 = 1 - 2*Dt*D0./(g_vs_xi(xi).*Dxi).^2;
end

function qm1 = q_left1(xi)
    global D0 Dxi Dt
    gx = gg_x_vs_xi(xi)./g_vs_xi(xi);
    qm1 = (Dt*D0./(g_vs_xi(xi)*Dxi).^2) .* (1 + 0.5*gx*Dxi);
end

function qp1 = q_right1(xi)
    global D0 Dxi Dt
    gx = gg_x_vs_xi(xi)./g_vs_xi(xi);
    qp1 = (Dt*D0./(g_vs_xi(xi)*Dxi).^2) .* (1 - 0.5*gx*Dxi);
end

function y=p_exact(x,t)
    global x_max D0 N_tot_disc
    y = (N_tot_disc*pi/(4*x_max))*cos(pi*x/(2*x_max)).*exp( - (pi/(2*x_max))^2*D0*t);
end
