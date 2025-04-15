function inverted_pendulum_simulation_collatz()
  close all;

  % --- System Parameters ---
  M = 1;      % kg, cart mass
  m = 0.1;    % kg, pendulum mass
  l = 1;      % m, pendulum length
  g = 9.81;   % m/s^2, gravity
  I = 0.006;  % kg*m^2, pendulum inertia

  % --- Simulation Parameters ---
  T = 10;     % seconds, simulation time
  N = 1000;   % Number of steps
  h = T / N;  % time step
  t = 0:h:T;  % time vector

  % --- Initial Conditions ---
  z0_no_control = [0; 0; pi / 6; 0]; % [x; x_dot; phi; phi_dot]
  z0_pid_std = [0; 0; pi / 36; 0];    % [x; x_dot; phi; phi_dot] for standard PID
  z0_pid_collatz = [0; 0; pi / 36; 0]; % [x; x_dot; phi; phi_dot] for Collatz PID

  % --- Standard PID Gains ---
  Kp_std = 100;
  Ki_std = 1;
  Kd_std = 20;

  % --- Collatz PID Parameters ---
  Kp0_collatz = 50;   % Base Kp
  alpha_collatz = 0.1; % Scaling factor for Collatz
  Ki_collatz = 1;
  Kd_collatz = 20;
  n0_collatz = 24000;  % Initial Collatz value

  % --- Function Definitions ---
  fz = @(t, z, M, m, l, g, I, F) [
    z(2);
    ((m * l) / ((M + m) * (I + m * l^2) - (m * l)^2)) * (-m * g * l * sin(z(3)) + m * l * z(4)^2 * sin(z(3)) * cos(z(3))) + ((M + m) / ((M + m) * (I + m * l^2) - (m * l)^2)) * F;
    z(4);
    ((I + m * l^2) / ((M + m) * (I + m * l^2) - (m * l)^2)) * (-m * g * l * sin(z(3)) + m * l * z(4)^2 * sin(z(3)) * cos(z(3))) - (m * l / ((M + m) * (I + m * l^2) - (m * l)^2)) * F
  ];

  % --- Simulation Without Control ---
  z_no_control = zeros(4, length(t));
  z_no_control(:, 1) = z0_no_control;

  for n = 1:N
    k1 = h * fz(t(n), z_no_control(:, n), M, m, l, g, I, 0);
    k2 = h * fz(t(n) + h / 2, z_no_control(:, n) + k1 / 2, M, m, l, g, I, 0);
    k3 = h * fz(t(n) + h / 2, z_no_control(:, n) + k2 / 2, M, m, l, g, I, 0);
    k4 = h * fz(t(n) + h, z_no_control(:, n) + k3, M, m, l, g, I, 0);
    z_no_control(:, n + 1) = z_no_control(:, n) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  end

  % --- Standard PID Control Simulation ---
  z_pid_std = zeros(4, length(t));
  z_pid_std(:, 1) = z0_pid_std;
  u_std = zeros(1, length(t)); % Control input
  int_e_std = 0;             % Integral of error

  for n = 1:N
    phi = z_pid_std(3, n);
    phi_dot = z_pid_std(4, n);
    e_std = phi - 0;          % Error (target phi = 0)
    int_e_std = int_e_std + e_std * h;
    der_e_std = phi_dot;

    u_std(n) = Kp_std * e_std + Ki_std * int_e_std + Kd_std * der_e_std;

    k1 = h * fz(t(n), z_pid_std(:, n), M, m, l, g, I, u_std(n));
    k2 = h * fz(t(n) + h / 2, z_pid_std(:, n) + k1 / 2, M, m, l, g, I, u_std(n));
    k3 = h * fz(t(n) + h / 2, z_pid_std(:, n) + k2 / 2, M, m, l, g, I, u_std(n));
    k4 = h * fz(t(n) + h, z_pid_std(:, n) + k3, M, m, l, g, I, u_std(n));
    z_pid_std(:, n + 1) = z_pid_std(:, n) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  end

  % --- Collatz PID Control Simulation ---
  z_pid_collatz = zeros(4, length(t));
  z_pid_collatz(:, 1) = z0_pid_collatz;
  u_collatz = zeros(1, length(t)); % Control input
  int_e_collatz = 0;             % Integral of error
  Kp_collatz_vec = zeros(1, length(t)); % Store Kp values
  n_collatz = zeros(1, length(t));
  n_collatz(1) = n0_collatz;

  for n = 1:N
    phi = z_pid_collatz(3, n);
    phi_dot = z_pid_collatz(4, n);
    e_collatz = phi - 0;          % Error (target phi = 0)
    int_e_collatz = int_e_collatz + e_collatz * h;
    der_e_collatz = phi_dot;

    % Collatz Sequence Update
    if mod(n_collatz(n), 2) == 0
      n_collatz(n + 1) = n_collatz(n) / 2;
    else
      n_collatz(n + 1) = 3 * n_collatz(n) + 1;
    end

    % Collatz-Modulated Kp
    Kp_collatz = Kp0_collatz + alpha_collatz * n_collatz(n);
    Kp_collatz_vec(n) = Kp_collatz; % Store Kp value

    u_collatz(n) = Kp_collatz * e_collatz + Ki_collatz * int_e_collatz + Kd_collatz * der_e_collatz;

    k1 = h * fz(t(n), z_pid_collatz(:, n), M, m, l, g, I, u_collatz(n));
    k2 = h * fz(t(n) + h / 2, z_pid_collatz(:, n) + k1 / 2, M, m, l, g, I, u_collatz(n));
    k3 = h * fz(t(n) + h / 2, z_pid_collatz(:, n) + k2 / 2, M, m, l, g, I, u_collatz(n));
    k4 = h * fz(t(n) + h, z_pid_collatz(:, n) + k3, M, m, l, g, I, u_collatz(n));
    z_pid_collatz(:, n + 1) = z_pid_collatz(:, n) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  end

  % --- Plotting ---
  figure;

  % No Control
  subplot(4, 2, 1);
  plot(t, z_no_control(1, :), 'b', 'LineWidth', 1.5);
  ylabel('x (m)');
  title('No Control - Cart Position');

  subplot(4, 2, 3);
  plot(t, z_no_control(3, :), 'r', 'LineWidth', 1.5);
  ylabel('\phi (rad)');
  xlabel('t (s)');
  title('No Control - Pendulum Angle');

  % Standard PID Control
  subplot(4, 2, 2);
  plot(t, z_pid_std(1, :), 'b');
  ylabel('x (m)');
  title('Standard PID - Cart Position');

  subplot(4, 2, 4);
  plot(t, z_pid_std(3, :), 'r');
  ylabel('\phi (rad)');

  subplot(4, 2, 6);
  plot(t, u_std, 'k');
  ylabel('F (N)');
  xlabel('t (s)');
  title('Standard PID - Force');

  % Collatz PID Control
  subplot(4, 2, 5);
  plot(t, z_pid_collatz(1, :), 'b');
  ylabel('x (m)');
  title('Collatz PID - Cart Position');

  subplot(4, 2, 7);
  plot(t, z_pid_collatz(3, :), 'r');
  ylabel('\phi (rad)');

  subplot(4, 2, 8);
  plot(t, u_collatz, 'k');
  ylabel('F (N)');
  xlabel('t (s)');
  title('Collatz PID - Force');

  % Collatz Kp
  subplot(4, 2, 6); % Overlapping with Standard PID Force
  plot(t, Kp_collatz_vec, 'm');
  ylabel('Kp');
  xlabel('t (s)');
  title('Collatz PID - Kp');

  % Overall Title
  sgtitle('Inverted Pendulum Simulation Comparison');
    % --- Calculate Performance Metrics ---
  % Overshoot
 % --- Calculate Performance Metrics ---
% Overshoot (using max absolute deviation)
max_deviation_std = max(abs(z_pid_std(3,:)));
max_deviation_collatz = max(abs(z_pid_collatz(3,:)));

% Steady-State Error (approximate) - Using the last value
steady_state_error_std = abs(0 - z_pid_std(3,end));
steady_state_error_collatz = abs(0 - z_pid_collatz(3,end));

% --- Display Results ---
disp('--- Performance Metrics ---');
disp('Max Deviation (Standard PID):');
disp(max_deviation_std);
disp('Max Deviation (Collatz PID):');
disp(max_deviation_collatz);
disp('Steady-State Error (Standard PID):');
disp(steady_state_error_std);
disp('Steady-State Error (Collatz PID):');
disp(steady_state_error_collatz);


end

% Run the simulation
inverted_pendulum_simulation_collatz();