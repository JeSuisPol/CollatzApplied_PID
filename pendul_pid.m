function inverted_pendulum_simulation()
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
  z0_pid = [0; 0; pi / 36; 0];       % [x; x_dot; phi; phi_dot]

  % --- PID Control Gains ---
  Kp = 100;
  Ki = 1;
  Kd = 20;

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

  % --- PID Control Simulation ---
  z_pid = zeros(4, length(t));
  z_pid(:, 1) = z0_pid;
  u = zeros(1, length(t)); % Control input
  int_e = 0;             % Integral of error

  for n = 1:N
    phi = z_pid(3, n);
    phi_dot = z_pid(4, n);
    e = phi - 0;          % Error (target phi = 0)
    int_e = int_e + e * h;
    der_e = phi_dot;

    u(n) = Kp * e + Ki * int_e + Kd * der_e;

    k1 = h * fz(t(n), z_pid(:, n), M, m, l, g, I, u(n));
    k2 = h * fz(t(n) + h / 2, z_pid(:, n) + k1 / 2, M, m, l, g, I, u(n));
    k3 = h * fz(t(n) + h / 2, z_pid(:, n) + k2 / 2, M, m, l, g, I, u(n));
    k4 = h * fz(t(n) + h, z_pid(:, n) + k3, M, m, l, g, I, u(n));
    z_pid(:, n + 1) = z_pid(:, n) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  end

  % --- Plotting ---
  figure;

  % No Control
  subplot(3, 2, 1);
  plot(t, z_no_control(1, :), 'b', 'LineWidth', 1.5);
  ylabel('x (m)');
  title('No Control - Cart Position');

  subplot(3, 2, 3);
  plot(t, z_no_control(3, :), 'r', 'LineWidth', 1.5);
  ylabel('\phi (rad)');
  xlabel('t (s)');
  title('No Control - Pendulum Angle');

  % PID Control
  subplot(3, 2, 2);
  plot(t, z_pid(1, :), 'b');
  ylabel('x (m)');
  title('PID Control - Cart Position');

  subplot(3, 2, 4);
  plot(t, z_pid(3, :), 'r');
  ylabel('\phi (rad)');

  subplot(3, 2, 6);
  plot(t, u, 'k');
  ylabel('F (N)');
  xlabel('t (s)');
  title('PID Control - Force');

  % Overall Title
  sgtitle('Inverted Pendulum Simulation Comparison');
end

% Run the simulation
inverted_pendulum_simulation();