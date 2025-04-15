function inverted_pendulum_collatz_iterations_plot_vertical()
  close all;

  % --- System Parameters ---
  M = 1;
  m = 0.1;
  l = 1;
  g = 9.81;
  I = 0.006;

  % --- Simulation Parameters ---
  T = 10;
  N = 10000;
  h = T / N;
  t = 0:h:T;

  % --- Collatz PID Parameters ---
  Kp0_collatz = 50;
  alpha_collatz = 0.1;
  Ki_collatz = 1;
  Kd_collatz = 20;

  % --- Function Definitions ---
  fz = @(t, z, M, m, l, g, I, F) [
    z(2);
    ((m * l) / ((M + m) * (I + m * l^2) - (m * l)^2)) * (-m * g * l * sin(z(3)) + m * l * z(4)^2 * sin(z(3)) * cos(z(3))) + ((M + m) / ((M + m) * (I + m * l^2) - (m * l)^2)) * F;
    z(4);
    ((I + m * l^2) / ((M + m) * (I + m * l^2) - (m * l)^2)) * (-m * g * l * sin(z(3)) + m * l * z(4)^2 * sin(z(3)) * cos(z(3))) - (m * l / ((M + m) * (I + m * l^2) - (m * l)^2)) * F
  ];

  % --- Iteration Setup ---
  num_iterations = 100;
  max_deviations = zeros(1, num_iterations);
  steady_state_errors = zeros(1, num_iterations);
  iterations = 1:num_iterations; % Vector of iteration numbers

  % --- Main Loop for Iterations ---
  for iteration = 1:num_iterations
    n0_collatz = iteration;

    % --- Initial Conditions (reset for each iteration) ---
    z0_pid_collatz = [0; 0; pi / 36; 0];
    z_pid_collatz = zeros(4, length(t));
    z_pid_collatz(:, 1) = z0_pid_collatz;
    u_collatz = zeros(1, length(t));
    int_e_collatz = 0;
    n_collatz = zeros(1, N + 1);
    n_collatz(1) = n0_collatz;

    % --- Collatz PID Control Simulation ---
    for n = 1:N
      phi = z_pid_collatz(3, n);
      phi_dot = z_pid_collatz(4, n);
      e_collatz = phi - 0;
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

      u_collatz(n) = Kp_collatz * e_collatz + Ki_collatz * int_e_collatz + Kd_collatz * der_e_collatz;

      k1 = h * fz(t(n), z_pid_collatz(:, n), M, m, l, g, I, u_collatz(n));
      k2 = h * fz(t(n) + h / 2, z_pid_collatz(:, n) + k1 / 2, M, m, l, g, I, u_collatz(n));
      k3 = h * fz(t(n) + h / 2, z_pid_collatz(:, n) + k2 / 2, M, m, l, g, I, u_collatz(n));
      k4 = h * fz(t(n) + h, z_pid_collatz(:, n) + k3, M, m, l, g, I, u_collatz(n));
      z_pid_collatz(:, n + 1) = z_pid_collatz(:, n) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end

    % --- Calculate and Store Results ---
    max_deviations(iteration) = max(abs(z_pid_collatz(3, :)));
    steady_state_errors(iteration) = abs(0 - z_pid_collatz(3, end));

    if mod(iteration, 1000) == 0
      fprintf('Iteration: %d\n', iteration);
    end
  end

  % --- Plotting ---
  figure;

  % Max Deviation vs. Iteration (Iterations on y-axis)
  subplot(2, 1, 1);
  plot(max_deviations, iterations, 'b'); 
  ylabel('Iteration (Initial Collatz Value)');
  xlabel('Max Deviation of Phi (rad)');
  title('Iteration vs. Max Deviation');
  grid on;
  

  % Steady-State Error vs. Iteration (Iterations on y-axis)
  subplot(2, 1, 2);
  plot(steady_state_errors, iterations, 'r');  
  ylabel('Iteration (Initial Collatz Value)');
  xlabel('Steady-State Error of Phi (rad)');
  title('Iteration vs. Steady-State Error');
  grid on;
 

  sgtitle('Collatz PID Performance: Iterations on Y-Axis');
end

% Run the simulation
