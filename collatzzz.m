function compare_pid_controllers_iterative()
  % COMPARE_PID_CONTROLLERS_ITERATIVE compares standard PID and Collatz-modulated PID
  % control over multiple iterations with varying initial Collatz values and saves
  % the results to a CSV file.

  % --- PID Controller Parameters (unchanged) ---
  Kp_std = 1;
  Ki_std = 1;
  Kd_std = 1;
  Kp0_collatz = 1;
  alpha_collatz = 0.1;
  Ki_collatz = 1;
  Kd_collatz = 1;

  % --- Simulation Parameters (mostly unchanged) ---
  sim_time = 50;
  dt = 0.01;
  t = 0:dt:sim_time;
  len_t = length(t);

  % --- Iteration Parameters ---
  num_iterations = 1000000;

  % --- Data Storage ---
  results = zeros(num_iterations, 5); % [Iteration, Initial_n, Overshoot_Collatz, SteadyStateError_Collatz, Overshoot_Std]

  % --- Main Loop ---
  for iteration = 1:num_iterations
    initial_n = iteration; % Start Collatz with the iteration number

    % --- Initialize Variables (reset for each iteration) ---
    y_std = zeros(len_t, 1);
    e_std = zeros(len_t, 1);
    u_std = zeros(len_t, 1);
    y_collatz = zeros(len_t, 1);
    e_collatz = zeros(len_t, 1);
    u_collatz = zeros(len_t, 1);
    Kp_collatz = zeros(len_t, 1);
    n = zeros(len_t, 1);
    n(1) = initial_n;

    % --- System Parameters (unchanged) ---
    zeta = 0.5;
    beta = 0.1;

    % --- Standard PID Simulation (unchanged, but in the loop) ---
    for i = 1:len_t - 1
      setpoint = 1;
      e_std(i) = setpoint - y_std(i);
      integral_term_std = sum(e_std(1:i)) * dt;
      derivative_term_std = (e_std(i) - e_std(max(1, i - 1))) / dt;
      u_std(i) = Kp_std * e_std(i) + Ki_std * integral_term_std + Kd_std * derivative_term_std;

      if i > 2
        dy_dt_std = (y_std(i) - y_std(i - 1)) / dt;
        dy2_dt2_std = (dy_dt_std - (y_std(i - 1) - y_std(i - 2)) / dt) / dt;
        next_dy_dt_std = dy_dt_std + dt * (u_std(i) - 2 * zeta * dy_dt_std - y_std(i) - beta * y_std(i)^3);
        y_std(i + 1) = y_std(i) + dt * next_dy_dt_std;
      elseif i == 2
        dy_dt_std = (y_std(i) - y_std(i - 1)) / dt;
        next_dy_dt_std = dy_dt_std + dt * (u_std(i) - 2 * zeta * dy_dt_std - y_std(i) - beta * y_std(i)^3);
        y_std(i + 1) = y_std(i) + dt * next_dy_dt_std;
      else
        y_std(2) = y_std(1) + dt * (u_std(1) - 2 * zeta * 0 - y_std(1) - beta * y_std(1)^3);
      end
    end

    % --- Collatz PID Simulation (unchanged, but in the loop) ---
    for i = 1:len_t - 1
      setpoint = 1;
      e_collatz(i) = setpoint - y_collatz(i);

      if mod(n(i), 2) == 0
        n(i + 1) = n(i) / 2;
      else
        n(i + 1) = 3 * n(i) + 1;
      end

      Kp_collatz(i) = Kp0_collatz + alpha_collatz * n(i);
      integral_term_collatz = sum(e_collatz(1:i)) * dt;
      derivative_term_collatz = (e_collatz(i) - e_collatz(max(1, i - 1))) / dt;
      u_collatz(i) = Kp_collatz(i) * e_collatz(i) + Ki_collatz * integral_term_collatz + Kd_collatz * derivative_term_collatz;

      if i > 2
        dy_dt_collatz = (y_collatz(i) - y_collatz(i - 1)) / dt;
        dy2_dt2_collatz = (dy_dt_collatz - (y_collatz(i - 1) - y_collatz(i - 2)) / dt) / dt;
        next_dy_dt_collatz = dy_dt_collatz + dt * (u_collatz(i) - 2 * zeta * dy_dt_collatz - y_collatz(i) - beta * y_collatz(i)^3);
        y_collatz(i + 1) = y_collatz(i) + dt * next_dy_dt_collatz;
      elseif i == 2
        dy_dt_collatz = (y_collatz(i) - y_collatz(i - 1)) / dt;
        next_dy_dt_collatz = dy_dt_collatz + dt * (u_collatz(i) - 2 * zeta * dy_dt_collatz - y_collatz(i) - beta * y_collatz(i)^3);
        y_collatz(i + 1) = y_collatz(i) + dt * next_dy_dt_collatz;
      else
        y_collatz(2) = y_collatz(1) + dt * (u_collatz(1) - 2 * zeta * 0 - y_collatz(1) - beta * y_collatz(1)^3);
      end
    end

    % --- Calculate Performance Metrics ---
    overshoot_std = (max(y_std) - 1) / 1 * 100;
    overshoot_collatz = (max(y_collatz) - 1) / 1 * 100;
    steady_state_error_collatz = abs(1 - y_collatz(end));

    % --- Store Results ---
    results(iteration, :) = [iteration, initial_n, overshoot_collatz, steady_state_error_collatz, overshoot_std];
  end

  % --- Save to CSV ---
  csv_filename = 'collatz_pid_results.csv';
  header = {'Iteration', 'Initial_n', 'Overshoot_Collatz', 'SteadyStateError_Collatz', 'Overshoot_Std'};
  fid = fopen(csv_filename, 'w');
  fprintf(fid, '%s,%s,%s,%s,%s\n', header{1}, header{2}, header{3}, header{4}, header{5}); % Write header

  for i = 1:num_iterations
    fprintf(fid, '%d,%f,%f,%f,%f\n', results(i, 1), results(i, 2), results(i, 3), results(i, 4), results(i, 5));
  end

  fclose(fid);

  disp(['Results saved to: ' csv_filename]);

  % --- Optional: Plotting (Consider if you want to plot aggregate results) ---
  % figure;
  % plot(results(:,1), results(:,3), 'r');
  % title('Collatz PID Overshoot vs. Iteration (Initial n)');
  % xlabel('Iteration (Initial n)');
  % ylabel('Overshoot (%)');
end

% --- Run the Simulation ---
compare_pid_controllers_iterative();