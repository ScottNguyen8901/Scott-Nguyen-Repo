clear
clc

% parse_ephemeris_data('C:\Users\scott\Documents\Folder\GTOC-1\body_info',...
%     'C:\Users\scott\Documents\Folder\GTOC-1\body_state')

% Constants
AU_to_km = 1.4959787066e8;
mu.sun = 1.32712428e11;

% Orbital elements for TW229 (Asteroid)
TW229_koe.a = 2.5897261 * AU_to_km;
TW229_koe.e = 0.2734625;
TW229_koe.i = deg2rad(6.40734);
TW229_koe.w = deg2rad(264.78691);
TW229_koe.W = deg2rad(128.34711);
TW229_koe.M = deg2rad(320.47955);

% Solve Kepler's equation to get the true anomaly
TW229_koe.F = kepler_solver(TW229_koe.M, TW229_koe.e);
TW229_koe.f = anom_to_true_anom(TW229_koe.F, TW229_koe.e);

% Initial state in orbital elements, converting to position and velocity
[r, v] = koe_to_rv(TW229_koe, mu.sun);

% Define time vector
time_vec = 1:(2462502.500000000 - 2455197.500000000);

% Function to calculate the acceleration due to gravity
orbital_equation = @(t, state) orbital_dynamics(t, state, mu.sun);

% Initial state vector (r and v) at time = 0 (beginning of the integration)
initial_state = [r; v];

% Time span (in days)
t_span = time_vec;

% Integrating the orbital equations using ode45
[t, solution] = ode45(orbital_equation, t_span, initial_state);

% Extract the positions of the asteroid from the solution
r_astroid = solution(:, 1:3);  % Asteroid position (X, Y, Z)

% Specify the folder containing the planet ephemeris files
folder_path = 'C:\Users\scott\Documents\Folder\GTOC-1\body_state'; % Replace with your folder path

% Create a new figure for plotting the orbits
figure;
hold on;

% Plot the asteroid orbit
plot3(r_astroid(:, 1), r_astroid(:, 2), r_astroid(:, 3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Asteroid TW229');

% Get a list of all .txt files in the folder
file_list = dir(fullfile(folder_path, '*.txt'));

% Loop through each file and extract data for planets
for k = 1:length(file_list)
    % Get the file name
    file_name = file_list(k).name;
    
    % Extract the planet name from the filename (e.g., mercury_ephemeris.txt -> mercury)
    planet_name = erase(file_name, '.txt'); % Removes the '.txt' extension
    
    % Construct the full file path
    file_path = fullfile(folder_path, file_name);
    
    % Read the data from the file (tab-delimited data)
    data = readtable(file_path, 'Delimiter', '\t');
    
    % Extract position (X, Y, Z) of the planet over time
    Position = [data{:, 'X'}, data{:, 'Y'}, data{:, 'Z'}];
    
    % Plot the planet's orbit separately
    plot3(Position(:, 1), Position(:, 2), Position(:, 3), 'LineWidth', 1.5, 'DisplayName', planet_name);
end

% Add a legend to distinguish between asteroid and planets
legend show;