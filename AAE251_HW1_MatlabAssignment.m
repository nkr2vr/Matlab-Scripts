% Rohrlich, Noah
% HW1 Matlab Assignment due Tuesday January 27 2015
%
% Recreate the standard atmosphere as in appendix B of Anderson, and
% calculate the speed of sound at each altitude (equation for speed of 
% sound can be found in Anderson, p.176).
%
% Produce correctly labeled plots for:
% - altitude vs temperature
% - altitude vs pressue
% - altitude vs density
% - altitude vs the speed of sound
% 
% Altitudes are geometric and must go from Sea Level (0 ft, 0 m) up to 100,000 ft
% or 30480 m.
% 
% Data must be shown in SI and English Engineering units.

%--CONSTANTS
g = 9.81;                    % Acceleration due to gravity in m/s^2

R = 287;                     % Gas constant for dry air in J/(kg 8 K)

r_e = 6.3781e6;              % Equatorial radius of Earth in m

alt = 0:200:30500;           % Geometric altitudes in increments of 200 m

GAMMA = 1.4;                 % Isentropic expansion factor

%--CALCULATIONS

% Geopotential altitudes in m
gpt_alt = (r_e ./ (r_e + alt)) .* alt;

% Create vectors for pressure, density, temperature, and speed of sound
pres = zeros(length(gpt_alt), 1);

dens = zeros(length(gpt_alt), 1);

temp = zeros(length(gpt_alt), 1);

speed = zeros(length(gpt_alt), 1);

% Sea level pressure, density, and temperature
pres(1) = 1.01325e5;                        % N/m^2

dens(1) = 1.2250;                           % kg/m^3

temp(1) = 288.16;                           % K

speed(1) = sqrt(GAMMA * R * temp(1));       % m/s

% Fill in rest of arrays using formulas from Anderson 3.4 (see attached write-up)
for y = 2:length(gpt_alt)
    % Check temperature lapse rate
    h = gpt_alt(y);
    a = 0;
    if and(h >= 0, h < 11000)
        a = -6.5e-3;
    else if and(h >= 25000, h < 47000)
        a = 3e-3;
        end
    end
    
    % Get change in altitude
    del_h = gpt_alt(y) - gpt_alt(y-1);
    % Calculate temperature at h
    temp(y) = temp(y-1) + a*(del_h);
    
    % Calculate change rate for pressure and density
    rate = exp(-g * del_h / (R * temp(y)));
    
    % Calculate pressure at h
    pres(y) = pres(y-1) * rate;
    % Calculate density at h
    dens(y) = dens(y-1) * rate;
    % Calculate speed of sound at h
    speed(y) = sqrt(GAMMA * R * temp(y));
end

% Convert data into engineering units
alt_ee = alt * 3.28;        % Altitude in ft
temp_ee = temp * 1.8;       % Temperature in Rankine
pres_ee = pres * 0.0209;    % Pressure in lb/ft^2
dens_ee = dens * 0.00194;   % Density in slug/ft^3
speed_ee = speed * 3.28;    % Speed in ft/s

%--OUTPUTS
% Temperature vs Altitude in SI Units
subplot(4, 2, 1)
plot(temp, alt)
xlabel('Temperature (K)')
ylabel('Altitude (m)')
axis([200 300 0 35000])

% Temperature vs Altitude in Engineering Units
subplot(4, 2, 2)
plot(temp_ee, alt_ee)
xlabel('Temperature (degrees R)')
ylabel('Altitude (ft)')
axis([350 550 0 110000])

% Pressure vs Altitude in SI Units
subplot(4, 2, 3)
plot(pres, alt)
xlabel('Pressure (Pa)')
ylabel('Altitude (m)')

% Pressure vs Altitude in Engineering Units
subplot(4, 2, 4)
plot(pres_ee, alt_ee)
xlabel('Pressure (lbf/ft^2)')
ylabel('Altitude (ft)')
axis([0 2500 0 110000])

% Density vs Altitude in SI Units
subplot(4, 2, 5)
plot(dens, alt)
xlabel('Air Density (kg/m^3)')
ylabel('Altitude (m)')

% Density vs Altitude in Engineering Units
subplot(4, 2, 6)
plot(dens_ee, alt_ee)
xlabel('Air Density (slug/ft^3)')
ylabel('Altitude (ft)')
axis([0 2.5e-3 0 110000])

% Speed of sound vs Altitude in SI Units
subplot(4, 2, 7)
plot(speed, alt)
xlabel('Speed of Sound (m/s)')
ylabel('Altitude (m)')
axis([290 341 0 35000])

% Speed of sound vs Altitude in Engineering Units
subplot(4, 2, 8)
plot(speed_ee, alt_ee)
xlabel('Speed of Sound (ft/s)')
ylabel('Altitude (ft)')
axis([960 1150 0 110000])