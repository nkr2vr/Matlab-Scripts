% Noah Rohrlich nkr2vr@virginia.edu
% Due date: 03/18/2016
% Transformer Design Script
% Purpose:
% Given basic theoretical requirements of transformer, output physical
% characteristics of transformer.
%
% Inputs:
% - None. Requirements are hard-coded.
%
% Final Outputs:
% - Number_of_laminations
% - Primary_wire_gauge
% - Primary_turns
% - Secondary_wire_gauge
% - Secondary_turns
% - R_Cu_final:           Total resistance of copper wire in transformer
% - jX_l_final:           Total leakage reactance in transformer
% - R_core_final:         Core resistance of transformer
% - jX_m_final:           Core reactance (magnetization) of transformer
% - Copper_losses
% - Core_losses
% - Vreg_final
% - Effcy_final
% - J_final
% - B_final
% - Wh_used:        Height of window used for windings
% - Ww_used:        Width of window used for windings
%% 
% Clear the screen and workspace
clear
clc
% DATA ------------------------------------------------------------------

% Requirements
VA_rated = 60;
Vp_rated = 120;
Vs_rated = 24;
Vreg_max = 0.04; %Vreg_act should be less than this
Effcy_min = 0.8; %Effcy_act should be greater than this

%% 
% Electromagnetic characteristics
f = 60; % Operating frequency in Hertz
% Convert m^2 to cm^2 for convenience during calculation
B_range = (1.3:0.01:1.6) / 100^2; % B values in Wb/cm^2
J_range = 25:0.1:300; % J values in A/cm^2
Effcy_plane = zeros(length(B_range), length(J_range));
Vreg_plane = zeros(length(B_range), length(J_range));

%% 
% Physical characteristics for EI-125H M6 lamination
t_lam = 0.036; % lamination thickness in cm
W_width = 1.5875; % Window height in cm
W_height = 4.7625; % Window width in cm
Aw_act = W_height * W_width; % Actual window area in cm^2
Core_width = 3.175; % Core height in cm
Core_height = ...
    W_height; % Creating a separate, identical variable for legibility
Weight_per_1k = 15.987; % kg per 1000 laminations
Stack_factor = 0.9; % Stacking factor for laminations
K_u = 0.3; % Wire utilization factor
K_f = 4.443; % Frequency excitation factor at 60Hz -> (2*PI)/Sqrt(2)

%% 
% Characteristics of windings
t_winding_tube = 0.2; % Thickness of winding tube wall in cm
window_clearance = ...
    0.2; % Height clearance inside window at top and bottom in cm
Wh_used = W_height - 2 * window_clearance;
Ww_used = ...
    W_width - window_clearance; % Wall only along core side of window

%%
% AWG table
gauge = ...
    csvread('AWG_Table.csv', 1, 0, [1,0,15,0]); % AWG rating
bare_area = ...
    csvread('AWG_Table.csv', 1, 1, [1,1,15,1]) / 10^3; % in cm^2
resistance = ...
    csvread('AWG_Table.csv', 1, 3, [1,3,15,3]) / 10^6; % in Ohm/cm @ 20degC
wire_area = ...
    csvread('AWG_Table.csv', 1, 4, [1,4,15,4]) / 10^3; % in cm^2
wire_diameter = ...
    csvread('AWG_Table.csv', 1, 6, [1,6,15,6]); % in cm
turns_per_length = ...
    csvread('AWG_Table.csv', 1, 8, [1,8,15,8]); % in turns/cm
turns_per_area = ...
    csvread('AWG_Table.csv', 1, 10, [1,10,15,10]); % in turns/cm^2
weight = ...
    csvread('AWG_Table.csv', 1, 12, [1,12,15,12]) / 10^3; % kg/cm

%% 
% Initialize outputs
Number_of_laminations = 0;
Primary_wire_gauge = 0;
Primary_turns = 0;
Secondary_wire_gauge = 0;
Secondary_turns = 0;
R_Cu_final = 0;
jX_l_final = 0;
R_core_final = 0;
jX_m_final = 0;
Copper_losses = 0;
Core_losses = 0;
Vreg_final = 0;
Effcy_final = 0;
J_final = 0;
B_final = 0;

%%
% CALCULATIONS -----------------------------------------------------------
% Required apparent power
P_t = VA_rated * (1 + 1/Effcy_min);
% Input and output current
I_in = VA_rated / (Vp_rated * Effcy_min);
I_out = VA_rated / Vs_rated;


%%
% Iterate over all combinations of B and J
for r = 1:length(B_range)
    for c = 1:length(J_range)
        B = B_range(r);
        J = J_range(c);
        % Core dimensions
        N_lam = P_t / (K_f * K_u * f * B * J * Aw_act * t_lam * Core_width);
        N_lam = ceil(N_lam); % Round up to lower flux density
        Core_depth = N_lam * t_lam / Stack_factor;
        % Cross-sectional area of core in cm^2
        A_core = Core_depth * Core_width;
        Core_weight = Weight_per_1k / 1000 * N_lam; % weight in kg
        
        % Must reconvert B to Wb/m^2 for use with this function
        W_per_kg = 0.0618 * exp(3.1555 * B * 100^2);
        VAR_per_kg = 0.0288 * exp(5.1482 * B * 100^2);
        
        % Core losses
        P_core = Core_weight * W_per_kg;
        Q_core = Core_weight * VAR_per_kg;
        
        R_core = Vp_rated^2 / P_core;
        jX_m = Vp_rated^2 / Q_core;
        %%
        % Windings
        %Primary winding
        % Number of turns on primary
        N_p = ceil(Vp_rated / (K_f * f * B * A_core));
        % Copper area required in primary winding at J value
        Awire_p = I_in / J;
        p_index = 2; % Index for searching primary wire in AWG table
        % Index must start at 2 so subtraction can occur later
        while (p_index <= 15) && (bare_area(p_index) > Awire_p)
            % Wire gauge unacceptable if bare area is insufficient for Awire
            % Stops after smallest possible gauge checked
            p_index = p_index + 1;
        end
        p_index = p_index - 1; % Use last acceptable gauge
        WG_p = gauge(p_index);
        % Round down because needs to fit
        p_trns_per_lyr = floor(Wh_used * turns_per_length(p_index));
        % Number of primary layers, rounded up because fraction of layer is
        % still a layer
        p_lyrs = ceil(N_p/p_trns_per_lyr);
        % Physical depth of coil in cm
        p_depth = p_lyrs * wire_diameter(p_index);
        
        %Secondary winding
        % Number of turns on secondary
        N_s = ceil(N_p * Vs_rated / Vp_rated * (1 + Vreg_max));
        % Copper area required in secondary winding at J value
        Awire_s = I_out / J;
        s_index = 2; % Index for searching secondary wire in AWG table
        while (bare_area(s_index) > Awire_s) || (s_index == 15)
            s_index = s_index + 1;
        end
        s_index = s_index - 1;
        WG_s = gauge(s_index);
        s_trns_per_lyr = floor(Wh_used * turns_per_length(s_index));
        s_lyrs = ceil(N_s/s_trns_per_lyr);
        s_depth = s_lyrs * wire_diameter(s_index);
        
        % Sanity check: ensuring the coils fit the window
        total_coil_depth = p_depth + s_depth;
        
        if(total_coil_depth > Ww_used)
            break % End this iteration if the coils won't fit
        end
        %%
        % Copper Losses
        % EXTRA CREDIT!!!!!
        % First set: assume primary is on inside, secondary outside
        %Mean Length per Turn
        MLT_p = 2 * (Core_depth + 2 * t_winding_tube)...
            + 2 * (Core_width + 2 * t_winding_tube) + pi * p_depth;
        MLT_s = 2 * (Core_depth + 2 * t_winding_tube) + 2 * ...
            (Core_width + 2 * t_winding_tube) + pi * (2 * p_depth + s_depth);
        R_Cu_p = MLT_p * N_p * resistance(p_index);
        R_Cu_s = MLT_s * N_s * resistance(s_index);
        %Copper Losses with primary on inside
        P_Cu = I_in^2 * R_Cu_p + I_out^2 * R_Cu_s;
        
        % Second set: assume primary is on outside, are the losses smaller?
        MLT_s2 = 2 * (Core_depth + 2 * t_winding_tube)...
            + 2 * (Core_width + 2 * t_winding_tube) + pi * s_depth;
        MLT_p2 = 2 * (Core_depth + 2 * t_winding_tube) + 2 * ...
            (Core_width + 2 * t_winding_tube) + pi * (2 * s_depth + p_depth);
        R_Cu_p2 = MLT_p2 * N_p * resistance(p_index);
        R_Cu_s2 = MLT_s2 * N_s * resistance(s_index);
        %Copper losses with primary on outside
        P_Cu2 = I_in^2 * R_Cu_p2 + I_out^2 * R_Cu_s2;
        
        % If putting the primary on the outside is more efficient (smaller
        % loss), then use that configuration
        if P_Cu2 < P_Cu
            P_Cu = P_Cu2;
            R_Cu_p = R_Cu_p2;
            R_Cu_s = R_Cu_s2;
            MLT_p = MLT_p2;
            MLT_s = MLT_s2;
        end
        
        %Total copper line resistance
        R_Cu = R_Cu_p + R_Cu_s * N_p^2 / N_s^2;
        Effcy_act = VA_rated / (VA_rated + P_Cu + P_core);
        % Assume voltage drop across leakage reactance is negligible, so
        % Vreg can be approximated
        Vreg_act = P_Cu / VA_rated;
        %%
        % Leakage Reactance
        % Leakage induction as seen by primary, in Henrys
        % Insulating layer between coils ~0
        L_p = 4 * pi * MLT_p * N_p^2 / Wh_used * (p_depth + s_depth)...
            / 3 * 10^-9;
        jX_l = L_p * f * 2 * pi;
        
        Effcy_plane(r,c) = Effcy_act;
        Vreg_plane(r,c) = Vreg_act;
        
        % Store all values if current Efficiency and Vreg are good
        if(Effcy_act > Effcy_min && Vreg_act < Vreg_max)
            % Store all current data if this efficiency is better than the
            % last one
            if(Effcy_act > Effcy_final)
                Number_of_laminations = N_lam;
                Primary_wire_gauge = WG_p;
                Primary_turns = N_p;
                Secondary_wire_gauge = WG_s;
                Secondary_turns = N_s;
                R_Cu_final = R_Cu;
                jX_l_final = jX_l;
                R_core_final = R_core;
                jX_m_final = jX_m;
                Copper_losses = P_Cu;
                Core_losses = P_core;
                Vreg_final = Vreg_act;
                Effcy_final = Effcy_act;
                J_final = J;
                B_final = B;
            end
        end
    end
end

subplot(1,2,1);
s1 = surf(J_range * 100^2, B_range * 100^2, Effcy_plane * 100);
set(s1, 'LineStyle', 'none')
xlabel('J (A/m^2)')
ylabel('B (Wb/m^2)')
zlabel('Power Efficiency (%)')
subplot(1,2,2);
s2 = surf(J_range * 100^2, B_range * 100^2, Vreg_plane * 100);
set(s2, 'LineStyle', 'none')
xlabel('J (A/m^2)')
ylabel('B (Wb/m^2)')
zlabel('Voltage Regulation (%)')

fprintf('\n-----------------------------');
fprintf('\nOptimized Transformer Data:');
fprintf('\n-----------------------------');
fprintf('\nNumber of laminations: %d', Number_of_laminations);
fprintf('\nPrimary Coil Wire Gauge (AWG): %d', Primary_wire_gauge);
fprintf('\nNumber of Turns on Primary Coil: %d', Primary_turns);
fprintf('\nSecondary Coil Wire Gauge (AWG): %d', Secondary_wire_gauge);
fprintf('\nNumber of Turns on Secondary Coil: %d', Secondary_turns);
fprintf('\nCircuit Model: ');
fprintf('\n\tCopper Resistance R_Cu: %1.3f Ohms', R_Cu_final);
fprintf('\n\tLeakage Reactance jX_l: %1.3f Ohms', jX_l_final);
fprintf('\n\tCore Resistance R_core: %4.0f Ohms', R_core_final);
fprintf('\n\tCore Magnetization jX_m: %4.0f Ohms', jX_m_final);
fprintf('\nVoltage Regulation: %1.2f%%', Vreg_final * 100);
fprintf('\nReal Power Efficiency: %2.1f%%', Effcy_final * 100);
fprintf('\nCore Loss P_core: %3.2f W', P_core);
fprintf('\nCopper Loss P_Cu: %3.2f W', P_Cu);
fprintf('\nOptimal Current Density J: %3.0f A/cm^2', J_final);
fprintf('\nOptimal Magnetic Flux Density B: %1.2f T', B_final * 100^2);
fprintf('\nWindow height used: %1.2f cm', Wh_used);
fprintf('\nWindow width used: %1.2f cm\n\n', Ww_used);
%----------------------------------------------------END OF PROGRAM