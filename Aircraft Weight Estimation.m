% Noah Rohrlich
% AAE251 HW04 Question 2 
% Due Tuesday February 24 2015
%
% Aircraft Sizing algorithm: This script will guide the user to enter the
% data needed to calculate the aircraft's estimated take-off weight.
%
% SAMPLE DATA:
%
% a) For range of 3500 naut mi and payload of 20000 lb and assuming
% anti-submarine warfare example from Raymer ch. 3: OUTPUT:
%   Will you be using the (much better) metric system? (y = 1, n = 0)0
%   Enter desired range (naut miles): 3500 Enter payload (lbm): 20000 
%   Enter total crew weight (lbm): 800
%   1. Turbojet 
%   2. Low-bypass turbofan 
%   3. High-bypass turbofan 
%   4. Piston prop 
%   5. Turboprop 
%   Please select engine type: 3 
%   Enter maximum Lift-over-Drag ratio (Raymer Fig. 3.6 may be of use): 16 
%   How many times will the craft loiter (durations will be entered next, 
%   one at a time)? 2 
%   Enter duration of one loiter (hrs): 3 
%   Enter duration of one loiter (hrs): 1/3 
%   1. Unpowered sailplane 
%   2. Powered sailplane 
%   3. Metal/Wood homebuilt 
%   4. Composite homebuilt 
%   5. Single engine general aviation 
%   6. Twin engine general aviation 
%   7. Agricultural aircraft
%   8. Twin turboprop 
%   9. Flying boat 
%   10. Jet trainer 
%   11. Jet fighter 
%   12. Cargo/Bomber 
%   13. Jet transport 
%   Which type of plane? 12
%
%   Based on the given design requirements, aircraft take-off weight has
%   been calculated to be 146775 lbm
%-----------------------------------------------------------------------
% b) For range of 2500 naut mi and payload of 25000 lb and assuming
% anti-submarine warfare example from Raymer ch. 3: OUTPUT:
%   Will you be using the (much better) metric system? (y = 1, n = 0)0
%   Enter desired range (naut miles): 2500
%   Enter payload (lbm): 25000
%   Enter total crew weight (lbm): 800
%   1. Turbojet
%   2. Low-bypass turbofan
%   3. High-bypass turbofan
%   4. Piston prop
%   5. Turboprop
%   Please select engine type: 3
%   Enter maximum Lift-over-Drag ratio (Raymer Fig. 3.6 may be of use): 16
%   How many times will the craft loiter (durations will be entered next, 
%   one at a time)? 2
%   Enter duration of one loiter (hrs): 3
%   Enter duration of one loiter (hrs): 1/3
%   1. Unpowered sailplane
%   2. Powered sailplane
%   3. Metal/Wood homebuilt
%   4. Composite homebuilt
%   5. Single engine general aviation
%   6. Twin engine general aviation
%   7. Agricultural aircraft
%   8. Twin turboprop
%   9. Flying boat
%   10. Jet trainer
%   11. Jet fighter
%   12. Cargo/Bomber
%   13. Jet transport
%   Which type of plane? 12
%   
%   Based on the given design requirements, aircraft take-off weight has
%   been calculated to be 242860 lbm

%--COLLECT DATA
% Ask user for all important design requirements of craft
metric = -1;
while metric < 0
    metric = input('\nWill you be using the (much better) metric system? (y = 1, n = 0)');
    if metric ~= 0 && metric ~= 1
        fprintf('Please enter 1 or 0')
        metric = -1;
    end
end

% This will be used when printing the answer at the end of the script
if metric
    w_unit = 'kg';
else
    w_unit = 'lbm';
end

% Get User input for design requirements
if metric
    range = input('Enter desired range (km): ');
    range = range * 1000; %Convert to meters
    payload = input('Enter payload (kg): ' );
    crew = input('Enter total crew weight (kg): ');
    
else
    range = input('Enter desired range (naut miles): ');
    range = range * 6076.12; %Convert to feet
    payload = input('Enter payload (lbm): ');
    crew = input('Enter total crew weight (lbm): ');
end
engine = input('1. Turbojet\n2. Low-bypass turbofan\n3. High-bypass turbofan\n4. Piston prop\n5. Turboprop\nPlease select engine type: ');
L_Dmax = input('Enter maximum Lift-over-Drag ratio (Raymer Fig. 3.6 may be of use): ');
numloit = input('How many times will the craft loiter (durations will be entered next, one at a time)? ');
loits = zeros(1,numloit);
for x = 1:numloit
    loits(x) = input('Enter duration of one loiter (hrs): ');
    % Convert to seconds from hours
    loits(x) = loits(x) * 3600;
end


% Assign specific fuel consumption, L/D estimate from tables
if engine == 1 %Turbojet
    C_cruise = 0.9;
    C_loiter = 0.8;
    L_Dcruise = 0.866*L_Dmax;
    L_Dloiter = L_Dmax;
else if engine == 2 %Low-BP Turbofan
        C_cruise = 0.8;
        C_loiter = 0.7;
        L_Dcruise = 0.866*L_Dmax;
        L_Dloiter = L_Dmax;
    else if engine == 3 %High-BP Turbofan
            C_cruise = 0.5;
            C_loiter = 0.4;
            L_Dcruise = 0.866*L_Dmax;
            L_Dloiter = L_Dmax;
        else if engine == 4 %Piston Prop
                C_cruise = 0.4;
                C_loiter = 0.5;
                L_Dcruise = L_Dmax;
                L_Dloiter = 0.866*L_Dmax;
            else if engine == 5 %Turboprop
                    C_cruise = 0.5;
                    C_loiter = 0.6;
                    L_Dcruise = L_Dmax;
                    L_Dloiter = 0.866*L_Dmax;
                end
            end
        end
    end
end
% C is in units of hr^-1, must be converted into units of s^-1
C_cruise = C_cruise / 3600;
C_loiter = C_loiter / 3600;


% Get constants for empty weight fraction equation
planetype = input('1. Unpowered sailplane\n2. Powered sailplane\n3. Metal/Wood homebuilt\n4. Composite homebuilt\n5. Single engine general aviation\n6. Twin engine general aviation\n7. Agricultural aircraft\n8. Twin turboprop\n9. Flying boat\n10. Jet trainer\n11. Jet fighter\n12. Cargo/Bomber\n13. Jet transport\nWhich type of plane? ');
if planetype == 1 %Unpowered sailplane
    A = 0.86;
    if metric
        A = 0.83;
    end
    c = -0.05;
else if planetype == 2 %Powered sailplane
        A = 0.91;
        if metric
            A = 0.88;
        end
        c = -0.05;
    else if planetype == 3 %Metal/wood homebuilt
            A = 1.19;
            if metric
                A = 1.11;
            end
            c = -0.09;
        else if planetype == 4 %Composite homebuilt
                A = 1.17;
                if metric
                    A = 1.07;
                end
                c = -0.09;
            else if planetype == 5 %Single engine general aviation
                    A = 2.36;
                    if metric
                        A = 2.05;
                    end
                    c = -0.18;
                else if planetype == 6  %Twin engine general aviation
                        A = 1.51;
                        if metric
                            A = 1.4;
                        end
                        c = -0.10;
                    else if planetype == 7  %Agricultural aircraft
                            A = 0.74;
                            if metric
                                A = 0.72;
                            end
                            c = -0.03;
                        else if planetype == 8  %Twin turboprop
                                A = 0.96;
                                if metric
                                    A = 0.92;
                                end
                                c = -0.05;
                            else if planetype == 9  %Flying boat
                                    A = 1.09;
                                    if metric
                                        A = 1.05;
                                    end
                                    c = -0.05;
                                else if planetype == 10 %Jet trainer
                                        A = 1.59;
                                        if metric
                                            A = 1.47;
                                        end
                                        c = -0.10;
                                    else if planetype == 11 %Jet fighter
                                            A = 2.34;
                                            if metric
                                                A = 2.11;
                                            end
                                            c = -0.13;
                                        else if planetype == 12 %Cargo/bomber
                                                A = 0.93;
                                                if metric
                                                    A = 0.88;
                                                end
                                                c = -0.07;
                                            else if planetype == 13 %Jet transport
                                                    A = 1.02;
                                                    if metric
                                                        A = 0.97;
                                                    end
                                                    c = -0.06;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%--SIZING ALGORITHM

% Initial weight
w_0 = 0;
w_g = 50000;

% Continue loop until weight is found, with small margin for error
%for debug = 1:2
while abs(1 - w_g/w_0) >= 0.0001 
    % 0.01% error addresses converging values that won't equal
    w_0 = w_g;
    
    % Calculate fuel weight fraction by finding mission fuel fraction From
    % historical data:     Warmup and take-off fraction:   0.970
    %                           Climb fraction:                 0.985
    %                           Landing fraction:               0.995
    
    % Warmup and takeoff
    misnfrac = 0.97;
    
    % Climb
    misnfrac = misnfrac * 0.985;
    
    % Cruises (assuming typical cruise altitude of 30000ft [~9144 m] and
    % velocity of 0.6 Mach)
    v = 0.6 * 994.8; %speed of sound at 30000ft, in ft/s
    if metric
        v =  0.6 * 303.215; %speed of sound at 9144m, in m/s
    end
    cruisefrac = exp(-range * C_cruise / (v * L_Dcruise));
    misnfrac = misnfrac * (cruisefrac^numloit);
    
    % Loiters (taking into account multiple loiters of different durations)
    for x = 1:numloit
        % loits(x) is the duration of the loiter
        loitfrac = exp(-loits(x) * C_loiter / L_Dloiter);
        misnfrac = misnfrac * loitfrac;
    end
    
    % Landing (w7/w0)
    misnfrac = misnfrac * 0.995;
    % Account for 6% trapped fuel
    fuelfrac = 1.06 * (1 - misnfrac);
    
    % Calculate empty weight fraction
    emptyfrac = A * w_0^c;
    
    % Recalculate gross weight
    if (fuelfrac + emptyfrac) < 1
        w_g = (crew + payload) / (1 - fuelfrac - emptyfrac);
    else
        % If denominator is negative, the following prevents a wrench from
        % getting thrown into the algorithm.
        w_g = (crew + payload) / (1 - 1.06*(misnfrac) - emptyfrac);
    end
end

% Round answer to whole number
w_g = round(w_g);

%--OUTPUT
fprintf('\nBased on the given design requirements, aircraft take-off weight has been calculated to be %d %s\n',w_g, w_unit)

