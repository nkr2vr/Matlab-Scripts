function [ dV_split, m_prop, m_i, m_f ] = deltaVsplit( mpay, deltaV, Isp, f_in )
%deltaVsplit Calculates optimal distribution of Delta V to reduce
%propellant mass for a 2-stage rocket.
%   Given the payload mass, desired deltaV, specific impulses for each
%   stage (can be input as a vector), and inert mass fraction, calculates
%   the required deltaV for each stage by minimizing the propellant mass.
%   Outputs the dV split as a vector, the minimized propellant mass, the
%   minimized initial mass before burns, and the final (empty) mass.

s = 0.3:0.01:0.7;
dV1 = s*deltaV;
stage1 = (exp(dV1./(9.8*Isp(1))) - 1) * (1-1/9) ./ (1-1/9*exp(dV1./(9.8*Isp(1))));
dV2 = deltaV - dV1;
stage2 = (exp(dV2./(9.8*Isp(2))) - 1) * (1-1/9) ./ (1-1/9*exp(dV2./(9.8*Isp(2))));

propmass = mpay .* stage1 .* stage2
[m_prop, dV_split(1)] = min(propmass)
dV_split(2) = deltaV - dV_split(1)
m_i = mpay + m_prop + f_in*(m_prop)/(1-f_in)
m_f = m_i - m_prop
end

