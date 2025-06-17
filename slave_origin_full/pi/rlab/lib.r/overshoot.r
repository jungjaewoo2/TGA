// RLaB diary file: Overshoot_vs__DampingRatio.diary.r. Opened Thu Mar  8 22:26:22 2018
// For a second order system, the damping ratio (zeta) determines the the percentage of overshoot in the time
//  domain
overshoot = function(zeta)
{
global(pi)
return exp(-zeta*pi.*(1-zeta.^2).^-0.5);
};
