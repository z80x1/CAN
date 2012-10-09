
atomsind

T=300;

N0 = 1
n = 1
R = 1
 

u_a = 10.^([-1:0.1:1]);

Evib = CC.k*T *u_a; % h*nu_a

H_vib = N0 * Evib ./ (exp(u_a)-1);

S_vib = n * R * (1./(u_a.*exp(u_a)-1)  - log(1-exp(-u_a)));

plot( u_a, H_vib - T*S_vib)