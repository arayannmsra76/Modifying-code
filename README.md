# Modifying-code
%% Parameters
tol = 1e-7;
dim = 3;
if dim == 3
    Ih = [1; 1; 1; 0; 0; 0];
    I = [1; 1; 1; 0; 0; 0];
    enToPhys = [1; 1; 1; 0.5; 0.5; 0.5];
    ncomp = 6;
    Ce = zeros(6); % You can replace this with your linear elastic material data
else
    Ih = [1; 1; 0; 0];
    I = [1; 1; 0; 0];
    enToPhys = [1; 1; 0.5; 0];
    ncomp = 4;
    Ce = zeros(4); % You can replace this with your linear elastic material data
end

IxI = I * I';
Is = 0.5 * (diag(I) + eye(ncomp));
Id = Is - (1/3) * IxI;

matnr = 1; % You need to set this to an appropriate material number
sig_y0 = 100; % Assign a specific value for sig_y0

G = 10;
K =50;
r = 5;
s = 1;
Rinf = 12;
gamma = 9;

K2 = 2 * K;
G2 = 2 * G;
G3 = 3 * G;
G6 = 6 * G;

%% Get Strain, Damage, Hardening (Assuming these variables are defined)
eps = zeros(ncomp, 1);
eps_pl = zeros(ncomp, 1);
Dam0 = 0; % Initial damage
Int0 = 1 - Dam0;
R0 = 0; % Initial hardening internal variable

%% Trial State
eps_e_tr = eps - eps_pl;
eps_e_hyd_tr = sum(eps_e_tr .* Ih);
sig_hyd_tr = K * eps_e_hyd_tr;
eps_e_dev_tr = eps_e_tr - (1/3) * eps_e_hyd_tr * Ih;

eps_e_dev_tr = eps_e_dev_tr .* enToPhys;
temp = eps_e_dev_tr.^2;
J_2 = G2^2 * (0.5 * sum(temp .* Ih) + sum(temp));
q_tr = sqrt(3 * J_2);
fsig0 = sig_y0 + Rinf * (1 - exp(-gamma * R0));
Phi = q_tr - fsig0;

if Phi >= 0
    plasticMult = Int0 * Phi / (3 * G);
    R = R0 + plasticMult;
    norm_F = 1;
    sig_hyd_tr2 = sig_hyd_tr^2;
    iter = 0;
    maxiter = 100;

    while norm_F >= tol && iter <= maxiter
        if iter == maxiter
            disp('Fatal: return mapping reached maximum iterations!');
            break;
        end

        fsig = sig_y0 + 3300 * (1 - exp(-0.4 * R));
        f1 = (G3 / (q_tr - fsig));
        Int = f1 * plasticMult;
        Y = -(fsig^2) / G6 - (sig_hyd_tr2) / K2;
        f2 = -Y / r;
        F = Int - Int0 + (-Y / r)^s / f1;
        norm_F = abs(F);
        dfsig = 1320 * exp(-0.4 * R);
        dY = -(fsig * dfsig) / G3;
        f = f1 + f1 * plasticMult * dfsig / (q_tr - fsig) - (dfsig / G3) * f2^s - (s * dY / (f1 * r)) * f2^(s-1);
        plasticMult = plasticMult - F / f;
        R = R0 + plasticMult;
        iter = iter + 1;
    end

    fsig = sig_y0 + Rinf * (1 - exp(-gamma * R));
    dfsig = gamma * Rinf * exp(-gamma * R);
    f1 = (G3 / (q_tr - fsig));
    Int = f1 * plasticMult;
    dInt = (G3 + Int * dfsig) / (q_tr - fsig);
    Y = -(fsig^2) / G6 - (sig_hyd_tr2) / K2;
    dY = -(fsig * dfsig) / G3;
    f2 = -Y / r;
    f = f1 + f1 * plasticMult * dfsig / (q_tr - fsig) - (dfsig / G3) * f2^s - (s * dY / (f1 * 1)) * f2^(s-1);

    if (Int < 1e-20)
        disp('GP integrity too small!');
    end

    Dam = 1 - Int;
    q = Int * fsig;
    sig_hyd = Int * sig_hyd_tr;
    sig_dev = G2 * (q / q_tr) * eps_e_dev_tr;
    sig = sig_dev + sig_hyd * Ih;
    plCor = G3 * plasticMult / (Int * q_tr);
    eps_e = (1 - plCor) * eps_e_dev_tr * enToPhys + (1/3) * eps_e_hyd_tr * Ih;
    deps_pl = eps_e_tr - eps_e;
    eldat.sigma(1i, :, ig) = sig;
    eldat.eps_pl_n(1i, :, ig) = eps_pl + deps_pl;
    eldat.damage_n(1i, ig) = Dam;
    eldat.R_n(1i, ig) = R;
    eldat.q_n(1i, ig) = q;
    eldat.triax_n(1i, ig) = sig_hyd / q;
    eldat.sigy_n(1i, ig) = fsig;

    snorm = sqrt(sum(enToPhys .* (sig - Ce * eps).^2));
    f3 = q_tr - fsig;
    a1 = (1 / f) * (Int / 3 - (1 / G3) * f2^s);
    a2 = -s * sig_hyd_tr * f3 / (G3 * r * K * f) * f2^(s-1);
    a3 = a2 * dInt;
    a4 = a1 * dInt - Int / f3;
    a = G2 * Int * fsig / q_tr;
    b = G2 * (a1 * dfsig * Int + a4 * fsig - Int * fsig / q_tr);
    b = b / (snorm^2);
    c = K * sqrt(2/3) * (a2 * dfsig * Int + a3 * fsig);
    c = c / snorm;
    d = G2 * sqrt(3/2) * sig_hyd_tr * a4;
    d = d / snorm;
    e = K * (Int + a3 * sig_hyd_tr);
    D = a * Id + b * (sig_dev * sig_dev') + c * (sig_dev * I') + d * (I * sig_dev') + e * (IxI);

else
    D = Int0 * Ce;
    sig_hyd = Int0 * sig_hyd_tr;
    sig_dev = G2 * Int0 * eps_e_dev_tr;
    sig = sig_dev + sig_hyd * Ih;
    eldat.sigma(1i, :, ig) = sig;
    eldat.eps_pl_n(1i, :, ig) = eps_pl;
    eldat.damage_n(1i, ig) = Dam0;
    eldat.R_n(1i, ig) = R0;
    eldat.q_n(1i, ig) = Int0 * q_tr;
    eldat.triax_n(1i, ig) = sig_hyd / (Int0 * q_tr);
    eldat.sigy_n(1i, ig) = fsig0;
end



