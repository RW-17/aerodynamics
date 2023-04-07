close all; clc; format long

foil = '5312';
P = 160;
alpha = 0;

[nodes, x, z_c] = nodes_ini(foil, P+1, 3);
[Cl, Cp, pos_c, Cmba, cli] = coefs(nodes, alpha);

% Airfoil plot with points, panels and control points
figure
plot(nodes(:,1), nodes(:,2), '-', 'Marker', '.')
hold on
plot(pos_c(:,1), pos_c(:,2), '.')
hold on
plot(x, z_c, '-')
hold off
axis equal

Cp_plot = zeros(P+1, 1); 
Cp_plot(1:P) = Cp; 
Cp_plot(P+1) = Cp(1);

figure (2)
plot(x, Cp_plot, '-', 'LineWidth', 0.8)
set(gca, 'YDir','reverse')

%% Function that creates the airfoil
function [nodes, x, z_c] = nodes_ini(foil, N, opt)

f = str2double(foil(1))*0.01;
x_fmax = str2double(foil(2))*0.1;
t_max = str2double(foil(3:4))*0.01;

x = zeros(N, 1);
z_c = zeros(N, 1);
nodes = zeros(N, 2);
P = N - 1;

% Equispaced distribution
if opt == 1
    for i = 0:P
        x(i+1) = abs(P - 2*i)/P;
    end
end

% More points at both edges 
if opt == 2
    for i = 0:P
        x(i+1) = 0.5*(1 - cos(2*i*pi/P + pi));
    end
end

% More points at the leading edge
if opt == 3
    for i = 0:P
        x(i+1) = (1 - cos(i*pi/P - pi/2));
    end
end

for i = 1:N
   n = [0.0, 0.0];
   if x(i) < x_fmax
       z_c(i) = f*(2*x_fmax*x(i) - x(i)^2)/(x_fmax^2);
       n = [-(f*(2*x_fmax-2*x(i)))/x_fmax^2, 1];
   end
   if x(i) >= x_fmax
       z_c(i) = f*((1-2*x_fmax) + 2*x_fmax*x(i) - x(i)^2)/(1-x_fmax)^2;
       n = [-(f*(2*x_fmax-2*x(i)))/(1-x_fmax)^2; 1];
   end

   n_hat = n/norm(n);

   z_esp = 5*t_max*(0.2969*sqrt(x(i)) - 0.1260*x(i)...
           - 0.3516*x(i)^2 + 0.2834*x(i)^3 - 0.1015*x(i)^4);
   if i <= N/2
       z_esp = -z_esp;       
   end
   nodes(i, :) = [x(i) + z_esp*n_hat(1); z_c(i) + z_esp*n_hat(2)];
end
end

%% Function that determines the lift coefficient
function [Cl, Cp, pos_c, Cmba, cli] = coefs(nodes, alpha)

alpha = deg2rad(alpha);
P = size(nodes, 1) - 1;

lp = zeros(P, 1);
pos_c = zeros(P, 2);
theta = zeros(P, 1);

[xp, zp, beta, a3, a4, a5, a6, a7, a8, a9, a10] = deal(zeros(P));
[B, C] = deal(zeros(P, P+1));
A = zeros(P+1);
b = zeros(P+1, 1);
[cli, vtg, Cp, delta, d, cmbai] = deal(zeros(P, 1));

for i = 1:P
    lp(i) = sqrt((nodes(i+1, 1) - nodes(i, 1))^2 + ...
        (nodes(i+1, 2) - nodes(i, 2))^2);
    pos_c(i,1) = 0.5*(nodes(i+1, 1) + nodes(i, 1));
    pos_c(i,2) = 0.5*(nodes(i+1, 2) + nodes(i, 2));
    theta(i) = atan2(nodes(i+1, 2) - nodes(i, 2), ...
                     nodes(i+1, 1) - nodes(i, 1));
end

for i = 1:P 
    for j = 1:P
        xp(i,j) = (pos_c(i, 1) - nodes(j, 1))*cos(theta(j))...
                + (pos_c(i, 2) - nodes(j, 2))*sin(theta(j));

        zp(i,j) = (pos_c(i, 2) - nodes(j, 2))*cos(theta(j))...
                - (pos_c(i, 1) - nodes(j, 1))*sin(theta(j));
        
        beta(i,j) = atan2(zp(i,j), xp(i,j) - lp(j)) ...
                                    - atan2(zp(i,j), xp(i,j));
        beta(i, i) = pi;
    end
end

for i = 1:P
    for j = 1:P
        r = sqrt((pos_c(i, 1) - nodes(j, 1))^2 ...
                 + (pos_c(i, 2) - nodes(j, 2))^2);

        rp1 = sqrt((pos_c(i, 1) - nodes(j+1, 1))^2 ...
                 + (pos_c(i, 2) - nodes(j+1, 2))^2);

        a6(i,j) = 1/(2*pi)*(log(rp1/r)...
        - 1/lp(j)*(log(rp1/r)*xp(i,j) + lp(j) - beta(i,j)*zp(i,j)));

        a4(i,j) = 1/(2*pi*lp(j))*(log(rp1/r)*xp(i,j) + lp(j) ...
                                                     - beta(i,j)*zp(i,j));
        a5(i,j) = 1/(2*pi)*(beta(i,j) ...
            - 1/lp(j)*(log(rp1/r)*zp(i,j) + beta(i,j)*xp(i,j)));

        a3(i,j) = 1/(2*pi*lp(j))*(log(rp1/r)*zp(i,j) + beta(i,j)*xp(i,j));

        a7(i,j) = a5(i,j)*cos(theta(j)) - a6(i,j)*sin(theta(j));
        a8(i,j) = a3(i,j)*cos(theta(j)) - a4(i,j)*sin(theta(j));

        a9(i,j) = a5(i,j)*sin(theta(j)) + a6(i,j)*cos(theta(j));
        a10(i,j) = a3(i,j)*sin(theta(j)) + a4(i,j)*cos(theta(j));
    end
end

for i = 1:P
    for j = 1:P
        if j == 1
            B(:,1) = a7(:,1);
            C(:,1) = a9(:,1);
        else
            B(i,j) = a8(i, j-1) + a7(i,j);
            C(i,j) = a10(i, j-1) + a9(i,j);
        end
        B(:, P+1) = a8(:, P);
        C(i, P+1) = a10(i, P);
        A(i,j) = - B(i,j)*sin(theta(i)) + C(i,j)*cos(theta(i));
    end
    A(i,P+1) = - B(i,P+1)*sin(theta(i)) + C(i,P+1)*cos(theta(i));
    b(i) = -sin(alpha - theta(i));
end

A(P+1, 1) = 1; A(P+1, P+1) = 1;
b(P+1) = 0;

gamma = A\b;
Cl = 0;
for i = 1:P
    cli(i) = gamma(i) + gamma(i+1);
    Cl = Cl + cli(i)*lp(i);
end

Cmba = 0;
for i = 1:P
    suma = 0;
    d(i) = sqrt(pos_c(i,1)^2 + pos_c(i,2)^2);
    delta(i) = atan2(pos_c(i,2), pos_c(i,1));
    for j = 1:P+1
        suma = suma + (B(i,j)*cos(theta(i)) ...
                         + C(i,j)*sin(theta(i)))*gamma(j);
    end
    vtg(i) = cos(theta(i) - alpha) + suma;
    Cp(i) = 1 - vtg(i).^2;
    cmbai(i) = cli(i)*d(i)*cos(delta(i) - alpha);
    Cmba = Cmba - cmbai(i)*lp(i);
end
end
