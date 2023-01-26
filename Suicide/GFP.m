ini = [998.14, 1, 1, 0, 0]; % 感觉DNA浓度越大越能让那个效率降低
tmesh = 0:1:420;
[t, y1] = ode45(@GFP37, tmesh, ini);
% figure;
% plot(t, y1(:, 5));
% xlabel("Time, min");
% ylabel("Fluorescence Intensity, a.u.");
% title('Fluorescence when 37');

tmp = y1(421, :);
[t, y2] = ode45(@GFP28, tmesh, tmp);
% figure;
% plot(y2(:, 5));
% xlabel("Time, min");
% ylabel("Fluorescence Intensity, a.u.");
% title('Fluorescence when changed to 28');

y = [y1; y2];
figure;
plot(y(:, 5));
xlabel("Time, min");
ylabel("Fluorescence Intensity, a.u.");
title('Fluorescence together');

function ydot = GFP37(t, y)
    D = 3.4; % 23~236
    phiT = 1; % T > 37 
    k = 0.49997147; % half life
    kts = 18.2;
    ktl = 16.1;
    kmat = 0.156;
    kcs = 1.1e-2;
    dmrna = 8.4e-2;
    dtlr = 4.5e-3;
    Ks = 4.2; % 0.8~4.2
    Kl = 30; % 14~30
    Ktlr = 6e-5;
%     dgfp = 0.00417;

    ydot = zeros(5, 1);
    % 函数命名 %
    % ydot(1) = d[DNA] / dt %
    % ydot(2) = d[TsR] / dt %
    % ydot(3) = d[TlR] / dt % 
    % ydot(4) = d[mRNA] / dt %
    % ydot(5) = d[GFP] / dt %
    % ydot(6) = d[GFP*] / dt %
    
%     ydot(2) = (-kcs*y(2)*y(1)/(Ks+y(1)))*(y(2) > 0 && y(2) > kcs*y(2)*y(1)/(Ks+y(1))); % 增加约束条件
%     ydot(3) = (-dtlr*y(3)/(Ktlr+y(3)))*(y(3) > 0 && y(3) > dtlr*y(3)/(Ktlr+y(3))); % 增加约束条件
%     ydot(2) = -kcs*y(2)*y(1)/(Ks+y(1));
%     ydot(3) = -dtlr*y(3)/(Ktlr+y(3));
    ydot(4) = kts*y(2)*y(1)/(Ks+y(1))*phiT-dmrna*y(4);
%     ydot(5) = ktl*y(3)*y(4)/(Kl+y(4))-kmat*y(5)-dgfp;
    ydot(5) = (ktl*y(3)*y(4)/(Kl+y(4))-kmat*y(5))*(y(5)+ktl*y(3)*y(4)/(Kl+y(4)) > kmat*y(5));
%     ydot(6) = kmat*y(5)-dgfp;

end

function ydot = GFP28(t, y)
    phiT = 0.132; % T <= 37 
    kts = 18.2;
    ktl = 16.1;
    kmat = 0.156;
    kcs = 1.1e-2;
    dmrna = 8.4e-2;
    dtlr = 4.5e-3;
    Ks = 4.2; % 0.8~4.2
    Kl = 30; % 14~30
    Ktlr = 6e-5;
%     dgfp = 0.00417;

    ydot = zeros(5, 1);
    % 函数命名 %
    % ydot(1) = d[DNA] / dt %
    % ydot(2) = d[TsR] / dt %
    % ydot(3) = d[TlR] / dt % 
    % ydot(4) = d[mRNA] / dt %
    % ydot(5) = d[GFP] / dt %
    % ydot(6) = d[GFP*] / dt %
    
%     ydot(2) = (-kcs*y(2)*y(1)/(Ks+y(1)))*(y(2) > 0 && y(2) > kcs*y(2)*y(1)/(Ks+y(1))); % 增加约束条件
%     ydot(3) = (-dtlr*y(3)/(Ktlr+y(3)))*(y(3) > 0 && y(3) > dtlr*y(3)/(Ktlr+y(3))); % 增加约束条件
%     ydot(2) = -kcs*y(2)*y(1)/(Ks+y(1));
%     ydot(3) = -dtlr*y(3)/(Ktlr+y(3));
    ydot(4) = kts*y(2)*y(1)/(Ks+y(1))*phiT-dmrna*y(4);
%     ydot(5) = ktl*y(3)*y(4)/(Kl+y(4))-kmat*y(5)-dgfp;
    ydot(5) = (ktl*y(3)*y(4)/(Kl+y(4))-kmat*y(5))*(y(5)+ktl*y(3)*y(4)/(Kl+y(4)) > kmat*y(5));
%     ydot(6) = kmat*y(5)-dgfp;

end