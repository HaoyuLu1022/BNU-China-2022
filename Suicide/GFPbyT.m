function ydot = GFPbyT(t, y)
    D = 3.4; % 23~236
%     phiT = 1; % T > 37 
    phiT = 0.172; % T <= 37 
%     phiT_MazF = 0.94;
    k = 0.49997147; % half life
    kts = 34.1;
    ktl = 16.1;
    kmat = 0.02;
    kcs = 1.1e-2;
    dmrna = 8.3e-2;
    dtlr = 4.5e-3;
    Ks = 4.2; % 0.8~4.2
    Kl = 30; % 14~30
    Ktlr = 6e-5;
    dgfp = 0.00417;

    ydot = zeros(6, 1);
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
    ydot(5) = (ktl*y(3)*y(4)/(Kl+y(4))-kmat*y(5)-dgfp)*(y(5)+ktl*y(3)*y(4)/(Kl+y(4)) > kmat*y(5)+dgfp);
%     ydot(6) = kmat*y(5)-dgfp;

end
