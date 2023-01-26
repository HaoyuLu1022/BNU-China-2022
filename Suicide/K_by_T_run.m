P = [0, 0, 0, 0, 0];
tmesh = 0:1:7200;
[t,y] = ode45(@K_by_T, tmesh, P);
figure;
plot(t, y(:, 3), 'r');
hold on;
plot(t, y(:, 4), 'b');
hold on;
plot(t, y(:, 5), 'g');
legend('mazF', 'mazE', 'mazEF');
xlabel('Time, s')
ylabel('Concentration, M')
title('< 37℃')
% title(">=37℃")
 
% plot(t, y(:, 3), 'r');
% title('mazF concentration (>37 degree Celsius)');
% xlabel('time, s');
% ylabel('concentration');