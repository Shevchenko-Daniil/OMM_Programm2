N_x = 300;
N_y = 300;
N_t = 100;
t_max = 1;
h_x = 1/(N_x-1);
h_y = 0.5*pi/(N_y-1);
tau = t_max/(N_t - 1);
x = 0:h_x:1;
y = 0:h_y:0.5*pi;
t = 0:tau:t_max;

t_point = 1;

u_t_point = zeros(N_x, N_y);

for i=1:N_x
    for j=1:N_y
        u_t_point(i, j) = sin(y(j)).*(T_0_func(t_point) + sum_func(x(i), t_point));
    end
end

figure(2)
    surf(x, y, u_t_point', EdgeColor = 'none')
    set(findall(figure(2),'type','axes'),'fontsize',13)
    xlabel('x', 'FontSize',18)
    ylabel('y', 'FontSize',18)
    zlabel('u', 'FontSize',18)
    colorbar
    grid on
%% этот кусок записывает avi-анимацию
v = VideoWriter('analit_resol.avi');
open(v)
for l=1:N_t
    t_point = t(l);
    for i=1:N_x
        for j=1:N_y
            u_t_point(i, j) = sin(y(j))*(T_0_func(t_point) + sum_func(x(i), t_point));
        end
    end

    imagesc(x, y, u_t_point', [0 1])
    title(['t = ', num2str(t_point), ' c'], 'FontSize', 36)
    colorbar
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%


function T_n = T_n_func(t, n)
    T_n = -(4*t.^2)./(pi^2 * n^2 * (pi^2 * n.^2 + 1)) + (8.*t)./(pi^2 * n.^2 .* (pi^2 * n.^2 + 1).^2) - 8.*(1 - (exp(-((pi^2)*(n^2).*t+t)) ) )./(pi^2 * n^2 * (pi^2 * n^2 + 1).^3);
end

function T_0 = T_0_func(t)
    T_0 = -exp(-t) + 0.5*t^2 - t + 1;
end
function sum = sum_func(x, t)
    sum = 0;
    for k=0:300
        sum = sum + cos(pi*(2*k+1)*x).*T_n_func(t, 2*k+1);
    end
end








