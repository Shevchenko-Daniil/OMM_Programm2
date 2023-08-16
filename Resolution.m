N_x = 300;
N_y = 300;
N_t = 400;
t_max = 1;
h_x = 1/(N_x-1);
h_y = 0.5*pi/(N_y-1);
tau = t_max/(N_t - 1);
x = 0:h_x:1;
y = 0:h_y:0.5*pi;
t = 0:tau:t_max;
gamma_x = tau/(h_x^2);
gamma_y = tau/(h_y^2);
u = zeros(N_x, N_y, 2*N_t+1);


u(:, :, 1) = 0;     %начальное условие
d_x = zeros(1 ,N_x);
sigma_x = zeros(1, N_x);
d_y = zeros(1, N_y);
sigma_y = zeros(1, N_y);




for j=2:2:2*N_t
    for i2=2:N_y-1
        d_x(2)= 1;          % из граничных условий в x=0
        sigma_x(2) = 0;     % из граничных условий в x=0
        A = 0.5*gamma_x;
        B = 1 + gamma_x;
        C= 0.5*gamma_x;
        for m=2:N_x-1
            F_m_x =  -(0.5.*gamma_y .* (u(m, i2-1, j-1) + u(m, i2+1, j-1))+ (1 - gamma_y) .* u(m, i2, j-1) + 0.5 .* tau .* x(m)*((tau.*(j + 1)./2).^2) .* sin(y(i2)));
            d_x(m+1) = C./(B-A.*d_x(m));
            sigma_x(m+1) = (F_m_x - A .* sigma_x(m))./(A .* d_x(m) - B);
        end
        d_x(2)= 1;
        sigma_x(2) = 0;
        u(N_x, i2, j) = sigma_x(N_x)./(1 - d_x(N_x));            % из граничных условий в x=1
        for m=N_x:-1:2
            u(m-1, i2, j) = d_x(m) .* u(m, i2, j) + sigma_x(m);
        end
    end
    
    
    for i1=2:N_x-1
        d_y(2) = 0;             % из граничных условий в y=0
        sigma_y(2) = 0;         % из граничных условий s y=0
        A =  0.5*gamma_y;
        B = 1 + gamma_y;
        C =  0.5*gamma_y;
        for m=2:N_y-1
            F_m_y =  -(0.5.*gamma_x.*(u(i1-1, m, j) + u(i1+1, m, j))+ (1 - gamma_x).*u(i1, m, j) + 0.5.*tau.*x(i1).* ((tau*(j-1)./2).^2) .* sin(y(m)));
            d_y(m+1) = C ./(B - A .* d_y(m));
            sigma_y(m+1) = (F_m_y - A .* sigma_y(m))./(A .* d_y(m) - B);
        end
    
        u(i1, N_y, j+1) = sigma_y(N_y)./(1 - d_y(N_y));            % из граничных условий в y=pi/2
        for m=N_y:-1:2
            u(i1, m-1, j+1) = d_y(m).*u(i1, m, j+1) + sigma_y(m);
        end
        
    end
end
%% этот кусок записывает avi-анимацию
v = VideoWriter('resol.avi');
open(v)

figure(10)

for l=2:2:2*N_t
    t_point = t(l./2);

    imagesc(x, y(1:N_y-1), u(:, 1:N_y-1, l)', [0 6])
    title(['t = ', num2str(t_point), ' c'], 'FontSize', 36)
    xlabel('x', 'FontSize',18)
    ylabel('y', 'FontSize',18)
    colorbar
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
%% построение решения
figure(1)
    surf(x, y(1:N_y-1), u(:, 1:N_y-1, 400)', EdgeColor='none')
    set(findall(figure(1),'type','axes'),'fontsize',13)
    xlabel('x', 'FontSize',18)
    ylabel('y', 'FontSize',18)
    zlabel('u', 'FontSize',18)
    colorbar
    grid on
    
    

    %% построение относительной ошибки, нужны данные из аналитического решения
    
    err = 100.*(u_t_point - u(:, :, 800))./u_t_point;
figure(32)
    surf(x(3:N_x-2), y(3:N_y-2), abs(err(3:N_x-2, 3:N_y-2)), EdgeColor='none')
    set(findall(figure(32),'type','axes'),'fontsize',13)
    xlabel('x', 'FontSize',18)
    ylabel('y', 'FontSize',18)
    zlabel('Отн. ошибка, %', 'FontSize',18)
    colorbar
    grid on







