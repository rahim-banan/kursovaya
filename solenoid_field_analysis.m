function solenoid_field_analysis()
    % Параметры соленоида
    Radius = 1.0;          % Радиус соленоида (м)
    Solenoid_Length = 2.0; % Длина соленоида (м)
    N_Coil = 100;          % Количество витков
    I = 1.0;               % Ток (А)
    
    % Параметры области расчёта
    Xmin = -2.0; Xmax = 2.0; % Диапазон по X (м)
    Ymin = -2.0; Ymax = 2.0; % Диапазон по Y (м)
    Zmin = -1.0; Zmax = 1.0; % Диапазон по Z (м)
    N_Knot = 20;             % Количество узлов для расчёта
    
    % Точка для анализа поля
    X0 = 0.5; Y0 = 0.5; Z0 = 0.0; % Координаты точки (м)
    
    % Параметры для анализа вдоль линии
    Mod_R = 0.5;   % Начальное расстояние от оси (м)
    Phi = pi/4;    % Угол направления (рад)
    dR = 0.1;      % Шаг по расстоянию (м)
    N_Step = 10;   % Количество шагов
    
    % Выбор режима визуализации
    % 1 - Плоскость YoZ (X = X0)
    % 2 - Плоскость XoZ (Y = Y0)
    % 3 - Плоскость XoY (Z = Z0)
    % 4 - 3D пространство
    visualization_mode = 1;
    
    % Расчёт и визуализация
    switch visualization_mode
        case 1
            visualize_YoZ_plane(Radius, Solenoid_Length, N_Coil, I, X0, Ymin, Ymax, Zmin, Zmax, N_Knot);
        case 2
            visualize_XoZ_plane(Radius, Solenoid_Length, N_Coil, I, Y0, Xmin, Xmax, Zmin, Zmax, N_Knot);
        case 3
            visualize_XoY_plane(Radius, N_Coil, I, Z0, Xmin, Xmax, Ymin, Ymax, N_Knot);
        case 4
            visualize_3D_space(Radius, Solenoid_Length, N_Coil, I, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, N_Knot);
    end
    
    % Дополнительный анализ вдоль линии
    analyze_along_line(Radius, Solenoid_Length, N_Coil, I, Mod_R, Phi, dR, N_Step, Z0);
end

% Функции для расчёта компонент магнитного поля
function z = Bx(x0, y0, z0, a, h, N, I)
    f = inline('I*((z0-psi*h/(2*pi))*a.*cos(psi)-h/(2*pi)*(y0-a*sin(psi)))./(x0^2+y0^2+z0^2+a^2+(h/(2*pi)*psi).^2-2*(x0*a*cos(psi)+a*y0*sin(psi)+h/(2*pi)*z0*psi)).^(3/2)', 'psi', 'x0', 'y0', 'z0', 'a', 'h', 'I');
    phi = 0:pi/50:2*pi*N;
    F = f(phi, x0, y0, z0, a, h, I);
    z = trapz(phi, F);
end

function z = By(x0, y0, z0, a, h, N, I)
    f = inline('I*((x0-a*cos(psi))*h/(2*pi)+a*sin(psi).*(z0-psi*h/(2*pi)))./(x0^2+y0^2+z0^2+a^2+(h/(2*pi)*psi).^2-2*(x0*a*cos(psi)+a*y0*sin(psi)+h/(2*pi)*z0*psi)).^(3/2)', 'psi', 'x0', 'y0', 'z0', 'a', 'h', 'I');
    phi = 0:pi/50:2*pi*N;
    F = f(phi, x0, y0, z0, a, h, I);
    z = trapz(phi, F);
end

function z = Bz(x0, y0, z0, a, h, N, I)
    f = inline('I*((y0-a*sin(psi))*a.*sin(psi)+a*cos(psi).*(x0-a*cos(psi)))./(x0^2+y0^2+z0^2+a^2+(h/(2*pi)*psi).^2-2*(x0*a*cos(psi)+a*y0*sin(psi)+h/(2*pi)*z0*psi)).^(3/2)', 'psi', 'x0', 'y0', 'z0', 'a', 'h', 'I');
    phi = 0:pi/50:2*pi*N;
    F = f(phi, x0, y0, z0, a, h, I);
    z = trapz(phi, F);
end

% Функции для визуализации в разных плоскостях
function visualize_YoZ_plane(Radius, Solenoid_Length, N_Coil, I, X0, Ymin, Ymax, Zmin, Zmax, N_Knot)
    h = Solenoid_Length / N_Coil;
    y = linspace(Ymin, Ymax, N_Knot);
    z = linspace(Zmin, Zmax, N_Knot);
    [Y, Z] = meshgrid(y, z);
    B = zeros(size(Y));
    
    for i = 1:numel(Y)
        B(i) = sqrt(Bx(X0, Y(i), Z(i), Radius, h, N_Coil, I)^2 + ...
                    By(X0, Y(i), Z(i), Radius, h, N_Coil, I)^2 + ...
                    Bz(X0, Y(i), Z(i), Radius, h, N_Coil, I)^2);
    end
    
    figure;
    contourf(Y, Z, reshape(B, size(Y)), 20);
    colorbar;
    title('Магнитное поле в плоскости YoZ');
    xlabel('Y (м)');
    ylabel('Z (м)');
end

function visualize_XoZ_plane(Radius, Solenoid_Length, N_Coil, I, Y0, Xmin, Xmax, Zmin, Zmax, N_Knot)
    h = Solenoid_Length / N_Coil;
    x = linspace(Xmin, Xmax, N_Knot);
    z = linspace(Zmin, Zmax, N_Knot);
    [X, Z] = meshgrid(x, z);
    B = zeros(size(X));
    
    for i = 1:numel(X)
        B(i) = sqrt(Bx(X(i), Y0, Z(i), Radius, h, N_Coil, I)^2 + ...
               By(X(i), Y0, Z(i), Radius, h, N_Coil, I)^2 + ...
               Bz(X(i), Y0, Z(i), Radius, h, N_Coil, I)^2);
    end
    
    figure;
    contourf(X, Z, reshape(B, size(X)), 20);
    colorbar;
    title('Магнитное поле в плоскости XoZ');
    xlabel('X (м)');
    ylabel('Z (м)');
end

function visualize_XoY_plane(Radius, N_Coil, I, Z0, Xmin, Xmax, Ymin, Ymax, N_Knot)
    h = 0; % Для плоского случая (2D)
    x = linspace(Xmin, Xmax, N_Knot);
    y = linspace(Ymin, Ymax, N_Knot);
    [X, Y] = meshgrid(x, y);
    B = zeros(size(X));
    
    for i = 1:numel(X)
        B(i) = sqrt(Bx(X(i), Y(i), Z0, Radius, h, N_Coil, I)^2 + ...
               By(X(i), Y(i), Z0, Radius, h, N_Coil, I)^2 + ...
               Bz(X(i), Y(i), Z0, Radius, h, N_Coil, I)^2);
    end
    
    figure;
    contourf(X, Y, reshape(B, size(X)), 20);
    colorbar;
    title('Магнитное поле в плоскости XoY');
    xlabel('X (м)');
    ylabel('Y (м)');
end

function visualize_3D_space(Radius, Solenoid_Length, N_Coil, I, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, N_Knot)
    h = Solenoid_Length / N_Coil;
    x = linspace(Xmin, Xmax, N_Knot);
    y = linspace(Ymin, Ymax, N_Knot);
    z = linspace(Zmin, Zmax, N_Knot);
    [X, Y, Z] = meshgrid(x, y, z);
    B = zeros(size(X));
    
    for i = 1:numel(X)
        B(i) = sqrt(Bx(X(i), Y(i), Z(i), Radius, h, N_Coil, I)^2 + ...
               By(X(i), Y(i), Z(i), Radius, h, N_Coil, I)^2 + ...
               Bz(X(i), Y(i), Z(i), Radius, h, N_Coil, I)^2);
    end
    
    figure;
    slice(X, Y, Z, B, [Xmax, Ymax, Zmin]);
    shading interp;
    colorbar;
    title('Магнитное поле в 3D пространстве');
    xlabel('X (м)');
    ylabel('Y (м)');
    zlabel('Z (м)');
end

function analyze_along_line(Radius, Solenoid_Length, N_Coil, I, Mod_R, Phi, dR, N_Step, Z0)
    h = Solenoid_Length / N_Coil;
    R = zeros(1, N_Step);
    Mh = zeros(1, N_Step);
    
    for i = 1:N_Step
        x = Mod_R * cos(Phi) + dR * (i-1) * cos(Phi);
        y = Mod_R * sin(Phi) + dR * (i-1) * sin(Phi);
        R(i) = sqrt(x^2 + y^2);
        Hx = Bx(x, y, Z0, Radius, h, N_Coil, I);
        Hy = By(x, y, Z0, Radius, h, N_Coil, I);
        Hz = Bz(x, y, Z0, Radius, h, N_Coil, I);
        Mh(i) = sqrt(Hx^2 + Hy^2) / abs(Hz);
    end
    
    figure;
    plot(R, Mh, '-o');
    title('Зависимость магнитного поля вдоль линии');
    xlabel('Расстояние от оси (м)');
    ylabel('|B| / |Bz|');
    grid on;
end