% Создание простого GUI с кнопкой и графиком
function startup(app)
    app.Button.Text = 'Построить график';
    app.Button.ButtonPushedFcn = @plotData;
end

function plotData(src, event)
    x = linspace(0, 10, 100);
    y = sin(x);
    plot(app.UIAxes, x, y);
end