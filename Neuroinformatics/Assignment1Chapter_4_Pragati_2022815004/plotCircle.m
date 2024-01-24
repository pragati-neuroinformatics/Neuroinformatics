function a =  plotCircle(x, y, radius, color)
    theta = linspace(0, 2*pi, 100);
    xc = x + radius * cos(theta);
    yc = y + radius * sin(theta);
    fill(xc, yc, color, 'EdgeColor', 'none');
end


