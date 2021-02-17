%{
    Max bend = pi
%}

function myArc (x, y, bend, scalings, reverse, varargin)

    if length(x) ~= 2
        error ("Vector x does not contain two values.");
    end
    
    if length(y) ~= 2
        error ("Vector y does not contain two values.");
    end

    if bend == 0
        error ("Just use 'plot' if you aren't going to give it any bend :)")
    else
        origDist = (sign(bend) * pi - bend) / bend;
    end
    
    if size(varargin) < 3
        lineStyle = '-';
    else
        lineStyle = varargin{3};
    end
    
    if size(varargin) < 2
        color = 'k';
    else
        color = varargin{2};
    end
    
    if size(varargin) < 1
        lineWidth = 1;
    else
        lineWidth = varargin{1};
    end
    
    centerLocX = (x(1)+x(2)) * 0.5;
    centerLocY = (y(1)+y(2)) * 0.5;

    angle = atan2(y(2) - y(1), x(2) - x(1));
    
    originLocX = centerLocX - sign(bend) * origDist * sin(angle);
    originLocY = sign(bend) * origDist * cos(angle) + centerLocY;
    radius = sqrt((originLocX - x(1))^2 + (originLocY - y(1))^2);
    
%     scatter(originLocX, originLocY)
%     hold on;
%     scatter(centerLocX, centerLocY)
%     scatter(x, y)    
    
    if bend < 0
        startAngle = (0.5 * pi + angle) + acos(origDist / radius);
        endAngle = (0.5 * pi + angle) - acos(origDist / radius);
    else
        startAngle = -(0.5 * pi - angle) - acos(origDist / radius);
        endAngle = -(0.5 * pi - angle) + acos(origDist / radius);
    end
    detail = 100;
    diff = endAngle - startAngle;
    arcX = zeros(detail, 1);    
    arcY = zeros(detail, 1);

    idx = 1;
    
%     if bend < 0
% %         tmpAngle = endAngle;
% %         startAngle = tmpAngle;
% %         endAngle = startAngle;
%         thetaRange = endAngle:-diff/(detail-1):startAngle;
%     else
        thetaRange = startAngle:diff/(detail-1):endAngle;
%     end
    for theta = thetaRange
        arcX(idx) = (radius * real(exp(1i*theta)) + originLocX);
        arcY(idx) = (radius * imag(exp(1i*theta)) + originLocY);
        idx = idx + 1;
    end
    
    plot(arcX, arcY, 'Linewidth', lineWidth, 'color', color, 'Linestyle', lineStyle)
    drawnow;
end