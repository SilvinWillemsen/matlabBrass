%{
    Arrow function that draws arrowhead based on ratio of xlim and ylim of
    current axis. It is thus also important to define these limits
    beforehand!
    
    Input arguments are 
        - vector of two x-locations
        - vector of two y-locations 
%}

function arrow (x, y, varargin)
    if length(x) ~= 2
        error ("Vector x does not contain two values.");
    end
    
    if length(y) ~= 2
        error ("Vector y does not contain two values.");
    end
    
    if size(varargin) < 4
        color = 'k';
    else
        color = varargin{4};
    end
    
    if size(varargin) < 3
        arrowHeadLength = 0.2;
    else
        arrowHeadLength = varargin{3};
    end
    
    if size(varargin) < 2
        arrowHeadWidth = 0.02;
    else
        arrowHeadWidth = varargin{2};
    end
    
    if size(varargin) < 1
        lineWidth = 1;
    else
        lineWidth = varargin{1};
    end
        
    % calculate magnitude and angle (with the y axis) of arrow
    mag = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
    angle = atan2(y(2) - y(1), x(2) - x(1)) - 0.5 * pi;
    
    % obtain figure limits and ratio
    % calculate width and height of arrowhead
    triWidth = arrowHeadWidth * mag;
    triHeight = arrowHeadLength * mag;

    % calculate location of base of arrowhead
    midPointX = x(2) + sin(angle) * triHeight;
    midPointY = y(2) - cos(angle) * triHeight;
    
    % plot line until 0.9 times inside arrowhead
    plotUntilX = x(2) + sin(angle) * triHeight * 0.9;
    plotUntilY = y(2) - cos(angle) * triHeight * 0.9;
    
    % plot arrow line
    plot([x(1), plotUntilX], [y(1), plotUntilY], color, 'Linewidth', lineWidth);
    hold on;
    
    xLimSave = xlim;
    yLimSave = ylim;
    fig = gcf;
    ax = gca;
    figWH = fig.Position(3:4);
    axWH = ax.Position(3:4);
    xyRatio = axWH(2) / axWH(1) * figWH(2) / figWH(1) * (xLimSave(2) - xLimSave(1)) / (yLimSave(2) - yLimSave(1));

    
    % build arrowhead
    xTri1 = midPointX + cos(angle) * triWidth * xyRatio;
    xTri2 = midPointX - cos(angle) * triWidth * xyRatio;

    yTri1 = midPointY + sin(angle) * triWidth / xyRatio;
    yTri2 = midPointY - sin(angle) * triWidth / xyRatio;

    xTri = [x(2), xTri1, xTri2]; % x coordinates of arrowhead vertices
    yTri = [y(2), yTri1, yTri2]; % y coordinates of arrowhead vertices
    
    patch(xTri, yTri, color, 'EdgeColor', 'none');
end