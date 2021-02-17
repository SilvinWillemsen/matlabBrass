%{
    Arrow function that draws arrowhead based on ratio of xlim and ylim of
    current axis. It is thus also important to prepare your figure and
    define these limits beforehand!
    
    Input arguments are 
        - vector of two x-locations
        - vector of two y-locations 
        - Linewidth
        - Arrowhead length (in px)
        - Arrowhead width (in px)
        - Color
        - Linestyle
        - Bend (0-pi rad)
%}

function myArrow (x, y, varargin)
    if length(x) ~= 2
        error ("Vector x does not contain two values.");
    end
    
    if length(y) ~= 2
        error ("Vector y does not contain two values.");
    end
    if size(varargin) < 7
        reverseBend = false;
    else
        reverseBend = varargin{7};
    end
    if size(varargin) < 6
        bend = 0;
    else
        bend = varargin{6};
    end
    if size(varargin) < 5
        lineStyle = '-';
    else
        lineStyle = varargin{5};
    end
    
    if size(varargin) < 4
        color = 'k';
    else
        color = varargin{4};
    end
    
    if size(varargin) < 3
        arrowHeadLength = 50;
    else
        arrowHeadLength = varargin{3};
    end
    
    if size(varargin) < 2
        arrowHeadWidth = 50;
    else
        arrowHeadWidth = varargin{2};
    end
    
    if size(varargin) < 1
        lineWidth = 1;
    else
        lineWidth = varargin{1};
    end
    hold on;
    
    % obtain figure limits and ratio
    xLimSave = xlim;
    yLimSave = ylim;
    fig = gcf;
    ax = gca;
    figWH = fig.Position(3:4);
    axWH = ax.Position(3:4);
    
    xDistInPx = (x(2) - x(1)) / (xLimSave(2) - xLimSave(1)) * (figWH(1) * axWH(1));
    yDistInPx = (y(2) - y(1)) / (yLimSave(2) - yLimSave(1)) * (figWH(2) * axWH(2));
    angle = atan2(yDistInPx, xDistInPx) + bend * 0.5;
    
    % calculate location of base of arrowhead
    pixelsTipToMidPoint = arrowHeadLength;
    
    pixelsTipToMidPointX = pixelsTipToMidPoint * cos(angle);
    pixelsTipToMidPointY = pixelsTipToMidPoint * sin(angle);

    xScaling = (xLimSave(2) - xLimSave(1)) / (figWH(1) * axWH(1));
    yScaling = (yLimSave(2) - yLimSave(1)) / (figWH(2) * axWH(2));

    midPointX = x(2) - xScaling * pixelsTipToMidPoint * cos(angle);    
    midPointY = y(2) - yScaling * pixelsTipToMidPoint * sin(angle);    

    % plot line until 0.9 times inside arrowhead
    plotUntilX = x(2) - xScaling * 0.9 * pixelsTipToMidPoint * cos(angle);    
    plotUntilY = y(2) - yScaling * 0.9 * pixelsTipToMidPoint * sin(angle);   %     xyRatio = axWH(2) / axWH(1) * figWH(2) / figWH(1) * (xLimSave(2) - xLimSave(1)) / (yLimSave(2) - yLimSave(1));
        
    % plot arrow line
    if bend == 0
        plot([x(1), plotUntilX], [y(1), plotUntilY], 'color', color, 'Linewidth', lineWidth, 'Linestyle', lineStyle);    
    else
        plotUntilX = x(2) - xScaling * pixelsTipToMidPoint * cos(angle);    
        plotUntilY = y(2) - yScaling * pixelsTipToMidPoint * sin(angle);   %     xyRatio = axWH(2) / axWH(1) * figWH(2) / figWH(1) * (xLimSave(2) - xLimSave(1)) / (yLimSave(2) - yLimSave(1));
%         myArc([x(1), x(2)], [y(1), y(2)], bend, [xScaling, yScaling], reverseBend, lineWidth, color, lineStyle);  
        myArc([x(1), midPointX], [y(1), midPointY], bend, [xScaling, yScaling], reverseBend, lineWidth, color, lineStyle);    

    end
    % build arrowhead
    pixelsMidPointToCorner = arrowHeadWidth * 0.5;
    xTri1 = midPointX - cos(0.5 * pi - angle) * xScaling * pixelsMidPointToCorner; 
    xTri2 = midPointX + cos(0.5 * pi - angle) * xScaling * pixelsMidPointToCorner; 
    
    yTri1 = midPointY + sin(0.5 * pi - angle) * yScaling * pixelsMidPointToCorner; 
    yTri2 = midPointY - sin(0.5 * pi - angle) * yScaling * pixelsMidPointToCorner; 

    xTri = [x(2), xTri1, xTri2]; % x coordinates of arrowhead vertices
    yTri = [y(2), yTri1, yTri2]; % y coordinates of arrowhead vertices

    patch(xTri, yTri, color, 'EdgeColor', 'none');
end