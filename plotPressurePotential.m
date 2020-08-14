%{
    Function to plot the pressure potential in a 1D air tube. Required are
    the current state of the system (state) and the geometry (S). The 
    vectors need to be the same length.
%}
function plotPressurePotential (varargin)
    state = varargin{1} ;
    S = varargin{2};
    if length(varargin) == 3
        offset = varargin{3};
    else
        offset = 0;
    end
    for i = 1:length(state)-1
        j = i + offset;
        % we want to use the average state of the current index and the next for the colour
        avgState = (state(i) + state(i+1)) * 0.5;
        patch([j j+1 j+1 j], ... x-positions corners
              [S(i) S(i+1) -S(i+1) -S(i)], ... y-positions corners
              [clamp(avgState + 1, 0, 1), clamp(1-abs(avgState), 0, 1), clamp(1-avgState, 0, 1)], ... colour mapping
              'EdgeColor', 'none'); % set edge to no colour
    end 
end

function [val] = clamp (input, min, max)
   if input < min
       val = min;
   elseif input > max
       val = max;
   else
       val = input;
   end
end