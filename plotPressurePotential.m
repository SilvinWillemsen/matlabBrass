%{
    Function to plot the pressure potential in a 1D air tube. Required are
    the 
%}

function plotPressurePotential (state, S)
    for i = 1:length(state)-1
        avgState = (state(i) + state(i+1)) * 0.25;
        patch([i i+1 i+1 i], ... x-positions corners
              [S(i) S(i+1) -S(i+1) -S(i)], ... y-positions corners
              [clamp(avgState + 1, 0, 1), clamp(1-abs(avgState), 0, 1), clamp(1-avgState, 0, 1)], ... colour mapping
              'EdgeColor', 'none'); % set edge to no colour
    end
