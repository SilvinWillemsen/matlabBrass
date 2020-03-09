function [val] = clamp (input, min, max)
   if input < min
       val = min;
   elseif input > max
       val = max;
   else
       val = input;
   end
end