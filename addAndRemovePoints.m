%% Add and remove points
if abs(N - NPrev) > 1
    disp('too fast')
end

useAvg = false;

% add point if N^n > N^{n-1}
if N > NPrev
 
    customIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
        2 * alf / (alf + 2); ...
        2 / (alf + 2); ...
        2 * alf / -((alf + 3) * (alf + 2))]';
    if alternatePV
        if ~connectedWithP % if the new N connects at v prepare the statevectors (by adding to v)
            upMp1Prev = customIp * [upPrev(end-1:end); wpPrev(1:2)];
            wpm1Prev = fliplr(customIp) * [upPrev(end-1:end); wpPrev(1:2)];

            if useAvg
                vNextAvg = (uvNextMph + wvNextmh) * 0.5;
                vAvg = (uvMph + wvmh) * 0.5;

                uvNext = [uvNext; vNextAvg];
                uv = [uv; vAvg];

                wvNext = [vNextAvg; wvNext];
                wv = [vAvg; wv];
            else
                uvNext = [uvNext; uvNextMph];
                uv = [uv; uvMph];

                wvNext = [wvNextmh; wvNext];
                wv = [wvmh; wv];
            end
        else % if the new N connects at p prepare the statevectors (by adding to p)
            % NOT SURE IF ALL OF THIS STILL WORKS NOW THAT VIRTUAL GRID POINTS ARE INCLUDED
            %% Here, also calculate p^n as that hasn't been done yet when switching to this setting
            uvMph = uv(end) * quadIp(3) + wv(2) * quadIp(2) + wv(3) * quadIp(1);
            wvmh = uv(end-1) * quadIp(1) + uv(end) * quadIp(2) + wv(2) * quadIp(3);

            %% Calculate p^n
            up(upRange) = upPrev(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uv(upRange) - SHalf(upRange-1) .* uv(upRange-1));
            upMp1 = upMp1Prev - rho * c * lambda / SBar(length(up)) * (SHalf(length(up)) .* uvMph - SHalf(length(up) - 1) .* uv(end));
            up(1) = upPrev(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uv(1));

            wp(wpRange-1) = wpPrev(wpRange-1) - rho * c * lambda ./ SBar(wpRange + length(up) - 2) .* (SHalf(wpRange + length(up) - 2) .* wv(wpRange) - SHalf(wpRange + length(up) - 3) .* wv(wpRange-1));
            wpm1 = wpm1Prev - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wv(1) - SHalf(length(up) - 1) .* wvmh);

            if radiation
                wp(end) = ((1 - rho * c * lambda * z3) * wpPrev(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wv(end))/SBar(end))) / (1 + rho * c * lambda * z3);
            end

            v1Next = v1 + k / (2 * Lr) * (wp(end) + wpPrev(end));
            p1Next = z1 / 2 * (wp(end) + wpPrev(end)) + z2 * p1;

            v1 = v1Next;
            p1 = p1Next;

            uvMph = customIp * [uv(end-1:end); wv(1:2)];
            wvmh = fliplr(customIp) * [uv(end-1:end); wv(1:2)];

            if useAvg
                pAvg = (upMp1 + wpm1) * 0.5;
                pPrevAvg = (upMp1 + wpm1) * 0.5;

                upNext = [upNext; pAvg];
                up = [up; pAvg];
                upPrev = [upPrev; pPrevAvg];

                wpNext = [pAvg; wpNext];
                wp = [pAvg; wp];
                wpPrev = [pPrevAvg; wpPrev];
            else
                upNext = [upNext; upMp1];
                up = [up; upMp1];
                upPrev = [upPrev; upMp1Prev];

                wpNext = [wpm1; wpNext];
                wp = [wpm1; wp];
                wpPrev = [wpm1Prev; wpPrev];
            end
        end
    else
        if mod(N,2) == 1
            vNextDiff = wvNext(1) - uvNext(end);
            vDiff = wv(1) - uv(end);
            vPrevDiff = wvPrev(1) - uvPrev(end);

            uvNext = [uvNext; customIp * [uvNext(end-1:end); wvNext(2:3)-vNextDiff]];
            uv = [uv; customIp * [uv(end-1:end); wv(2:3)-vDiff]];
            uvPrev = [uvPrev; customIp * [uvPrev(end-1:end); wvPrev(2:3)-vPrevDiff]];
            
%             hold off;
%             plot([uvNext(end-1:end); wvNext(2:3)-vNextDiff])
%             hold on;
%             plot([uv(end-1:end); wv(2:3)-vDiff])
%             plot([uvPrev(end-1:end); wvPrev(2:3)-vPrevDiff])
%             title("added to u")
%             drawnow;
%             pause(0.5)
            
            upNext = [upNext; customIp * [upNext(end-1:end); wpNext(1:2)]];
            up = [up; customIp * [up(end-1:end); wp(1:2)]];
            upPrev = [upPrev; customIp * [upPrev(end-1:end); wpPrev(1:2)]];
            
        else 
            vNextDiff = wvNext(1) - uvNext(end);
            vDiff = wv(1) - uv(end);
            vPrevDiff = wvPrev(1) - uvPrev(end);
            
            wvNext = [fliplr(customIp) * [uvNext(end-2:end-1)+vNextDiff; wvNext(1:2)]; wvNext];
            wv = [fliplr(customIp) * [uv(end-2:end-1)+vDiff; wv(1:2)]; wv];
            wvPrev = [fliplr(customIp) * [uvPrev(end-2:end-1)+vPrevDiff; wvPrev(1:2)]; wvPrev];
          
%             hold off;
%             plot([uvNext(end-2:end-1)+vNextDiff; wvNext(1:2)])
%             hold on;
%             plot([uv(end-2:end-1)+vDiff; wv(1:2)])
%             plot([uvPrev(end-2:end-1)+vPrevDiff; wvPrev(1:2)])
%             title("added to w")
%             pause(0.5)
%             drawnow;
            
            wpNext = [fliplr(customIp) * [upNext(end-1:end); wpNext(1:2)]; wpNext];
            wp = [fliplr(customIp) * [up(end-1:end); wp(1:2)]; wp];
            wpPrev = [fliplr(customIp) * [upPrev(end-1:end); wpPrev(1:2)]; wpPrev];

        end
    end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n, setToOnes);
    
    disp("point added " + num2str(alf) + " " + num2str(N))
%     plotStatesAfterAddingPoint;
end

% remove point if N^n < N^{n-1}
if N < NPrev
    if alternatePV
        if connectedWithP % remove from v to be connected at p
            % NOT SURE IF ALL OF THIS STILL WORKS NOW THAT VIRTUAL GRID POINTS ARE INCLUDED

            %% Here, also calculate p^n as that hasn't been done yet when switching to this setting
            uv(end) = uv(end-1) * quadIp(3) + wv(1) * quadIp(2) + wv(2) * quadIp(1);
            wv(1) = uv(end-2) * quadIp(1) + uv(end-1) * quadIp(2) + wv(1) * quadIp(3);

            %% Calculate p^n
            up(upRange) = upPrev(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uv(upRange) - SHalf(upRange-1) .* uv(upRange-1));
            upMp1 = upMp1Prev - rho * c * lambda / SBar(length(up)) * (SHalf(length(up)) .* uvMph - SHalf(length(up) - 1) .* uv(end));
            up(1) = upPrev(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uv(1));

            wp(wpRange-1) = wpPrev(wpRange-1) - rho * c * lambda ./ SBar(wpRange + length(up) - 2) .* (SHalf(wpRange + length(up) - 2) .* wv(wpRange) - SHalf(wpRange + length(up) - 3) .* wv(wpRange-1));
            wpm1 = wpm1Prev - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wv(1) - SHalf(length(up) - 1) .* wvmh);

            if radiation
                wp(end) = ((1 - rho * c * lambda * z3) * wpPrev(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wv(end))/SBar(end))) / (1 + rho * c * lambda * z3);
            end

            v1Next = v1 + k / (2 * Lr) * (wp(end) + wpPrev(end));
            p1Next = z1 / 2 * (wp(end) + wpPrev(end)) + z2 * p1;

            v1 = v1Next;
            p1 = p1Next;
            
            uvNextMph = uvNext(end);
            uvMph = uv(end);
            
            uvNext = uvNext(1:end-1);
            uv = uv(1:end-1);
            uvPrev = uvPrev(1:end-1);
            
            wvNextmh = wvNext(1);
            wvmh = wv(1);
            
            wvNext = wvNext(2:end);
            wv = wv(2:end);
            wvPrev = wvPrev(2:end);

        else  % remove from p to be connected at v
            upNext = upNext(1:end-1);
            up = up(1:end-1);
            upPrev = upPrev(1:end-1);

            wpNext = wpNext(2:end);
            wp = wp(2:end);
            wpPrev = wpPrev(2:end);

        end
    else
        if mod(N,2) == 0
            % useavg is gone here now (using virtual gridpoints within the
            % vectors)
            
            uvNext = uvNext(1:end-1);
            uv = uv(1:end-1);
            uvPrev = uvPrev(1:end-1);
            
            upNext = upNext(1:end-1);
            up = up(1:end-1);
            upPrev = upPrev(1:end-1);
        else 
           % useavg is gone here now (using virtual gridpoints within the
            % vectors)

            wvNext = wvNext(2:end);
            wv = wv(2:end);
            wvPrev = wvPrev(2:end);

            wpNext = wpNext(2:end);
            wp = wp(2:end);
            wpPrev = wpPrev(2:end);
        end
    end
%     if flag
    disp("point removed " + num2str(alf) + " " + num2str(N))
%     end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n, setToOnes);
        
end
if connectedWithP
    upRange = 2:length(up);         
    wpRange = 1:length(wp)-1;
else
    upRange = 2:length(up);         
    wpRange = 2:length(wp);
end


