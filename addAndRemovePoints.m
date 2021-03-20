%% Add and remove points
if abs(N - NPrev) > 1
    disp('too fast')
end

% add point if N^n > N^{n-1}
if N > NPrev
    customIpNorm = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
        2 * alf / (alf + 2); ...
        2 / (alf + 2); ...
        2 * alf / -((alf + 3) * (alf + 2))]';
     customIp = [-((alf + 1)*(alf + 2))/((alf + 3)*(alf + 4)), ...
                      (2*(alf + 1))/(alf + 3), ...
                                  2/(alf + 3), ...
         -(2*(alf + 1))/((alf + 3)*(alf + 4))];
    if ~connectedWithP % if the new N connects at v prepare the statevectors (by adding to v)
        %             upMp1 = customIp * [up(end-1:end); wp(1:2)]; % unnecessary?
        upMp1Prev = customIpNorm * [upPrev(end-1:end); wpPrev(1:2)];
        
        %             wpm1 = fliplr(customIp) * [up(end-1:end); wp(1:2)]; % unnecessary?
        wpm1Prev = fliplr(customIpNorm) * [upPrev(end-1:end); wpPrev(1:2)];
        
%         uvNextInnerBoundarySave = [uvNext(end-1:end); wvNextmh; wvNext(1)];
%         uvInnerBoundarySave = [uv(end-1:end); wvmh; wv(1)];
%         
%         wvNextInnerBoundarySave = [uvNext(end); uvNextMph; wvNext(1:2)];
%         wvInnerBoundarySave = [uv(end); uvMph; wv(1:2)];
        innerBoundarySaveNext = [uvNext(end-1:end); wvNext(1:2)];
        innerBoundarySave = [uv(end-1:end); wv(1:2)];

        pointToAddNext = (uvNextMph + wvNextmh) * 0.5;
        pointToAdd = (uvMph + wvmh) * 0.5;

%         uvNext = [uvNext; customIp * innerBoundarySaveNext];
%         uv = [uv; customIp * innerBoundarySave];
%         
%         wvNext = [fliplr(customIp) * innerBoundarySaveNext; wvNext];
%         wv = [fliplr(customIp) * innerBoundarySave; wv];
        uvNext = [uvNext; pointToAddNext];
        uv = [uv; pointToAddNext];
        
        wvNext = [pointToAddNext; wvNext];
        wv = [pointToAdd; wv];

    else % if the new N connects at p prepare the statevectors (by adding to p)
        %% Calculate interpolated velocities
        uvMph = uv(end) * quadIp(3) + wv(1) * quadIp(2) + wv(2) * quadIp(1);
        wvmh = uv(end-1) * quadIp(1) + uv(end) * quadIp(2) + wv(1) * quadIp(3);

        %% Calculate p^n
        up(upRange) = upPrev(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uv(upRange) - SHalf(upRange-1) .* uv(upRange-1));
        upMp1 = upMp1Prev - rho * c * lambda / SBar(length(up)) * (SHalf(length(up)) .* uvMph - SHalf(length(up) - 1) .* uv(end));
%         if connectedToLip
            up(1) = upPrev(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uv(1));
%         end

        wp(wpRange-1) = wpPrev(wpRange-1) - rho * c * lambda ./ SBar(wpRange + length(up) - 1) .* (SHalf(wpRange + length(up) - 2) .* wv(wpRange) - SHalf(wpRange + length(up) - 3) .* wv(wpRange-1));
        wpm1 = wpm1Prev - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wv(1) - SHalf(length(up) - 1) .* wvmh);
        
        if radiation
            wp(end) = ((1 - rho * c * lambda * z3) * wpPrev(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wv(end))/SBar(end))) / (1 + rho * c * lambda * z3);
        end
        
        v1Next = v1 + k / (2 * Lr) * (wp(end) + wpPrev(end));
        p1Next = z1 / 2 * (wp(end) + wpPrev(end)) + z2 * p1;
        
        v1 = v1Next;
        p1 = p1Next;
        
        innerBoundarySaveNext = [upNext(end-1:end); wpNext(1:2)];
        innerBoundarySave = [up(end-1:end); wp(1:2)];
        innerBoundarySavePrev = [upPrev(end-1:end); wpPrev(1:2)];

%         upNextInnerBoundarySave = [upNext(end); upMp1; wpNext(1:2)];
%         upInnerBoundarySave = [up(end); upMp1; wp(1:2)];
%         upPrevInnerBoundarySave = [upPrev(end); upMp1Prev; wpPrev(1:2)];
% 
%         wpNextInnerBoundarySave = [upNext(end-1:end); wpm1; wpNext(1)];
%         wpInnerBoundarySave = [up(end-1:end); wpm1; wp(1)];
%         wpPrevInnerBoundarySave = [upPrev(end-1:end); wpm1Prev; wpPrev(1)];

        uvNextMph = customIpNorm * [uvNext(end-1:end); wvNext(1:2)]; % unnecessary?
        uvMph = customIpNorm * [uv(end-1:end); wv(1:2)];
        
        wvNextmh = fliplr(customIpNorm) * [uvNext(end-1:end); wvNext(1:2)]; % unnecessary?
        wvmh = fliplr(customIpNorm) * [uv(end-1:end); wv(1:2)];
        
        pointToAdd = (upMp1 + wpm1) * 0.5;
        pointToAddPrev = (upMp1Prev + wpm1Prev) * 0.5;

%         upNext = [upNext; customIp * innerBoundarySaveNext];
%         up = [up; customIp * innerBoundarySave];
%         upPrev = [upPrev; customIp * innerBoundarySavePrev];
%         
%         wpNext = [fliplr(customIp) * innerBoundarySaveNext; wpNext];
%         wp = [fliplr(customIp) * innerBoundarySave; wp];
%         wpPrev = [fliplr(customIp) * innerBoundarySavePrev; wpPrev];

        upNext = [upNext; pointToAdd];
        up = [up; pointToAdd];
        upPrev = [upPrev; pointToAddPrev];
        
        wpNext = [pointToAdd; wpNext];
        wp = [pointToAdd; wp];
        wpPrev = [pointToAddPrev; wpPrev];

    end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n, setToOnes);
    % insert matrix creation here
    
    statesSave = [statesSave; [up(end-1), up(end), wp(1), wp(2), uv(end-1), uv(end), wv(1), wv(2), uvMph, wvmh] ];
%     disp("point added")
%     plotStatesAfterAddingPoint;
end

% remove point if N^n < N^{n-1}
if N < NPrev
    if connectedWithP % remove from v to be connected at p
        uvNext = uvNext(1:end-1);
        uv = uv(1:end-1);
        wvNext = wvNext(2:end);
        wv = wv(2:end);
    else  % remove from p to be connected at v
        upNext = upNext(1:end-1);
        up = up(1:end-1);
        upPrev = upPrev(1:end-1);
        
        wpNext = wpNext(2:end);
        wp = wp(2:end);
        wpPrev = wpPrev(2:end);
        
    end
    if flag
        disp("point removed")
    end
    [S, SHalf, SBar] = setTube(N+1, NnonExtended, n);
    
    statesSave = [statesSave; [up(end-1), up(end), wp(1), wp(2), uv(end-1), uv(end), wv(1), wv(2), uvMph, wvmh] ];
    
end
if connectedWithP
    upRange = 2:length(up)-1;         % range without boundaries
    wpRange = 2:length(wp)-1;
else
    upRange = 2:length(up);         % range without boundaries
    wpRange = 2:length(wp);
end


