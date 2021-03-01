function [s, h_z, h_r] = geom_thickness(zw, rw, zt, rt, ztm, rtm)

    % Assuming that at the input we get an array of points of the section 
    % of the workpiece with an equal step in z, we calculate this step
    hzw = zw(2) - zw(1);
    
    % We equate arrays of cutting lengths and thicknesses to zero values 
    % if the cutting edge is to the left of the left end face or to the 
    % right of the right end face of the workpiece
    if( (zw(1) - zt(length(zt)) >= 0) || (zt(1) - zw(length(zw)) >= 0) )
        s = 0;
        h_z = 0;
        h_r = 0;
    else
        % Select from arrays of points of the workpiece surface
        % zw, rw only those that are below (or above) the tool

        % the index of the first point in the array zw located under the tool
        qmin = max((floor(zt(1) / hzw) + 1), 1);

        % the index of the last point in the array zw located under the tool
        qmax = min((ceil(zt(length(zt)) / hzw) + 1), length(zw));


        % Create new geometry arrays for points below (or above) the cutting edge
        zwc = zeros(1, qmax - qmin + 1);
        rwc = zeros(1, qmax - qmin + 1);
        q = qmin;
        for i=1:(length(zwc))
            zwc(i) = zw(q);
            rwc(i) = rw(q);
            q = q+1;
        end

        % Further, we work only with this array.
        % Check if the cutting edge is much higher than
        % the material surface and there is no cutting.
        if(min(rt) >= max(rwc))
            s = 0;
            h_z = 0;
            h_r = 0;
        else
            % Otherwise we are looking for points of entry and exit from the material
            j = 1;
            enter = 0;
            ext = 0;

            % Cycle walking through the elements of cutting edge
            while((j < length(zt)) && (ext == 0))
                % We are looking for the first point of cutting edge that dives into the material
                q = max((floor(zt(j) / hzw) - qmin + 2), 1);
                q_max = min((ceil(zt(j+1) / hzw) - qmin + 2), length(zwc));
                
                % Cycle walking through the elements of the workpiece
                while((q < q_max) && (enter == 0))
                    [zp, rp] = geom_cross1(zt(j), rt(j), zt(j+1), rt(j+1), zwc(q), rwc(q), zwc(q+1), rwc(q+1));
                    % If this point lies within both elements, then the entry point
                    % to the material has been found.
                    if((zt(j) <= zp) && (zp <= zt(j+1)) && (zwc(q) <= zp) && (zp <= zwc(q+1)))
                        z_ent = zp;
                        r_ent = rp;
                        enter = 1;
                        % Memorize the index of the tool element entering into material
                        j_ent = j;
                        q_ent = q;
                    else
                        % Otherwise continue iterations through the elements of the workpiece.
                        q = q + 1;
                    end
                end
                
                % Then, we are looking for exit point of cutting edge from the material
                while((q < q_max) && (enter == 1) && (ext == 0))
                    [zp, rp] = geom_cross1(zt(j), rt(j), zt(j+1), rt(j+1), zwc(q), rwc(q), zwc(q+1), rwc(q+1));
                    if((zt(j) <= zp) && (zp <= zt(j+1)) && (zwc(q) <= zp) && (zp <= zwc(q+1)) && (zp > z_ent))
                        % if find intersection point
                        z_ext = zp;
                        r_ext = rp;
                        ext = 1;
                        % Memorize the index of the tool element exiting from material
                        j_ext = j;
                        q_ext = q;
                    else
                        % Otherwise continue iterations through the elements of the workpiece.
                        q = q + 1;
                    end
                end
                    j = j + 1;
            end
            
            % If enter point to material has been found, but exit point hasn't
            % we assume that the exit point is the last point of the workpirce section
            if((enter == 1) && (ext == 0))
                 z_ext = zt(length(zt));
                 r_ext = rt(length(rt));
                 j_ext = length(zt);
                 q_ext = min((ceil(zt(length(zt)) / hzw) - qmin + 2), length(zwc)) - 1;
            end
            
            % If no enter and no exit were found
            % So, there is no cutting process on this step
            if((enter == 0) && (ext == 0))
                s = 0;
                h_z = 0;
                h_r = 0;
            else
                % But if both enter and exit point have been found,
                % we have to calculate s and h.
                % But firstly, write the tool points in the cutting zone into an array
                ztc = zeros(1, j_ext - j_ent + 2);
                rtc = zeros(1, j_ext - j_ent + 2);
                ztc(1) = z_ent;
                rtc(1) = r_ent;
                ztc(length(ztc)) = z_ext;
                rtc(length(rtc)) = r_ext;
                i = 2;
                j = j_ent + 1;

                while((j_ent ~= j_ext) && (i < length(ztc)))
                    ztc(i) = zt(j);
                    rtc(i) = rt(j);
                    i = i + 1;
                    j = j + 1;
                end

                % Create arrays to store s and h components
                s = zeros(1, j_ext - j_ent + 1);
                h_z = zeros(1, length(s));
                h_r = zeros(1, length(s));

                % Cycle through the elements of the cutting edge immersed into workpiece material
                for j=1:(length(s))
                    % Here we calculate the length of the element 
                    s(j) = sqrt(((rtc(j+1) - rtc(j))^2) + ((ztc(j+1) - ztc(j))^2));

                    % Here we calculate the intersection of the tool element normal ray
                    % with the workpiece section in order to calculate cutting deep
                    q = q_ent;
                    normal = 0;
                    % In this cycle, we look for the intersection
                    % of the normal with the workpiece section
                    while((q <= q_ext) && (normal == 0))
                        % Calculate intersection point
                        [zp, rp] = geom_cross2(ztc(j), rtc(j), ztc(j+1), rtc(j+1), zwc(q), rwc(q), zwc(q+1), rwc(q+1));

                        % If the intersection point lies on the workpiece edge and 
                        % this point is above the cutting edge, then we found the intersection, 
                        % remember this and calculate the projection of the cut depth segment
                        if( (zwc(q) <= zp) && (zp <= zwc(q+1)) && (rp > ( (rtc(j) + rtc(j+1)) / 2 )) )
                            h_z(j) = zp - ((ztc(j+1) + ztc(j)) / 2);
                            h_r(j) = rp - ((rtc(j+1) + rtc(j)) / 2);
                            normal = 1;
                        else
                            % Otherwise continue iterations
                            q = q + 1;
                        end
                    end
                    % Calculate alternative cutting depths
                    % (from the element of the cutting edge to the middle section of the cutting edge)
                    % These values will be used to avoid over-estimation of forces for cutting depths
                    % which overlap other depths
                    [zp, rp] = geom_cross2(ztc(j), rtc(j), ztc(j+1), rtc(j+1), ztm(1), rtm(1), ztm(2), rtm(2));
                    h_z2 = zp - ((ztc(j+1) + ztc(j)) / 2);
                    h_r2 = rp - ((rtc(j+1) + rtc(j)) / 2);
                    if ((rp > ( (rtc(j) + rtc(j+1)) / 2 )) && (h_z2 < h_z(j)) && (h_r2 < h_r(j)))
                        h_z(j) = h_z2;
                        h_r(j) = h_r2;
                   end
                end
            end
        end
    end