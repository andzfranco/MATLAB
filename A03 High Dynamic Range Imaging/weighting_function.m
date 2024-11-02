function w = weighting_function(z)
    %disp(z);    
    z=double(z)/255;
    %disp(z);
    z_min = 0;
    z_max = 1;
   
    %w = double(z <= (z_min + z_max) / 2) * (z - z_min) + double(z > (z_min + z_max) / 2) * (z_max - z);
    if z <= 1/2 * (z_min + z_max)
        w = z;
    else
        w = z_max - z;
    end
    %disp(w);
end