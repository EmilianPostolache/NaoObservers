function diff = angdiff( theta_d, theta )
    % how to properly compute orientation errors

    % build the corresponding unit vectors

    n_d=[cos(theta_d) sin(theta_d)];
    n=[cos(theta) sin(theta)];

    % the scalar product of the two vectors gives the cosine of the
    % angle alpha between them

    alpha_abs=acos(n_d*n');

    % the sign of alpha can be derived from the cross product n x n_d
    % if its z component is positive, then sign is +; else, it is -

    % add zero third components to define the vectors in R^3
    n3=[n(1) n(2) 0];
    n_d3=[n_d(1) n_d(2) 0];

    n_c3=cross(n3,n_d3);

    % now look at the third component of n_c3
    if n_c3(3)>0
        alpha=alpha_abs;
    else
        alpha=-alpha_abs;
    end

    % here is the final orientation error
    diff= alpha;
end

