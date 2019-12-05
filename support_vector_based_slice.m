function polytope_2D = support_vector_based_slice(polytope, ...
    direction_numbers, slice_values)

    directions = [cos(linspace(0, 2*pi, direction_numbers + 1));
                  sin(linspace(0, 2*pi, direction_numbers + 1))];
    n_vertices = size(polytope.V,1);
    polytope_vertices_2D = zeros(direction_numbers, 2);
    
    for dir_indx = 1:direction_numbers    % Skip 2*pi at the end
        fprintf('Analyzing direction: %d/%d\n', dir_indx, direction_numbers)
        cvx_begin quiet
            variable cvx_comb(n_vertices, 1) nonnegative;
            variable support_vector(polytope.Dim, 1);
            
            minimize (directions(:, dir_indx)'*support_vector(1:2))
            
            subject to
                sum(cvx_comb) == 1;
                support_vector == polytope.V'*cvx_comb;
                support_vector(3:end) == slice_values;
        cvx_end
        switch cvx_status
            case {'Solved', 'Solved/Inaccurate'}
                polytope_vertices_2D(dir_indx,:) = support_vector(1:2);
            otherwise
                error('CVX failed');
        end
    end    
    polytope_2D = Polyhedron('V', polytope_vertices_2D);
end