% load('workspace-vars-10-27-2020.mat');
representation_points = Representative_points(:,1:2);
n_points = sqrt(size(representation_points, 1));
x = representation_points(end-n_points+1:end, 1);
Z = reshape(Problem_Solution, n_points, n_points);
[C_DP]=contourc(x, x, Z, [prob_thresh prob_thresh]);

% Parse the contour matrix
col_indx = 1;
poly_array_vertices = [];
while col_indx <= length(C_DP)
    % Contour matrix has a specific structure. See
    % https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.contour-properties.html#d119e148367_panel-group
    current_level = C_DP(1, col_indx);
    current_level_indx = find(abs(prob_thresh - current_level)<eps);
    no_of_points = C_DP(2, col_indx);
    poly_array_vertices = [poly_array_vertices, ...
        C_DP(:,col_indx+1:col_indx+no_of_points)];
    col_indx = col_indx + no_of_points + 1;        
end

% Construct the polyhedrons
originalPolytope = Polyhedron('V', poly_array_vertices');
originalPolytope.minVRep();
originalVertices = originalPolytope.V;
n_orig_vertices  = size(originalVertices,1);   
corners = safe_set.outerApprox.V;
if n_orig_vertices == 0 && max(Problem_Solution) < prob_thresh
    % Skip the check if there are no vertices and empty set
else
    for corner_indx = 1:size(corners,1)
        tempPolytope = Polyhedron('V',[originalVertices;
                                       corners(corner_indx,:)]);
        tempPolytope.minVRep();
        if (size(tempPolytope.V,1) == n_orig_vertices + 1) ||...
                (size(tempPolytope.V,1) == n_orig_vertices)
            poly_array_vertices = [poly_array_vertices,corners(corner_indx,:)'];
        end
    end
end
poly_faust2 = Polyhedron('V', poly_array_vertices');        

figure(1);
clf;
surf(x,x,Z);
xlabel('x');
ylabel('y');
zlabel('Value function (FAUST2)');

figure(2);
clf;
plot(poly_faust2);
xlabel('x');
ylabel('y');