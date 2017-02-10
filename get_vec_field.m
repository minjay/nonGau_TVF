function [vec_field_curl, vec_field_div, vec_field_sum,...
    mag_curl_grid, mag_div_grid, mag_sum_grid] = get_vec_field(N, N_grid, gamma,...
    A_2_1, A_2_2, A_3_1, A_3_2, A_2_1_grid, A_2_2_grid, A_3_1_grid, A_3_2_grid,...
    theta_grid)

% compute vec_field_curl and vec_field_div
vec_field_curl = zeros(N, 2);
vec_field_curl(:, 1) = A_2_1*gamma(:, 1); 
vec_field_curl(:, 2) = A_2_2*gamma(:, 1);
vec_field_div = zeros(N, 2);
vec_field_div(:, 1) = A_3_1*gamma(:, 2); 
vec_field_div(:, 2) = A_3_2*gamma(:, 2);

vec_field_curl_grid = zeros(N_grid, 2);
vec_field_curl_grid(:, 1) = A_2_1_grid*gamma(:, 1);
vec_field_curl_grid(:, 2) = A_2_2_grid*gamma(:, 1);
mag_curl_grid = reshape(sqrt(vec_field_curl_grid(:, 1).^2+vec_field_curl_grid(:, 2).^2), size(theta_grid));

vec_field_div_grid = zeros(N_grid, 2);
vec_field_div_grid(:, 1) = A_3_1_grid*gamma(:, 2);
vec_field_div_grid(:, 2) = A_3_2_grid*gamma(:, 2);
mag_div_grid = reshape(sqrt(vec_field_div_grid(:, 1).^2+vec_field_div_grid(:, 2).^2), size(theta_grid));

% compute sum
vec_field_sum = vec_field_curl+vec_field_div;
vec_field_sum_grid = vec_field_curl_grid+vec_field_div_grid;
mag_sum_grid = reshape(sqrt(vec_field_sum_grid(:, 1).^2+vec_field_sum_grid(:, 2).^2), size(theta_grid));

end