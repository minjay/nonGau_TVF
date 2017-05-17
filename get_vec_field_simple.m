function [vec_field_curl, vec_field_div, vec_field_sum] = get_vec_field_simple(N, gamma, A_2_1, A_2_2, A_3_1, A_3_2)

% compute vec_field_curl and vec_field_div
vec_field_curl = zeros(N, 2);
vec_field_curl(:, 1) = A_2_1*gamma(:, 1); 
vec_field_curl(:, 2) = A_2_2*gamma(:, 1);
vec_field_div = zeros(N, 2);
vec_field_div(:, 1) = A_3_1*gamma(:, 2); 
vec_field_div(:, 2) = A_3_2*gamma(:, 2);

% compute sum
vec_field_sum = vec_field_curl+vec_field_div;

end