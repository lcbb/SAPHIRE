%%%%%%%%%%%%%%%%%%%%
% Orders the states from lowest to largest distnace from origin
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function [states_ordered, params_ordered] = order_states(states,params)

temp = find(sum(params.mu_emit,1)~=0);
mag = sum(params.mu_emit(:,temp).*params.mu_emit(:,temp),1);
[~,idx] = sort(mag);
final_order = temp(idx);

% Re-order the states
if iscell(states)
    states_ordered = cell(size(states));
    for i=1:length(states)
        states_ordered{i} = states{i} + max(states{i});
        for j=1:max(states{i})
            states_ordered{i}(states_ordered{i}==final_order(j)+max(states{i})) = j;
        end
    end
else
    states_ordered = states + max(states);
    for i=1:max(states)
        states_ordered(states_ordered==final_order(i)+max(states)) = i;
    end
end

% Re-order the parameters
params_ordered = params;
params_ordered.p_start = params.p_start(final_order);
params_ordered.p_trans = params.p_trans(final_order,final_order);
params_ordered.mu_emit = params.mu_emit(:,final_order);
params_ordered.sigma_emit = params.sigma_emit(final_order);

end
