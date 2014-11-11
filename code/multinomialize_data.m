new_data = function multinomialize_data(data_cell_array);
% Yes, I know this code is really dumb - I just want this done by tomorrow if possible
new_data = cell(length(data_cell_array), 1);

for i = 1:length(data_cell_array)

	old_array = data_cell_array{i};
	new_array = zeros(sum(old_array));

	for j = 1:length(old_array)

		while old_array(j) > 0
			new_array(j) = new_array(j) + 1;
		end % while

	end % for - array index

	new_data{i} = new_array;

end
