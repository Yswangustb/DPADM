load_data <- function(file_path) {
	csv_data = read.csv(file_path, 
		skip=1,              # ignore first row
		header=FALSE,        # do not request the header
		encoding="utf-8",    # encoding type
		sep='\t');           # split character
	# delete header
	csv_data_rows = nrow(csv_data); csv_data_cols = length(csv_data);
	csv_data = csv_data[1: csv_data_rows, 2: csv_data_cols];

	# record the row between 10 and 15
	sample_data = data.matrix(csv_data[c(11: 15),]);
	# delete the row between 10 and 15
	csv_data = data.matrix(csv_data[-c(11: 15),]);
	return(list(csv_data, sample_data));
};


# calculate the variance of two vector
calculate_variance <- function(vector_1, vector_2) {
	if (length(vector_1) != length(vector_2)) {
		return(0);
	};

	vector_1.mean = mean(vector_1); vector_2.mean = mean(vector_2);
	vector_1 = vector_1 - vector_1.mean;
	vector_2 = vector_2 - vector_2.mean;
	sigma = sum((vector_1 * vector_2)) / (length(vector_2) - 1);
	return(sigma);
};

# single line as a vector, calculate the variance
calculate_covariance <- function(matrix_data) {
	rows = nrow(matrix_data);
	result = matrix(data=0, nrow=rows, ncol=rows);
	for (row_index in 1: rows) {
		for (col_index in 1: rows) {
			result[row_index,col_index] =
				calculate_variance(matrix_data[row_index,],matrix_data[col_index,]);
		};
	};
	return(result);
};


calculate_p_value <- function(csv_data, sample_data) {
	csv_data.rows_number = nrow(csv_data); csv_data.cols_number = length(data);
	p_value = matrix(data=0, nrow=csv_data.rows_number, ncol=csv_data.cols_number);

	for (row_index in c(1: csv_data.rows_number)) {
		row_data = csv_data[row_index, ];
		stack_sample_data = rbind(sample_data, row_data);
		stack_sample_data.cov = calculate_covariance(stack_sample_data);
		stack_sample_data.row_mean = apply(stack_sample_data, 1, mean);
		sigma = solve(stack_sample_data.cov);
		sigma_0 = sigma;

		# set the last row and the last col to 0
		sigma_0.rows = nrow(sigma_0);
		last_value = sigma_0[sigma_0.rows, sigma_0.rows]
		sigma_0[sigma_0.rows,] = 0; sigma_0[,sigma_0.rows] = 0;
		sigma_0[sigma_0.rows, sigma_0.rows] = last_value;
		
		# hyposis
		rows_number = nrow(stack_sample_data);
		data_number = length(stack_sample_data) / nrow(stack_sample_data);
		accept_sum_rate = 0; refuse_sum_rate = 0;
		for (col_index in c(1: data_number)) {
			col_data = stack_sample_data[,col_index];
			accept_rate = t((col_data - stack_sample_data.row_mean)) %*% 
						  solve(sigma) %*%
						  (col_data - stack_sample_data.row_mean)
			refuse_rate = t((col_data - stack_sample_data.row_mean)) %*% 
						  solve(sigma_0) %*%
						  (col_data - stack_sample_data.row_mean)
			accept_sum_rate = accept_sum_rate + sum(accept_rate);
			refuse_sum_rate = refuse_sum_rate + sum(refuse_rate);
		}


		accept = ((data_number * log(((2*pi)^(rows_number/2))*det(sigma)))/(-2)) - accept_sum_rate;
		refuse = ((data_number * log(((2*pi)^(rows_number/2))*det(sigma_0)))/(-2)) - refuse_sum_rate;
		score = matrix(data=accept-refuse,nrow=nrow(sample_data),ncol=1);

		# TODO
	
##chi-square value
		
		p_value[row_index] = pchisq(q=score, df=5, lower.tail=FALSE)
		
			# p_value[row_index] = chisq.test(abs(score))[2];
	};

	return(unlist(p_value));
};

write_to_file <- function(p_value, out_filename) {
	decrease_value = as.matrix(sort(p_value));
	result = as.matrix((decrease_value < 0.05));
	result = cbind(decrease_value, result);
	write.table(result, out_filename, sep='\t', eol='\n');
};


data <- load_data('D:/YS-Wang/思路文章/Pathway and drug sensitivity.R/cmap.origdata.txt');
csv_data <- data[[1]]; sample_data = data[[2]]; rm(data);
p_value <- calculate_p_value(csv_data, sample_data);
#write_to_result = write_to_file(p_value, "d:/dev/extras/data/result.txt");


#source('D:/YS-Wang/思路文章/Pathway and drug sensitivity.R/chi2test.R')
