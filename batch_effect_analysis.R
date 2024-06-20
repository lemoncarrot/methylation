

GSE72308 <- read.csv("GSE72308_beta.csv")

#create a 200 column by x row matrix
#fill every other column with beta values and the column directly after with how many elements have been processed
#algorithm for filling the matrix
#for every row
#for every column
#subset the input df for the age of interest, and extract indices
#fill in value of matrix with average of all indices, then input the length of the indices to the next column
#for filling in, if next col value is 0, proceed with normal averaging
#if next col value is >0, multiply currect average beta with next col value, add, then average accordingly
#update next col value

