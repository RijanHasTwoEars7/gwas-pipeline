#ncol in blues object
numMetabolites <- 4701

indexes = seq(1, numMetabolites, 100)

#if(numMetabolites %% 31 != 0) {
indexes = c(indexes, numMetabolites)
#}
print(indexes)

#inputs to the parallelize script will be:
#1. The full metabolite table
#2. The starting indices of each chunk
#3. the stop indices of each chunk
tablePath = "/home/lconnelly/GWASparallel/most_recent/sorg_blues_july_14_23.csv"

for(index in 2:length(indexes)) {
#  if(index == 3) {
#break
#}
  start = indexes[index-1]+1
  stop = indexes[index]


  
  chunkToRunString = paste(tablePath, start, stop, sep = ",")
  #print(chunkToRunString)
  theString = paste("bash /home/lconnelly/GWASparallel/GWASparallelSorg.sh", index-1, chunkToRunString)
  print(theString)
  #this is what actually executes the R script with the arguments provided in the above lines
  system(theString)
}



