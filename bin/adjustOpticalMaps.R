c# Lee Mendelowitz
# LMendelo@umiacs.umd.edu
# File: adjustOpticalMaps.R
# 
# Apply loess regression to remove bias in reported optical map fragment sizes
###############################################################################

############################################
# Convert a map to a string for output
mapToStr <- function(map)
{
  sout = paste(map$mapId, map$L, map$numFrags, sep="\t")
  fragStr = paste(map$frags, collapse="\t")
  sout = paste(sout, fragStr, sep="\t")
  return(sout)
}

############################################
# Convert a string (line) into a map
strToMap <- function(line)
{
  fields = strsplit(line, "\t")[[1]]
  n = length(fields)
  mapId = fields[1]
  L = fields[2]
  numFrags = fields[3]
  frags = as.numeric(fields[4:n])
  bp = sum(frags)
  stopifnot(length(frags)==numFrags)
  stopifnot(L==bp)
  map = list(mapId=mapId, L=L, numFrags=numFrags,frags=frags)
  class(map) <- "SOMAMap"
  return(map)
}

############################################
# Write the list of maps to the file handle
writeMaps <- function(maps, fname)
{
  handle = file(description=fname, open="w")
  mapStrs = lapply(maps, mapToStr)
  writeMap <- function(myStr) {write(myStr, file = handle, append=TRUE, sep="\n")}
  for(mapStr in mapStrs)
  {
    writeMap(mapStr)
  }
  close(handle)
}

############################################
# Read maps from a file, return a list of maps
readMaps <- function(fname)
{
  lines = readFile(fname)
  maps = lapply(lines, strToMap)
}

############################################
# Read lines from a file
readFile <- function(fname)
{
    lines <- scan(file=fname, what = character(), sep="\n", quiet=TRUE)
    return(lines)
}

############################################
runLoess <- function(alignedChunkFile)
{
  t <- read.table(alignedChunkFile, header=TRUE, sep="\t")
  fit <- loess(t$contig ~ t$optical)
  return(fit)
}

###########################################
# Adjust the map fragments by applying new frags
adjustMap <- function(map, fit)
{
  newMap = map
  newMap$frags <- predict(fit, map$frags)
  # If the fit gave NA, use the original value
  naInd = which(is.na(newMap$frags))
  if(length(naInd) > 0)
  {
    cat(sprintf("\n\nWARNING: Loess adjustment gave NA for these fragments in map %s:\n", map$mapId))
    origVals = map$frags[naInd]
    ind2Str <- function(ind) {s = sprintf('[%d]=%d',ind,map$frags[ind])}
    cat(paste(lapply(naInd, ind2Str)))
    cat ("\nKeeping these fragments as original values.\n\n")
    newMap$frags[naInd] = map$frags[naInd]    
  }
  newMap$frags = round(newMap$frags)
  newMap$L = sum(newMap$frags)
  return(newMap)
}

############################################
run <- function(mapFileIn, alignedChunkFileIn, adjustMapFileOut)
{
  maps <- readMaps(mapFileIn)
  cat(sprintf("Read %d maps from file %s.\n", length(maps), mapFileIn))
  
  fit = runLoess(alignedChunkFileIn)
  cat(sprintf("Ran loess fit using aligned chunk file %s\n", alignedChunkFileIn))
  
  newMaps = lapply(maps, function(map){adjustMap(map,fit)} )
  cat("Done adjusting maps.\n")
  
  writeMaps(newMaps, adjustMapFileOut)
  cat(sprintf("Wrote adjusted maps to output file %s.\n", adjustMapFileOut))
}

############################################
main <- function()
{
  args = parseArgs()
  mapFileInput = args[1]
  matchedChunksFile = args[2]
  mapFileOutput = args[3]
  run(mapFileInput, matchedChunksFile, mapFileOutput)
}

############################################
printUsage <- function()
{
  msg = paste("To run this script from commandline:\n\n",
          "R --quiet --slave --args MAP_FILE_INPUT MATCHED_CHUNKS_FILE MAP_FILE_OUTPUT < adjustOpticalMaps.R\n\n\n",
          "DESCRIPTION:\n",
          "This will adjust the input maps specified in MAP_FILE_INPUT by performing loess\n",
          "regression on the matched chunks in MATCHED_CHUNKS_FILE.\n",
          "The adjusted maps will be written to MAP_FILE_OUTPUT\n")
 cat(msg)  
}

############################################
parseArgs <- function()
{
  args = commandArgs(trailingOnly=TRUE)
  
  if (length(args) < 3)
  {
    cat('ERROR: Too few arguments!\n')
    printUsage()
    q()
  }
  
  if (length(args) > 3)
  {
    error('ERROR: Too many arguments!\n')
    printUsage()
    q()
  }
  
  mapFileInput = args[1]
  matchedChunksFile = args[2]
  mapFileOutput = args[3]
  
  msg = sprintf('MapFileInput: %s\nMatchedChunkInput: %s\nMapFileOutput: %s\n',
                mapFileInput, matchedChunksFile, mapFileOutput)
  cat(msg)
  return(c(mapFileInput, matchedChunksFile, mapFileOutput))
}

############################################n

main()
