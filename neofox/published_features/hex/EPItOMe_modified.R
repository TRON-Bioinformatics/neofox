#' Get position specific weights, higher scores for the mid bases compared to tails.
#'
#' Compute the position specific scores for the amino acid bases with heavier weight for the middle bases 
#' and lower score for bases towards the edges. Scores are a function of amino acid bases vector length
#' and the magnitude of score elevation by base distance from edges, bases on the third quartile of edges
#' are simplified to magnitude of 1 beginning from 1.
#'
#' @param x A character vector of Amino Acid bases.
#' @param mag Numeric value to use for seq() function to create distant sequence of weights
#' @return A numeric vector of scores equal to the length of the input vector of Amino Acids
#' @examples
#' \dontrun{
#' get_pos_weights(x=AA.vec, mag=4)
#' }
#' @keywords internal
#' @export
get_pos_weights <- function(x, mag=4){
	posScore <- c()
	xLen <- length(x)
	xMid <- xLen/2
	if(xLen %% 2){
		#print("ODD?")
		xMidCeil <- ceiling(xMid)
		xMidScore <- xMidCeil*mag
		posScore <- seq(1,xMidScore,mag)
		posScore <- c(posScore, rev(posScore[1:length(posScore)-1]))
		xTopFloor <- floor(xLen/3)
		posScore[1:xTopFloor] <- 1:xTopFloor
		xTail <- xLen-xTopFloor+1
		posScore[xTail:xLen] <- xTopFloor:1
		
	}else{
		#print("EVEN?")
		xMidScore <- xMid*mag
		posScore <- seq(1,xMidScore,mag)
		posScore <- c(posScore, rev(posScore))
		xTopFloor <- floor(xLen/3)
		posScore[1:xTopFloor] <- 1:xTopFloor
		xTail <- xLen-xTopFloor+1
		posScore[xTail:xLen] <- xTopFloor:1
	}
	return(posScore)
}

#' Alignment of singular epitope vector against a reference matrix of epitope, rows as epitopes with names and bases in 
#' columns by position.
#'
#' Align the epitope sequence to the reference set of epitope sequences. Sequences are alignment globally
#' from end to end with matches assigned a postitive score as the position specific weight and mis-matches
#' are assigned a negative score as the negative of the position specific weight. Final alignment score is
#' the sum of all the positive and neagative scores from each base. This ensures that the alignment is matches
#' in the middle are scored higher than the alignments with mis-matches in the middle.
#'
#' @param x A character vector of Amino Acid bases.
#' @param ref A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @param mag Numeric value to use for seq() function to create distant sequence of weights
#' @return A numeric vecdtor of scores for x alignment to the each epitope in ref
#' @examples
#' \dontrun{
#' align_to_ref_epitopes(x, ref, mag=4)
#' }
#' @keywords internal
#' @export
align_to_ref_epitopes <-
  function(x,
           ref,
           mag = 4,
           aln_matrix = "BLOSUM62") {
    if (!exists(aln_matrix)) {
      #data(aln_matrix)
      data(aln_matrix, envir = environment())
    }
    aln_matrix <- get(aln_matrix)
    xPosWt <- get_pos_weights(x, mag)
    #print(xPosWt)
    xLen <- length(x)
    
    xAlnScoreVec <- c()
    
    if (is.matrix(ref)) {
      refSeqCount <- nrow(ref)
    } else if (is.character(ref)) {
      refSeqCount <- 1
    }
    else{
      stop("Unsupported object as ref. Must be either matrix or character")
    }
    
    xAlnScoreVec <- apply(ref, 1, function(refSeq){
      xAlnMat <- diag(aln_matrix[x, refSeq])
      xAlnAmp <- xPosWt * xAlnMat
      xAlnScore <- sum(xAlnAmp)
      return(xAlnScore)
    })
    names(xAlnScoreVec) <- rownames(ref)
    
    return(xAlnScoreVec)
  }

#' Read the epitope amino acid sequences from the fasta file
#'
#' Read epitope sequences from a fasta file, provided with the correct path to the file.
#' The read sequences are provided as a matrix with sequence in rows and each base of the
#' sequence as a column.
#'
#' @importFrom Biostrings readAAStringSet as.matrix
#'
#' @param iFile Path to a file containing epitope sequences in fasta format.
#' @return A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @examples
#' \dontrun{
#' read_epitope_seq(iFile=input.file)
#' }
#' @keywords internal
#' @export
read_epitope_seq <- function(iFile){
	epitopes.AASet <- Biostrings::readAAStringSet(iFile)
	epitopes.mat <- Biostrings::as.matrix(epitopes.AASet, use.names=T)
	return(epitopes.mat)
}



#' Convert the epitopes data frame to amino acid matrix with bases separated into columns
#'
#' Data frame of predict epitopes is converted to a matrix of epitope sequences by rows and their
#' bases in columns.
#'
#' @importFrom stats setNames
#' @importFrom Biostrings AAStringSet as.matrix
#'
#' @param sequence.name Name of the input protein sequence for predicting epitopes.
#' @param epitope.DF Data frame of the epitopes predicted by the IEDB server.
#' @return A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns
#' @examples
#' \dontrun{
#' epitopeDF_to_AAmat(sequence.name, epitope.DF)
#' }
#' @keywords internal
#' @export
epitopeDF_to_AAmat <- function(sequence.name, epitope.DF){
        if(length(grep("#", sequence.name))>0){
		sequence.name <- gsub("#", "|", sequence.name)
	}
	res.vec <- stats::setNames(epitope.DF$peptide, paste0(sequence.name, "#", epitope.DF$start, "#", epitope.DF$end, "#", epitope.DF$percentile_rank, "#", epitope.DF$ann_ic50))
	res.AASet <- Biostrings::AAStringSet(res.vec)
	res.mat <- Biostrings::as.matrix(res.AASet, use.names=T)
	return(res.mat)
}

#' Align set of query epitope sequences against set of reference epitope sequences
#'
#' Aligning set of query epitopes by multi-threaded parallel execution of align_to_red_epitopes()
#' function. The alignment score between all pairs of query and ref epitopes is reported as a matrix.
#'
#' @importFrom parallel makeCluster clusterExport parApply stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @param query.set A matrix of Amino Acid epitope sequences to align against the reference epitope sequences.
#' @param ref.set A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns.
#' @return A matrix of alignment scores for query epitopes against reference epitopes
#' @examples
#' \dontrun{
#' align_sets(query.set, ref.set, env=NULL, num.cores=2)
#' }
#' @keywords internal
#' @export
# paralellisation was removed
align_sets <- function(query.set, ref.set, env=NULL, num.cores=2, aln_matrix="BLOSUM62"){
        if(!exists(aln_matrix)){
                print(paste0("Getting data:", aln_matrix))
                #data(aln_matrix)
                data(aln_matrix, envir=environment())
        }
	ref.set.name <- deparse(substitute(ref.set))
	ref.set.name <- strsplit(ref.set.name, "\\[")[[1]][1]
        if(is.null(env)){
            env <- environment()
        }

	# align the neoepitope sequence against a set of viral epitopes
	alnRes <- align_to_ref_epitopes(query.set, ref=ref.set, aln_matrix=aln_matrix)

	return(alnRes)
}

#' Get the pairs of query and reference epitopes with the best alignment score or above the 
#' specified alignment score threshold
#'
#' Aligning set of query epitopes by multi-threaded parallel execution of align_to_red_epitopes()
#' function. The alignment score between all pairs of query and ref epitopes is reported as a matrix.
#'
#' @importFrom parallel makeCluster clusterExport parApply stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @param alnRes A matrix of alignment scores for query epitopes against reference epitopes.
#' @param query.set A matrix of Amino Acid epitope sequences to align against the reference epitope sequences.
#' @param ref.set A matrix of Amino Acid epitope sequences, epitopes with names in rows and bases by position in columns.
#' @param best A boolean flag to switch between selecting the pairs of epitopes with best alignment score or use the 
#' input score value as a threshold.
#' @param score A numerical value to use as alignment score threshold.
#' @return A data frame of alignment result for the selected pairs of aligned epitopes
#' @examples
#' \dontrun{
#' get_best_alignment(alnRes, query.set, ref.set, best=TRUE, score=NULL)
#' }
#' @keywords internal
#' @export
get_best_alignment <- function(alnRes, query.set, ref.set, best=TRUE, score=NULL){
        if(best==TRUE){
                bestAln <- which(alnRes==max(alnRes), arr.ind = TRUE)
        }else if(is.null(score)){
                bestAln <- which(alnRes==max(alnRes), arr.ind = TRUE)
        }else{
               bestAln <- which(alnRes>=score, arr.ind = TRUE) 
        }

        if(length(bestAln)==0){
                return(NULL)
        }

	bestAlnCount <- nrow(bestAln)
	resDF <- data.frame(query=character(), q_start=numeric(), q_end=numeric(), q_epitope_rank=numeric(), q_epitope_ic50=numeric(), q_epitope_seq=character(), ref=character(), r_start=numeric(), r_end=numeric(), r_epitope_rank=numeric(), r_epitope_ic50=numeric(), r_epitope_seq=character(), score=numeric(), stringsAsFactors=F)
	for(i in 1:bestAlnCount){
		bestAlnRowName <- rownames(alnRes)[bestAln[i, "row"]]
		bestAlnColName <- colnames(alnRes)[bestAln[i, "col"]]
		bestAlnScore <- as.numeric(alnRes[bestAln[i,1],bestAln[i,2]])
		qSeq <- paste0(query.set[bestAlnColName,], collapse="")
		rSeq <- paste0(ref.set[bestAlnRowName,], collapse="")

		bestAlnRowName.split <- strsplit(bestAlnRowName, "#")[[1]]
		refName <- bestAlnRowName.split[1]
		rStart <- as.numeric(bestAlnRowName.split[2])
		rEnd <- as.numeric(bestAlnRowName.split[3])
		rRank <- as.numeric(bestAlnRowName.split[4])
		rIC50 <- as.numeric(bestAlnRowName.split[5])

		bestAlnColName.split <- strsplit(bestAlnColName, "#")[[1]]
		queryName <- bestAlnColName.split[1]
		qStart <- as.numeric(bestAlnColName.split[2])
		qEnd <- as.numeric(bestAlnColName.split[3])
		qRank <- as.numeric(bestAlnColName.split[4])
		qIC50 <- as.numeric(bestAlnColName.split[5])

		tmpDF <- data.frame(query=queryName, q_start=qStart, q_end=qEnd, q_epitope_rank=qRank, q_epitope_ic50=qIC50, q_epitope_seq=qSeq, ref=refName, r_start=rStart, r_end=rEnd, r_epitope_rank=rRank, r_epitope_ic50=rIC50, r_epitope_seq=rSeq, score=bestAlnScore, stringsAsFactors=F)
		resDF <- rbind(resDF, tmpDF)
	}
	return(resDF)
}


# B-score/ log-likelihood score 
# this code was not working and was fixed
align_to_cpl <- function(x, ref){
  xlen <- length(x)
  reflen <- nrow(ref)
  
  cpl <- matrix(nrow=xlen, ncol=20)
  aacols <- Biostrings::AAString("ACDEFGHIKLMNPQRSTVWY")
  
  cpl <- sapply(1:20, function(j) {
    sapply(1:xlen, function(xl) {
      if (x[xl] == as.character(aacols[j])) {
        100
      } else{
        1
      }
    })
  })

  colnames(cpl) <- c("A","C","D","E","F","G","H","I","K","L",
                     "M","N","P","Q","R","S","T","V","W","Y")
  
  score_vec <- c()
  # scoring matrix
  cpl <- log(cpl / rowSums(cpl))
  
  score_vec <- apply(ref, 1, function(peptide){
    # this has been changed in comparison to the original version
    # the loop over i  was not over the length of the peptide 
    sum(sapply(1:xlen, function(l) {
      sum(score <-
            sapply(1:20, function(j) {
              if (peptide[l] == as.character(aacols[j])) {
                # if(peptide[i] == aacols[j]){
                cpl[l, j]
              } else{
                NA
              }
            }), na.rm = T)
      return(score)
    }), na.rm = T)
  })
  names(score_vec) <- rownames(ref)
  return(score_vec)
}

align_cpl_sets <- function(query.set, ref.set, env=NULL, num.cores=2, aln_matrix="BLOSUM62"){
 
  ref.set.name <- deparse(substitute(ref.set))
  ref.set.name <- strsplit(ref.set.name, "\\[")[[1]][1]
  if(is.null(env)){
    env <- environment()
  }
  cplRes <- align_to_cpl(query.set, ref.set)

  return(cplRes)
}

