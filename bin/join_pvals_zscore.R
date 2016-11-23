library(evd)

args <- commandArgs(trailingOnly = TRUE)

data <- as.matrix(read.table(args[1], sep="\t", header=FALSE)) 

data.pvalues <- data[,c(16,18,26,28)]
pvalues <- matrix(as.numeric(unlist(data.pvalues)),nrow=nrow(data.pvalues)) 

# read weights
weightdata <- matrix(, nrow(pvalues), 4)

for(x in 1:nrow(pvalues)){
    weightdata[x,] = c(5,1,5,1)
}

#weightdata <- rep(1, nrow(data))#read.table("zscore.weight", sep=";", header=FALSE)
rownames(weightdata)<-NULL
#refids <- as.character(weightdata$V1)

# intervals seq(0,1,by=0.1), seq(prelim.rho-0.1,prelim.rho+0.1,by=0.01)               ### first call                           ### second call
compute_rho <- function(interval, pvalues, weightdata) {  # left:0 right:1 data:tags.clusterd // left:prelim.rho-0.1 right:prelim.rho+0.1 data:tags.clusterd
    error <- 1000000000 # initialize
    optimal_rho <- interval[1]
    final.out_lines <- c()
    # for every rho
    for (rho in interval) {
        all_pvalues <- c()
        final.out_lines_temp <- c()
        # for every cluster
        for (i in 1:nrow(pvalues)) {

            p <- pvalues[i,]
           
            # fix for problem with p-values == 1 // qnrom(1) returns Inf
            p[which(p==1)] <- 0.9999999999999999 
 
            # transform to normal distribution
            p <- qnorm(p, lower.tail = FALSE)
        
            weight <- weightdata[i,]

            cols <- c(1)

            if(pvalues[i,2]!=1){
                cols <- c(cols,2)
            }
            if(pvalues[i,3]!=1){
                cols <- c(cols,3)
            }
            if(pvalues[i,4]!=1){
                cols <- c(cols,4)
            }

            p <- p[cols]
            weight <- weight[cols]
       
            # formula 2.2 from "A Note on Combning Dependent Tests of Significance" Joachim Hartung
            cp <- pnorm(sum(weight * p)/sqrt((1-rho)*sum(weight^2)+rho*(sum(weight)^2)), lower.tail = FALSE)
            if(is.nan(cp) == TRUE) { cp <- 0 } #############################  workaround for the NaN problem (usually for artifact predictions)
            all_pvalues <- c(all_pvalues,cp)

            final.out_lines_temp <- c(final.out_lines_temp,cp)
            
        } 

        h <- hist(all_pvalues, plot=FALSE, breaks=100)
        temp_err <- sum((1-h$intensities)**2)
        # check for error reduction
        if(temp_err < error) {
            error <- temp_err
            optimal_rho <- rho
            final.out_lines <- final.out_lines_temp
            output <- list(optimal_rho,final.out_lines)
        }
    } 
    return(output)
}

prelim_rho_output<-compute_rho(seq(0,1,by=0.1), pvalues, weightdata)
prelim_rho <- prelim_rho_output[[1]]

if (prelim_rho == 0) {
    prelim_rho <- 0.1
}

final_rho_output<-compute_rho(seq((prelim_rho-0.1),(prelim_rho+0.1),by=0.01), pvalues, weightdata)
final_rho <- final_rho_output[[1]]

final_out_lines <- final_rho_output[[2]]


final.matrix <- matrix(nrow = nrow(data), ncol = (ncol(data)+1)) 
final.matrix[,1:ncol(data)] <- data[,1:ncol(data)]
final.matrix[,ncol(final.matrix)] <- final_out_lines

data.table <- as.table(final.matrix)

colnames(data.table) <- c("Transcript_ID","Gene_ID","Gene_Symbol","Target_Start","Target_End","Transcript_ID","Gene_ID","Gene_Symbol","lncRNA_Start","lncRNA_End","Accessibility_E","Hybrid_E","Ratio","Free_Eneregy","LT_P","LT_FDR","LT_Pcor","LT_Pcor_P","LT_Dir_P","LT_Dir","PT_Start","PT_End","Protein_Gene_ID","Protein_Gene_Symbol","PT_Score","PT_P","PT_Pcor","PT_Pcor_P","PT_Dir_P","PT_Dir","Mechanism","Joint_P")

write.table(file=paste(args[1],".pvalues.csv", sep = ""), data.table[order(data.table[,"Joint_P"]),], sep = ",", row.names=FALSE, col.names=TRUE, quote=FALSE)
