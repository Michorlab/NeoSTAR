library(PureCN)

path="/michorlab/jacobg/Ellisen/purecn/results"
curated_path="/michorlab/jacobg/Ellisen/purecn/curation/"

# Create empty dataframe to store results
chromosome_instability <- data.frame(
    sample = character(),
    dominant_cin = numeric(),
    normal_cin = numeric(),
    stringsAsFactors = FALSE
)

files = list.files(path, full.names=TRUE, pattern = "*.rds")
for (f in files) {
    print(f)
    samplename <- gsub(".*/(.+)\\.rds$", "\\1", f)
    
    # Handle curated files
    if (samplename %in% c("BIDMC5_pre", "BIDMC8_pre", "MGH13_pre","MGH9_pre")){
        print("curated!!!")
        ret = readCurationFile(paste0(curated_path, samplename, ".rds"))
    } else {
        ret <- readRDS(f)
    }
    
    print(samplename)
    print(names(ret))
    
    # Calculate CIN scores for both reference states
    dominant_cin <- callCIN(ret, allele.specific = TRUE, reference.state = 'dominant')
    normal_cin <- callCIN(ret, allele.specific = TRUE, reference.state = 'normal')
    
    # Add results to dataframe
    chromosome_instability <- rbind(chromosome_instability, 
                                  data.frame(
                                      sample = samplename,
                                      dominant_cin = dominant_cin,
                                      normal_cin = normal_cin
                                  ))
}

# Write results to CSV
write.csv(chromosome_instability, 
          file = "chromosome_instability.csv", 
          row.names = FALSE)