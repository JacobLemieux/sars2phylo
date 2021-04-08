### Functions for regional data analysis
# February 28 2021
# lemieux@broadinstitute.org

#' read sample set and merge metadata
#'
#' Read a sample set from a combination of nextclade and gisaid metadata
#' @param nextclade_csv nextclade csv file
#' @param ss_metadata gisaid metadata
#' @return dataframe with merged values
#' @examples 
#' ss <- readSampleSet("nextclade.csv", "metadata_2021_03_29.csv.gz");
#' @export

readSampleSet <- function(nextclade_csv, ss_metadata){
  SS_meta <- read_tsv(ss_metadata)
  SS_meta <- SS_meta[,c("strain", "date", "pango_lineage")]
  names(SS_meta) <- c("seqName", "Date", "pango_lineage")
  SS_meta$Date <- ymd(SS_meta$Date)
  SS <- readNextClade("~/Dropbox/COVID/regional/data/", nextclade_csv)
  #SS$seqName <- sapply(str_split(SS$seqName, "202[01]-"), function(x) x[[1]][1])
  SS_meta$seqName <- gsub("hCoV-19/", "", SS_meta$seqName)
  SS <- left_join(SS, SS_meta, by = "seqName")
  SS$state <- substr(SS$seqName, 5,6)
  SS$Sgeno <- sapply(SS$aaSubstitutions, extractS) # parse nextclade format to extract S genotype
  SS$epiweek <- epiweek(SS$Date)
  SS <- SS %>% mutate(month = month(Date)) %>%
    mutate(year = year(Date)) %>%
    mutate(week = week(Date)) %>%
    filter(!is.na(year) & !is.na(month)) 
  SS <- convertMonth(SS)
  SS$epidate <- get_date(week = SS$week, year = SS$year)
  SS$beast <- paste(SS$seqName, "|", SS$Date, sep="")
  SS <- encodeLineage(SS)
  SS
}

#' read terra assemblies tsv
#'
#' @param terra_csv terra assemblies tsv file
#' @return dataframe with merged values
#' @examples 
#' ss <- readTerra("assemblies.tsv");
#' @export
#' 
readTerra <- function(terra_csv){
  SS_meta <- read_tsv(terra_csv)
  names(SS_meta)[c(1,12,27, 32)] <- c("seqName", "date", "pango_lineage", "purpose_of_sequencing")
  SS_meta <- SS_meta[,c("seqName", "date", "pango_lineage", "purpose_of_sequencing")]
  names(SS_meta) <- c("seqName", "Date", "pango_lineage", "purpose_of_sequencing")
  SS_meta$ctl <- ifelse(substr(SS_meta$seqName, 5, 7) == "H2O" | substr(SS_meta$seqName, 1, 3) == "H2O", TRUE, FALSE)
  SS_meta <- SS_meta %>% filter(ctl == FALSE)
  SS_meta$Date <- ymd(SS_meta$Date)
  SS_meta$state <- substr(SS_meta$seqName, 5,6)
  SS <- SS_meta %>% mutate(month = month(Date)) %>%
    mutate(year = year(Date)) %>%
    mutate(week = week(Date)) %>%
    filter(!is.na(year) & !is.na(month)) %>% 
    filter(!is.na(pango_lineage))
  SS <- convertMonth(SS)
  SS$epidate <- get_date(week = SS$week, year = SS$year)
  encodeLineage(SS)
}


#' read and parse nextclade results
#'
#' @param fpath string specifying file path of directory containing one or more nextclade csv files
#' @param exten prefix for matching nextclade csv files (in the case of multiple) or filename (if single)
#' @return dataframe with merged values
#' @examples 
#' ss <- readNextClade("path","nextclade.csv");
#' @export
#' 

readNextClade <- function(fpath, exten){
  fnames <- list.files(fpath)[grep(exten, list.files(fpath))]
  sstmp <- read_delim(paste(fpath, fnames[1], sep=""), ";")
  if(length(fnames) > 1){
    for (i in 2:length(fnames)){
      sstmp <- rbind(sstmp, read_delim(paste(fpath, fnames[i], sep=""), ";") )
    }
  }
  sstmp
}

#' extract S haplotypes from a list of mutations (e.g. called by nextclade)
#'
#' @param mutlist comma-separated vector
#' @return S haplotype
#' @export
#' 


extractS <- function(mutlist){
  splitlist <- strsplit(mutlist, ",")[[1]]
  paste(splitlist[grep("S:", splitlist)], collapse=",")
}

#' convert month to sequential month (jmonth) and week to sequential week (jweek):
#'
#' @param SS dataframe of sample set, including SS$week and SS$month as columns
#' @return dataframe augment with jweek and jmonth columns
#' @examples 
#' ss <- convertMonth(ss);
#' @export
#' 

convertMonth <- function(SS){
  SS$jmonth <- SS$month + (12*(SS$year - 2020))
  SS$jweek <- SS$week + (53*(SS$year - 2020))
  SS
}



#' relabel sequences for BEAST
#'
#' @param DNA_string_set Biostrings DNAstringset
#' @param ss_metadata sample set dataframe with a column ss$beast (with | delimiters) and a column ss$seqName
#' @return DNA_string_set with names formatted for a beast analysis
#' @examples 
#' DNAss <- relabel_DNAss(DNAss);
#' @export
#' 

relabel_DNAss <- function(DNA_string_set, ss_metadata){
  for(i in 1:length(names(DNA_string_set))){
    if(names(DNA_string_set)[i] %in% ss_metadata$seqName){
      names(DNA_string_set)[i] <- ss_metadata$beast[ss_metadata$seqName == names(DNA_string_set)[i]]
    }
  }
  DNA_string_set
}



#' construct a table of haplotypes by month
#'
#' @param SS dataframe of sample set, including SS$jmonth and SS$Sgeno as columns
#' @return table of haplotypes by month
#' @examples 
#' MA_haplo <- constructHaploTable(ss);
#' @export
#' 

constructHaploTable <- function(SS){
  ct <- as.data.frame.matrix(table(SS[,c("jmonth", "Sgeno")]))
  ct <- ct[,-c(1)]
  ct$seqs_by_month <- apply(ct, 1, sum)
  ct$month <- rownames(ct)
  
  for(i in 1:(ncol(ct) - 2)){
    ct[,i] <- ct[,i] / ct$seqs_by_month
  }
  ct <- as_tibble(ct)
  ct <- select(ct, -c(seqs_by_month, month))
  ct 
}

#' encode genotypes as columns using one-hot encoding
#'
#' @param SS dataframe of sample set
#' @param spike restrict analysis only to Spike gene
#' @param smallversion if true, limit number of genotypes to N
#' @param N limits number of genotypes encoded 
#' @return dataframe augmented with columns as one-hot encoded genotypes
#' @examples 
#' ss <- encodeGenotype(ss);
#' @export
#' 

encodeGenotype <- function(SS, spike = FALSE, smallversion=TRUE, N = 100){
  genos <- Reduce(union, str_split(SS$aaSubstitutions, ","))
  genos <- genos[!is.na(genos)]
  if(spike == TRUE){
    genos <- genos[grep("S:", genos)]
  }
  for(i in 1:length(genos)){
    new_col_index <- ncol(SS) + 1
    SS[,new_col_index] <- 0
    SS[grep(genos[i], SS$aaSubstitutions),new_col_index] <- 1
    colnames(SS)[new_col_index] <- genos[i]
    if(smallversion == TRUE){
      if(i > N){
        break
      }
    }
    print(i)
  }
  SS
}



#' one-hot encoding lineage as binary outcome, as augmented columns
#'
#' @param SS dataframe of sample set
#' @return dataframe augmented with one-hot encoded vectors of each lineage
#' @examples 
#' ss <- endodeLineage(ss);
#' @export
#' 

encodeLineage <- function(SS){
  one_hot_encoding <- one_hot(as.data.table(as.factor(SS$pango_lineage)))
  names(one_hot_encoding) <- gsub("V1_", "", names(one_hot_encoding))
  cbind(SS, one_hot_encoding)
}

#' Plot lineage as successes/failures over time, and fit logistic regression as smoother
#'
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param LIN character string of particular lineage
#' @param bystate subset by state; SS must have state as a column
#' @param xlower lower bound of x range
#' @param xupper upper bound of x range, will extrapolate the logistic regression smoother
#' @return creates a plot of successes (samples belong to lineage) vs failures (not lineage) over time
#' @examples 
#' plotLineage(ss, "B.1.1.7", bystate = TRUE);
#' @export
#' 

plotLineage <- function(SS, LIN, bystate=FALSE, xlower = "2020-05-01", xupper = "2021-05-01"){
  xlims <- ymd(c(xlower, xupper))
  if(bystate == TRUE){
    alphapoint = 0.5
    alphaline = 0.8
    alphashade = 0.05
    ggplot(SS, aes(x = Date, y = .data[[LIN]], color = state)) + 
      geom_jitter(width = 0.1, height = 0.05, alpha = alphapoint) + 
      theme_bw()+ 
      scale_x_date(limits=xlims) +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") + 
      stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0) 
  }else{
    alphapoint = 0.15
    alphaline = 0.7
    alphashade = 0.3
    ggplot(SS, aes(x = Date, y = .data[[LIN]], color = NULL)) + 
      geom_jitter(width = 0.1, height = 0.05, alpha = alphapoint) + 
      theme_bw()+ 
      scale_x_date(limits=xlims) +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") + 
      stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0)
  }
}


#' Create plotting function for lineage with x bins as week, not day
#'
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param LIN character string of particular lineage
#' @param bystate subset by state; SS must have state as a column
#' @param xlower lower bound of x range
#' @param xupper upper bound of x range, will extrapolate the logistic regression smoother
#' @return creates a plot of successes (samples belong to lineage) vs failures (not lineage) over time
#' @examples 
#' plotVarWeekly(ss, "B.1.1.7", bystate = FALSE, xlower = 25, xupper = 75);
#' @export
#' 

plotVarWeekly <- function(SS, LIN, bystate=FALSE, xlower, xupper){
  xlims <- c(xlower, xupper)
  if(bystate == TRUE){
    alphapoint = 0.05
    alphaline = 0.5
    alphashade = 0.15
    ggplot(SS, aes(x = Week, y = .data[[LIN]], color = state)) + 
      geom_jitter(width = 0.1, height = 0.05, alpha = alphapoint) + 
      theme_bw()+ 
      scale_x_continuous(limits=xlims) +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") + 
      stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0)
  }else{
    alphapoint = 0.15
    alphaline = 0.7
    alphashade = 0.3
    ggplot(SS, aes(x = Week, y = .data[[LIN]], color = NULL)) + 
      geom_jitter(width = 0.1, height = 0.05, alpha = alphapoint) + 
      theme_bw()+ 
      scale_x_continuous(limits=xlims) +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") + 
      stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0)
  }
}

#' fit logistic regression and return coefficients and p values
#' @param SS sample set
#' @param START column at which the one-hot encoding of the lineages begins
#' @return list length (ncol(SS) - START) of model fits
#' @examples 
#' modFit <- logRegFit(ss, 50)
#' @export
#' 

logRegFit <- function(SS, START){
  var_model_list <- vector("list", length = ncol(SS[,START:ncol(SS)]))
  for(i in 1:length(var_model_list)){
    indexvar <- START + i - 1
    var_model_list[[i]] <- glm(unlist(SS[,indexvar]) ~ SS$Date, family = "binomial")
    print(i)
  }
  var_model_coefs <- sapply(var_model_list, function(x) coefficients(summary(x))[2,1])
  var_coef_std_err <- sapply(var_model_list, function(x) coefficients(summary(x))[2,2])
  var_model_p <- sapply(var_model_list, function(x) coefficients(summary(x))[2,4])
  var_model_df <- data.frame(lineage = names(SS)[START:ncol(SS)],coefs = var_model_coefs,
                             stderrs = var_coef_std_err,
                             pvals = var_model_p, 
                             logp = -log10(var_model_p), 
                             transmissibility = exp(var_model_coefs))
  var_model_df$logp[is.infinite(var_model_df$log)] <- 250
  var_model_df
}

#' plot proportion by time interval with 95 percent CI
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param lineage character string of particular lineage
#' @param timeInt time interval, usually weekly ("epidate")
#' @param smoother display loess smoother
#' @param logreg display logistic regression smoother
#' @return creates a plot of proportion of lineage by time interval
#' @examples 
#' plotProp(ss, "B.1.1.7", "epidate");
#' @export
#' 

plotProp <- function(SS, lineage, timeInt, smoother = FALSE, logreg = FALSE){
  l_int <- SS %>% 
    filter(Date > "2020-10-31") %>% 
    group_by(.data[[timeInt]]) %>% 
    summarise(Proportion = mean(.data[[lineage]], na.rm=T), K= sum(.data[[lineage]], na.rm=T), N = n())
  l_int$LowerCI <- BinomCI(l_int$K, l_int$N)[,2]
  l_int$UpperCI <- BinomCI(l_int$K, l_int$N)[,3]
  pp <- ggplot(l_int, aes(x = .data[[timeInt]], y = Proportion)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.5) + 
    theme_bw() + 
    ggtitle(lineage) + 
    xlab("Date")  
  if(smoother == TRUE){
    pp <- pp + geom_smooth(method = "loess", span = 2)
  }
  if(logreg == TRUE){
    alphapoint = 0.5
    alphaline = 0.5
    alphashade = 0.15
    pp <- pp +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") #+ 
    #stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0)
  }
  pp
}

#' plot absolute counts by timeInt
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param lineage character string of particular lineage
#' @param timeInt time interval, usually weekly ("epidate")
#' @param smoother display loess smoother
#' @return creates a plot of proportion of lineage by time interval
#' @examples 
#' plotCounts(ss, "B.1.1.7", "epidate");
#' @export
#' 

plotCounts <- function(SS, lineage, timeInt, smoother = FALSE){
  xlims <- ymd(c("2020-10-25", "2021-04-05"))
  l_int <- SS %>% 
    filter(Date > "2020-10-31") %>% 
    group_by(.data[[timeInt]]) %>% 
    summarise(Variant = sum(.data[[lineage]], na.rm=T), All = n())
  print(names(l_int))
  names(l_int)[2] <- lineage
  l_long <- pivot_longer(l_int, cols = c(All, lineage), names_to = "Type", values_to = "Counts")
  pp <- ggplot(l_long, aes(x = .data[[timeInt]], y = Counts, color = Type, fill = Type)) + 
    scale_x_datetime(limits=xlims) +
    geom_bar(stat = "identity") + 
    theme_bw() + 
    ggtitle(lineage) + 
    xlab("Date")
  if(smoother == TRUE){
    pp <- pp + geom_smooth(method = "loess", span = 2)
  }
  pp
}

#' plot variants of concern (B.1.1.7, B.1.351, P.1)
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param timeInt time interval, usually weekly ("epidate")
#' @param smoother display loess smoother
#' @param ggtit title for graph
#' @return creates a plot of absolute counts by lineage by week
#' @examples 
#' plotVOC(ss, "B.1.1.7", "epidate", "VOC by Week");
#' @export
#' 

plotVOC <- function(SS, timeInt, smoother = FALSE, ggtit){
  l_int <- SS %>% 
    filter(Date > "2020-10-31") %>% 
    group_by(.data[[timeInt]]) %>% 
    summarise(Lin1 = sum(.data[["B.1.1.7"]], na.rm=T),Lin2 = sum(.data[["B.1.351"]], na.rm=T), Lin3 = sum(.data[["P.1"]], na.rm=T), 
              All_Other = n() - ((sum(.data[["B.1.1.7"]], na.rm=T) + sum(.data[["B.1.351"]], na.rm=T) + sum(.data[["P.1"]], na.rm=T))))
  print(names(l_int))
  names(l_int)[2:4] <- c("B.1.1.7", "B.1.351", "P.1")
  l_long <- pivot_longer(l_int, cols = c(All_Other, B.1.1.7, B.1.351, P.1), names_to = "Type", values_to = "Counts")
  l_long$Counts[l_long$Counts == 0] <- NA
  pp <- ggplot(l_long, aes(x = .data[[timeInt]], y = Counts, color = Type, fill = Type)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    ggtitle(ggtit) + 
    xlab("Date")
  if(smoother == TRUE){
    pp <- pp + geom_smooth(method = "loess", span = 2)
  }
  pp
}

#' create count of a given lineage by given time interval
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param timeInt time interval, usually weekly ("epidate")
#' @return datafram of absolute counts by lineage by week
#' @examples 
#' countByIntrval(ss, "B.1.1.7", "epidate");
#' @export
#' 

countByInterval <- function(SS, lineage, timeInt){
  l_int <- SS %>% 
    filter(Date > "2020-10-31") %>% 
    group_by(.data[[timeInt]]) %>% 
    summarise(Variant = sum(.data[[lineage]], na.rm=T), All = n())
  l_int
}


#' plot lineages as successes by week (for accurate uncertainty in logistic regression smoother)
#' @param SS dataframe of sample set, with one-hot encoding of each lineage as a column
#' @param lineage character string specifying lineage of interest
#' @param timeInt time interval, usually weekly ("epidate")
#' @param smoother display loess smoother
#' @param logreg display logistic regression smoother
#' @return creates a plot of successes by week for a given lineage
#' @examples 
#' plotWeekly(ss, "B.1.1.7", "epidate", logreg=TRUE);
#' @export
#' 


plotWeekly <- function(SS, lineage, timeInt, smoother = FALSE, logreg = FALSE){
  pp <- ggplot(SS, aes(x = .data[[timeInt]], y = .data[[lineage]])) + 
    geom_jitter(width = 0.1, height = 0.05, alpha = 0.1) + 
    theme_bw() + 
    ggtitle(lineage)
  if(smoother == TRUE){
    pp <- pp + geom_smooth(method = "loess", span = 2)
  }
  if(logreg == TRUE){
    alphapoint = 0.5
    alphaline = 0.5
    alphashade = 0.15
    pp <- pp +
      geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") + 
      stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0)
  }
  pp
}

#' plot top scoring lineages
#' @param mn_sum summary table of multinomial regression
#' @param K the number of top lineages to plot
#' @return creates a plot (using plotProp) of top lineages by week, with CIs and best fit
#' @examples
#' plotTop(multinom_sum, 6)
#' @export

plotTop <- function(mn_sum, K){
  p_list <- vector("list", length = K)
  topK <- mn_sum$lineage[order(mn_sum$coefs, decreasing=TRUE)]
  topK <- topK[1:K]
  for(i in 1:length(topK)){
    p_list[[i]] <- plotProp(ss, topK[i], "epidate", logreg=T)
  }
  plot_grid(plotlist = p_list, ncol = K/2)
}
