# -------------- #
#### Functions ####
# -------------- #


#' Run fast-dm
#'
#' Estimate a diffusion model using fast-dm.
#' @param data Dataframe containing all relevant columns (Subject, Conditions, TIME and RESPONSE).
#' @param Subject Name of the column containing the subject id.
#' @param Conditions Name(s) of the column(s) containing the levels of conditions across which parameters might vary.
#' @param TIME Name of the column containing response times (in seconds).
#' @param RESPONSE Name of the column containing responses.
#' @param method Method of parameter estimation. Can be "ml", "ks" or "cs" for maximum likelihood, kolmogorov smirnov or chi squared.
#' @param precision Precision of calculation. VALUE corresponds roughly to the number of decimals of the predicted CDFs that are calculated accuratly.
#' @param fix_to List containing fixed parameter values. Names are parameters, elements are values.
#' @param depend_on_condition List containing parameters which may vary across conditions. Names are parameters, elements are conditions.
#' @param invariant Vector of names of parameters, which are not fixed but invariant across conditions.
#' @param delete Should temporary files be deleted?
#' @return Returns a list of the input data, individual and aggregated parameter estimates as well as individual and aggregated cdf values per condition.
#' @seealso For further information regarding fast-dm see \url{https://www.psychologie.uni-heidelberg.de/ae/meth/fast-dm/}.
#' @examples
#'# simulate data for 2 subjects
#'data = data.frame(sub = rep(c(1,2), each = 100),
#'                  cnd = rep(c(1,2), times = 100),
#'                  RESPONSE = sample (c(0,1), 200, p=c(0.1,0.9), replace = TRUE),
#'                  TIME = round((rnorm(200,400,30) + rexp(200,0.01))/1000, 2))
#'
#'# load package
#'library(FastDMinR)
#'
# run analysis
#'results <- fast_dm(data,
#'                   Subject = "sub",
#'                   Conditions = "cnd",
#'                   TIME = "TIME",
#'                   RESPONSE = "RESPONSE",
#'                   precision = 5.0,
#'                   method = "ks",
#'                   fix_to = list(p = 0, d = 0, sv = 0, st0 = 0, szr = 0),
#'                   depend_on_condition = list(a = "cnd"),
#'                   invariant = c("zr", "v", "t0"))
#'
#'# plot cdf values with ggplot2
#'library(ggplot2)
#'ggplot(results$cdf$aggr_cdf, aes(x = RT, y = CDF)) +
#'  geom_line(aes(lty = cdf_Type), lwd = 1) +
#'  facet_grid(. ~ cnd)
#'
#' @keywords fast-dm
#' @keywords diffusion modeling
#' @import tidyr
#' @import sys
#' @import stringr
#' @export

fast_dm <- function(data,
                    Subject,
                    Conditions,
                    TIME,
                    RESPONSE,
                    fast_dm_path,
                    method,
                    precision,
                    fix_to,
                    depend_on_condition,
                    invariant,
                    delete = TRUE){

  data = IAT_dataframe_S1
  Subject = "id"
  Conditions = "compatible"
  TIME = "rt"
  RESPONSE = "correct"
  method = "ml"
  precision = "5.0"
  fix_to = list("zr" = .5, "d" = 0, "szr" = 0, "sv" = 0, "st0" = 0, "p" = 0)
  depend_on_condition = list("a" = "compatible", "v" = "compatible", "t0" = "compatible")
  invariant = NULL
  delete = TRUE

  require(tidyr)
  require(sys)
  require(stringr)

  fast_dm_path <- paste0(find.package("FastDMinR"),"/fast-dm")

  start <- Sys.time()

  # checks
  if(!is.data.frame(data)){stop("data must be a data.frame!")}
  if(all(c(Subject,Conditions,TIME,RESPONSE) %in% names(data)) == FALSE  ){stop('Check column names!')}
  if(!file.exists(paste0(fast_dm_path,"/fast-dm.exe"))){stop("Could not find fast-dm.exe!")}
  if(!file.exists(paste0(fast_dm_path,"/plot-cdf.exe"))){stop("Could not find plot-cdf.exe!")}
  if(!(method %in% c("ml","ks","cs")))stop("Check method!")
  if(!setequal(c(names(fix_to), names(depend_on_condition), invariant), c("a","v","zr","d","t0","szr","sv","st0","p")))stop("All parameters (a, v, t0, d, zr, p, sv, st0 and szr) need to be specified as fixed, depending on a condition or invariant!")

  # results list
  Diffusion_results = list()

  # create clean  data.frame
  Diffusion_data <- data[,c(Subject, TIME, RESPONSE, Conditions)]
  colnames(Diffusion_data) <- c("Subject", "TIME", "RESPONSE", Conditions)

  # factor to numeric
  Diffusion_data[Conditions][,sapply(Diffusion_data[Conditions], is.factor)] <- as.numeric(Diffusion_data[Conditions][,sapply(Diffusion_data[Conditions], is.factor)])

  # ms to s
  if(max(Diffusion_data$TIME) > 100 ){Diffusion_data$TIME <- Diffusion_data$TIME / 1000
  cat("\nRTs seem to be measured in ms, converting to seconds ...\n")}


  # create tmp directory
  dir.create("Diffusion", showWarnings = FALSE)

  # Count outliers
  Diffusion_estimates <- data.frame(Subject = unique(Diffusion_data$Subject))
  Diffusion_estimates$n_outlier <- numeric(nrow(Diffusion_estimates))

  # Conditions grid
  Conditions_grid <- unique(Diffusion_data[,Conditions])
  Conditions_grid <-as.data.frame(Conditions_grid)
  names(Conditions_grid) <- Conditions

  # Write data and compute ecdf
  cat("\nPreparing data ... \n")

  for (i in unique(Diffusion_data$Subject)){

    df <- Diffusion_data[Diffusion_data$Subject == i, c("TIME","RESPONSE", Conditions)]

    NAMES <- names(df)
    NAMES[1] <- paste0("#", NAMES[1])

    outlier <- which(df$TIME %in% boxplot.stats(df$TIME)$out | df$TIME < .1 | df$TIME > 5)

    Diffusion_estimates[Diffusion_estimates$Subject == i, "n_outlier"] <- length(outlier)

    if ( length(outlier) != 0){df <- df[-outlier,]}

    write.table(df, file = paste0("Diffusion//",i,".csv"), row.names = FALSE, quote = FALSE, col.names =  NAMES)


    cdf_list <- list()

    for (j in 1:nrow(Conditions_grid)) {

      cond_comb <- paste(Conditions_grid[j,], collapse = "-")
      cond_comb_indicator <- apply(df[Conditions],1, paste0, collapse = "-") == cond_comb

      cdf_data <- data.frame(TIME = ifelse(df$RESPONSE[cond_comb_indicator] == 0,  -1* df$TIME[cond_comb_indicator], df$TIME[cond_comb_indicator]))
      obs_cdf <- data.frame(Subject = i,
                            RT = as.numeric(names(cumsum(prop.table(table(cdf_data$TIME))))),
                            CDF = cumsum(prop.table(table(cdf_data$TIME))),
                            cdf_Type = "empirical",
                            cond_comb = cond_comb, row.names = NULL)

      cdf_list[[cond_comb]] <- obs_cdf

    }

    Diffusion_results$cdf_indiv[[i]]<- do.call("rbind", cdf_list)

  }

  obs_cdf_indiv <- do.call("rbind", Diffusion_results$cdf_indiv)
  obs_cdf_indiv <- separate(obs_cdf_indiv,"cond_comb",Conditions)
  rownames(obs_cdf_indiv) <- NULL

  # aggregated analysis
  df <- Diffusion_data[, c("TIME","RESPONSE",Conditions)]
  outlier <- which(df$TIME %in% boxplot.stats(df$TIME)$out | df$TIME < .1 | df$TIME > 5)

  if ( length(outlier) != 0){df <- df[-outlier,]}

  NAMES <- names(df)
  NAMES[1] <- paste0("#", NAMES[1])

  write.table(df, file = paste0("Diffusion//all.csv"), row.names = FALSE, quote = FALSE, col.names =  NAMES)



  # delete files
  if(delete == FALSE){
    cat(paste("wrote", length(unique(Diffusion_data$Subject)) + 1, "files to", paste0(getwd(),"/Diffusion")  ))
  }


  # write ctl
  set_lines <- paste("set", names(fix_to), unlist(fix_to), sep = " ")
  depends_lines <- paste("depends", names(depend_on_condition), sapply(depend_on_condition, paste, collapse=" "), sep = " ")
  format_line <- paste("format", paste(names(Diffusion_data)[-1], collapse=" "))

  fileConn <- file(paste0(fast_dm_path,"/experiment.ctl"))
  writeLines(c(paste0("method ", method),
               paste0("precision ", precision),
               set_lines,
               depends_lines,
               format_line,
               paste0("load \"",normalizePath(getwd()), "\\Diffusion\\*.csv\""),
               paste0("log \"", normalizePath(getwd()),"\\Diffusion\\estimates.txt\"")), fileConn
  )
  close(fileConn)

  # for archiving
  fileConn <- file(paste0("Diffusion//session_",Sys.Date(),".ctl"))
  writeLines(c(paste0("method ", method),
               paste0("precision ", precision),
               set_lines,
               depends_lines,
               format_line,
               paste0("load \"", normalizePath(getwd()), "\\Diffusion\\*.csv\""),
               paste0("log \"", normalizePath(getwd()), "\\Diffusion\\estimates.txt\"")), fileConn
  )
  close(fileConn)

  # write bat
  fileConn <- file(paste0(fast_dm_path,"/fast-dm.bat"))
  writeLines(c(paste0(sub("^([[:alpha:]]*).*", "\\1", fast_dm_path),":"),
               paste0("cd ", fast_dm_path),
               "fast-dm.exe"), fileConn)
  close(fileConn)

  # launch fast-dm
  cat("\nRunning fast-dm ... \n")
  system(shQuote(paste0(fast_dm_path,"/fast-dm.bat"), type = "cmd"), wait = TRUE, intern = TRUE)


  # read results
  tmp <- read.table("Diffusion//estimates.txt", sep = "", header = TRUE)
  aggr_estimates <- tmp[tmp$dataset == "all",]
  rownames(aggr_estimates) <- NULL
  tmp <- tmp[tmp$dataset != "all",]
  # tmp$Subject <- as.numeric(as.character(tmp$dataset))
  tmp$Subject <- as.character(tmp$dataset)
  tmp$dataset <- NULL
  indiv_estimates <- merge(tmp, Diffusion_estimates)


  # launch plot
  cat("\nCalculating cdf's ... \n")

  pb = txtProgressBar(min = 0, max = length(unique(Diffusion_data$Subject)), initial = 0)
  stepi <- 1
  for (i in unique(Diffusion_data$Subject)){

    cdf_list <- list()

    for (j in 1:nrow(Conditions_grid)) {

      a <- capture.output(sys::exec_wait(cmd = paste0(fast_dm_path,"/plot-cdf.exe"),
                                         args = c(paste0("-a ", ifelse("a" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("a", paste(Conditions_grid[j, depend_on_condition$a], collapse = "_"), sep = "_")],
                                                                       ifelse("a" %in% names(fix_to), fix_to$a,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "a"]))),
                                                  paste0("-z ", ifelse("zr" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("zr", paste(Conditions_grid[j, depend_on_condition$z], collapse = "_"), sep = "_")],
                                                                       ifelse("zr" %in% names(fix_to), fix_to$zr,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "zr"]))),
                                                  paste0("-v ", ifelse("v" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("v", paste(Conditions_grid[j, depend_on_condition$v], collapse = "_"), sep = "_")],
                                                                       ifelse("v" %in% names(fix_to), fix_to$v,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "v"]))),
                                                  paste0("-t ", ifelse("t0" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("t0", paste(Conditions_grid[j, depend_on_condition$t0], collapse = "_"), sep = "_")],
                                                                       ifelse("t0" %in% names(fix_to), fix_to$t0,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "t0"]))),
                                                  paste0("-d ", ifelse("d" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("d", paste(Conditions_grid[j, depend_on_condition$d], collapse = "_"), sep = "_")],
                                                                       ifelse("d" %in% names(fix_to), fix_to$d,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "d"]))),
                                                  paste0("-Z ", ifelse("szr" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("szr", paste(Conditions_grid[j, depend_on_condition$szr], collapse = "_"), sep = "_")],
                                                                       ifelse("szr" %in% names(fix_to), fix_to$szr,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "szr"]))),
                                                  paste0("-V ", ifelse("sv" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("sv", paste(Conditions_grid[j, depend_on_condition$sv], collapse = "_"), sep = "_")],
                                                                       ifelse("sv" %in% names(fix_to), fix_to$sv,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "sv"]))),
                                                  paste0("-T ", ifelse("st0" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("st0", paste(Conditions_grid[j, depend_on_condition$st0], collapse = "_"), sep = "_")],
                                                                       ifelse("st0" %in% names(fix_to), fix_to$st0,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "st0"]))),
                                                  paste0("-p ", ifelse("p" %in% names(depend_on_condition),
                                                                       indiv_estimates[indiv_estimates$Subject == i, paste("p", paste(Conditions_grid[j, depend_on_condition$p], collapse = "_"), sep = "_")],
                                                                       ifelse("p" %in% names(fix_to), fix_to$p,
                                                                              indiv_estimates[indiv_estimates$Subject == i, "p"])))
                                         ), std_err = FALSE))
      a <- str_remove(a[-length(a)],"\r")
      cdf<- separate(data.frame(a), col = 1, sep = " ", into = c("RT","CDF"))
      cdf$RT <- as.numeric(cdf$RT)
      cdf$CDF <- as.numeric(cdf$CDF)
      cdf <- cdf[cdf$RT <= 3 & cdf$RT >= -3,]
      cdf$Subject <- i
      cdf$cond_comb <- paste(Conditions_grid[j,], collapse = "-")
      cdf$cdf_Type <- "predicted"


      cdf_list[[j]] <- cdf

    }

    Diffusion_results$cdf_indiv2[[i]]<- do.call("rbind", cdf_list)

    setTxtProgressBar(pb,stepi)
    stepi <- stepi + 1
  }


  pred_cdf_indiv <- do.call("rbind", Diffusion_results$cdf_indiv2)
  pred_cdf_indiv <- separate(pred_cdf_indiv, "cond_comb", Conditions)

  for (j in 1:nrow(Conditions_grid)) {

    a <- capture.output(sys::exec_wait(cmd = paste0(fast_dm_path,"/plot-cdf.exe"),
                                       args = c(paste0("-a ", ifelse("a" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("a", paste(Conditions_grid[j, depend_on_condition$a], collapse = "_"), sep = "_")],
                                                                     ifelse("a" %in% names(fix_to), fix_to$a,
                                                                            aggr_estimates[, "a"]))),
                                                paste0("-z ", ifelse("zr" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("zr", paste(Conditions_grid[j, depend_on_condition$z], collapse = "_"), sep = "_")],
                                                                     ifelse("zr" %in% names(fix_to), fix_to$zr,
                                                                            aggr_estimates[, "zr"]))),
                                                paste0("-v ", ifelse("v" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("v", paste(Conditions_grid[j, depend_on_condition$v], collapse = "_"), sep = "_")],
                                                                     ifelse("v" %in% names(fix_to), fix_to$v,
                                                                            aggr_estimates[, "v"]))),
                                                paste0("-t ", ifelse("t0" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("t0", paste(Conditions_grid[j, depend_on_condition$t0], collapse = "_"), sep = "_")],
                                                                     ifelse("t0" %in% names(fix_to), fix_to$t0,
                                                                            aggr_estimates[, "t0"]))),
                                                paste0("-d ", ifelse("d" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("d", paste(Conditions_grid[j, depend_on_condition$d], collapse = "_"), sep = "_")],
                                                                     ifelse("d" %in% names(fix_to), fix_to$d,
                                                                            aggr_estimates[, "d"]))),
                                                paste0("-Z ", ifelse("szr" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("szr", paste(Conditions_grid[j, depend_on_condition$szr], collapse = "_"), sep = "_")],
                                                                     ifelse("szr" %in% names(fix_to), fix_to$szr,
                                                                            aggr_estimates[, "szr"]))),
                                                paste0("-V ", ifelse("sv" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("sv", paste(Conditions_grid[j, depend_on_condition$sv], collapse = "_"), sep = "_")],
                                                                     ifelse("sv" %in% names(fix_to), fix_to$sv,
                                                                            aggr_estimates[, "sv"]))),
                                                paste0("-T ", ifelse("st0" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("st0", paste(Conditions_grid[j, depend_on_condition$st0], collapse = "_"), sep = "_")],
                                                                     ifelse("st0" %in% names(fix_to), fix_to$st0,
                                                                            aggr_estimates[, "st0"]))),
                                                paste0("-p ", ifelse("p" %in% names(depend_on_condition),
                                                                     aggr_estimates[, paste("p", paste(Conditions_grid[j, depend_on_condition$p], collapse = "_"), sep = "_")],
                                                                     ifelse("p" %in% names(fix_to), fix_to$p,
                                                                            aggr_estimates[, "p"])))
                                       ), std_err = FALSE))
    a <- str_remove(a[-length(a)],"\r")
    cdf<- separate(data.frame(a), col = 1, sep = " ", into = c("RT","CDF"))
    cdf$RT <- as.numeric(cdf$RT)
    cdf$CDF <- as.numeric(cdf$CDF)
    cdf <- cdf[cdf$RT <= 3 & cdf$RT >= -3,]
    cdf$cond_comb <- paste(Conditions_grid[j,], collapse = "-")
    cdf$cdf_Type <- "predicted"


    cdf_list[[j]] <- cdf}

  pred_cdf_aggr <- do.call("rbind", cdf_list)
  pred_cdf_aggr <- separate(pred_cdf_aggr, "cond_comb", Conditions)
  rownames(pred_cdf_aggr) <- NULL

  # - - - - -  - -  - - -  - - -  - - - -  - - - - - - - -  - - -  -


  df <- Diffusion_data[, c("TIME","RESPONSE", Conditions)]

  outlier <- which(df$TIME %in% boxplot.stats(df$TIME)$out | df$TIME < .1 | df$TIME > 5)

  if ( length(outlier) != 0){df <- df[-outlier,]}

  cdf_list <- list()


  for (j in 1:nrow(Conditions_grid)) {

    cond_comb <- paste(Conditions_grid[j,], collapse = "-")
    cond_comb_indicator <- apply(df[Conditions],1, paste0, collapse = "-") == cond_comb

    cdf_data <- data.frame(TIME = ifelse(df$RESPONSE[cond_comb_indicator] == 0,  -1* df$TIME[cond_comb_indicator], df$TIME[cond_comb_indicator]))
    obs_cdf <- data.frame(RT = as.numeric(names(cumsum(prop.table(table(cdf_data$TIME))))),
                          CDF = cumsum(prop.table(table(cdf_data$TIME))),
                          cdf_Type = "empirical",
                          cond_comb = cond_comb, row.names = NULL)

    cdf_list[[cond_comb]] <- obs_cdf}


  obs_cdf_aggr <- do.call("rbind", cdf_list)
  obs_cdf_aggr <- separate(obs_cdf_aggr, "cond_comb", Conditions)
  rownames(obs_cdf_aggr) <- NULL

  indiv_cdf <- rbind(obs_cdf_indiv, pred_cdf_indiv)
  aggr_cdf <- rbind(obs_cdf_aggr, pred_cdf_aggr)




  # delete files
  cat("\n\nDeleting temporary files ... \n")
  if(delete == TRUE){
    unlink(paste0(getwd(),"/Diffusion"), recursive = TRUE)}

  unlink(paste0(fast_dm_path,"/fast-dm.bat"), recursive = TRUE)
  unlink(paste0(fast_dm_path,"/experiment.ctl"), recursive = TRUE)

  stop <- Sys.time()

  Duration <- capture.output(stop - start)
  Duration <- str_remove(Duration, "Time difference of ")

  cat(paste0("\nDone! \nDuration: ", Duration))

  return(list(raw_data = Diffusion_data,
              indiv_estimates = indiv_estimates,
              aggr_estimates = aggr_estimates,
              cdf = list(indiv_cdf = indiv_cdf,
                         aggr_cdf = aggr_cdf)))
}




