#' Truncated Monte Carlo Shapley
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#' @param T test dataset
#' @param tol trancation tolerance
#' @param convTol convergence tolerance
#'
#' @return Shapley value of training points
#' @export
#'
dataShapley<-function(D,A,V,T,tol=0.05,convTol=tol){
  N<-dim(D)[1]
  phi<-list()
  val<-list()
  sd<-list()
  permL<-list()
  model<-A(D)
  vTot<-V(model,T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[1]]<-rep(0.0,N)
  sd[[1]]<-rep(0.0,N)
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
    cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
    }else{
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
    }
    perm<-makePerm(N)
    model<-A(D[FALSE,])
    vNull<-V(model,T)
    v<-c()
    v[1]<-vNull
    for (j in (1:N)){
      if(abs(vTot-v[j])< perfTolerance){
        v[j+1]<-v[j]
      }else{
        model<-A(D[perm[1:j],])
        v[j+1]<-V(model,T)
      }
    }
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]*(t-1)/t+(v[2:N+1]-v[1:N])/t
    sd[[t]]<-rep(0.0,N)
    sd[[t]][perm]<-sd[[t]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm])
    val[[t]]<-v
    permL[[t]]<-perm
  }
  return(list(phi=phi,val=val,perm=permL))
}

#' Data Shapley with different truncation rule
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#' @param T test dataset
#' @param tol trancation tolerance
#' @param convTol convergence tolerance
#' @param log.file name of file to be used in `cat` calls
#' @param log.append logical parameter indicating if `cat` should append its output to log-file
#' @param rdata.name name of RData file should be used to save intermediate calculation results
#'
#' @return Shapley value of training points
#' @export
#'
dataShapleyI5<-function(D,A,V,T,tol=0.01,convTol=tol*5, log.file="", log.append=F, rdata.name="tmpShapley", .continue = TRUE){
  conv_check_step <- 100
  rdata.directory <- file.path(dirname(rdata.name), "temp_data")
  N <- dim(D)[1]
  phi <- list()
  sd <- list()
  val <- list()
  alph <- c(0.01, 0.05, 0.1)
  Z <- qnorm(alph, lower.tail = FALSE)
  m2 <- list()
  permL <- list()
  model <- A(D)
  tolMS <- tolMeanScore(model, V, T)
  vTot <- tolMS$mean
  v <- rep(0.0, N)
  vNull <- V(NULL, T)
  perfTolerance <- tol * vTot
  t <- 1
  phi[[t]] <- rep(0.0, N)
  val[[t]] <- rep(0.0, N)
  sd[[t]] <- rep(0.0, N)
  m2[[t]] <- rep(0.0, N)
  cat(rdata.directory, "\n", file = log.file, append = log.append)
  if (!dir.exists(rdata.directory)) {
    dir.create(rdata.directory, recursive = TRUE)
  } else {
    if (.continue) {
      prev_files_list <- list.files(rdata.directory)
      if (length(prev_files_list) > 1) {
        cat(
          format(Sys.time(), "%b %d %X"),
          "Try loaded data from previous calculation",
          "data files count =", length(prev_files_list), "\n", file = log.file, append = log.append
        )
        file_numbers <- as.integer(
          stringi::stri_replace_all_fixed(
            prev_files_list,
            paste0("_", basename(rdata.name), ".RData"), ""
          )
        )
        if (max(file_numbers) / conv_check_step == length(file_numbers) + 1) {
          last_rdata <- paste0(
            max(file_numbers), "_",
            basename(rdata.name), ".RData"
          )
          load(file.path(rdata.directory, last_rdata))
        } else {
          last_rdata <- paste0(
            max(file_numbers) - conv_check_step, "_",
            basename(rdata.name), ".RData"
          )
          load(file.path(rdata.directory, last_rdata))
        }
        phi[[t]] <- phiLast
        phi[((t - conv_check_step): (t - 1))] <- lapply(
          ((t - conv_check_step): (t - 1)),
          function(x) {
            rep_len(0, N)
          }
        )
        m2[[t]] <- m2Last
        permL[[t]] <- permLLast
        val[[t]] <- valLast
        cat(
          format(Sys.time(), "%b %d %X"),
          "Loaded data from previous calculation",
          "t=", t, "\n", file = log.file, append = log.append
        )
      }
    }
  }
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
      cat(
        format(
          Sys.time(), "%b %d %X"),
          't=',t,'\n', file = log.file, append = log.append
      )
    }else if(t%%100==0){
      rdata.file.name <- file.path(rdata.directory, paste0(t, '_', basename(rdata.name), '.RData'))
      sd<-m2[[t-1]]/(t-2)
      e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t-1)))
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      phiLast <- phi[[t - 1]]
      valLast <- val[[t - 1]]
      permLLast <- permL[[t - 1]]
      m2Last <- m2[[t - 1]]
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV, '\n', file = log.file, append = log.append)
      save(
        phiLast,t,N,vTot,v,valLast,permLLast,sd,perfTolerance,vNull,tolMS,
        m2Last,e,file = rdata.file.name
      )
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'Save is completed','\n', file = log.file, append = log.append)
    }
    perm<-makePerm(N)
    v<-rep(0.0,N)
    newRes<-vNull
    belowIdx<-0
    for (j in (1:N)){
      oldRes<-newRes
      model<-A(D[perm[1:j],])
      if (is.null(model)) {
        newRes <- vNull
        cat(format(Sys.time(), "%b %d %X"), t, "Model is null, j =", j, "\n", file = log.file, append = log.append)
      } else {
        newRes<-V(model,T)
      }
      if(abs(vTot-newRes)< perfTolerance){
        belowIdx<-belowIdx+1
      }else{
        belowIdx<-0
      }
      if(belowIdx>5){
        v[j:N]<-0
        cat(format(Sys.time(), "%b %d %X"),t,"Tolerance break:",j, length(phi),"\n", file = log.file, append = log.append)
        break()
      }
      v[j]<-newRes-oldRes
    }
    val[[t]]<-rep(0.0,N)
    val[[t]][perm]<-v
    permL[[t]]<-perm
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]+(v-phi[[t-1]][perm])/t
    m2[[t]]<-rep(0.0,N)
    m2[[t]][perm]<-m2[[t-1]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm])
  }
  sd<-m2[[t]]/(t-1)
  e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t)))
  tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
  cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n', file = log.file, append = log.append)
  save(phi,t,N,vTot,v,val,permL,perfTolerance,vNull,tolMS,m2,e,file = paste0(rdata.name, '.RData'))
  return(list(phi=phi,
              val=val,
              perm=permL,
              sd=sd,
              err01=e[,1],
              err05=e[,2],
              err10=e[,3]))
}



#' Data Shapley with different truncation rule (multithreading)
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#' @param T test dataset
#' @param tol truncation tolerance
#' @param convTol convergence tolerance
#' @param log.file name of file to be used in `cat` calls
#' @param log.append logical parameter indicating if `cat` should append its output to log-file
#' @param rdata.name name of RData file should be used to save intermediate calculation results
#' @param cluster.size number of CPU cores to be used for multiparallel computing
#' @param conv_check_step convergence is to be computed and checked every this number of iterations
#' @param base.seed Seed value to be used as basic value in different workers
#'
#' @return Shapley value of training points
#' @export
#'
dataShapleyI5.MT <- function(
  D, A, V, T, tol = 0.01, conv_tol = tol * 5, log.file = "",
  log.append = FALSE, rdata_name = "tmpShapleyML", cluster_size = 4,
  conv_check_step = 100, base_seed = as.numeric(Sys.time()),
  .continue = TRUE, .packages = c()
) {
  library(foreach)
  library(doParallel)

  cl <- makeCluster(cluster_size, outfile = log.file)
  registerDoParallel(cl)
  N <- dim(D)[1]
  phi <- list()
  sd <- list()
  val <- list()
  alph <- c(0.01, 0.05, 0.1)
  Z <- qnorm(alph, lower.tail = FALSE)
  m2 <- list()
  permL <- list()
  model <- A(D)
  cat("default model calc\n")
  tolMS <- tolMeanScore(model, V, T)
  vTot <- tolMS$mean # V(D,model,T)
  v <- rep(0.0, N)
  vNull <- V(NULL, T)
  perfTolerance <- tol * vTot
  t <- 1

  rdata.directory <- file.path(dirname(rdata_name), "temp_data")
  cat(rdata.directory, "\n")
  if (!dir.exists(rdata.directory)) {
    dir.create(rdata.directory, recursive = TRUE)
  } else {
    if (.continue) {
      prev_files_list <- list.files(rdata.directory)
      if (length(prev_files_list) > 1) {
        file_numbers <- as.integer(stringi::stri_replace_all_fixed(prev_files_list, paste0("_", basename(rdata_name), ".RData"), ""))
        if (max(file_numbers) / conv_check_step == length(file_numbers)) {
          last_rdata <- paste0(max(file_numbers), "_", basename(rdata_name), ".RData")
          load(file.path(rdata.directory, last_rdata))
          t <- ind_to_save
        } else {
          last_rdata <- paste0(max(file_numbers) - conv_check_step, "_", basename(rdata_name), ".RData")
          load(file.path(rdata.directory, last_rdata))
          t <- ind_to_save - conv_check_step
        }
      }
    }
  }
  while (!convCriteria(phi, conv_tol)) {
    t <- t + conv_check_step
    if (t <= 101 + conv_check_step) {
      cat(format(Sys.time(), "%b %d %X"), "t=", t, "\n", file = log.file, append = log.append)
    } else if ((t - conv_check_step - 1) %% 100 == 0) {
      ind_to_save <- t - conv_check_step
      rdata.file.name <- file.path(rdata.directory, paste0(ind_to_save, "_", basename(rdata_name), ".RData"))
      sd <- m2[[conv_check_step]] / (conv_check_step - 1)
      e <- sapply(Z, function(.x) sqrt((.x^2 * sd) / conv_check_step))
      tolV <- sum(abs(phi[[conv_check_step]] - phi[[conv_check_step - 100]]) / (1e-5 + abs(phi[[conv_check_step]])))
      cat(format(Sys.time(), "%b %d %X"), "ind_to_save =", ind_to_save, "tol=", tolV, "\n", file = log.file, append = log.append)
      save(phi, ind_to_save, N, vTot, v, val, permL, sd, perfTolerance, vNull, tolMS, m2, e, file = rdata.file.name)
      cat(format(Sys.time(), "%b %d %X"), "ind_to_save =", ind_to_save, "Save is completed", "\n", file = log.file, append = log.append)
    }
    perm_lists <- lapply(1:conv_check_step, function(x) makePerm(N))
    resV <- foreach(i = 1:conv_check_step, .combine = combResults, .init = list(val = val, permL = permL), .packages = .packages) %dopar% {
      set.seed(base_seed + i + t - conv_check_step)
      perm <- perm_lists[[i]]
      newRes <- vNull
      belowIdx <- 0
      v <- rep(0.0, N)

      for (j in (1:N)) {
        oldRes <- newRes
        model <- A(D[perm[1:j], ])
        if (is.null(model)) {
          newRes <- vNull
        } else {
          newRes <- V(model, T)
        }
        if (abs(vTot - newRes) < perfTolerance) {
          belowIdx <- belowIdx + 1
        } else {
          belowIdx <- 0
        }
        if (belowIdx > 5) {
          v[j:N] <- 0
          cat(format(Sys.time(), "%b %d %X"), "Worker #", i + t - conv_check_step, "Tolerance break:", j, "\n")
          break()
        }
        v[j] <- newRes - oldRes
      }
      list(i = i, perm = perm, v = v)
    }
    phi_old <- phi
    m2_old <- m2
    phi <- list()
    sd <- list()
    val <- list()
    m2 <- list()
    permL <- list()
    phi[[1]] <- rep(0.0, N)
    val[[1]] <- rep(0.0, N)
    sd[[1]] <- rep(0.0, N)
    m2[[1]] <- rep(0.0, N)

    if (t - conv_check_step > 1) {
      perm <- resV$permL[[1]]
      v <- resV$val[[1]]
      val[[1]][perm] <- v
      phi[[1]][perm] <- phi_old[[conv_check_step]][perm] + (v - phi_old[[conv_check_step]][perm]) / (t - conv_check_step)
      m2[[1]][perm] <- m2_old[[conv_check_step]][perm] + (v - phi_old[[conv_check_step]][perm]) * (v - phi[[1]][perm])
      permL[[1]] <- perm
    } else {
      permL[[1]] <- rep(0.0, N)
    }
    for (i in 2:conv_check_step) {
      perm <- resV$permL[[i]]
      v <- resV$val[[i]]
      val[[i]] <- rep(0.0, N)
      phi[[i]] <- rep(0.0, N)
      m2[[i]] <- rep(0.0, N)
      val[[i]][perm] <- v
      phi[[i]][perm] <- phi[[i - 1]][perm] + (v - phi[[i - 1]][perm]) / (t - conv_check_step + i - 1)
      m2[[i]][perm] <- m2[[i - 1]][perm] + (v - phi[[i - 1]][perm]) * (v - phi[[i]][perm])
      permL[[i]] <- perm
      if (convCriteria(phi,conv_tol)) {
        cat(format(Sys.time(), "%b %d %X"),
            "Convergency criteria has been met at",
            t - conv_check_step + i,
            "Phi length:",
            length(phi), "\n", file = log.file, append = log.append)
        break()
      }
    }
  }
  stopCluster(cl)
  phi_count <- length(phi)
  sd <- m2[[phi_count]] / (t - conv_check_step + phi_count - 1)
  e <- sapply(Z, function(.x) sqrt((.x^2 * sd) / (phi_count)))
  tolV <- sum(abs(phi[[phi_count]] - phi[[phi_count - 100]]) / (1e-5 + abs(phi[[phi_count]])))
  cat(format(Sys.time(), "%b %d %X"), "t =", t, "tol =", tolV, "\n", file = log.file, append = log.append)
  save(phi, t, N, vTot, v, val, permL, perfTolerance, vNull, tolMS, m2, e, file = paste0(rdata_name, ".RData"))
  return(list(
    phi = phi,
    val = val,
    perm = permL,
    sd = sd,
    err01 = e[, 1],
    err05 = e[, 2],
    err10 = e[, 3],
    temp_rdata_dir = rdata.directory
  ))
}
