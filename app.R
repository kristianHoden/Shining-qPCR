library(shiny)
library(dplyr) #for %>%
library(ggplot2)
library(ggsignif)
library(multcompView)
library(agricolae)
library(shinyjqui)
library(shinydashboard)



# Define UI for the application 
ui <- shinyUI(dashboardPage(
  
  # Application title
  dashboardHeader(title = "Shining qPCR"),  
  
  # Application sidebar
  dashboardSidebar(
    # Application sidebar menu
    sidebarMenu(id = "tabs",
                menuItem("Settings", tabName = "Settings", icon = icon("th")),
                menuItem("Evaluate", tabName = "Evaluate", icon = icon("dashboard")),
                menuItem("Results", tabName = "Results", icon = icon("bar-chart-o"))
    )
  ),
  
  # Application body
  dashboardBody(
    shinyjs::useShinyjs(),
    # Backgroud color
    tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
      }
    '
    ))),
    
    # Application sidebar menu
    tabItems(
      
      # Settings tab content
      tabItem(tabName = "Settings",
              sidebarPanel(
                fileInput('file1', 'Choose excel File',
                          accept=c('xlsx', 'xlsx', '.xlsx')),
                
              ),
              
              sidebarPanel(
                tabsetPanel(id = "set1",
                            tabPanel("Sheet", uiOutput("sheet")),
                            tabPanel("Separator",
                                     textInput("separator", "Sample separator", "-"),
                                     verbatimTextOutput("sampSep"))
                )
              ),
              
              sidebarPanel(
                tabsetPanel(id = "set2",
                            tabPanel("Reference", uiOutput("hk")),
                            tabPanel("Control", uiOutput("ctrl"))
                )
              ),
              
              sidebarPanel(
                tabsetPanel(id = "set3",
                            tabPanel("Samples", uiOutput("checkboxSamp")),
                            tabPanel("Targets", uiOutput("checkboxTarg"),
                                     actionButton("buttonDdct", strong("ddCt"))
                            )
                )
              ),
              
              sidebarPanel(
                tabsetPanel(id = "set4",
                            tabPanel("Others", 
                                     checkboxInput('color', 'Color', TRUE),
                                     checkboxInput('anova', 'Anova', T),
                                     checkboxInput('sem', 'SEM', TRUE),
                                     #checkboxInput('remBadMelt', 'Remove bad melt', F),
                                     checkboxInput('multi', 'MultiPlot', F),
                                     checkboxInput('skipLev', 'Skip Levene', TRUE)),
                            tabPanel("Export", downloadButton('export'))
                )
              ),
              
              sidebarPanel(
                tabsetPanel(id = "set5",
                            tabPanel("Order", uiOutput("order")),
                            tabPanel("Rename", uiOutput("rename"),
                                     actionButton("buttonRename", strong("Rename")))
                )
              ),
      ),
      
      # Evaluate tab content
      tabItem(tabName = "Evaluate",
              fluidRow(
                box(title = "Evaluate samples", tableOutput("Eval"),
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    status = "primary",),
                
                box(title = "Exclude outliers", uiOutput("outliers"),
                    actionButton("buttonExclude", strong("Exclude")),
                    status = "warning", 
                    solidHeader = TRUE,
                    collapsible = TRUE
                )
              ),
              
      ),
      
      # Results tab content
      tabItem(tabName = "Results",
              tabsetPanel(id = "main",
                          tabPanel("Plot", jqui_resizable(plotOutput("plot2",width = '400px', height = '400px'))), 
                          tabPanel("Table", tableOutput("table")),
                          tabPanel("Statistics", textOutput("statistics", container = pre))
                          #pre(id = "console")
              )
      )
    )
  )
)
)


df2processing <- function(df2,checkMelt = F,samplSep = "-"){
  # Processes the from excel imported dataframe to fit to further analysis 
  df2 <- data.frame(lapply(df2, as.character), stringsAsFactors=FALSE)
  allmisscols <- sapply(df2, function(x) all(is.na(x) | x == '' ))
  if(any(allmisscols)){
    df2 <- df2[-which(allmisscols)]
  }
  
  if(ncol(df2) == 2 ){ #&& checkMelt == F
    df2 <- na.omit(df2)
    df2new <- data.frame()
    targSam <- strsplit(df2[,1], split = samplSep)
    df2new <- cbind(Target = sapply(targSam, "[", 2),
                    Sample = paste0(sapply(targSam, "[", 1),samplSep,sapply(targSam, "[", 3)),
                    Cq = df2$Cq
    )
    df2 <- data.frame(df2new)
  } else if(ncol(df2) == 3 && checkMelt == T){
    df2 <- na.omit(df2)
    df2new <- data.frame()
    targSam <- strsplit(df2[,1], split = samplSep)
    df2new <- cbind(Target = sapply(targSam, "[", 2),
                    Sample = paste0(sapply(targSam, "[", 1),samplSep,sapply(targSam, "[", 3)),
                    Cq = df2$Cq,
                    Melt.Temperature = df2$Melt.Temperature 
    )
    df2 <- data.frame(df2new)
  }else if(ncol(df2) == 3 && checkMelt == F){
    if(colnames(df2)[3] == "Melt.Temperature"){
      df2 <- df2[,-3]
      df2 <- na.omit(df2)
      df2new <- data.frame()
      targSam <- strsplit(df2[,1], split = samplSep)
      df2new <- cbind(Target = sapply(targSam, "[", 2),
                      Sample = paste0(sapply(targSam, "[", 1),samplSep,sapply(targSam, "[", 3)),
                      Cq = df2$Cq
      )
      df2 <- data.frame(df2new)
    }
    
  }
  return(df2)
}

printP <- function(...){
  # Pastes and prints
  print(paste(...))
}

whichIsNotNaMatch <- function(arg1,arg2){
  # Checks which occurences in a match that are not NA 
  return(which(!is.na(match(as.character(arg1), as.character(arg2)))))
}

substrFromEnd <- function(x, startFromEnd,endFromEnd=0){
  # substr but counted from the end
  
  if(endFromEnd==0){
    endFromEnd <- nchar(x)
  }
  return(IRanges::reverse(substr(IRanges::reverse(x),startFromEnd,endFromEnd))) 
}
substrFromStartToXFromEnd <- function(word1,endFromEnd=0){
  # substr from start but to X char from the end
  startFromEnd = nchar(word1)
  endFromEnd = endFromEnd+1  
  if(endFromEnd==0){
    endFromEnd <- nchar(word1)
  }
  return(IRanges::reverse(substr(IRanges::reverse(word1),endFromEnd,startFromEnd))) 
}
replNrev <- function(x,fn,rp,n) {
  # Replaces the from the end nth occurrence of fn with rp of x.
  
  x <- IRanges::reverse(x)
  sel <- rep(fn, n*length(rp))
  sel[seq_along(rp)*n] <- rp
  regmatches(x, gregexpr(fn, x)) <- list(sel)
  return(IRanges::reverse(x))
}

replN <- function(x,fn,rp,n) {
  # Replaces the nth occurence of fn with rp of x.
  sel <- rep(fn, n*length(rp))
  sel[seq_along(rp)*n] <- rp
  regmatches(x, gregexpr(fn, x)) <- list(sel)
  x
}

noNAqPCR2 <- function(df3){
  # Removes al NAs
  for(x in seq_along(colnames(df3))){
    if(any(is.na(df3[,x]))){
      df3 <- df3[-which(is.na(df3[,x])),]
    }
    if(any(df3[,x] == "NaN")){
      df3 <- df3[-which(df3[,x] == "NaN" ),]
    }
  }
  return(df3)
}

noNAqPCR <- function(df2){
  # Removes al NAs
  if(any(is.na(df2$Cq))){
    df2 <- df2[-which(is.na(df2$Cq)),]
  }
  if(any(is.na(df2$Sample))){
    df2 <- df2[-which(is.na(df2$Sample)),]
  }
  if(any(is.na(df2$Target))){
    df2 <- df2[-which(is.na(df2$Target)),]
  }
  if(any(is.na(df2$TechMean))){
    df2 <- df2[-which(is.na(df2$Target)),]
  }
  if(any(is.na(as.numeric(df2$TechMean)))){
    df2 <- df2[-which(is.na(as.numeric(df2$TechMean))),]
  }
  return(df2)
}

gsubMany <- function(argMatrix1 = argMatrix$Sample ,sampLorder,sampLnames){
  #gsub for many occurences
  for(x in(1:length(sampLorder))){
    argMatrix1 <- gsub(argMatrix1,pattern = sampLorder[x], replacement = sampLnames[x], fixed = T)
  }
  return(argMatrix1)
}

gen_data <- function(toPlot,repli = 3){
  # Genereates representative samples to mean sd and number of replicates   
  colnames(toPlot) <- tolower(colnames(toPlot))
  means<-toPlot$mean
  sds <- toPlot$sd
  n.grp <- length(means)
  samplesizes <- rep(repli,n.grp)
  grps <- factor(rep(1:n.grp, samplesizes))
  namesSamples <- rownames(toPlot)[grps]
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*as.numeric(sds[i]) + as.numeric(means[i])})
  y <- do.call(rbind, dat)
  out1 <- data.frame(group = grps, y = y, names = namesSamples)
  return(out1)
}

maxLimsQpcr <- function(spp4h){
  # Defines the max limits for the qPCR plot
  max(as.numeric(spp4h$pow2posSd))+((max(as.numeric(spp4h$pow2posSd)))/5)
}

statTest <- function(chosen_data, skipLevTest = F, LevTest = NULL, decTukey = T){
  #Performs the statistical tests
  if(is.null(LevTest)){
    LevTest <- car::leveneTest(y ~ groups2, data=chosen_data)
    Print("Warning, producing approximated Levene's test")
  }
  ret <- list()
  if(LevTest$`Pr(>F)`[1] < 0.05 & skipLevTest == F){
    
    library("onewaytests")
    bfTest <- bf.test(y ~ groups2, data=chosen_data)
    if(bfTest$p.value > 0.05){
      stop('Violating Null Hyp Levene test! And failing Brown-Forsythe test \n\n') 
    }
    wTest <- welch.test(y ~ groups2, data=chosen_data, alpha = 0.05)
    if(wTest$p.value > 0.05){
      stop('Violating Null Hyp Levene test! And failing Welch test \n\n') 
    }

    res <- userfriendlyscience::oneway(chosen_data$groups2, y = chosen_data$y, posthoc = 'games-howell');
    ret[[1]] <- res
    ### Extract dataframe with post hoc test results,
    res2 <- res$intermediate$posthoc;
    ### Extract p-values and comparison 'names'
    pValues <- res2$p;
    ### Create logical vector, assuming alpha of .05
    dif3 <- pValues < .05;
    ### Assign names (row names of post hoc test dataframe)
    names(dif3) <- row.names(res2);
    ### convert this vector to the letters to compare
    ### the group means (see `?multcompView` for the
    ### references for the algorithm):
    letters1 <- multcompView::multcompLetters(dif3)
    letters2 <- letters1$Letters[order(names(letters1$Letters))]
    message2 <- "Welch test -> Games-Howell"
  }else{
    message2 <- "Tukey"
    av <- aov(y ~ groups2, data = chosen_data)
    ret[[1]] <- av
    tuk <- multcomp::glht(av, linfct = multcomp::mcp(groups2 = "Tukey"))
    ret[[2]] <- tuk
    if(decTukey == T){
      tuk.cld <- multcomp::cld(tuk,decreasing = T)   
    }else{
      tuk.cld <- multcomp::cld(tuk,decreasing = F)
    }
    letters1 <- as.factor(tuk.cld$mcletters$Letters)
    letters2 <- letters1
  }
  return(list(letters2,message2,ret))
}

substrFromStartToChar <- function(text1 = names(fas), char = "_"){
  # substr from start to the first occurrence of char 
  text2 <- text1[grep(text1,pattern = char)]
  text3 <- substr(text2,1, regexpr(text2,pattern = char,fixed= T)-1)
  text1[grep(text1,pattern = char)] <- text3
  return(text1)
}

toNumericDfIfPos <- function(df1){
  #Converts columns of a df to numeric if possible
  for(x in 1:ncol(df1)){
    if(!any(is.na(suppressWarnings(tryCatch(as.numeric(as.character(df1[,x]))))))){
      df1[,x] <- suppressWarnings(as.numeric(as.character(df1[,x])))
    }
  }
  return(df1)
}

toNumericAndRoundIfposDf <- function(df1,digs = 2){
  #Converts columns of a df to numeric and rounds the decimals if possible
  
  df1 <- toNumericDfIfPos(df1)
  df1 <- df1 %>% dplyr::mutate_if(is.numeric, round, digits=digs)
  return(df1)
}

reorderAnova<-function(inV){
  # Reorder Anova letters
  
  inV <- as.character(inV)
  collapsed <- paste(inV,sep="",collapse = "")
  u <- unique(strsplit(collapsed,"")[[1]])
  if(length(u)<2){
    return(inV)
  }
  u <- u[order(u)]
  m <- matrix(nrow=NROW(inV),ncol=length(u))
  m[]<-F
  for(i in 1:length(inV)){
    s <- strsplit(inV[i],"")[[1]]
    index <- match(s,u)
    m[i,index] <- T
  }
  for(i in 1:(length(u)-1)){
    firstColT <- match(T,m[,i])[1] 
    firstT <- match(T,rowSums(m[,i:length(u)] > 0))[1] 
    if(firstT < firstColT){
      colT <- match(T,m[firstT,i:length(u)])[1]
      colT <- colT + i - 1 
      tmp <- m[,colT]
      m[,colT] <- m[,i]
      m[,i] <- tmp
    }
  }
  res <- vector(mode = "character", length=length(inV))
  for(i in 1:length(inV)){
    l <- u[m[i,]]
    res[i] <- paste(l,sep="",collapse = "")
  }
  return(res)
}

qPCR_ddct <- function(hk = "Actin",
                      ctrl ="mock", 
                      path1 = "/home/kristianh/p/ravi/fromZhen/R.xlsx",
                      sheetndex1 = 1,
                      samplSep = "-",
                      changeStart = F, #changeStart alternatives "Pi-", "1st_"
                      sem = F,
                      lengthBio = 1, #bio2_1 bio2_2 -> bio2 = 4
                      checkMelt = T,
                      checkSampleDiffFromMean = T,
                      removeBio = "",
                      noTechnicalReps = F){
  # Reads an excel file with qPCR data and analyzes the data with the ddCt method.    
  
  df2 <- xlsx::read.xlsx(file = path1, 
                         header = TRUE,sheetName = sheetndex1) 
  
  df2 <- df2processing(df2,checkMelt = checkMelt)
  
  if(changeStart == "Pi-" ){
    df2$Sample[which(startsWith(df2$Sample, "Pi-"))] <- 
      replN(df2$Sample[which(startsWith(df2$Sample, "Pi-"))], 
            fn="-", rp=c("_"), n=2)
    df2$Sample[-which(startsWith(df2$Sample, "Pi-"))] <- 
      replN(df2$Sample[-which(startsWith(df2$Sample, "Pi-"))], 
            fn="-", rp=c("_"), n=1)
    df2$Sample[-which(startsWith(df2$Sample, "Pi-"))] <- 
      replN(df2$Sample[-which(startsWith(df2$Sample, "Pi-"))], 
            fn="_", rp=c("-"), n=2)
  }
  
  if(changeStart == "1st_" ){
    df2$Sample <-  replNrev(x = df2$Sample, 
                            fn="_", rp=c(samplSep), n=2)
  }
  
  if (checkMelt == T){
    
    if(ncol(df2) > 3){
      if(colnames(df2[4]) != "Melt.Temperature"){
        stop("4th col is not called Melt. Temperature")
      }
      if(any(is.na(df2[4]))){
        remNone <- as.integer(rownames(plyr::match_df(df2[which(is.na(df2[4])),c("Target","Sample")],df2[which(!is.na(df2[4])),c("Target","Sample")])))
        if(any(remNone)) {
          df2 <- df2[-remNone,]
        }
        
      }
      
      if(any(df2[4] == "None")){
        remNone <- as.integer(rownames(plyr::match_df(df2[which(df2[4] == "None"),c("Target","Sample")],df2[which(df2[4] != "None"),c("Target","Sample")])))
        if(any(remNone)) {
          df2 <- df2[-remNone,]
        }
        
      }
      if(any(df2[4] == "")){
        remNone <- as.integer(rownames(plyr::match_df(df2[which(df2[4] == ""),c("Target","Sample")],df2[which(df2[4] != ""),c("Target","Sample")])))
        if(any(remNone)) {
          df2 <- df2[-remNone,]
        }
      }
    }
  }else{
    if(ncol(df2) > 3){
      if(colnames(df2[4]) != "Melt.Temperature"){
        stop("4th col is not called Melt Temperature")
      }
      df2 <- df2[,-4]
    }
  }
  
  ##If deleting samples with bad melting type these rows removes the corresponding hk or gene samples   
  df2hk <- df2[which(df2$Target == hk),]
  df2a <- df2[which(df2$Target != hk),]
  df2hkSampleMissing <- df2hk$Sample[which(is.na(match(df2hk$Sample, df2a$Sample)))]  
  if(length(df2hkSampleMissing)> 0){
    df2 <- df2[-whichIsNotNaMatch(df2$Sample,df2hkSampleMissing),]
  }
  df2aSampleMissing <- df2a$Sample[which(is.na(match(df2a$Sample, df2hk$Sample)))]  
  if(length(df2aSampleMissing)> 0){
    df2 <- df2[-whichIsNotNaMatch(df2$Sample,df2aSampleMissing),]
  }

  df2 <- noNAqPCR(df2)
  df2 <- df2[complete.cases(df2),]
  df2[] <- lapply(df2, as.character)
  if (any(df2$Target== "")){
    df2<- df2[-whichIsNotNaMatch(df2$Target, ""),] 
  }
  
  #Check if some row contains empty info, if so remove the row
  df2emt <- df2 == ""
  if(any(apply(df2emt, 1, any)) == T){
    df2 <- df2[-which(apply(df2emt, 1, any)),]
  }
  
  df3 <- data.frame()
  sampl <- unique(df2$Sample)
  targ <- unique(df2$Target)
  for(x in seq_along(targ)){
    for(y in seq_along(sampl)){
      sampl1 <- as.character(sampl[y])
      targ1 <- as.character(targ[x])
      targs <- which(df2$Target == targ1)
      samps <- which(df2$Sample == sampl1)
      both1 <- samps[na.omit(match(targs,samps))]
      hit <- as.numeric(df2$Cq[both1])
      sd2 <- sd(hit)
      if(is.na(sd2)) {sd2 <- 0} 
      if(!is.null(df2$Melt.Temperature)){
        if(df2$Melt.Temperature[both1][1] != "None"){
          mT <- mean(as.numeric(df2$Melt.Temperature[both1]))
        }else{
          mT <- "None"
        }
        df3 <- rbind(df3, cbind(as.character(targ1),as.character(sampl1),mean(hit),sd2,mT))
      }else{
        df3 <- rbind(df3, cbind(as.character(targ1),as.character(sampl1),mean(hit),sd2))
      }
    }
  }
  if(length(colnames(df3)) == 4){
    colnames(df3) <- c("Target","Sample", "TechMean", "TechSd")
  }else if(length(colnames(df3)) == 5){
    colnames(df3) <- c("Target","Sample", "TechMean", "TechSd", "Melt")
  }
  if(noTechnicalReps == T){
    df3 <-df3[,-which(colnames(df3) == "TechSd")]
  }else{
    if(all(is.na(df3$TechSd))){
      print("seems like no tech reps? Maybe put noTechnicalReps = T")
    }
  }
  df3 <- noNAqPCR2(df3 = df3)
  df3[] <- lapply(df3, as.character)
  if(removeBio != ""){
    for(x in seq_along(removeBio)){
      df3 <- df3[-grep(x = df3$Sample, pattern = removeBio[x]),]
      printP("Removed ",removeBio[x])
    }
  }
  #Investigate which samples deviate most from the mean
  if(checkSampleDiffFromMean == T){
    df3hk <- df3[which(df3$Target == hk),]
    df3a <- df3[which(df3$Target != hk),]
    df3a$diffRef <- NA
    for(x in seq_along(rownames(df3a))){
      hkNr <- which(df3a$Sample[x] == df3hk$Sample)
      df3a$diffRef[x] <- as.numeric(df3a$TechMean[x]) - as.numeric(df3hk$TechMean[hkNr])
    }
    x<-5
    df3a$diffMean <- NA
    for(x in seq_along(unique(df3a$Target))){
      targUniHit <- which(df3a$Target == unique(df3a$Target)[x])
      sampUni <- unique(substrFromStartToXFromEnd(word1 = df3a$Sample,endFromEnd = 3))
      for(y in seq_along(sampUni)){
        sampUniHit <- which(startsWith(df3a$Sample,sampUni[y]))
        sampUniHitMatch <- sampUniHit[na.omit(match(targUniHit,sampUniHit))]
        df3a$diffMean[sampUniHitMatch] <- df3a$diffRef[sampUniHitMatch] - mean(df3a$diffRef[sampUniHitMatch])
      }
    }  
    returnList <- list() 
    returnList[[1]] <- df3a[order(abs(df3a$diffMean), decreasing = T),]
    biol <- substrFromEnd(x = df3a$Sample,startFromEnd = 1, endFromEnd = 2)
    biolU <- unique(biol)
    biolUa <- cbind(biolU,NA)
    for(x in seq_along(biolU)){
      biolUa[x,2] <- sum(df3a$diffMean[which(biol == biolU[x])])
    }
    returnList[[2]] <- biolUa 
  }
  
  #Performs Levene's test on the cq data. Separating all samples and targets genes. 
  #If wanting to do all targets together, chage to df3stat2 
  df3stat <- cbind(toNumericDfIfPos(df3),lev = paste0(df3$Target, substrFromStartToChar(df3$Sample,samplSep)))
  #df3stat2 <- cbind(toNumericDfIfPos(df3),lev =  substrFromStartToChar(df3$Sample,samplSep))
  returnList[[4]] <- suppressWarnings(car::leveneTest(TechMean ~ lev, data=df3stat))
  
  df4 <- data.frame()
  if(any(sampl == "")){
    sampl <- sampl[-which(sampl == "")]
  }
  for(x in seq_along(targ)){
    for(y in seq_along(sampl)){
      sampl1 <- sampl[y]
      sampl2 <- substr(sampl1,1, regexpr(text = sampl1,pattern = samplSep, fixed = T)-1) 
      targ1 <- targ[x]
      targs <- which(df3$Target == targ1)
      samps <- which(startsWith(as.character(df3$Sample), prefix =sampl2))
      both1 <- samps[na.omit(match(targs,samps))]
      hit2 <- sapply(base::strsplit(as.character(df3$Sample)[both1],  split = samplSep),"[[",2)  
      check <- substr(hit2,1,lengthBio)
      check2 <- names(which(table(check) > 1))
      if(length(check2) > 0){
        for (z in(seq_along(check2))){
          check3 <- df3[both1,][whichIsNotNaMatch(check,check2[z]),]
          meanOther <- mean(as.numeric(df3[both1,][-whichIsNotNaMatch(check,check2[z]),][,"TechMean"]))
          
          keepNot <- which(abs(as.numeric(check3$TechMean) - meanOther) != min(abs(as.numeric(check3$TechMean) - meanOther)))
          both1 <- both1[-whichIsNotNaMatch(check,check2[z])[keepNot]]
        }
      }
      hit <- as.numeric(df3$TechMean[both1])
      # #computation of the standard error of the mean
      # sem<-sd(x)/sqrt(length(x))
      # #95% confidence intervals of the mean
      # c(mean(x)-2*sem,mean(x)+2*sem)
      sd1 <- sd(hit)
      mean1 <- mean(hit)
      sem1 <- sd(hit)/sqrt(length(hit))
      ci <- c(mean(hit)-2*sem1,mean(hit)+2*sem1)  
      ciD <- ci[1]
      ciU <- ci[2]
      hits <- cbind(targ1,sampl2,mean1,sd1,sem1 ,ciU,ciD, length(hit))
      if(any(is.na(hits))){
        #print(paste("df4 NAs changed to 0", colnames(is.na(hits))[is.na(hits)], "set to 0"))
        hits[is.na(hits)] <- 0
        
      }
      df4 <- rbind(df4,hits )
    }
  }
  colnames(df4) <- c("Target","Sample", "Mean", "Sd","SEM","CI_U","CI_D", "Samples")
  df4 <- dplyr::distinct(df4)
  df4[] <- lapply(df4, as.character)
  df4hk <- df4[which(df4$Target == hk),]
  df4a <- df4[which(df4$Target != hk),]
  df4a[] <- lapply(df4a, as.character)
  df4hk[] <- lapply(df4hk, as.character)
  samplB <- unique(df4a$Sample)
  targB <- unique(df4a$Target)
  df5 <- data.frame(stringsAsFactors = F) 
  
  for(x in seq_along(targB)){
    for(y in seq_along(samplB)){
      hit1 <- df4a[df4a$Target == targB[x] & df4a$Sample == samplB[y],]
      akt <- df4hk[df4hk$Sample == samplB[y],]
      df5 <- rbind(df5, cbind(targB[x],
                              samplB[y],
                              as.numeric(hit1$Mean)-as.numeric(akt$Mean),
                              sqrt((as.numeric(hit1$Sd)^2)+(as.numeric(akt$Sd)^2)),
                              sqrt((as.numeric(hit1$SEM)^2)+(as.numeric(akt$SEM)^2))))
    }
  }
  colnames(df5) <- c("Target","Sample", "dCt", "dCt_Sd","dCt_SEM")
  df5[] <- lapply(df5, as.character)
  for( k in(seq_along(df5$Sample))){
    df5$Samples[k] <- min(as.numeric(df4$Samples[which(df4$Sample == df5$Sample[k] & df4$Target == df5$Target[k])]),
                          as.numeric(df4hk$Samples[df4hk$Sample == df5$Sample[k]]))
  }
  dCtDf <- df5
  samplC <- unique(df5$Sample)
  targC <- unique(df5$Target)
  df6 <- data.frame()
  for(x in seq_along(targC)){
    work <- df5[which(df5$Target == targC[x]),]
    work$strCt <- as.numeric(work$dCt)/as.numeric(work$dCt[which(work$Sample == ctrl)])
    work$ddCt <- as.numeric(work$dCt)-as.numeric(work$dCt[which(work$Sample == ctrl)])
    work$pow2 <- 2^-as.numeric(work$ddCt)
    if(sem == T){
      work$pow2posSd <- 2^-(as.numeric(work$ddCt)-as.numeric(work$dCt_SEM))
      work$pow2negSd <- 2^-(as.numeric(work$ddCt)+as.numeric(work$dCt_SEM))
    }else{
      work$pow2posSd <- 2^-(as.numeric(work$ddCt)-as.numeric(work$dCt_Sd))
      work$pow2negSd <- 2^-(as.numeric(work$ddCt)+as.numeric(work$dCt_Sd))
    }
    work$log2 <- log(work$pow2,2) # = -work$ddCt
    df6 <- rbind(df6, work)
  }
  returnList[[3]] <- df6
  return(returnList)
}

barPlotQpcr <- function(argMatrix = qPCRmiR8788mimic2,
                        yLab = "Oomycete DNA",
                        wBars = 0.6,
                        blackWhite = F,
                        byType = F, # names with _, sets after _ to Type
                        shiftLevels = F,
                        byTypeSample = F,
                        scaleToFirst = T,
                        anova1 = T, # anova1 = F
                        cols = NA,
                        fontSize =0.3,
                        goi_ref = T,
                        maxIs0dpi = F,
                        byTickInterval = 1,
                        titL = "",
                        sampLorder = c("mock","A.S","NL"),
                        sampLnames = c("H2O","As","Pi"),
                        multPlot = F,
                        decTukey = F,
                        skipLevTest = F,
                        targetOrder = NULL, #c("StAgo10a", "StAgo1b", "StAgo4a") 
                        LevTest = NULL
                        
){
  # Computes the staistical analyses and and the plots 
  
  message2 <- ""
  argMatrix <- argMatrix[order(match(argMatrix$Sample,sampLorder)),]
  groups <- match(argMatrix$Target,unique(argMatrix$Target))
  argMatrix <- argMatrix[order(groups),]
  groups <- match(argMatrix$Target,unique(argMatrix$Target))
  argMatrix$Sample <- gsubMany(argMatrix1 = argMatrix$Sample ,sampLorder,sampLnames)
  argMatrix$origSample <- argMatrix$Sample
  
  if(length(argMatrix$Sample) > length(unique(argMatrix$Sample))){
    argMatrix$Sample <- paste(argMatrix$Sample,argMatrix$Target,sep = "_")
  }
  
  toPlot <- argMatrix[,c("log2","dCt_Sd","Sample")]
  rownames(toPlot) <- argMatrix$Sample
  colnames(toPlot) <- c("mean","sd","Sample")
  h2oC <- data.frame()
  
  for(k in 1: nrow(toPlot)){
    h2oCa <- gen_data(toPlot = toPlot[k,],repli = argMatrix$Samples[k])
    h2oCa$group <- k
    h2oC <- rbind(h2oC, cbind(h2oCa, Target = argMatrix$Target[k]))
  }
  
  if(scaleToFirst == T){
    argMatrix$pow2 <- argMatrix$pow2/ argMatrix$pow2[1]
    argMatrix$pow2posSd <- argMatrix$pow2posSd/ argMatrix$pow2[1]
    argMatrix$pow2negSd<- argMatrix$pow2negSd/ argMatrix$pow2[1]
  }
  
  h2oC$groups2 <- as.factor(h2oC[, "group"])
  library("ggplot2")
  
  #Statistics
  chosen_data <- as.data.frame(h2oC)
  x <- 3
  let2 <- lapply(1:length(unique(chosen_data$Target)), function(x){
    chosen_data1 <- chosen_data[which(chosen_data$Target == unique(chosen_data$Target)[x]),]
    let1 <- statTest(chosen_data = chosen_data1, skipLevTest = skipLevTest, LevTest = LevTest)
    return(let1)
  })
  ret <- lapply(let2,"[[",3)
  message2 <- sapply(let2,"[[",2)
  let2 <- lapply(let2,"[[",1)
  ret[["chosen_Data"]] <- chosen_data
  
  
  
  for(x in seq_along(let2)){
    let2x <- let2[[x]]
    let2[[x]] <- factor(reorderAnova(inV = let2[[x]]))
    # many <- which(nchar(as.character(let2x)) == 2)
    # if(length(many) > 0){
    #   for(y in seq_along(let2x[many]){
    #     let2x[many]
    #   }
    # }
    # let2x
    # let2[[1]]
    # order(let2[[x]])
    # levels(let2[[x]])
  }
  
  
  
  
  
  #unlist(let2)
  
  #toPlot <- as.data.frame(h2oDf)
  argMatrix <- cbind(argMatrix, Tukey = unlist(let2))# letters2)
  
  #levels(toPlot$Sample) <- levels(as.factor(toPlot$Sample))[c(1:length(toPlot$Sample))]
  
  theTable <- within(argMatrix, 
                     Sample <- factor(Sample, 
                                      levels=unique(Sample)))
  
  toPlot <- argMatrix[,c("log2","dCt_Sd", "Sample")]
  
  toPlot2 <- argMatrix[,c("pow2","pow2posSd", "Sample")]
  
  
  colnames(toPlot2) <- c("Mean", "Sd","Sample")
  
  colnames(theTable)
  #maxLim <- maxLims(spp4h = toPlot2)
  maxLim <- maxLimsQpcr(spp4h = argMatrix)
  
  #maxIs0dpi = T
  #maxIs0dpi = F
  if(maxIs0dpi == T){
    maxLim = maxLim*(1+(0.065* (nrow(argMatrix)-2)))
  }
  
  # if(byType == T){
  #   if(!all(grepl("_",theTable$Sample))){
  #     print("No _ found, if more than 2 samples, and T-test is wanted, please separate them with a _")
  #     type1 <- theTable$Sample
  #   }else{
  #     type1 <- sapply(strsplit(as.character(theTable$Sample), "_"),"[[",2)
  #   }
  #   theTable <- cbind(theTable, Type = as.factor(type1))
  #   if(shiftLevels == T){
  #     theTable$Type <- factor(theTable$Type,levels(theTable$Type)[c(2:1)])
  #     #levels() <- levels(theTable$Type)[c(2,1)]
  #   }
  # }else(theTable <- cbind(theTable, Type = theTable$Sample))
  # #byTypeSample = F
  # if(byTypeSample == T){
  #   theTable$Sample2 <- sapply(strsplit(as.character(theTable$Sample), "_"),"[[",1)
  # }else(theTable$Sample2 <- theTable$Sample)
  
  #yLab<- "test"
  #cols <- NA
  
  
  
  theTable$origSample <- factor(theTable$origSample, levels = theTable$origSample[1:length(unique(theTable$origSample))])
  
  #p4d
  #multPlot = T
  if(multPlot == T){
    #theTable$Target <- factor(theTable$Target)
    if(!is.null(targetOrder)){
      theTable$Target <- factor(theTable$Target, levels = targetOrder)
    }
    p4d <-  ggplot(dplyr::arrange(theTable, origSample), aes(x=Target,
                                                             y=pow2,
                                                             fill=origSample,
                                                             width= wBars ))+ #order=origSample,
      labs( x = "X-AXIS TITLE", y = "Y-AXIS TITLE")
  }else{
    p4d <-  ggplot(theTable, aes(x=Sample, y=pow2,fill=origSample, width= wBars ))+
      labs( x = "X-AXIS TITLE", y = "Y-AXIS TITLE")
  }
  
  if(!is.na(cols)){
    p4d <- p4d + geom_bar(position=position_dodge() ,stat = "identity",fill = cols)
  }else{
    p4d <- p4d + geom_bar(position=position_dodge() ,stat = "identity")
    #str(p4d)
    #aes(fill=Sample)) 
  }
  
  
  p4d <-  p4d + geom_errorbar(aes(ymin=pmax(0,pow2negSd), ymax=pow2posSd),
                              size=.3,
                              # Thinner lines
                              width=.1,
                              position=position_dodge(wBars)) +
    scale_x_discrete(name="")+   #labels=c(namesVar2) , limits=c(0, 30)
    scale_y_continuous(expand = c(0,0),
                       limits = c(0, maxLim),
                       name=yLab,
                       breaks = seq.int(from = 0, to = maxLim, by = byTickInterval))+
    guides(fill=FALSE)+ 
    theme_classic() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=fontSize,face = "plain", colour = "Black"), #family = "Arial"
          axis.title=element_text(size=fontSize) )+ #family = "Arial"
    ggtitle(titL)
  
  
  p4d
  
  if(blackWhite == T){
    p4d <- p4d + scale_fill_grey() #scale_color_grey()
  }
  
  
  
  
  # #devtools::install_github("s6juncheng/ggpval")
  # library("ggpval")
  # 
  #  
  # add_pval(p4d, pairs = list(c(1, 2)))
  # 
  # library(ggsignif)
  # 
  # p4d + geom_text(aes(label=Tukey, y=(Mean+Sd)), position=position_dodge(width=wBars), vjust= -0.5)
  # 
  # p4d + geom_signif(comparisons = list(c(0.5,1)), 
  #             map_signif_level=TRUE)
  # 
  # p4d + geom_signif(comparisons = list("NL","NL"), 
  #                   map_signif_level=TRUE)
  # 
  if(anova1 == T){
    p4d <- p4d + geom_text(aes(label=Tukey, y=(pow2posSd)), position=position_dodge(width=wBars), vjust= -0.5)
  }else{
    #x <- 2
    #x <- 3
    
    
    chosen_data <- chosen_data[order(chosen_data$Target),]
    
    
    for(y in 1:length(unique(chosen_data$Target))){
      chosen_data1 <- chosen_data[which(chosen_data$Target == unique(chosen_data$Target)[y]),]
      minCD <- min(chosen_data1$group)
      maxCD <- max(chosen_data1$group)
      highest <- 1
      count1 <- 0
      for(x in (minCD+1): maxCD){
        
        ttestList <- data.frame()
        
        pos4 <- na.omit(match(chosen_data1$group, c(minCD,x)))
        if(is.null(LevTest)){
          lev1 <- car::leveneTest(y ~ groups2, data=chosen_data1[seq_along(pos4),])      
        }else{
          lev1 <- LevTest
        }
        
        g1 <- chosen_data1$y[which(chosen_data1$group == minCD)]
        g2 <- chosen_data1$y[which(chosen_data1$group == x)]
        print(g2)
        print(g2)
        print(lev1$`Pr(>F)`[1])
        if(lev1$`Pr(>F)`[1] > 0.05){      
          
          ttest1 <- t.test(g1,g2, var.equal = T)$p.value
          message2 <- "Students's t-test"
          
        }else{
          # print(paste0(x,"st variable failed levene's test"))
          print("Doing Welch test instead of students t-test")
          message2 <- "Welch's t-test"
          ttest1 <- t.test(g1,g2, var.equal=F)$p.value
          
        }
        
        if(dplyr::between(ttest1, 0.05, 1)){
          ttestSign <- "" #"NS"
        }else if(dplyr::between(ttest1, 0.01,0.05)){
          ttestSign <- "*" #✱
        }else if(dplyr::between(ttest1, 0.001,0.01)){
          ttestSign <- "**"
        }else if(ttest1 < 0.001){
          ttestSign <- "***"
        }
        
        ttestRow <- cbind(Sample = argMatrix$Sample[x], ttest1, ttestSign)#, heightTtset 
        ttestList <- rbind(ttestList, ttestRow)
        
        #if(nrow(ttestList)== 1){
        if(multPlot == F){
          xmin1 <- minCD
          xmax1 <- x
        }else{
          n <- maxCD - minCD + 1 #3
          adj <- (n-1)/(n*2) #första wBars/ (maxCD - minCD + 1) #* wBars * 0.17 # 1/.136  wBars = 0.4
          #noIG/(noIG-1)  #(maxCD - minCD) * wBars * 0.17
          xmin1 <- y - (adj * wBars)
          xmax1 <- xmin1 + ((x - minCD)/n * wBars)   #((adj*2/noIG)*(x - minCD))  #(noIG-1) *
        }
        
        # xmin1 <- as.character(minCD)
        # xmax1 <- as.character(x)
        
        
        #maxes <- argMatrix[minCD:maxCD,]$pow2posSd
        extOver <- max(argMatrix[,"pow2posSd"]) * 0.06
        
        addHeight <- max(argMatrix[minCD:x,"pow2posSd"]) + extOver 
        #(max(argMatrix[minCD:x,"pow2posSd"]) * 0.05)
        
        maxMore <- which(argMatrix[minCD:x,"pow2posSd"] == max(argMatrix[minCD:x,"pow2posSd"]))
        if(maxMore != 1){
          if(highest != maxMore){
            count1 <- 0
          }
          highest <- maxMore
        }
        # if(maxMore == 1){
        #   maxMore <- 2
        # }
        
        #maxMore2 <- (x - minCD + 1) - maxMore
        #addHeight <- addHeight + maxMore2 * extOver
        addHeight <- addHeight + count1 * extOver
        
        
        if(ttestList$ttestSign != ""){
          count1 <- count1 + 1
          p4d <- p4d + 
            geom_signif(y_position=addHeight,
                        xmin=xmin1,
                        xmax=xmax1,
                        annotation=ttestList$ttestSign,
                        tip_length=0.01,
                        size = 0.3,#0.6,
                        textsize = 4)
          # geom_signif(
          #   comparisons = list(minCD, x),
          #   y_position=addHeight,
          #   #map_signif_level= ttestSign
          #   annotations=ttestList$ttestSign,
          #   #margin_top = 0.01
          #   )
        } 
        
        
        
        ret[["ttestList"]][[x-1]] <- ttestList
        
        #}
      }
    }
  }
  
  #print(p4d)
  if(exists("message2")){
    p4d <- list(p4d, message2)
  }
  return(list(p4d,ret))
}

server <- shinyServer(function(input, output, session){
  # Server that analyzes the qPCR data in an excel file and generates a barplot of the data  
  
  # Defines the reactive variables
  vals <- reactiveValues(bp=NULL,
                         qTable=NULL,
                         df3=NULL, 
                         order=NULL, 
                         evalTable=NULL,
                         #evalTable2 = NULL,
                         #qTable2=NULL, 
                         #outL =NULL,
                         outliers = NULL,
                         rename = NULL,
                         #melt = NULL,
                         stat = NULL,
                         mess = NULL)
  
  # Resets these reactive variables when a new file is loaded
  observeEvent(input$file1, {
    vals$evalTable <- NULL
    vals$outliers <- NULL
    vals$order <- NULL
    vals$rename = NULL
  })
  
  # UIs
  # Sheet UI
  output$sheet <- renderUI({
    inFile <- input$file1
    if(is.null(inFile))
      return(NULL)
    wb <- xlsx::loadWorkbook(inFile$datapath)  #inFile <- "KD1_5_6_Mir8878 qpcr data.xlsx"
    sheets <- names(xlsx::getSheets(wb))
    updateTabsetPanel(session, "set2", selected = "Reference")
    selectInput("sheet", "Sheet name:", sheets, selected = sheets[length(sheets)])
  }) 
  
  # Reference gene UI 
  output$hk <- renderUI({
    if(is.null(input$sheet))
      return(NULL)
    df2 <- xlsx::read.xlsx( input$file1$datapath, header = TRUE,sheetName = input$sheet)
    vals$df3 <- df2processing(df2)
    selectInput("hk", "Ref gene:", unique(vals$df3$Target))
  })
  
  # Control UI 
  output$ctrl <- renderUI({
    if(is.null(input$sheet))
      return(NULL)
    samp <- unique(vals$df3$Sample)
    samp2 <- unique(substrFromStartToChar(samp,input$separator))
    selectInput("ctrl", "ctrl:", samp2)
  })
  
  # Sample checkbox UI
  output$checkboxSamp <- renderUI({
    if(is.null(input$sheet))
      return(NULL)
    samp <- unique(vals$df3$Sample)
    samp2 <- unique(substrFromStartToChar(samp,input$separator))
    choice <- samp2[-which(samp2 == input$ctrl)]
    checkboxGroupInput("checkboxSamp","Select samples", choices = choice, selected = choice)
  })
  
  # Target checkbox UI
  output$checkboxTarg <- renderUI({
    if(is.null(input$sheet))
      return(NULL)
    targ <- unique(vals$df3$Target)
    targ <- targ[-which(targ == input$hk)]
    checkboxGroupInput(inputId = "checkboxTarg","Select targets", choices = targ, selected = targ)
  })
  
  # Rename UI
  output$rename <- renderUI({
    if(is.null(input$order))
      return(NULL)
    textInput("rename", "Relabel",value = paste(unique(input$order),collapse = " ")) #simulation code can go here
    
  })
  # Order UI
  output$order <- renderUI({
    if(!is.null(vals$rename)){
      rename1 <- unlist(strsplit(vals$rename, split = " ",fixed = T))
      vals$evalTable$Sample <- gsubMany(vals$evalTable$Sample,
                                        sampLorder = unique(vals$qTable$Sample),
                                        sampLnames = rename1
      ) 
      vals$qTable$Sample <-
        rename1[match(vals$qTable$Sample, unique(vals$qTable$Sample))]
      rename1[match(vals$evalTable$Sample, unique(vals$evalTable$Sample))]
    }
    
    if(is.null(vals$order)){
      orderInput(inputId = "order",label =  "Decide order",items = unique(vals$qTable$Sample)) #samps
    }else{
      orderInput(inputId = "order",label =  "Decide order",items = unique(vals$order)) #samps
    }
  })
  
  # Outlier UI
  output$outliers <- renderUI({
    if(is.null(input$sheet)) #checkboxTarg
      return(NULL)
    # vals$qTable2
    # targ <- unique(vals$df3$Target)
    # targ <- targ[-which(targ == input$hk)]
    print("vals$evalTable$Samples")
    print(vals$evalTable$Samples)
    checkboxGroupInput(inputId = "outliers","Select outliers", 
                       choices = unique(vals$evalTable$Sample),
                       selected = vals$outliers
    ) #paste(vals$evalTable$Target,vals$evalTable$Sample,sep=input$separator)
  })  

  # Tables
  # ddCT summary table
  output$table <- renderTable({
    toNumericAndRoundIfposDf(vals$qTable)
  })
  
  # Evaluation table for all samples 
  output$Eval <- renderTable({
    if(!is.null(input$outliers)){
      vals$evalTable[-whichIsNotNaMatch(vals$evalTable$Sample,input$outliers),]
    }else{
      vals$evalTable
    }
  })

  # Plots the barplot
  output$plot2<-renderPlot({
    if(is.null(input$order))
      return(NULL)
    
    plotOmess <- barPlotQpcr(argMatrix = vals$qTable,
                             yLab = "Relative expression",
                             wBars = 0.4,
                             blackWhite = input$color == F,
                             byType = F, # names with _, sets after _ to Type
                             shiftLevels = F,
                             byTypeSample = F,
                             scaleToFirst = F,
                             anova1 = input$anova,
                             cols = NA,
                             fontSize =10,
                             goi_ref = T,
                             sampLorder = input$order,
                             sampLnames = input$order,#  if(is.null(relab)){input$order} else{unlist(strsplit(relab,split = " "))},
                             multPlot = input$multi,
                             decTukey = T,
                             skipLevTest = input$skipLev,
                             LevTest = vals$levTest
                             
    )
    vals$bp <- plotOmess[[1]][[1]] 
    vals$mess <- plotOmess[[1]][[2]]
    vals$stat <- plotOmess[[2]]
    vals$bp
  })
  
  # Prints the statistics
  output$statistics <- renderPrint({
    if(is.null(vals$stat[[2]])){
      av <- vals$stat[[1]]
      tt <- vals$stat[["ttestList"]]
      if(is.null(vals$levTest)){
        return("No value for Levene")
      }
      printP("Performed Levene's: P =", round(vals$levTest$`Pr(>F)`[1],digits = 3))
      printP("Performed test:", vals$mess)
      print(av)
    }else{
      if(is.null(vals$levTest)){
        return("No value for Levene")
      }
      printP("Performed Levene's: P =", round(vals$levTest$`Pr(>F)`[1],digits = 3))
      printP("Performed test:", unique(vals$mess))
      cat("\n")
      print("In the Anova calculations further down, sample names are")
      print("represented by the numbers in the following table:") 
      chosen_Data <- vals$stat[["chosen_Data"]]
      chosen_Data_Info <- chosen_Data[-which(duplicated(chosen_Data[,"names"])),c("names","groups2")]
      colnames(chosen_Data_Info) <- c("Sample", "Number")
      print(chosen_Data_Info)
      for(x in 1:length(vals$stat)){
        if(!is.data.frame(vals$stat[[x]][[1]]) && !is.integer(vals$stat[[x]][[1]])){
          print("Anova")
          av <- vals$stat[[x]][[2]][[1]]
          print(av)
          summary(multcomp::glht(vals$stat[[1]][[2]][[1]], linfct = multcomp::mcp(groups2 = "Tukey")))
          print(summary(multcomp::glht(av, linfct = multcomp::mcp(groups2 = "Tukey"))))
        }else if (is.data.frame(vals$stat[[x]][[1]])){
          printP(unique(vals$mess))
          ttests <- do.call(rbind,vals$stat[[x]])
          ttests[,2] <- signif(as.numeric(as.character(ttests[,2])),digits = 3)
          if(any(is.na(ttests[,1]))){
            ttests <- ttests[which(!is.na(ttests[,1])),]
          }
          colnames(ttests) <- c("Sample", "p-val", "Sign")
          print(ttests)
        }
        cat("\n")
      }
    }
  })
  
  # Exports the plot
  output$export = downloadHandler(
    filename = function() {"plots.pdf"},
    content = function(file) {
      Cairo::CairoPDF(file, width = input$plot2_size$width/100,height = input$plot2_size$height/100)#onefile = TRUE
      gridExtra::grid.arrange(vals$bp)
      dev.off()
    }
  )

  # Switches set2 to view the Reference tab when a sheet is set
  observeEvent(input$sheet, {
    updateTabsetPanel(session, "set2", selected = "Reference")
  })
  
  # Switches set2 to view the Control tab when the reference tab is set 
  observeEvent(input$hk, {
    updateTabsetPanel(session, "set2", selected = "Control")
  })
  
  # Switches set3 to view the Targets tab when any Sample checkbox is changed 
  observeEvent(input$checkboxSamp, {
    updateTabsetPanel(session, "set3", selected = "Targets")
  })
  
  # Switches set5 to view the Order tab when any Target checkbox is changed 
  observeEvent(input$checkboxTarg, {
    updateTabsetPanel(session, "set5", selected = "Order")
  })

  # Switches set5 to view the Order tab when the rename button is clicked 
  observeEvent(input$buttonRename, {
    vals$rename <- input$rename
    updateTabsetPanel(session, "set5", selected = "Order")
  })
    
  # Viewes the Settings main tab and automatically clicks the ddCt button in set3 when the Exclude button in clicked
  observeEvent(input$buttonExclude, {
    updateTabsetPanel(session, "set3", selected = "buttonDdct")
    updateTabsetPanel(session, inputId = "tabs",selected = "Settings")
    shinyjs::click("buttonDdct",asis = T)
  })
  
  #Selects the Plot tab among the main tabs (updating the plot data) every time sidebarMenu tab is changed  
  observeEvent(input$tabs, {
    updateTabsetPanel(session, "main", selected = "Plot")
  })
  
  # Runs the ddCt calculations upon clicking the ddCt button 
  observeEvent(input$buttonDdct, {
    if(is.null(input$checkboxTarg))
      return(NULL)
    returnList <- qPCR_ddct(hk = input$hk,#EF actin
                            ctrl = input$ctrl,
                            path1 = input$file1$datapath, #"KD1_5_6_Mir8878 qpcr data.xlsx",
                            sheetndex1 = input$sheet,
                            samplSep = input$separator,
                            changeStart = F,
                            sem = input$sem,
                            lengthBio = 1,
                            checkSampleDiffFromMean = T,
                            removeBio = if(is.null(vals$outliers)){""}else{vals$outliers}  #c("B7")c("B5","B6","B7") c("B5")
    )
    qPCRacE <- returnList[[3]]
    if(is.null(vals$evalTable)){
      vals$evalTable <- toNumericAndRoundIfposDf(returnList[[1]])
    }
    vals$levTest <- returnList[[4]]
    qPCRacE2 <- qPCRacE[which(!is.na(match(qPCRacE$Target,input$checkboxTarg))),]
    qPCRacE2 <- qPCRacE2[which(!is.na(match(qPCRacE2$Sample,c(input$ctrl,input$checkboxSamp)))),]
    vals$qTable <- qPCRacE2
  })
  
  # Updates outliers reactive variable at changes in the outlier input
  observeEvent(input$outliers,{
    vals$outliers <- input$outliers
  })
}
)

# Create Shiny object
shinyApp(ui = ui, server = server)