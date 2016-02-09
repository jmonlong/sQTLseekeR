##' Find the splicing event(s) that differenciate pairs of transcripts.
##'
##' The input transcript structure is a data.frame with information about the CDSs and UTRs start/end positions as well as the DNA strand. The columns should be: 'transId' for the transcript ID; 'strand' for the DNA strand; 'cdsStart' for the CDS start positions; 'cdsEnds' for the end positions; 'utrStarts' and 'utrEnds' for the UTRs. The position in 'cdsStarts'/'cdsEnds'/'utrStart'/'utrEnds' must be a concatenation of the position with ',' as separator. See example.
##'
##' The classification code follows mostly the one defined by AStalavista (\url{http://genome.crg.es/astalavista/FAQ.html}). Numbers represent exon boundaries and are ordered by genomic position, '-' represents exon body, '^' a splice junction and ')' the end of the transcript. The state of the two transcript is separated by a ','. For example, '1-2^,3-4^' represent mutually exclusive exons, ',1^2-' means intron retention. In addition, we added an extra formatting: '<>' means that the event involves UTRs. Hence '<>,1^),2^)' means alternative 3' UTR.
##' @title Classify splicing events
##' @param df a data.frame that includes pairs of transcript IDs in columns 'tr.first' and 'tr.second'.
##' @param trans.struct a data.frame with the transcript structure, i.e. the location of its exons and eventually its UTRs. See Details for format.
##' @return a list with
##' \item{res}{input data.frame with two new columns: 'classCode' and 'classEvent'. See Details for interpretation.}
##' \item{stats}{a data.frame with the occurence of each event in the data.}
##' @author Jean Monlong
##' @export
##' @examples
##'
##' ## Creating a fake transcript structure
##' tr.str = data.frame(transId=c("t1","t2","t3"),strand="+",
##' cdsStarts=c("10,40,100","10,20,100","10,40,100"),
##' cdsEnds=c("15,55,130","15,30,130","15,55,130"),
##' utrStarts=c("5,130","5,130","5,130"),
##' utrEnds=c("10,135","10,135","10,150"))
##' tr.str
##'
##' ## Creating the data.frame with the transcript pairs
##' tr.df = data.frame(tr.first=c("t1","t1"), tr.second=c("t2","t3"))
##'
##' ## Calling the function
##' classify.events(tr.df, tr.str)
##'
classify.events <- function(df, trans.struct){
  ## Convert to character in case
  for(col in colnames(trans.struct)){
    trans.struct[,col] = as.character(trans.struct[,col])
  }
  rownames(trans.struct) = trans.struct$transId

  ## Translate a splicing code into a splicing event names
  translate.event <- function(events){
    ev.tr = c(
      ",1-2^"="Skipped exon",
      ",1^2-"="Intron retention",
      "1-2^,3-4^"="Mutually exclusive exon",
      "1-,2-"="Alternative 5' splice site",
      "1^,2^"="Alternative 3' splice site",
      "<>,(1-2^,(3-4^"="Alternative 5' UTR",
      "<>,1-2^),3-4^)"="Alternative 3' UTR",
      "<>,1^),2^)"="Tandem 3' UTR",
      "<>,(1-,(2-"="Tandem 5' UTR",
      "(1-2^,(3-4^"="Alternative first exon",
      "1-2^),3-4^)"="Alternative last exon",
      "1-2^,"="Skipped exon",
      "1^2-,"="Intron retention",
      "3-4^,1-2^"="Mutually exclusive exon",
      "2-,1-"="Alternative 5' splice site",
      "2^,1^"="Alternative 3' splice site",
      "<>,(3-4^,(1-2^"="Alternative 5' UTR",
      "<>,3-4^)1-2^)"="Alternative 3' UTR",
      "<>,2^),1^)"="Tandem 3' UTR",
      "<>,(2-,(1-"="Tandem 5' UTR",
      "(3-4^(1-2^"="Alternative first exon",
      "3-4^),1-2^)"="Alternative last exon")
    ev.res = rep("Complex splicing event",length(events))
    ev.res[grepl("\\(",events)] = "Complex 5' event"
    ev.res[grepl("\\)",events)] = "Complex 3' event"
    ev.res[events %in% names(ev.tr)] = ev.tr[events[events %in% names(ev.tr)]]
    ev.res
  }

  ## Retrieve splicing code and event names from two transcript ids
  classify.pair <- function(tr1 = "ENST00000451283.1",tr2 = "ENST00000215882.5"){
    ## Compute event code from two sets of positions.
    find.events <- function(tr1.pos,tr2.pos,check.5p=FALSE,check.3p=FALSE,posStrand=TRUE){
      pos.df = data.frame(pos=sort(unique(c(tr1.pos,tr2.pos)),decreasing=!posStrand))
      pos.df$tr1 = pos.df$pos %in% tr1.pos
      pos.df$tr2 = pos.df$pos %in% tr2.pos
      ii.3p = c(max(c(which(pos.df$tr1),1)),max(c(which(pos.df$tr2),1)))
      ii.5p = c(min(c(which(pos.df$tr1),1)),min(c(which(pos.df$tr2)),1))
      sep.c = c("-","^")
      s.1 = s.2 = 1
      c.cpt = 1
      c.tr1 = c.tr2 = ""
      ev.l = NULL
      for(ii in 1:nrow(pos.df)){
        if(pos.df[ii,"tr1"]!=pos.df[ii,"tr2"]){
          if(pos.df[ii,"tr1"]){
            c.tr1 = paste(c.tr1,c.cpt,sep.c[s.1],sep="")
          } else {
            c.tr2 = paste(c.tr2,c.cpt,sep.c[s.2],sep="")
          }
          c.cpt = c.cpt + 1
          if(check.5p & ii.5p[1]==ii) c.tr1 = paste("(",c.tr1,sep="")
          if(check.5p & ii.5p[2]==ii) c.tr2 = paste("(",c.tr2,sep="")
          if(check.3p & ii.3p[1]==ii) c.tr1 = paste(c.tr1,")",sep="")
          if(check.3p & ii.3p[2]==ii) c.tr2 = paste(c.tr2,")",sep="")
        } else if(c.cpt != 1){
          ev.l = c(ev.l, paste(c.tr1,c.tr2,sep=","))
          c.cpt = 1
          c.tr1 = c.tr2 = ""
        }
        if(pos.df[ii,"tr1"]) s.1 = 3 - s.1
        if(pos.df[ii,"tr2"]) s.2 = 3 - s.2
      }
      if(c.cpt != 1) ev.l = c(ev.l, paste(c.tr1,c.tr2,sep=","))
      ev.l
    }
    if(!any(tr1 %in% rownames(trans.struct)) | !any(tr2 %in% rownames(trans.struct)))
      return(list(class.df=data.frame(),class.stats=character(0)))
    tr1.cds.pos = as.integer(unlist(strsplit(c(trans.struct[tr1,"cdsStarts"],trans.struct[tr1,"cdsEnds"]),",")))
    tr2.cds.pos = as.integer(unlist(strsplit(c(trans.struct[tr2,"cdsStarts"],trans.struct[tr2,"cdsEnds"]),",")))
    mean.mm <- function(e)mean(c(min(e),max(e)))
    tr.center = c(mean.mm(tr1.cds.pos),mean.mm(tr2.cds.pos))
    posStrand = trans.struct[tr1,"strand"]=="+"
    tr1.utr.pos = as.integer(unlist(strsplit(c(trans.struct[tr1,"utrStarts"],trans.struct[tr1,"utrEnds"]),",")))
    tr2.utr.pos = as.integer(unlist(strsplit(c(trans.struct[tr2,"utrStarts"],trans.struct[tr2,"utrEnds"]),",")))
    tr1.utr5.pos = tr1.utr.pos[ifelse(posStrand,1,-1)*(tr1.utr.pos-tr.center[1])<0]
    tr2.utr5.pos = tr2.utr.pos[ifelse(posStrand,1,-1)*(tr2.utr.pos-tr.center[2])<0]
    tr1.utr3.pos = tr1.utr.pos[ifelse(posStrand,1,-1)*(tr1.utr.pos-tr.center[1])>0]
    tr2.utr3.pos = tr2.utr.pos[ifelse(posStrand,1,-1)*(tr2.utr.pos-tr.center[2])>0]
    ## Find events involving exons
    ev = find.events(tr1.cds.pos,tr2.cds.pos,check.5p=TRUE,check.3p=TRUE,posStrand=posStrand)
    ## Find events involving UTRs
    if(length(c(tr1.utr5.pos,tr2.utr5.pos))>0){
      utr.ev = find.events(tr1.utr5.pos,tr2.utr5.pos,check.5p=TRUE,posStrand=posStrand)
      if(length(utr.ev)>0){
        ev = c(ev,paste("<>",utr.ev,sep=","))
      }
    }
    if(length(c(tr1.utr3.pos,tr2.utr3.pos))>0){
      utr.ev = find.events(tr1.utr3.pos,tr2.utr3.pos,check.3p=TRUE,posStrand=posStrand)
      if(length(utr.ev)>0){
        ev = c(ev,paste("<>",utr.ev,sep=","))
      }
    }
    ## Tranlate events
    ev.t = unique(translate.event(ev))
    list(class.df = data.frame(trPair=paste(tr1,tr2,sep="-"),classCode=paste(unique(ev),collapse=";"),classEvent=paste(sort(ev.t),collapse=";")),
         class.stats = ev.t, anySpl=!all(grepl("\\(",ev) | grepl("\\)",ev)))
  }

  ## Find unique pairs of transcripts to compare
  tr.to.class = unique(t(apply(df[,c("tr.first","tr.second")],1,sort)))
  ## Find splicing events
  class.res = apply(tr.to.class,1,function(trs)classify.pair(trs[1],trs[2]))
  class.df = do.call(rbind, lapply(class.res,function(l)l$class.df))
  class.stats = do.call(rbind, lapply(class.res,function(l){
    if(length(l$class.stats)>0){
      data.frame(count=1,event=l$class.stats)
    } else {
      data.frame()
    }
  }))
  class.stats = aggregate(count~event,data=class.stats,sum)
  class.stats$prop = class.stats$count / sum(class.stats$count)
  class.stats$prop.sqtl = class.stats$count / nrow(tr.to.class)
  spl.count = sum(unlist(lapply(class.res,function(l)if(length(l$class.stats)>0) l$anySpl)))
  class.stats = rbind(class.stats, data.frame(event="any splicing",count=spl.count,prop=NA,prop.sqtl=spl.count/nrow(tr.to.class)))

  df$trPair = apply(df[,c("tr.first","tr.second")],1,function(trs)paste(sort(trs),collapse="-"))
  df = merge(df,class.df,all.x=TRUE)
  df$trPair = NULL

  return(list(res=df, stats=class.stats))
}
