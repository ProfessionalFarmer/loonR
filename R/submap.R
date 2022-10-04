#############################################################
# submap.R                             Oct.14, 2008
#
#   Subclass Mapping (SubMap)
#     Hoshida Y, Brunet JP, Tamayo P, Golub TR, Mesirov JP
#     Broad Institute of MIT and Harvard University
#           
#############################################################

### submap (main) ###



message <- function (..., domain = NULL, appendLF = TRUE) {
  
}

submap.main <- function(
  input.data.A,
  input.data.B,
  input.cls.A,
  input.cls.B,
  output.filename="SubMap",
  ntag=100,
  nperm=50,
  nperm.fisher=1000,
  weighted.score.type=1,
  null.dist="pool",
  p.corr="FDR",
  clust.row=1,
  clust.col=1,
  nom.p.mat="F",
  create.legend="T",
  rnd.seed=47365321
)
{
  # for GenePattern
  
  ntag<-as.integer(ntag)
  
  nperm<-as.integer(nperm)
  
  nperm.fisher<-as.integer(nperm.fisher)
  
  weighted.score.type<-as.integer(weighted.score.type)
  
  rnd.seed<-as.numeric(rnd.seed)
  
  if (clust.row=="NA")
  {
    clust.row<-NA
  }
  else
  {
    clust.row<-as.integer(clust.row)
  }
  
  if (clust.col=="NA")
  {
    clust.col<-NA
  }
  else
  {
    clust.col<-as.integer(clust.col)
  }
  
  set.seed(rnd.seed)
  
  null.dist.fig <- "F"   # null distribution for Fisher statistic
  
  # Read input data & clustids
  original.data.A<-read.gct(filename=input.data.A)
  
  original.data.B<-read.gct(filename=input.data.B)
  
  data.A<-take.intersection(original.data.A,original.data.B)
  
  data.B<-take.intersection(original.data.B,original.data.A)
  
  num.intersection<-length(data.A[,1])
  
  if (num.intersection<=ntag)
  {
    stop(paste("### Number of intersection genes is ",num.intersection," (too small) ###",sep=""))
  }
  
  # Read clustid
  cand.cls.label.A<-read.cls(filename=input.cls.A)
  
  cand.cls.label.B<-read.cls(filename=input.cls.B)
  
  # Mutual mapping to generate ES matrices
  ES.1<-submap.perm(data.A,data.B,cand.cls.label.A,cand.cls.label.B,ntag,nperm,weighted.score.type)
  
  ES.2<-submap.perm(data.B,data.A,cand.cls.label.B,cand.cls.label.A,ntag,nperm,weighted.score.type)
  
  n.row<-length(ES.1[,1,1])
  
  n.col<-length(ES.1[1,,1])
  
  # Convert ES matrices into nominal p-values
  nom.p.matrix.1<-array(9,c(n.row,n.col,(nperm+1)))
  
  nom.p.matrix.2<-array(9,c(n.col,n.row,(nperm+1)))
  
  if (null.dist=="each")
  {
    for (r in 1:n.row)
    {
      for (c in 1:n.col)
      {
        each.ES.1<-ES.1[r,c,]
        
        each.ES.2<-ES.2[c,r,]
        
        ES.rank.1<-rank(each.ES.1)
        
        ES.rank.2<-rank(each.ES.2)
        
        nom.p.matrix.1[r,c,]<-((nperm+2-ES.rank.1)/(nperm+1))
        
        nom.p.matrix.2[c,r,]<-((nperm+2-ES.rank.2)/(nperm+1))
      }
    }
  }
  
  if (null.dist=="pool")
  {
    perm.ES.vector.1<-as.vector(ES.1[,,2:(nperm+1)])
    
    perm.ES.vector.2<-as.vector(ES.2[,,2:(nperm+1)])
    
    num.element<-length(perm.ES.vector.1)
    
    nom.p.perm.ES.vector.1<-((num.element+2-rank(perm.ES.vector.1))/(num.element+1))
    
    nom.p.perm.ES.vector.2<-((num.element+2-rank(perm.ES.vector.2))/(num.element+1))
    
    nom.p.perm.ES.1<-array(nom.p.perm.ES.vector.1,c(n.row,n.col,nperm))
    
    nom.p.perm.ES.2<-array(nom.p.perm.ES.vector.2,c(n.col,n.row,nperm))
    
    for (r in 1:n.row)
    {
      for (c in 1:n.col)
      {
        rank.comb.ES.vector.1<-rank(c(ES.1[r,c,1],perm.ES.vector.1))
        
        rank.comb.ES.vector.2<-rank(c(ES.2[c,r,1],perm.ES.vector.2))
        
        nom.p.matrix.1[r,c,1]<-((num.element+2-rank.comb.ES.vector.1[1])/(num.element+1))
        
        nom.p.matrix.2[c,r,1]<-((num.element+2-rank.comb.ES.vector.2[1])/(num.element+1))
      }
    }
    
    nom.p.matrix.1[,,2:(nperm+1)]<-nom.p.perm.ES.1
    
    nom.p.matrix.2[,,2:(nperm+1)]<-nom.p.perm.ES.2
  }
  
  
  # compute nominal-p of Fisher's statistics from "each cell" or "pool of cells"
  nom.fisher.p.matrix<-matrix(9,nrow=n.row,ncol=n.col)
  
  each.perm.fisher.stat<-vector(length=nperm.fisher,mode="numeric")
  
  if (null.dist=="each")
  {
    for (r in 1:n.row)
    {
      for (c in 1:n.col)
      {
        obs.fisher<-(-log(nom.p.matrix.1[r,c,1])-log(nom.p.matrix.2[c,r,1]))
        
        perm.p.pool.1<-nom.p.matrix.1[r,c,2:(nperm+1)]
        
        perm.p.pool.2<-nom.p.matrix.2[c,r,2:(nperm+1)]
        
        for (f in 1:nperm.fisher)
        {
          rand.nom.p.1<-sample(perm.p.pool.1,1,replace=T)
          
          rand.nom.p.2<-sample(perm.p.pool.2,1,replace=T)
          
          each.perm.fisher.stat[f]<-(-log(rand.nom.p.1)-log(rand.nom.p.2))
        }
        
        comb.vector<-c(obs.fisher,each.perm.fisher.stat)
        
        fisher.rank<-rank(comb.vector)
        
        nom.fisher.p.matrix[r,c]<-((nperm.fisher+2-fisher.rank[1])/(nperm.fisher+1))
      }
    }
  }
  
  if (null.dist=="pool")
  {
    perm.p.pool.1<-as.vector(nom.p.matrix.1[,,2:(nperm+1)])
    
    perm.p.pool.2<-as.vector(nom.p.matrix.2[,,2:(nperm+1)])
    
    for (f in 1:nperm.fisher)
    {
      rand.nom.p.1<-sample(perm.p.pool.1,1,replace=T)
      
      rand.nom.p.2<-sample(perm.p.pool.2,1,replace=T)
      
      each.perm.fisher.stat[f]<-(-log(rand.nom.p.1)-log(rand.nom.p.2))
    }
    
    for (r in 1:n.row)
    {
      for (c in 1:n.col)
      {
        obs.fisher<-(-log(nom.p.matrix.1[r,c,1])-log(nom.p.matrix.2[c,r,1]))
        
        comb.vector<-c(obs.fisher,each.perm.fisher.stat)
        
        fisher.rank<-rank(comb.vector)
        
        nom.fisher.p.matrix[r,c]<-((nperm.fisher+2-fisher.rank[1])/(nperm.fisher+1))
      }
    }
  }
  
  
  # p-value correction (SA matrix)
  bonf.SA.matrix<-fdr.SA.matrix<-matrix(9,nrow=n.row,ncol=n.col)
  
  bonf.SA.matrix<-nom.fisher.p.matrix*(n.row*n.col)
  
  bonf.SA.matrix[bonf.SA.matrix>1]<-1
  
  rank.matrix<-matrix(rank(nom.fisher.p.matrix),nc=length(nom.fisher.p.matrix[1,]))
  
  fdr.SA.matrix<-nom.fisher.p.matrix*(n.row*n.col)/rank.matrix
  
  fdr.SA.matrix[fdr.SA.matrix>1]<-1
  
  
  # output: text
  header=rbind("*** Subbclass Mappping Results ***","",
               
               paste("input.data.A:  ",as.character(input.data.A),sep=""),
               
               paste("input.data.B:  ",as.character(input.data.B),sep=""),
               
               paste("# of intersection genes:  ",num.intersection,sep=""),
               
               paste("input.cls.A:  ",as.character(input.cls.A),sep=""),
               
               paste("input.cls.B:  ",as.character(input.cls.B),sep=""),
               
               paste("# of marker genes:  ",ntag,sep=""),
               
               paste("# of permutation for nominal-p of ES:  ",nperm,sep=""),
               
               paste("# of permutation for nominal-p of Fisher's statistics:  ",nperm.fisher,sep=""),
               
               paste("SNR weight for ES (1=yes, 0=no):  ",weighted.score.type,sep=""),
               
               paste("choice of null distribution:  ",null.dist,sep=""),
               
               paste("p-value correction method (for Fisher's statistics):  ",p.corr,sep=""),
               
               ""
  )
  
  outputTable.file <- paste(output.filename,"_SubMapResult.txt",sep="")
  
  write.table(header,outputTable.file,quote=F,sep="\t",row.names=F,col.names=F)
  
  row.label<-paste("A",as.character(seq(1:n.row)),sep="")
  
  col.label<-paste("B",as.character(seq(1:n.col)),sep="")
  
  rownames(bonf.SA.matrix)<-rownames(fdr.SA.matrix)<-rownames(nom.fisher.p.matrix)<-rownames(nom.p.matrix.1)<-colnames(nom.p.matrix.2)<-row.label
  
  colnames(bonf.SA.matrix)<-colnames(fdr.SA.matrix)<-colnames(nom.fisher.p.matrix)<-colnames(nom.p.matrix.1)<-rownames(nom.p.matrix.2)<-col.label
  
  if (p.corr=="Bonferroni")
  {
    
    results<-list(SA.matrix=as.data.frame(bonf.SA.matrix),nominal.p.matrix.Fisher=as.data.frame(nom.fisher.p.matrix),nominal.p.matrix.ES.A.on.B=as.data.frame(nom.p.matrix.1[,,1]),nominal.p.matrix.ES.B.on.A=as.data.frame(t(nom.p.matrix.2[,,1])))
  }
  
  if (p.corr=="FDR")
  {
    results<-list(SA.matrix=as.data.frame(fdr.SA.matrix),nominal.p.matrix.Fisher=as.data.frame(nom.fisher.p.matrix),nominal.p.matrix.ES.A.on.B=as.data.frame(nom.p.matrix.1[,,1]),nominal.p.matrix.ES.B.on.A=as.data.frame(t(nom.p.matrix.2[,,1])))
  }
  
  if (p.corr=="both")
  {
    results<-list(Bonferroni.SA.matrix=as.data.frame(bonf.SA.matrix),FDR.SA.matrix=as.data.frame(fdr.SA.matrix),nominal.p.matrix.Fisher=as.data.frame(nom.fisher.p.matrix),nominal.p.matrix.ES.A.on.B=as.data.frame(nom.p.matrix.1[,,1]),nominal.p.matrix.ES.B.on.A=as.data.frame(t(nom.p.matrix.2[,,1])))
  }
  
  write.table(capture.output(results),outputTable.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
  
  # output: gct for SA matrix
  if (p.corr=="Bonferroni")
  {
    Name <- Description <- rownames(bonf.SA.matrix)
    
    output.bonf.SA.matrix <- cbind(Name,Description,bonf.SA.matrix)
    
    colnames.output.bonf.SA.matrix <- colnames(output.bonf.SA.matrix)
    
    output.bonf.SA.matrix <- t(cbind(colnames.output.bonf.SA.matrix,t(output.bonf.SA.matrix)))
    
    write.table("#1.2",paste(output.filename,"_Bonferroni_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_Bonferroni_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.bonf.SA.matrix,paste(output.filename,"_Bonferroni_SAmatrix.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
  }
  
  if (p.corr=="FDR")
  {
    Name <- Description <- rownames(fdr.SA.matrix)
    
    output.fdr.SA.matrix <- cbind(Name,Description,fdr.SA.matrix)
    
    colnames.output.fdr.SA.matrix <- colnames(output.fdr.SA.matrix)
    
    output.fdr.SA.matrix <- t(cbind(colnames.output.fdr.SA.matrix,t(output.fdr.SA.matrix)))
    
    write.table("#1.2",paste(output.filename,"_FDR_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_FDR_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.fdr.SA.matrix,paste(output.filename,"_FDR_SAmatrix.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
  }
  
  if (p.corr=="both")
  {
    Name <- Description <- rownames(bonf.SA.matrix)
    
    output.bonf.SA.matrix <- cbind(Name,Description,bonf.SA.matrix)
    
    colnames.output.bonf.SA.matrix <- colnames(output.bonf.SA.matrix)
    
    output.bonf.SA.matrix <- t(cbind(colnames.output.bonf.SA.matrix,t(output.bonf.SA.matrix)))
    
    write.table("#1.2",paste(output.filename,"_Bonferroni_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_Bonferroni_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.bonf.SA.matrix,paste(output.filename,"_Bonferroni_SAmatrix.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    Name <- Description <- rownames(fdr.SA.matrix)
    
    output.fdr.SA.matrix <- cbind(Name,Description,fdr.SA.matrix)
    
    colnames.output.fdr.SA.matrix <- colnames(output.fdr.SA.matrix)
    
    output.fdr.SA.matrix <- t(cbind(colnames.output.fdr.SA.matrix,t(output.fdr.SA.matrix)))
    
    write.table("#1.2",paste(output.filename,"_FDR_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_FDR_SAmatrix.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.fdr.SA.matrix,paste(output.filename,"_FDR_SAmatrix.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
  }
  
  
  # output: gct for each nominal-p matrices of ES
  if (nom.p.mat=="T")
  {
    output.nom.p.matrix.A.on.B <- cbind(Name,Description,nom.p.matrix.1[,,1])
    
    colnames.output.nom.p.matrix <- colnames(output.nom.p.matrix.A.on.B)
    
    output.nom.p.matrix.A.on.B <- t(cbind(colnames.output.nom.p.matrix,t(output.nom.p.matrix.A.on.B)))
    
    write.table("#1.2",paste(output.filename,"_nominal_p_matrix_AonB.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_nominal_p_matrix_AonB.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.nom.p.matrix.A.on.B,paste(output.filename,"_nominal_p_matrix_AonB.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    output.nom.p.matrix.B.on.A <- cbind(Name,Description,t(nom.p.matrix.2[,,1]))
    
    colnames.output.nom.p.matrix <- colnames(output.nom.p.matrix.B.on.A)
    
    output.nom.p.matrix.B.on.A <- t(cbind(colnames.output.nom.p.matrix,t(output.nom.p.matrix.B.on.A)))
    
    write.table("#1.2",paste(output.filename,"_nominal_p_matrix_BonA.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    
    write.table(paste(n.row,n.col,sep="\t"),paste(output.filename,"_nominal_p_matrix_BonA.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    write.table(output.nom.p.matrix.B.on.A,paste(output.filename,"_nominal_p_matrix_BonA.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
  }
  
  
  
  # output: heatmap of SA matrix
  min.matrix<-max.matrix<-9
  
  if (p.corr=="Bonferroni")
  {
    if (min(bonf.SA.matrix)==max(bonf.SA.matrix))
    {
      min.matrix<-.9999999999
      max.matrix<-1
    }
    else
    {
      min.matrix<-min(bonf.SA.matrix)
      
      max.matrix<-max(bonf.SA.matrix)
    }
    
    if (capabilities("png"))
    {
      png(paste(output.filename,"_Bonferroni_SAmatrix.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_Bonferroni_SAmatrix.pdf",sep=""))
    }
    
    heatmap(bonf.SA.matrix,Rowv=clust.row,Colv=clust.col,scale="none",col=rainbow(1000,start=(0+.7*min.matrix),end=.7*max.matrix))
    
    dev.off()
  }
  
  if (p.corr=="FDR")
  {
    if (min(fdr.SA.matrix)==max(fdr.SA.matrix))
    {
      min.matrix<-.9999999999
      
      max.matrix<-1
    }
    else
    {
      min.matrix<-min(fdr.SA.matrix)
      
      max.matrix<-max(fdr.SA.matrix)
    }
    
    if (capabilities("png"))
    {
      png(paste(output.filename,"_FDR_SAmatrix.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_FDR_SAmatrix.pdf",sep=""), bg="white")
    }
    
    heatmap(fdr.SA.matrix,Rowv=clust.row,Colv=clust.col,scale="none",col=rainbow(1000,start=(0+.7*min.matrix),end=.7*max.matrix))
    
    dev.off()
  }
  
  if (p.corr=="both")
  {
    if (min(bonf.SA.matrix)==max(bonf.SA.matrix))
    {
      min.matrix<-.9999999999
      
      max.matrix<-1
    }
    else
    {
      min.matrix<-min(bonf.SA.matrix)
      
      max.matrix<-max(bonf.SA.matrix)
    }
    
    if (capabilities("png"))
    {
      png(paste(output.filename,"_Bonferroni_SAmatrix.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_Bonferroni_SAmatrix.pdf",sep=""))
      
    }
    
    heatmap(bonf.SA.matrix,Rowv=clust.row,Colv=clust.col,scale="none",col=rainbow(1000,start=(0+.7*min.matrix),end=.7*max.matrix))
    
    dev.off()
    
    if (min(fdr.SA.matrix)==max(fdr.SA.matrix))
    {
      min.matrix<-.9999999999
      
      max.matrix<-1
    }
    else
    {
      min.matrix<-min(fdr.SA.matrix)
      
      max.matrix<-max(fdr.SA.matrix)
    }
    
    if (capabilities("png"))
    {
      png(paste(output.filename,"_FDR_SAmatrix.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_FDR_SAmatrix.pdf",sep=""))
    }
    
    heatmap(fdr.SA.matrix,Rowv=clust.row,Colv=clust.col,scale="none",col=rainbow(1000,start=(0+.7*min.matrix),end=.7*max.matrix))
    
    dev.off()
  }
  
  # output: heatmap of each nominal-p matrices of ES
  if (nom.p.mat=="T")
  {
    if (capabilities("png"))
    {
      png(paste(output.filename,"_nominal_p_matrix_AonB.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_nominal_p_matrix_AonB.pdf",sep=""))
      
    }
    
    heatmap(nom.p.matrix.1[,,1],Rowv=NA,Colv=NA,scale="none",col=rainbow(1000,start=(0+.7*min(nom.p.matrix.1[,,1])),end=.7*max(nom.p.matrix.1[,,1])))
    
    dev.off()
    
    t.nom.p.matrix.2<-t(nom.p.matrix.2[,,1])
    
    if (capabilities("png"))
    {
      png(paste(output.filename,"_nominal_p_matrix_BonA.png",sep=""))
    }
    else
    {
      pdf(paste(output.filename,"_nominal_p_matrix_BonA.pdf",sep=""))
    }
    
    heatmap(t.nom.p.matrix.2,Rowv=NA,Colv=NA,scale="none",col=rainbow(1000,start=(0+.7*min(t.nom.p.matrix.2)),end=.7*max(t.nom.p.matrix.2)))
    
    dev.off()
  }
  
  # legend of heatmap
  if (create.legend=="T")
  {
    if (capabilities("png"))
    {
      png("legend.png")
    }
    else
    {
      pdf("legend.pdf")
      
    }
    
    par(plt=c(.1,.9,.45,.5))
    
    a=matrix(seq(1:1000),nc=1)
    
    image(a,col=rainbow(1000,start=0,end=.7),xlim=c(0,1),yaxt="n")
    
    box()
    
    dev.off()
  }
  
  # histgram of pooled Fisher's statistics
  if (null.dist.fig=="T")
  {
    if (capabilities("png"))
    {
      png("null.distribution.of.fisher.png")
    }
    else
    {
      pdf("null.distribution.of.fisher.pdf")
    }
    
    hist(each.perm.fisher.stat,br=30,probability=T,main="Null distribution of Fisher's statistics (pooled)")
    
    lines(density(each.perm.fisher.stat),col="red",lwd=2)
    
    dev.off()
  }
  
  #  return(nom.fisher.p.matrix)
}




### submap.perm ###

submap.perm <- function(data.mkr,data.rnk,cand.cls.label.mkr,cand.cls.label.rnk,ntag=100,nperm=1000,weighted.score.type=1)
{
  # cls label, box
  cand.cls.box.mkr<-cand.cls.box(cand.cls.label.mkr)
  
  cand.cls.box.rnk<-cand.cls.box(cand.cls.label.rnk)
  
  
  num.subcls.mkr<-length(cand.cls.box.mkr)
  
  num.subcls.rnk<-length(cand.cls.box.rnk)
  
  # initialize
  perm.ES.matrices<-array(0,c(num.subcls.mkr,num.subcls.rnk,(nperm+1)))
  
  
  # observed ES matrix
  ES.matrix<-matrix(0,nrow=num.subcls.mkr,ncol=num.subcls.rnk)
  
  perm.ES.matrices[,,1]<-submap(data.mkr,data.rnk,cand.cls.label.mkr,cand.cls.label.rnk,ES.matrix,num.subcls.mkr,num.subcls.rnk,ntag,weighted.score.type)
  
  
  # permutated ES matrices
  for (i in 2:(nperm+1))
  {
    print(paste("start perm ",i,"-1",sep=""))
    
    perm.cand.cls.label.rnk<-sample(cand.cls.label.rnk)
    
    perm.ES.matrices[,,i]<-submap(data.mkr,data.rnk,cand.cls.label.mkr,perm.cand.cls.label.rnk,ES.matrix,num.subcls.mkr,num.subcls.rnk,ntag,weighted.score.type)
  }
  
  return(perm.ES.matrices)
}    ## end of submap.perm (Main)



### submap ###

submap <- function(data.mkr,data.rnk,cand.cls.label.mkr,cand.cls.label.rnk,ES.matrix,num.subcls.mkr,num.subcls.rnk,ntag=100,weighted.score.type=1)
{
  for (rnk in 1:num.subcls.rnk)
  {
    # convert rnk cls to binary
    
    cls.bin.rnk<-cand.cls.label.rnk
    
    cls.bin.rnk[cls.bin.rnk!=rnk]<-0
    
    cls.bin.rnk[cls.bin.rnk==rnk]<-1
    
    
    cls1<-cls0<-cls.bin.rnk
    
    cls1[cls1==0]<-NA
    
    cls0[cls0==1]<-NA
    
    cls0[cls0==0]<-1
    
    # order gene in data.rnk
    
    data<-as.matrix(data.rnk)
    
    data1.tr <-apply(t(data)*cls1, 2, na.omit)
    data1.tr <- as.matrix(data1.tr)
    
    data0.tr <- apply(t(data)*cls0, 2, na.omit)
    data0.tr <- as.matrix(data0.tr)
    
    
    mean1<-apply(data1.tr, 2, mean)
    
    mean0<-apply(data0.tr, 2, mean)
    
    sd1<-apply(data1.tr,2,sd)
    
    sd0<-apply(data0.tr,2,sd)
    
    #    rm(data1.tr)
    
    #    rm(data0.tr)
    
    #    gc()
    
    s2n<-(mean1-mean0)/(sd1+sd0)
    
    #    rm(mean1)
    
    #    rm(mean0)
    
    #    rm(sd1)
    
    #    rm(sd0)
    
    #    gc()
    
    order.gene.s2n.rnk<-cbind(order(s2n,decreasing=T),sort(s2n,decreasing=T))
    
    #    rm(s2n)
    
    #    gc()
    
    for (mkr in 1:num.subcls.mkr)
    {
      print(paste("### rank: ",rnk,"/",num.subcls.rnk,", mkr: ",mkr,"/",num.subcls.mkr," ###",sep=""))
      
      # convert mkr cls to binary
      
      cls.bin.mkr<-cand.cls.label.mkr
      
      cls.bin.mkr[cls.bin.mkr!=mkr]<-0
      
      cls.bin.mkr[cls.bin.mkr==mkr]<-1
      
      cls1<-cls0<-cls.bin.mkr
      
      cls1[cls1==0]<-NA
      
      cls0[cls0==1]<-NA
      
      cls0[cls0==0]<-1
      
      #        rm(cls.bin)
      
      #        gc()
      
      # order gene, take marker
      
      data<-as.matrix(data.mkr)
      
      data1.tr<-apply(t(data)*cls1,2,na.omit)
      
      data0.tr<-apply(t(data)*cls0,2,na.omit)
      
      
      data1.tr <- as.matrix(data1.tr)
      data0.tr <- as.matrix(data0.tr)
      
      mean1<-apply(data1.tr,2,mean)
      
      mean0<-apply(data0.tr,2,mean)
      
      sd1<-apply(data1.tr,2,sd)
      
      sd0<-apply(data0.tr,2,sd)
      
      #        rm(data1.tr)
      
      #        rm(data0.tr)
      
      #        gc()
      
      s2n<-(mean1-mean0)/(sd1+sd0)
      
      #        rm(mean1)
      
      #        rm(mean0)
      
      #        rm(sd1)
      
      #        rm(sd0)
      
      #        gc()
      
      order.gene.mkr<-order(s2n,decreasing=T)
      
      marker<-order.gene.mkr[1:ntag]
      
      # ES
      tag.indicator <- sign(match(order.gene.s2n.rnk[,1], marker, nomatch=0))
      
      no.tag.indicator <- 1 - tag.indicator
      
      N <- length(order.gene.s2n.rnk[,1])
      
      Nh <- length(marker)
      
      Nm <-  N - Nh
      
      if (weighted.score.type == 0)
      {
        correl.vector <- rep(1, N)
      }
      else
      {
        alpha <- weighted.score.type
        
        correl.vector <- abs(order.gene.s2n.rnk[,2]**alpha)
      }
      
      sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
      
      norm.tag    <- 1.0/sum.correl.tag
      
      norm.no.tag <- 1.0/Nm
      
      RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
      
      max.ES <- max(RES)
      
      min.ES <- min(RES)
      
      #        rm(RES)
      
      #        gc()
      
      if (max.ES > - min.ES)
      {
        ES.matrix[mkr,rnk] <- max.ES
      }
      else
      {
        ES.matrix[mkr,rnk] <- min.ES
      }
    }      # loop end num.subcls.rnk
  }        # loop end num.subcls.mkr
  
  return(ES.matrix)
}  # end of submap



# read clustid & generate candidate cls labels

read.cls<-function(filename="NULL")
{
  if (regexpr(".cls$",filename)==-1)
  {
    stop("### class data should be .cls file! ###")
  }
  
  
  line1 <- scan(filename, nlines=1, what="character", quiet=TRUE)
  
  numberOfSamples <- as.integer(line1[1])
  numberOfClasses <- as.integer(line1[2])
  
  
  cls<-as.vector(as.matrix(read.delim(filename,header=F,sep=" ",skip=2)))
  
  
  if (is.na(as.numeric(cls[1])))
  {
    stop("### 3rd line of cls file should be numeric! ###")
  }
  
  if (min(as.numeric(cls))==0)
  {
    cls<-as.numeric(cls)+1
  }
  
  if(numberOfSamples!=length(cls)) {
    stop(paste("\nFound ", length(cls), "labels on line 3 of", filename,  "but expected", numberOfSamples, sep=' '))
  }
  
  if(numberOfClasses!=length(unique(cls))) {
    stop(paste("\nFound ", length(unique(cls)), "classes on line 3 of", filename,  "but expected", numberOfClasses, sep=' '))
  }
  
  return(cls)
}



cand.cls.box<-function(clustid)
{
  table.clustid<-table(clustid)
  
  cls.box<-as.numeric(rownames(table.clustid))
  
  return(cls.box)
}



# read .gct file

read.gct<-function(filename="NULL")
{
  if (regexpr(".gct$",filename)==-1)
  {
    stop("### input data should be .gct file! ###")
  }
  
  # line 2 dimensions
  dimensions <- scan(filename, what=list("integer", "integer"), nmax=1, skip=1, quiet=TRUE)
  rows <- dimensions[[1]]
  columns <- dimensions[[2]]
  
  cat("\ndimensions: ")
  print(dimensions)
  
  data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
  
  data<-data[order(data[,1]),]
  
  data<-data[-1]
  
  if(ncol(data) != columns)
  {
    stop(paste("\nFound", ncol(data), "samples in", filename, "but expected", columns, "samples", sep=' ' ))
    
  }
  
  if(nrow(data) != rows)
  {
    stop(paste("\nFound", nrow(data), "features in", filename, "but expected", rows, "features", sep=' ' ))
  }
  
  return(data)
}



take.intersection<-function(data.A,data.B)
{
  feature.name<-rownames(data.A)
  
  data.A<-cbind(feature.name,data.A)
  
  feature.name.B<-as.data.frame(rownames(data.B))
  
  colnames(feature.name.B)<-"feature.name"
  
  intersected<-merge(feature.name.B,data.A)
  
  rownames(intersected)<-intersected[,1]
  
  intersected<-intersected[-1]
  
  return(intersected)
}




