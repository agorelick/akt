library(data.table)
library(ggplot2)
library(cowplot)
library(Rpdb)
library(ggrepel)
library(ggsignif)
library(here)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3a: heatmap of residue displacement from WT in E17K and P68-77Dup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d.wt <- fread(here('data/rr_distance/average_wt_rr_distance_map_3o96.tsv'),skip=4)
names(d.wt) <- paste0('pos',1:ncol(d.wt))
d.wt <- cbind(residue=paste0('pos',1:nrow(d.wt)),d.wt)

d.dup <- fread(here('data/rr_distance/average_dup_rr_distance_map_3o96.tsv'),skip=4)
names(d.dup) <- paste0('pos',1:ncol(d.dup))
d.dup <- cbind(residue=paste0('pos',1:nrow(d.dup)),d.dup)

d.e17k <- fread(here('data/rr_distance/average_e17k_rr_distance_map_3o96.tsv'),skip=4)
names(d.e17k) <- paste0('pos',1:ncol(d.e17k))
d.e17k <- cbind(residue=paste0('pos',1:nrow(d.e17k)),d.e17k)

m.wt <- melt(d.wt,id.var='residue')
m.wt$residue <- factor(gsub('pos','',m.wt$residue),levels=c(1:480))
m.wt$variable <- factor(gsub('pos','',m.wt$variable),levels=c(1:480))
names(m.wt) <- c('res1','res2','ca_distance')
m.wt <- m.wt[!is.na(m.wt$ca_distance),]

m.dup <- melt(d.dup,id.var='residue')
m.dup$residue <- factor(gsub('pos','',m.dup$residue),levels=c(1:480))
m.dup$variable <- factor(gsub('pos','',m.dup$variable),levels=c(1:480))
names(m.dup) <- c('res1','res2','ca_distance')
m.dup <- m.dup[!is.na(m.dup$ca_distance),]

m.e17k <- melt(d.e17k,id.var='residue')
m.e17k$residue <- factor(gsub('pos','',m.e17k$residue),levels=c(1:480))
m.e17k$variable <- factor(gsub('pos','',m.e17k$variable),levels=c(1:480))
names(m.e17k) <- c('res1','res2','ca_distance')
m.e17k <- m.e17k[!is.na(m.e17k$ca_distance),]

delta <- m.wt
setnames(delta,'ca_distance','wt_ca_distance')
delta$dup_ca_distance <- m.dup$ca_distance
delta$e17k_ca_distance <- m.e17k$ca_distance
delta$dup_delta <- delta$dup_ca_distance - delta$wt_ca_distance
delta$e17k_delta <- delta$e17k_ca_distance - delta$wt_ca_distance

m.dup <- reshape(delta[,c('res1','res2','dup_delta'),with=F],idvar=c('res1'),timevar=c('res2'),direction='wide')
m.e17k <- reshape(delta[,c('res1','res2','e17k_delta'),with=F],idvar=c('res1'),timevar=c('res2'),direction='wide')
key_residues <- c(17,308)

hm1 <- as.matrix(m.dup[,c(2:ncol(m.dup)),with=F])
rownames(hm1) <- 1:480
colnames(hm1) <- 1:480
hm1 <- hm1[rownames(hm1) %in% 5:108,colnames(hm1) %in% 150:408]
hm1 <- melt(hm1)
hm1 <- as.data.table(hm1); names(hm1) <- c('PH','Kinase','distance')
hm1$AKT1 <- '68-77 Duplication'
hm1$distance[abs(hm1$distance) < 3.1] <- 0

hm2 <- as.matrix(m.e17k[,c(2:ncol(m.e17k)),with=F])
rownames(hm2) <- 1:480
colnames(hm2) <- 1:480
hm2 <- hm2[rownames(hm2) %in% 5:108,colnames(hm2) %in% 150:408]
hm2 <- melt(hm2)
hm2 <- as.data.table(hm2); names(hm2) <- c('PH','Kinase','distance')
hm2$distance[abs(hm2$distance) < 3.1] <- 0
hm2$AKT1 <- 'E17K'

hm <- rbind(hm1,hm2)
hm$AKT1 <- factor(hm$AKT1,levels=c('68-77 Duplication','E17K'))
setnames(hm,'distance','change in distance from wt')
lims <- range(hm1$distance)

p.dup <- ggplot() +
    scale_y_continuous(breaks=17) + scale_x_continuous(breaks=308) + 
    geom_tile(data=hm1,aes(x=Kinase,y=PH,fill=distance)) +
    scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=lims,name='Change in distance [A] ') +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position='bottom') +
    labs(x='Kinase domain residues',y='PH domain residues',title='68-77 duplication')

p.e17k <- ggplot() +
    scale_y_continuous(breaks=17) + scale_x_continuous(breaks=308) + 
    geom_tile(data=hm2,aes(x=Kinase,y=PH,fill=distance)) +
    scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=lims,name='Change in distance [A] ') +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position='bottom') +
    labs(x='Kinase domain residues',y='PH domain residues',title='E17K mutation')

p <- plot_grid(p.e17k,p.dup,ncol=2,nrow=1,rel_widths=c(1,1))
ggsave(here('fig_3a.pdf'),width=7,height=4)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig S3a: heatmap of residue displacement from WT in E17K and P68-77Dup for 4ekk
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d.wt <- fread(here('data/rr_distance/average_wt_rr_distance_map_4ekk.tsv'),skip=4)
names(d.wt) <- paste0('pos',1:ncol(d.wt))
d.wt <- cbind(residue=paste0('pos',1:nrow(d.wt)),d.wt)

d.dup <- fread(here('data/rr_distance/average_dup_rr_distance_map_4ekk.tsv'),skip=4)
names(d.dup) <- paste0('pos',1:ncol(d.dup))
d.dup <- cbind(residue=paste0('pos',1:nrow(d.dup)),d.dup)

d.e17k <- fread(here('data/rr_distance/average_e17k_rr_distance_map_4ekk.tsv'),skip=4)
names(d.e17k) <- paste0('pos',1:ncol(d.e17k))
d.e17k <- cbind(residue=paste0('pos',1:nrow(d.e17k)),d.e17k)

m.wt <- melt(d.wt,id.var='residue')
m.wt$residue <- factor(gsub('pos','',m.wt$residue),levels=c(1:480))
m.wt$variable <- factor(gsub('pos','',m.wt$variable),levels=c(1:480))
names(m.wt) <- c('res1','res2','ca_distance')
m.wt <- m.wt[!is.na(m.wt$ca_distance),]

m.dup <- melt(d.dup,id.var='residue')
m.dup$residue <- factor(gsub('pos','',m.dup$residue),levels=c(1:480))
m.dup$variable <- factor(gsub('pos','',m.dup$variable),levels=c(1:480))
names(m.dup) <- c('res1','res2','ca_distance')
m.dup <- m.dup[!is.na(m.dup$ca_distance) & !is.na(res1) & !is.na(res2),]

m.e17k <- melt(d.e17k,id.var='residue')
m.e17k$residue <- factor(gsub('pos','',m.e17k$residue),levels=c(1:480))
m.e17k$variable <- factor(gsub('pos','',m.e17k$variable),levels=c(1:480))
names(m.e17k) <- c('res1','res2','ca_distance')
m.e17k <- m.e17k[!is.na(m.e17k$ca_distance),]

delta <- m.wt
setnames(delta,'ca_distance','wt_ca_distance')
delta$dup_ca_distance <- m.dup$ca_distance
delta$e17k_ca_distance <- m.e17k$ca_distance
delta$dup_delta <- delta$dup_ca_distance - delta$wt_ca_distance
delta$e17k_delta <- delta$e17k_ca_distance - delta$wt_ca_distance

m.dup <- reshape(delta[,c('res1','res2','dup_delta'),with=F],idvar=c('res1'),timevar=c('res2'),direction='wide')
m.e17k <- reshape(delta[,c('res1','res2','e17k_delta'),with=F],idvar=c('res1'),timevar=c('res2'),direction='wide')
key_residues <- c(17,308)

hm1 <- as.matrix(m.dup[,c(2:ncol(m.dup)),with=F])
rownames(hm1) <- 1:480
colnames(hm1) <- 1:480
hm1 <- hm1[rownames(hm1) %in% 5:108,colnames(hm1) %in% 150:408]
hm1 <- melt(hm1)
hm1 <- as.data.table(hm1); names(hm1) <- c('PH','Kinase','distance')
hm1$AKT1 <- '68-77 Duplication'

hm2 <- as.matrix(m.e17k[,c(2:ncol(m.e17k)),with=F])
rownames(hm2) <- 1:480
colnames(hm2) <- 1:480
hm2 <- hm2[rownames(hm2) %in% 5:108,colnames(hm2) %in% 150:408]
hm2 <- melt(hm2)
hm2 <- as.data.table(hm2); names(hm2) <- c('PH','Kinase','distance')
hm2$AKT1 <- 'E17K'

hm <- rbind(hm1,hm2)
hm$AKT1 <- factor(hm$AKT1,levels=c('68-77 Duplication','E17K'))
setnames(hm,'distance','change in distance from wt')

lims <- c(-15,15)

p.dup <- ggplot() +
    scale_y_continuous(breaks=17) + scale_x_continuous(breaks=308) + 
    geom_tile(data=hm1,aes(x=Kinase,y=PH,fill=distance)) +
    scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=lims,name='Change in distance [A] ') +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position='bottom') +
    labs(x='Kinase domain residues',y='PH domain residues',title='68-77 duplication')

p.e17k <- ggplot() +
    scale_y_continuous(breaks=17) + scale_x_continuous(breaks=308) + 
    geom_tile(data=hm2,aes(x=Kinase,y=PH,fill=distance)) +
    scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=lims,name='Change in distance [A] ') +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position='bottom') +
    labs(x='Kinase domain residues',y='PH domain residues',title='E17K mutation')

p <- plot_grid(p.e17k,p.dup,ncol=2,nrow=1,rel_widths=c(1,1))
ggsave(here('fig_S3a.pdf'),width=7,height=4)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we find pairs of hydrophobic-interacting residues in the wild-type structure (3o96, time-averaged) and then check how the RMSD between these pairs of residues change in the different alleles/conformtions. This indicates that when starting from closed conformation (3o96), the DUP causes a significant increase in distance between these interacting residues, more closely resembling the open state. In the open state, all alleles and wild-type show similarly large distances. 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define helper functions
theme_std <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_classic(base_size = base_size, base_family = 'ArialMT', base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme( line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"), text = element_text(family = 'ArialMT', face = "plain", colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F), axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)), axis.ticks = element_line(colour = "black", size=rel(1)), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = rel(1)), legend.key = element_blank(), strip.background = element_blank())
}


table.freq <- function(value) {
    if (is.null(value) | length(value) == 0) {
        tbl <- data.table(value = NA, N = NA)
    }
    else {
        tbl <- as.data.table(table(value))
        tbl <- tbl[order(tbl$N, decreasing = T), ]
    }
    tbl
}

adt <- function(d) as.data.table(d)


## Residues in PH/Kinase domains (via UNIPROT): 
## PH domain: 5-108
## Kinase domain: 150-408

## these are the hydrophobic interaction residues from WT 3o96 (time-averaged)
intx <- fread(here('data/ph_kinase_hydrophobic_interactions.txt'))
field_order <- names(intx)

## replace 3-letter amino acids with 1-letter
dict <- fread(here('data/amino_acid_table.tsv'))
intx <- merge(intx, dict, by.x='aa_ph', by.y='abbrev', all.x=T)
intx$aa_ph <- intx$amino_acid
intx[,amino_acid:=NULL]
intx <- merge(intx, dict, by.x='aa_kinase', by.y='abbrev', all.x=T)
intx$aa_kinase <- intx$amino_acid
intx[,amino_acid:=NULL]
intx <- intx[,(field_order),with=F]

## format table
intx$pair <- paste(intx$resid_ph,intx$resid_kinase,sep='::')
intx$resid_ph <- as.integer(gsub('[.]A','',intx$resid_ph))
intx$resid_kinase <- as.integer(gsub('[.]A','',intx$resid_kinase))
intx$label_ph <- paste0(intx$aa_ph,intx$resid_ph)
intx$label_kinase <- paste0(intx$aa_kinase,intx$resid_kinase)
intx <- intx[order(resid_ph,resid_kinase),]
ph_levels <- unique(intx$label_ph[order(intx$resid_ph)])
kinase_levels <- unique(intx$label_kinase[order(intx$resid_kinase)])

## subset the interactions for pairs of residues with 10+ interactions
tbl <- table.freq(intx$pair)
intx <- intx[pair %in% tbl$value[tbl$N >= 10],]
intx <- intx[!duplicated(pair),]
intx <- intx[,c('pair','resid_ph','resid_kinase'),with=F]

load_pdb <- function(pdbfile, intx, allele, conformation, dup=F) {
	message(pdbfile)
	d <- read.pdb(pdbfile)$atom
	d <- as.data.table(d)
	d <- d[elename=='CA',]
	d <- d[,c('resid','x1','x2','x3'),with=F]
	if(dup==T) {
		d <- d[!resid %in% 78:87]
		d[resid >=78,resid:=resid-10]
	}
	intx_out <- merge(intx, d, by.x='resid_ph', by.y='resid', all.x=T)
	intx_out <- merge(intx_out, d, by.x='resid_kinase', by.y='resid', all.x=T)
	intx_out$rmsd <- sqrt((intx_out$x1.x-intx_out$x1.y)^2 + (intx_out$x2.x-intx_out$x2.y)^2 + (intx_out$x3.x-intx_out$x3.y)^2)
	intx_out$allele <- allele
	intx_out$conformation <- conformation
	intx_out <- intx_out[,c('allele','conformation','pair','rmsd'),with=F]
	intx_out
}

intx_WT_3o96 <- load_pdb('data/structures/WT_3o96_average.pdb', intx, 'WT', 'Closed')
intx_E17K_3o96 <- load_pdb('data/structures/E17K_3o96_average.pdb', intx, 'E17K', 'Closed')
intx_P68_C77_3o96 <- load_pdb('data/structures/P68-C77dup_3o96_average.pdb',intx, 'P68-C77dup','Closed', dup=T)
intx_3o96 <- rbind(intx_WT_3o96, intx_E17K_3o96, intx_P68_C77_3o96)

intx_WT_4ekk <- load_pdb('data/structures/WT_4ekk_average.pdb', intx, 'WT', 'Open')
intx_E17K_4ekk <- load_pdb('data/structures/E17K_4ekk_average.pdb', intx, 'E17K', 'Open')
intx_P68_C77_4ekk <- load_pdb('data/structures/P68-C77dup_4ekk_average.pdb',intx, 'P68-C77dup','Open', dup=T)
intx_4ekk <- rbind(intx_WT_4ekk, intx_E17K_4ekk, intx_P68_C77_4ekk)

res <- rbind(intx_3o96, intx_4ekk)
res$conformation <- factor(res$conformation, levels=c('Closed','Open'))
res$allele <- factor(res$allele, levels=c('WT','E17K','P68-C77dup'))

get_median <- function(res) list(mid=median(res$rmsd))
info <- res[,get_median(.SD),by=c('allele','conformation')]

getp <- function(test.allele1, test.allele2, test.conformation1, test.conformation2=NA) { 
    if(is.na(test.conformation2)) test.conformation2 <- test.conformation1
    tmp1 <- res[allele==test.allele1 & conformation==test.conformation1,c('pair','rmsd'),with=F]
    tmp2 <- res[allele==test.allele2 & conformation==test.conformation2,c('pair','rmsd'),with=F]
    tmp <- merge(tmp1,tmp2,by='pair')
    p=wilcox.test(tmp$rmsd.x,tmp$rmsd.y)$p.value
    p=paste0('P=',prettyNum(p,digits=1))
    p
}

p_Open <- c(
  getp('WT','E17K','Open'),
  getp('WT','P68-C77dup','Open'),
  getp('E17K','P68-C77dup','Open'))

p_Closed <- c(
  getp('WT','E17K','Closed'),
  getp('WT','P68-C77dup','Closed'),
  getp('E17K','P68-C77dup','Closed'))

comparisons <- list(c(1,2),c(1,3),c(2,3))
y_pos <- 34 + (0:5)*3

p1 <- ggplot(res[conformation=='Open'], aes(x=allele,y=rmsd)) +
    scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30)) +
    geom_point(position=position_jitter(width=0.2,seed=42),color='#BFBFBF') +
    geom_boxplot(fill=NA,width=0.3,outlier.shape=NA) +
    geom_signif(comparisons=comparisons, annotations=p_Open, y_position=y_pos,tip_length=NA) +
    theme_std(base_size=14) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    labs(x=NULL,y='RMSD (hydrophobic-\ninteracting residues) [Å]',subtitle='Open') 

p2 <- ggplot(res[conformation=='Closed'], aes(x=allele,y=rmsd)) +
    scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30)) +
    geom_point(position=position_jitter(width=0.2,seed=42),color='#BFBFBF') +
    geom_boxplot(fill=NA,width=0.3,outlier.shape=NA) +
    geom_signif(comparisons=comparisons, annotations=p_Closed, y_position=y_pos,tip_length=NA) +
    theme_std(base_size=14) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    labs(x=NULL,y='RMSD (hydrophobic-\ninteracting residues) [Å]',subtitle='Closed')

p <- plot_grid(p2,p1,align='h',ncol=2)
ggsave(here('fig_S3b.pdf'),width=6,height=4)






