library(data.table)
library(ggplot2)
library(cowplot)
library(here)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3a: heatmap of residue displacement from WT in E17K and P68-77Dup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d.wt <- fread(here('data/average_wt_rr_distance_map.tsv'),skip=4)
names(d.wt) <- paste0('pos',1:ncol(d.wt))
d.wt <- cbind(residue=paste0('pos',1:nrow(d.wt)),d.wt)

d.dup <- fread(here('data/average_dup_rr_distance_map.tsv'),skip=4)
names(d.dup) <- paste0('pos',1:ncol(d.dup))
d.dup <- cbind(residue=paste0('pos',1:nrow(d.dup)),d.dup)

d.e17k <- fread(here('data/average_e17k_rr_distance_map.tsv'),skip=4)
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



