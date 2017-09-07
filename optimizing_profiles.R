rm(list = ls())
#whole genome

library(ape)
library(seqinr)
library(reldna)
U00096.2 <- as.character(read.fasta('/home/jane/Документы/Misha/Genome_to_dataset_pro_U00096.2.txt', as.string = T)[[1]])
mpots_U00096.2 <- reldna::lseqspline1D.BP(s = U00096.2, width = 1, bound = c(1, nchar(U00096.2)), ref = 1)


library(zoo)

rollmean <- rollapply(mpots_U00096.2$mpot, 50, mean)
int <- (which.min(rollmean)-100):(which.min(rollmean)+100)

plot(mpots_U00096.2$mpot[int], type = 'l')
lines(rollmean[int], lty = 2 )

int_seq <- substr(U00096.2, min(int), max(int))

#all possible snps
int_seq_copy <- int_seq

snps <- c()
for (i in 1:nchar(int_seq)) {
  nuc <- substr(int_seq, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps <- c(snps, int_seq_copy)
    int_seq_copy <- int_seq
    
  }
}


snps_mpots <- c()
for (i in seq_along(snps)) {
  res <- lseqspline1D(snps[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps_mpots <- cbind(snps_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, sum))
plot(apply(snps_mpots, 2, min))
minimorum1 <- which.min(apply(snps_mpots, 2, min))
min_rollmean60 <- apply(snps_mpots, 2, function(x){min(rollapply(x, 60, mean))})
plot(min_rollmean60)
matplot((snps_mpots), type = 'l', lty = 1, col = topo.colors(ncol(snps_mpots)))

par(mfrow = c(2,4))
plot(dend, xlab = NULL)
matplot(snps_mpots[,which(cut ==1)], type = 'l')
matplot(snps_mpots[,which(cut ==2)], type = 'l')
matplot(snps_mpots[,which(cut ==3)], type = 'l')
matplot(snps_mpots[,which(cut ==4)], type = 'l')

res <- lseqspline1D(int_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(res$x%in%-300:300)
mpot <- res$mpot[inds]
x <- -300:300
dir.create('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps')

# #unlink('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps',recursive = T)
for (i in 1:ncol(snps_mpots)) {
  png(paste0('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/',  i, '.png'),width = 1000, height = 750, res =200)
  plot(x, snps_mpots[,i], type = 'l', xlab = 'Sequence (angstrom)', ylab ='EP vaue', ylim = range(snps_mpots), main = i)
  lines(x, mpot, col = 'orange')
  
  dev.off()
}
#and bendability
snps_bends <- c()
for (i in seq_along(snps)) {
  res <- bendability(snps[i], bound = c(1,201), width = 1)
  snps_bends <- cbind(snps_bends, res)
}


plot(apply(snps_bends, 2, sum))
plot(apply(snps_bends, 2, min))
minimorum1 <- which.min(apply(snps_bends, 2, min))
min_rollmean60 <- apply(snps_bends, 2, function(x){min(rollapply(x, 60, mean))})
plot(min_rollmean60)
matplot(snps_bends, type = 'l', lty = 1, col = topo.colors(201))

init_bend <- bendability(substr(mpots_U00096.2$seq, min(int), max(int)),  bound = c(1,201), width = 1)
x <- -100:100
dir.create('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps')
# #unlink('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps', recursive = T)
for (i in 1:ncol(snps_mpots)) {
  png(paste0('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps/', i, '.png'),width = 1000, height = 750, res =200)
  plot(x, snps_bends[,i], type = 'l', xlab = 'Sequence (angstrom)', ylab ='Bendability', ylim = range(snps_bends, na.rm = T), main = i)
  lines(x, init_bend, col = 'orange')
  
  dev.off()
}
library(fastcluster)
library(dendextend)

colfunc <- colorRampPalette(c('tomato', 'royalblue'))
colnames(snps_mpots) <- 1:603
hclusted <- hclust.vector(snps_mpots, method = 'ward')
hclusted %>% as.dendrogram %>% color_branches(k = length(hclusted$order), col = colfunc(length(hclusted$order))[hclusted$order])%>%raise.dendrogram(2)-> dend
plot(dend, xlab = NULL)


cut <- cutree(hclusted, k=4)




#system('cd /home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/
#convert -delay 10 *.png ani.gif')
system('cd /home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/

convert -delay 10 $(for i in $(seq 0 1 603); do echo to_gif${i}.png; done) -loop 0 animated.gif')


#next iteration with minimorum 375
snps1[minimorum1]
#all possible snps1
int_seq_copy <- snps1[minimorum1]

snps1 <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps1 <- c(snps1, int_seq_copy)
    
  }
}


snps1_mpots <- c()
for (i in seq_along(snps1)) {
  res <- lseqspline1D(snps1[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps1_mpots <- cbind(snps1_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)

minimorum2 <- which.min(apply(snps1_mpots, 2, min))

matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps1_mpots))
matlines(snps1_mpots, type = 'l', col =2, lty =1, ylim = range(snps1_mpots))

apply(snps_mpots, 2, min)[minimorum]
apply(snps1_mpots, 2, min)[minimorum1]


apply(snps1_mpots, 2, min)[minimorum1]

#next iteration with minimorum 221
snps1[minimorum2]
#all possible snps
int_seq_copy <- snps1[minimorum2]

snps2 <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps2 <- c(snps2, int_seq_copy)
    
  }
}


snps2_mpots <- c()
for (i in seq_along(snps2)) {
  res <- lseqspline1D(snps2[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps2_mpots <- cbind(snps2_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)
points(apply(snps2_mpots, 2, min), col =3)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)
points(apply(snps2_mpots, 2, sum), col =3)

minimorum2 <- which.min(apply(snps2_mpots, 2, min))

matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps2_mpots))
matlines(snps1_mpots, type = 'l', col =2, lty =1, ylim = range(snps2_mpots))
matlines(snps2_mpots, type = 'l', col =3, lty =1, ylim = range(snps2_mpots))

apply(snps_mpots, 2, min)[minimorum]
apply(snps2_mpots, 2, min)[minimorum1]


apply(snps2_mpots, 2, min)[minimorum1]



#bulding a giant loop

for (iter in 2:7) {

#next iteration with minimorum 375
assign(int_seq_copy, 
       get(paste0('snps', iter-1, '_mpots'))[get(paste0('minimorum', iter-1))])
#all possible snps1
#int_seq_copy <- get(paste0('snps', i+1))

tmp <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    tmp<- c(tmp, int_seq_copy)
    
  }
}

assign(paste0('snps', iter), tmp)

tmp1 <- c()
for (i in seq_along(tmp)) {
  res <- lseqspline1D(tmp[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  tmp1<- cbind(tmp1, res$mpot[inds])
}

assign(paste0('snps', iter, '_mpots'), tmp1)

assign(paste0('minimorum', iter), which.min(apply(tmp1, 2, min)))
}##
plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)

minimorum2 <- which.min(apply(snps1_mpots, 2, min))

matplot(snps_mpots, type = 'n', col =1, lty =1, ylim = range(snps1_mpots))
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)
for (i in 1:7) {
matlines(get(paste0('snps', i, '_mpots')), type = 'l', col =i+1, lty =1, ylim = range(snps1_mpots))
print(mean(apply(get(paste0('snps', i, '_mpots')), 2, min)))
}

#matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps1_mpots))


library(RColorBrewer)
colfunc <- brewer.pal(n = 7, name = 'Dark2')
svg('/home/jane/Документы/Misha/search_for_lowest_ep_snp.svg', width = 8, height = 10)
par(mfrow = c (4,2),
    mar = rep(0.5, 4))
matplot(get(paste0('snps', '_mpots')), type = 'l', col = colfunc[i], lty =1, ylim = range(snps1_mpots), main = 'E. coli, lowest EP region')
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)

for (i in 1:7) {
  matplot(get(paste0('snps', i, '_mpots')), type = 'l', col = colfunc[i], lty =1, ylim = range(snps1_mpots), main = paste0('Lowest EP SNP ', i,' iteration'))
  abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)
}
dev.off()

apply(snps_mpots, 2, min)[minimorum]
apply(snps1_mpots, 2, min)[minimorum1]


apply(snps1_mpots, 2, min)[minimorum1]

#to look at initial lowest EP
library(ape)
library(seqinr)
as.DNAbin(strsplit(int_seq, ''))
GC(strsplit(int_seq, '')[[1]])


library(zoo)

rollmean <- rollapply(mpots_U00096.2$mpot, 50, mean)


rollmean200 <- rollapply(mpots_U00096.2$mpot, 200, mean)

minN <- function(x, N=2){
  x <- -x
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

mins <- -minN(rollmean200, N = 1:200)

mins_pos <- sapply(mins, function(x) {which(rollmean200 == x)})

#EP profile itself
mins <- -minN(mpots_U00096.2$mpot, N = 1:200)

mins_pos <- sapply(mins, function(x) {which(mpots_U00096.2$mpot == x)})

plot(rollmean200, type ='l')
abline(v = mins_pos, col =2)

plot(mpots_U00096.2$mpot[(mins_pos[1]-50):(mins_pos[1]+50)], type = 'n', ylim = range(mpots_U00096.2$mpot))
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 2)

nos <- which(diff(sort(mins_pos))!=1)
sapply(mins_pos[nos], function(x) {
  lines(mpots_U00096.2$mpot[(x-50):(x+50)], type = 'l')
  }
)

#plot(mpots_U00096.2$mpot[(mins_pos-50):(mins_pos+50)], type = 'l')
#substr(U00096.2, start = mins_pos-50, stop = mins_pos+50)

substrings_lowest_ep <- sapply(mins_pos, function(x) {substr(U00096.2, start = x-50, stop = x+50)})
as.DNAbin(strsplit(substrings_lowest_ep, ''))

library(dplyr)
library(janeaustenr)
library(tidyr)
library(ngram)
d <- data_frame(txt = substrings_lowest_ep)
d
#  unnest_tokens(substrings_lowest_ep, text, token = "ngrams", n = 2) 
d %>%
  unnest_tokens(output = bigram, input = txt, token = "ngrams", n = 2) -> d_bigrams

d_bigrams %>%
  separate(bigram, c("word1", "word2"), sep = " ")
d_bigrams %>%
  separate(bigram, c("word1", "word2"), sep = " ") ->
  tidy_chekhov_bigrams

sa <- unlist(strsplit(substrings_lowest_ep, split = ''))
sa1 <- paste(sa, collapse = ' ')
ng2 <- ngram(sa1, n=2)
ng3 <- ngram(sa1, n=3)
ng4 <- ngram(sa1, n=4)
ng5 <- ngram(sa1, n=5)
ng6 <- ngram(sa1, n=6)

print(ng, output=  "full")
rbind(get.phrasetable(ng2), get.phrasetable(ng3), get.phrasetable(ng4),get.phrasetable(ng5),get.phrasetable(ng6))  

phr_ng6 <- get.phrasetable(ng6)
#artificial_seq <- paste0(unlist(strsplit(paste0(phr_ng6$ngrams[1:16], collapse = ''), split = ' ')), collapse = '')
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =10), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
artificial_output <- lseqspline1D(artificial_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(artificial_output$x%in%-300:300)
artificial_mpots <- artificial_output$mpot[inds]

plot(artificial_mpots, type = 'l', ylim = range(mpots_U00096.2$mpot))
#artificial against natural EP-lowest
sapply(mins_pos[nos], function(x) {
  lines(mpots_U00096.2$mpot[(x-300):(x+300)], type = 'l')
}
)

abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 'seagreen')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col = 'royalblue')
lines()

#searching through generated seqs

artificial_mpots <- c()
for (i in 1:15000) {
  
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =i), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
artificial_output <- lseqspline1D(artificial_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(artificial_output$x%in%-300:300)
tmp <- artificial_output$mpot[inds]
artificial_mpots <- cbind(artificial_mpots, tmp)
}

boxplot(apply(artificial_mpots, 2, min))
boxplot(apply(artificial_mpots, 2, mean))

which_min_min <- which.min(apply(artificial_mpots, 2, min))
which_min_mean <- which.min(apply(artificial_mpots, 2, mean))

plot(artificial_mpots[,which_min_min], type ='l')
lines(artificial_mpots[,which_min_mean], type ='l', col =2)
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 'seagreen')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col = 'royalblue')


#dynamical properties
dynchars<-function(seq, interval_size) {
  if (missing(seq))
    stop("Need to specify sequence (as a vector of chars)")
  
  if (missing(interval_size))
    stop("Need to specify interval size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  
  seq<-toupper(seq)
  seq<-c(seq, seq[2:(interval_size)])
  
  a<-3.4*10^(-10)
  I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)
  tau<-c(127, 99, 140, 84)
  
  csA<-cumsum(seq=='A') 
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')
  
  countA = csA[interval_size:length(csA)]-c(0, csA[1:(length(csA)-interval_size)])
  countT = csT[interval_size:length(csT)]-c(0, csT[1:(length(csT)-interval_size)])
  countG = csG[interval_size:length(csG)]-c(0, csG[1:(length(csG)-interval_size)])
  countC = csC[interval_size:length(csC)]-c(0, csC[1:(length(csC)-interval_size)])
  
  M<-cbind(countA, countT, countG, countC)/interval_size
  M_comp<-cbind(countT, countA, countC, countG)/interval_size
  M_comp<-apply(t(M_comp),1,rev) 
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  E01<-(8*(Ks*Vs)^0.5)* 6E23 / 4184
  d1<-((Ks*a^2)/Vs)^(0.5)/a;
  c1<-(Ks*a^2/Is)^0.5
  m1<-E01/c1/6.011E-26
  taus1<-as.numeric(M%*%tau) #!as.numeric conversion
  gc1 = M[,3] + M[,4]
  
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  
  E02<- 8*(Ks*Vs)^0.5  * 6E23 / 4184;
  d2<-((Ks*a^2)/Vs)^(0.5)/a;
  c2<-(Ks*a^2/Is)^0.5;
  m2<-E02/c2/6.011E-26;
  taus2<-as.numeric(M_comp%*%tau)
  gc2 = M_comp[,3] + M_comp[,4]
  
  dynchars_return<-list(E01=E01, d1=d1, c1=c1, m1=m1, taus1=taus1, gc1=gc1, E02=E02, d2=d2, c2=c2, m2=m2, taus2=taus2, gc2=gc2)
  
  return(dynchars_return)
  
}
