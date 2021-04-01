library(lhs)
library(scatterplot3d)
#randomLHS(n,k):建立拉丁方格設計(LHD)抽樣,n=樣本大小,k=變數(factor)個數or維度
set.seed(109225004)
b <- randomLHS(n=4,k = 3,preserveDraw = T) #3維拉丁方格取4點
layout(matrix(c(1,2,3,4),2,2))
scatterplot3d(b,color = 1:4,xlab = "x",ylab = "y",zlab = "z",main="3D")
plot(b[,c(1,2)],ylim=c(0,1),xlim=c(0,1),xlab = "x",ylab="y",col = 1:4,main = "x-y")
abline(h=seq(0,1,0.25),v=seq(0,1,0.25),lty=2)
text(b[,1],b[,2]+0.1,c(1:4))
plot(b[,c(1,3)],ylim=c(0,1),xlim=c(0,1),xlab = "x",ylab="z",col = 1:4,main = "x-z")
abline(h=seq(0,1,0.25),v=seq(0,1,0.25),lty=2)
text(b[,1],b[,3]+0.1,c(1:4))
plot(b[,c(2,3)],ylim=c(0,1),xlim=c(0,1),xlab = "y",ylab="z",col = 1:4,main = "y-z")
abline(h=seq(0,1,0.25),v=seq(0,1,0.25),lty=2)
text(b[,2],b[,3]+0.1,c(1:4))
#augmentLHS(lhs,m):在現有的LHD(利用randomLHS生成)以符合拉丁方格屬性方式添加新的樣本點,m:新增樣本數
set.seed(109225004)
c <- augmentLHS(b,1) #新增1樣本點
layout(matrix(c(1,2,3,4),2,2))
scatterplot3d(c,color = 1:nrow(c),xlab = "x",ylab = "y",zlab = "z",main="3D")
plot(c[,c(1,2)],ylim=c(0,1),xlim=c(0,1),xlab = "x",ylab="y",col = 1:nrow(c),main = "x-y")
abline(h=seq(0,1,0.2),v=seq(0,1,0.2),lty=2)
text(c[,1],c[,2]+0.1,c(1:nrow(c)))
plot(c[,c(1,3)],ylim=c(0,1),xlim=c(0,1),xlab = "x",ylab="z",col = 1:nrow(c),main = "x-z")
abline(h=seq(0,1,0.2),v=seq(0,1,0.2),lty=2)
text(c[,1],c[,3]+0.1,c(1:nrow(c)))
plot(c[,c(2,3)],ylim=c(0,1),xlim=c(0,1),xlab = "y",ylab="z",col = 1:nrow(c),main = "y-z")
abline(h=seq(0,1,0.2),v=seq(0,1,0.2),lty=2)
text(c[,2],c[,3]+0.1,c(1:nrow(c)))
#create_oalhs(n,k,bChooseLargerDesign,bverbose):OA&LHD結合
#bChooseLargerDesign:是否選擇比n,k大的OA
set.seed(109225004)
d <- create_oalhs(9, 4, TRUE, F)
#geneticLHS(n=10,k=2,pop=100,gen=4,pMut=0.1,criterium="S",verbose)="F"):LHD&Genetic Algorithm結合.
set.seed(109225004)
e <- geneticLHS(4, 3, 50, 5, .25)
#improvedLHS(n,k,dup=1):Improved Latin Hypercube Sample
set.seed(109225004)
f <- improvedLHS(4, 3, 2)
#maximinLHS(n,k,method="build",dup=1,eps=0.05,maxIter=100,optimize.on="grid",debug=F):Maximin Latin Hypercube Sample
#method:build,iterative
#optimize.on:grid,
set.seed(109225004)
g <- maximinLHS(4, 3, dup=2)
#oa_to_oalhs(n,k,oa):Latin hypercube from an orthogonal array
set.seed(109225004)
oa <- createBose(3, 4, TRUE)
h <- oa_to_oalhs(9, 4, oa) #(9,4)與oa的維度相同
#optAugmentLHS(lhs,m=1,mult=2):Optimal Augmented Latin Hypercube Sample
set.seed(109225004)
i <- optAugmentLHS(b, 2, 3)
#optimumLHS(n=10,k=2,maxSweeps=2,eps=0.1,verbose=F):Optimum Latin Hypercube Sample
set.seed(109225004)
j <- optimumLHS(4, 3, 5, .05)
#optSeededLHS(seed,m=0,maxSweeps=2,eps=0.1,verbose=FALSE):Optimum Seeded Latin Hypercube Sample
set.seed(109225004)
k <- optSeededLHS(b, 2, 2, .1)

#正交陣列(orthgonal array)性質:每行中不同數字出現次數相同,每兩行數字一定會有全部組合。即每個factor的level與另一個factor的level必有完全組合
#OA(n,k,q,t):n列,k行,q個symbol,t強度(任t行具完全組合)。對k具有限制
#createAddelKemp(q,ncol,bRandom):以Addelman-Kempthorne algorithm建立OA(a*q^2,k,q,2)
#q:level數,最好是奇數&質數,ncol:factor
set.seed(109225004)
OA1 <- createAddelKemp(3,3) #18列
#createAddelKemp3(q,ncol,bRandom):建立OA(2*q^3,k,q,2)
set.seed(109225004)
OA2 <- createAddelKemp3(3,3) #54列
#createAddelKempN(q,ncol,exponent,bRandom):建立OA(2*q^n,k,q,2)
set.seed(109225004)
OA3 <- createAddelKempN(3,3,4) #162列
#createBose(q,ncol,bRandom):以Bose algorithm建立OA(q^2,k,q,2)
set.seed(109225004)
OA4 <- createBose(3,3) #9列
#createBoseBush(q,ncol,bRandom):以Bose-Bush algorithm建立OA(2*q^2,k,q,2),q為2的次方數
set.seed(109225004)
OA5 <- createBoseBush(4,3) #32列
#createBoseBushl(q,ncol,lambda,bRandom):Bose-Bush algorithm with alternate strength>=3,OA(lambda*q^2,k,q,2)
set.seed(109225004)
OA6 <- createBoseBushl(3, 3, 3, TRUE)
#createBush(q,ncol,bRandom):Bush algorithm.OA(q^3,k,q,3)
set.seed(109225004)
OA7 <- createBush(3,3)
#createBusht(q,ncol,strength,bRandom):Bush algorithm.OA(q^t,k,q,t),t>=3
set.seed(109225004)
OA8 <- createBusht(5, 4, 4)

#runifint(n,min_int,max_int):與sample(seq(min_int,max_int),n,replace=T)相同,從離散均勻分配抽n個點
set.seed(109225004)
runifint(10,0,10) #7 4 8 8 7 1 5 2 4 1
set.seed(109225004)
sample(1:10,10,replace = T) #4 9 9 3 9 7 5 8 2 5
