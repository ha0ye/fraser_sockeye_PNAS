normalize <- function(block)
{
    if(NCOL(block) > 1)
    {
        n <- NROW(block)
        means <- sapply(block, mean, na.rm = TRUE)
        sds <- sapply(block, sd, na.rm = TRUE)
        return((block - matrix(rep(means, each = n), nrow = n)) / 
                   matrix(rep(sds, each = n), nrow = n))
    }
    else
        return((block - mean(block, na.rm = TRUE)) / sd(block, na.rm = TRUE))
}

normalize_by_cycle_line <- function(ts)
{
    n <- length(ts)
    means <- rep.int(NA, times = 4)
    sds <- rep.int(NA, times = 4)
    mu <- rep.int(NA, times = n)
    sigma <- rep.int(NA, times = n)
    for(k in 1:4)
    {
        index <- seq(from = k, to = n, by = 4)
        means[k] <- mean(ts[index], na.rm = TRUE)
        sds[k] <- sd(ts[index], na.rm = TRUE)
        mu[index] <- means[k]
        sigma[index] <- sds[k]
    }
    ts <- (ts - mu) / sigma
    df <- data.frame(cbind(ts, mu, sigma))
    return(df)
}

compute_stats <- function(obs, pred)
{
    # computes performance metrics for how well predictions match observations
    # obs = vector of observations
    # pred = vector of prediction
    
    N = sum(is.finite(obs) & is.finite(pred))
    rho = cor(obs, pred, use = "pairwise.complete.obs")
    mae = mean(abs(obs-pred), na.rm = TRUE)
    return(data.frame(N = N, rho = rho, mae = mae))
}

rho_comp <- function(x1, x2, y)
    # computes p-value for cor(x1, y) > cor(x2, y) using 
    # t-test with df = length(y) - 2
{
    if(identical(x1, x2))
        return(0.5)
    n <- sum(is.finite(x1) & is.finite(x2) & is.finite(y))
    x1y <- cor(x1, y, use = "pairwise")
    x2y <- cor(x2, y, use = "pairwise")
    err <- TWOpov_err(as.matrix(cbind(x1, x2)), y)
    p_value <- 1 - pt((x1y - x2y) / err, df = n-2, lower.tail = TRUE)
    return(data.frame(df = n-2, statistic = (x1y - x2y) / err, p.value = p_value))
}

# from Wilcox' Robust Statistics R package
TWOpov<-function(x,y,alpha=.05,CN=F)
{
    #
    # Comparing two dependent correlations: Overlapping case
    #
    # x is assumed to be a matrix with 2 columns
    #
    #  Compare correlation of x[,1] with y to x[,2] with y
    #
    if(!is.matrix(x))stop("x should be a matrix")
    if(ncol(x)!=2)stop("x should be a matrix with two columns")
    xy=elimna(cbind(x,y))
    x1=xy[,1]
    x2=xy[,2]
    y=xy[,3]
    r12=cor(x1,y)
    r13=cor(x2,y)
    r23=cor(x1,x2)
    ci12=pcorhc4(x1,y,alpha=alpha,CN=CN)$ci
    ci13=pcorhc4(x2,y,alpha=alpha,CN=CN)$ci
    corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
    term1=2*corhat*(r12-ci12[1])*(ci13[2]-r13)
    term2=2*corhat*(r12-ci12[2])*(ci13[1]-r13)
    L=r12-r13-sqrt((r12-ci12[1])^2+(ci13[2]-r13)^2-term1)
    U=r12-r13+sqrt((r12-ci12[2])^2+(ci13[1]-r13)^2-term2)
    c(L,U)
}

TWOpov_err <- function(x,y,CN=F)
{
    #
    # Comparing two dependent correlations: Overlapping case
    #
    # x is assumed to be a matrix with 2 columns
    #
    #  Compare correlation of x[,1] with y to x[,2] with y
    #
    # returns p-value
    if(!is.matrix(x))stop("x should be a matrix")
    if(ncol(x)!=2)stop("x should be a matrix with two columns")
    xy=elimna(cbind(x,y))
    x1=xy[,1]
    x2=xy[,2]
    y=xy[,3]
    r12=cor(x1,y)
    r13=cor(x2,y)
    r23=cor(x1,x2)
    err12 <- pcorhc4_err(x1,y,CN=CN)
    err13 <- pcorhc4_err(x2,y,CN=CN)
    corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
    err_correction_term = 2*corhat*(err12)*(err13)
    err_diff <- sqrt(err12^2 + err13^2 - err_correction_term)
    return(err_diff)
}

# from Wilcox' Robust Statistics R package
elimna<-function(m)
{
    #
    # remove any rows of data having missing values
    #
    if(is.null(dim(m)))m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    elimna<-m[ikeep[ikeep>=1],]
    elimna
}

pcorhc4_err <- function(x, y, CN = FALSE)
{
    z1 <- (x - mean(x)) / sd(x)
    z2 <- (y - mean(y)) / sd(y)
    ans <- olshc4(z1, z2, alpha = 0.05, CN = CN)
    return(ans$ci[2, 6])
}

# from Wilcox' Robust Statistics R package
pcorhc4<-function(x,y,alpha=.05,CN=F)
{
    #
    #   Compute a .95 confidence interval for Pearson's correlation coefficient.
    #   using the HC4 method
    #
    # CN=F, degrees of freedom are n-p; seems better for general use.
    # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
    #
    xy<-elimna(cbind(x,y))
    x<-xy[,1]
    y<-xy[,2]
    z1=(x-mean(x))/sqrt(var(x))
    z2=(y-mean(y))/sqrt(var(y))
    ans=olshc4(z1,z2,alpha=alpha,CN=CN)
    list(r=ans$r,ci=ans$ci[2,3:4],p.value=ans$ci[2,5])
}

# from Wilcox' Robust Statistics R package
olshc4<-function(x,y,alpha=.05,CN=FALSE,xout=FALSE,outfun=outpro,HC3=FALSE,...)
{
    #
    # Compute confidence for least squares
    # regression using heteroscedastic method
    # recommended by Cribari-Neto (2004).
    # CN=F, degrees of freedom are n-p
    # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
    # All indications are that CN=F is best for general use.
    #
    #  HC3=TRUE, will replace the HC4 estimator with the HC3 estimator.
    #
    x<-as.matrix(x)
    if(nrow(x) != length(y))stop("Length of y does not match number of x values")
    m<-cbind(x,y)
    m<-elimna(m)
    y<-m[,ncol(x)+1]
    x=m[,1:ncol(x)]
    n=length(y)
    nrem=n
    n.keep=length(y)
    x<-as.matrix(x)
    if(xout){
        flag<-outfun(x,...)$keep
        x<-as.matrix(x)
        x<-x[flag,]
        y<-y[flag]
        n.keep=length(y)
        x<-as.matrix(x)
    }
    temp<-lsfit(x,y)
    x<-cbind(rep(1,nrow(x)),x)
    xtx<-solve(t(x)%*%x)
    h<-diag(x%*%xtx%*%t(x))
    n<-length(h)
    d<-(n*h)/sum(h)
    for(i in 1:length(d)){
        d[i]<-min(4, d[i])
    }
    if(HC3)d=2
    hc4<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^d)%*%x%*%xtx
    df<-nrow(x)-ncol(x)
    crit<-qt(1-alpha/2,df)
    if(CN)crit=qnorm(1-alpha/2)
    al<-ncol(x)
    p=al-1
    ci<-matrix(NA,nrow=al,ncol=6)
    lab.out=rep("Slope",p)
    dimnames(ci)<-list(c("(Intercept)",lab.out),c("Coef.","Estimates",
                                                  "ci.lower","ci.upper","p-value","Std.Error"))
    for(j in 1:al){
        ci[j,1]<-j-1
        ci[j,2]<-temp$coef[j]
        ci[j,3]<-temp$coef[j]-crit*sqrt(hc4[j,j])
        ci[j,4]<-temp$coef[j]+crit*sqrt(hc4[j,j])
        test<-temp$coef[j]/sqrt(hc4[j,j])
        ci[j,5]<-2*(1-pt(abs(test),df))
        if(CN)ci[j,5]<-2*(1-pnorm(abs(test),df))
    }
    ci[,6]=sqrt(diag(hc4))
    list(n=nrem,n.keep=n.keep,ci=ci, cov=hc4)
}