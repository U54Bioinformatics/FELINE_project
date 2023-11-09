# Clear previous environment
rm(list = ls())
# time units=seconds
# Load some useful packages
require(deSolve);require(data.table);require(ggplot2)
require(rootSolve)
# Initial conditions 
inits <- c( C=1000, M=1000)

#Parameters: here just the carrying Capacity
par_vals <- c(r_C= 0.001138915, r_M=1e-5, mu_C=1e-3,#1e-10, 
              mu_M=1e-3, 
              beta_C=0.0004188661, beta_M=0.2067994,X=1.04) 




# Write a function defining rate equations of ODE model : here exponential growth
ode_fun <- function(t, y, pars) {
  # specify vectors of parameters (pars) and state variables (y) and make available in function  
  with( as.list( c(pars, y) ), {
    # rate of change
    dC = r_C*((1+beta_C*M)/(1+X))*C*(1-C/1e5) -mu_C*C
    dM = r_M*(1 + beta_M*C)*M*(1-M/1e5) -mu_M*M
    
    #dC = r_C*(1-X)*C*(1-C) - mu_C*C/(1 + beta_C*M)
    #dM = r_M*M*(1-M) - mu_M*M/(1 + beta_M*C)
    
    # return output
    list(c( dC, dM))
  })
}


# Time points to evaluate the ODE
timepoints <- seq(0, 100, length= 100)

#Simulate ODE model using numerical integration to evaluate the evolution of the state at different timepoints
sim <- data.table( 
  ode( y= inits, parms= par_vals, t= timepoints, func= ode_fun )
)

#plot state dynamics over time
ggplot(data=sim , aes(y= C, x= time) ) + geom_path() + theme_classic()
ggplot(data=sim , aes(y= M, x= time) ) + geom_path() + theme_classic()


#steady(y= inits, func=ode_fun, parms= par_vals,time=c(0,Inf),method="runsteady")$y
# steadyStates<-function(p){
#   with( as.list( c(p) ), {
#     Cstar <- (-beta_C*r_M*mu_C + r_M*r_C + beta_C*r_M*r_C + mu_M*beta_M*r_C - r_M*beta_M*r_C - beta_C*r_M*beta_M*r_C + sqrt(4*(1 + beta_C)*r_M*beta_M*r_C*((-mu_M + r_M)*r_C + beta_C*r_M*(-mu_C + r_C)) + ((r_M*(-1 + beta_M) - mu_M*beta_M)*r_C + beta_C*r_M*(mu_C + (-1 + beta_M)*r_C))^2))/(2*(1 + beta_C)*r_M*beta_M*r_C)
#     Mstar <- (-beta_M*r_C*mu_M + r_C*r_M + beta_M*r_C*r_M + mu_C*beta_C*r_M - r_C*beta_C*r_M - beta_M*r_C*beta_C*r_M + sqrt(4*(1 + beta_M)*r_C*beta_C*r_M*((-mu_C + r_C)*r_M + beta_M*r_C*(-mu_M + r_M)) + ((r_C*(-1 + beta_C) - mu_C*beta_M)*r_M + beta_M*r_C*(mu_M + (-1 + beta_C)*r_M))^2))/(2*(1 + beta_M)*r_C*beta_C*r_M)
# return(c("Cstar"=Cstar,"Mstar"=Mstar))
#   })
# }
# steadyStates(par_vals)
# modulate beta_C and evaluate steady state of C and M

resA <- rbindlist(lapply( seq(0, par_vals["beta_C"] , length=100), function(i){
  par_tt<-par_vals
  par_tt["X"]<-0.5
  par_tt["beta_C"]<-i
  
  
  return(data.table(beta_C=i,
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(resA, aes(x=beta_C/max(beta_C), y=abs(C), col=(abs(M)))) + 
  geom_path(size=2.5)+
  geom_point(size=1.7)+ 
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  labs(y="Cancer \n abundance",x="Cancer growth facilitation by fibroblasts \n (relative to maximum)")+
  #scale_color_continuous(name="Cancer abundance")+
  scale_color_continuous(name="Fibroblast \n abundance",
                         #breaks=c(0,0.1,0.2,0.3),
                         high=ggsci::pal_aaas()(2)[1], 
                         low=ggsci::pal_aaas()(2)[2])#
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/Cancer equilibrium abundance across fibroblast facilitation levels.pdf", width=8, height=8,dpi=320)
ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/Finalize Cancer equilibrium abundance across fibroblast facilitation levels.pdf", width=8, height=8,dpi=320)


res <- rbindlist(lapply(seq(0,par_vals["beta_C"],length=30), function(i){
  par_tt<-par_vals
  par_tt["X"]<-i
  return(data.table(X=i,
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(x=X,y=C,col=M)) + geom_path() +theme_classic()


lu <- data.table(expand.grid(beta_C=seq(0, par_vals["beta_C"], length=100),
                           X =seq(0, 1.04, length=100) ))

res<-rbindlist(lapply(1:nrow(lu), function(i){
  cat(i)
  par_tt<-par_vals
  par_tt["X"]<-lu[i]$X
  par_tt["beta_C"]<-lu[i]$beta_C
  return(data.table(lu[i],
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))

ggplot(res[], aes(x=(X), y=(beta_C)/max(beta_C), fill=round(abs(C),5))) + 
  geom_raster(interpolate=T) +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  labs(x="Endocrine treatment",y="Cancer growth facilitation \n by fibroblasts \n (relative to maximum)")+
  #scale_color_continuous(name="Cancer abundance")+
  scale_fill_gradient(name="Cancer \n abundance",
                      #breaks=c(0,0.2,0.4,0.6,0.8),
                      high=ggsci::pal_aaas()(2)[2], 
                      low="white")
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/coculture/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/FINALIZE Cancer equilibrium abundance across fulv doses and fibroblast facilitation levels.pdf", width=8, height=8,dpi=320)

ggplot(res[#X<1][beta_C<2.5
          ], aes(x=(X), y=(beta_C), fill=1+round(abs(C),5))) + 
  geom_raster() +
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  labs(x="Endocrine treatment",y="Cancer growth facilitation \n by fibroblasts")+
  #scale_color_continuous(name="Cancer abundance")+
  scale_fill_gradient(name="Cancer \n abundance",
                      breaks=c(0,0.2,0.4,0.6,0.8),
                      high=ggsci::pal_aaas()(2)[2], 
                      low="floralwhite")#white")#
#low=ggsci::pal_aaas()(2)[1])
#ggsave(file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/R Analysis/Cancer_Fibroblast_Experiments/20230305-13_Feng  Fulvestrant dose CAMA1 MCF7 T47D/Modeling_JG/Cancer equilibrium abundance across fulv doses and fibroblast facilitation levels.pdf", width=8, height=8,dpi=320)

write.csv( resA, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5e/Cancer facilitation ODE leftpanel.csv")

write.csv( res, file="/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/FELINE2/Cancer GF analyses/SourceData/Figure5/Fig5e/Cancer facilitation ODE rightpanel.csv")








# Clear previous environment
rm(list = ls())
# time units=seconds
# Load some useful packages
require(deSolve);require(data.table);require(ggplot2)
require(rootSolve)
# Initial conditions 
inits <- c( C=0.01,M=0.01)

#Parameters: here just the carrying Capacity
par_vals <- c(r_C=0.0011, r_M=-0.0062, mu_C=0.01, mu_M=0.01, 
              beta_C=10, beta_M=10,X=1) 




# Write a function defining rate equations of ODE model : here exponential growth
ode_fun <- function(t, y, pars) {
  # specify vectors of parameters (pars) and state variables (y) and make available in function  
  with( as.list( c(pars, y) ), {
    
    # rate of change
    #dC = r_C*(1/(1+X))*C*(1-C) - mu_C*C/(1+beta_C*M)
    dC = r_C*(1-X)*C*(1-C) - mu_C*C/(1+beta_C*M)
    dM = r_M*M*(1-M) - mu_M*M/(1+beta_M*C)
    
    # return output
    list(c( dC, dM))
  })
}


# Time points to evaluate the ODE
timepoints <- seq(0, 1000, length= 100)

#Simulate ODE model using numerical integration to evaluate the evolution of the state at different timepoints
sim <- data.table( 
  ode( y= inits, parms= par_vals, t= timepoints, func= ode_fun )
)

#plot state dynamics over time
ggplot(data=sim , aes(y= C, x= time) ) + geom_path() + theme_classic()
ggplot(data=sim , aes(y= M, x= time) ) + geom_path() + theme_classic()


steady(y= inits, func=ode_fun, parms= par_vals,time=c(0,Inf),method="runsteady")$y
# steadyStates<-function(p){
#   with( as.list( c(p) ), {
#     Cstar <- (-beta_C*r_M*mu_C + r_M*r_C + beta_C*r_M*r_C + mu_M*beta_M*r_C - r_M*beta_M*r_C - beta_C*r_M*beta_M*r_C + sqrt(4*(1 + beta_C)*r_M*beta_M*r_C*((-mu_M + r_M)*r_C + beta_C*r_M*(-mu_C + r_C)) + ((r_M*(-1 + beta_M) - mu_M*beta_M)*r_C + beta_C*r_M*(mu_C + (-1 + beta_M)*r_C))^2))/(2*(1 + beta_C)*r_M*beta_M*r_C)
#     Mstar <- (-beta_M*r_C*mu_M + r_C*r_M + beta_M*r_C*r_M + mu_C*beta_C*r_M - r_C*beta_C*r_M - beta_M*r_C*beta_C*r_M + sqrt(4*(1 + beta_M)*r_C*beta_C*r_M*((-mu_C + r_C)*r_M + beta_M*r_C*(-mu_M + r_M)) + ((r_C*(-1 + beta_C) - mu_C*beta_M)*r_M + beta_M*r_C*(mu_M + (-1 + beta_C)*r_M))^2))/(2*(1 + beta_M)*r_C*beta_C*r_M)
# return(c("Cstar"=Cstar,"Mstar"=Mstar))
#   })
# }
# steadyStates(par_vals)
# modulate beta_C and evaluate steady state of C and M
res<-rbindlist(lapply(seq(0,0.5,by=0.01), function(i){
  par_tt<-par_vals
  par_tt["beta_C"]<-i
  return(data.table(beta_C=i,
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(C,M)) + geom_point() 
ggplot(res, aes(x=beta_C,y=C,col=M)) + geom_path() +theme_classic()


res <- rbindlist(lapply(seq(0,0.1,by=0.01), function(i){
  par_tt<-par_vals
  par_tt["r_C"]<-i
  return(data.table(beta_C=i,
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(x=beta_C,y=C,col=M)) + geom_path() +theme_classic()

res <- rbindlist(lapply(seq(0,0.02,by=0.001), function(i){
  par_tt<-par_vals
  par_tt["mu_C"]<-i
  return(data.table(beta_C=i,
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(x=beta_C,y=C,col=M)) + geom_path() +theme_classic()


lu<-data.table(expand.grid(beta_C=seq(0,0.5,length=10),
                           mu_C =seq(0,0.02,length=30) ))

res<-rbindlist(lapply(1:nrow(lu), function(i){
  par_tt<-par_vals
  par_tt["mu_C"]<-lu[i]$mu_C
  par_tt["beta_C"]<-lu[i]$beta_C
  return(data.table(lu[i],
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(C,M)) + geom_point() 
ggplot(res, aes(x=mu_C, y=(C), col=beta_C, group=interaction(beta_C))) + 
  geom_path() +theme_classic()+
  theme(aspect.ratio=1)+
  labs(x="Chemotherapy treatment",y="Cancer abundance")+
  scale_color_continuous(name="Cancer growth facilitation \n by fibroblasts")

lu<-data.table(expand.grid(beta_C=seq(0,1,length=10),
                           X =seq(0,0.1,length=30) ))

res<-rbindlist(lapply(1:nrow(lu), function(i){
  par_tt<-par_vals
  par_tt["X"]<-lu[i]$X
  par_tt["beta_C"]<-lu[i]$beta_C
  return(data.table(lu[i],
                    t(steady(y= inits, func=ode_fun, parms= par_tt,time=c(0,Inf),method="runsteady")$y )))
}))
ggplot(res, aes(x=(X), y=(C), col=beta_C, group=interaction(beta_C))) + 
  geom_path() +theme_classic()+
  theme(aspect.ratio=1)+
  labs(x="Endocrine treatment",y="Cancer abundance")+
  scale_color_continuous(name="Cancer growth facilitation \n by fibroblasts")

ggplot(res[beta_C==max(beta_C)], aes(x=(X), y=(C), col=beta_C, group=interaction(beta_C))) + 
  geom_path() +theme_classic()+
  theme(aspect.ratio=1)+
  labs(x="Endocrine treatment",y="Cancer abundance")+
  scale_color_continuous(name="Cancer growth facilitation \n by fibroblasts")

ggplot(res, aes(x=(X), y=(beta_C), fill=sqrt(abs(C)))) + 
  geom_tile() +theme_classic()+
  theme(aspect.ratio=1)+
  labs(x="Endocrine treatment",y="Cancer abundance")+
  scale_color_continuous(name="Cancer growth facilitation \n by fibroblasts")











ggplot(res, aes(x=beta_C,y=mu_C,fill=C)) + geom_raster(interpolate=T) +theme_classic()

ggplot(res, aes(x=beta_C,y=mu_C,fill=C)) + geom_raster(interpolate=T) +theme_classic()











par_tt<-par_vals
par_tt
Cstar <- with( as.list( c(par_vals, inits) ), {
  (-beta_C*r_M*mu_C + r_M*r_C + beta_C*r_M*r_C + mu_M*beta_M*r_C - r_M*beta_M*r_C - beta_C*r_M*beta_M*r_C + sqrt(4*(1 + beta_C)*r_M*beta_M*r_C*((-mu_M + r_M)*r_C + beta_C*r_M*(-mu_C + r_C)) + ((r_M*(-1 + beta_M) - mu_M*beta_M)*r_C + beta_C*r_M*(mu_C + (-1 + beta_M)*r_C))^2))/(2*(1 + beta_C)*r_M*beta_M*r_C)
})
Mstar <- with( as.list( c(par_vals, inits) ), {
  (-beta_M*r_C*mu_M + r_C*r_M + beta_M*r_C*r_M + mu_C*beta_C*r_M - r_C*beta_C*r_M - beta_M*r_C*beta_C*r_M + sqrt(4*(1 + beta_M)*r_C*beta_C*r_M*((-mu_C + r_C)*r_M + beta_M*r_C*(-mu_M + r_M)) + ((r_C*(-1 + beta_C) - mu_C*beta_M)*r_M + beta_M*r_C*(mu_M + (-1 + beta_C)*r_M))^2))/(2*(1 + beta_M)*r_C*beta_C*r_M)
})


Cstar= 1- ( mu_C/(r_C*(1+beta_C*Mstar)) )
Mstar= 1- ( mu_M/(r_M*(1+beta_M*Cstar)) )

})

Cstar= 1- ( mu_C/(r_C*(1+beta_C* (1- ( mu_M/(r_M*(1+beta_M*Cstar)) ))    )) )

b=beta_C
β=beta_M

r=r_M
ρ=r_C

μ=mu_C
m=mu_M



