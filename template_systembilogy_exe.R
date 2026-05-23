rm(list=ls()) # Clear the memory
library(deSolve)
############# ODE function #########
my_ode<-function(t,state,parms){
  with(as.list(state),{
    dndt=rep(0,length(state))
    #-------My Equations----------------
    dndt[1] = -a*M
    dndt[2] = a*M-k*N
    #-------------------------------------
    return(list(dndt)) # Return
  })
}
############ END of function ##########
M=400
N=0 # Initial value = 100 μM
a=0.100
init=c(M=M,N=N) # Create a vector with initial values 
k=0.01 # Exponential decay constant (/h)
t=seq(0,300,1) # Run for 80 time steps
out <- ode(y =init, times = t, func = my_ode, parms = NULL)
tail(out) # Prints the first time points
plot(out,type="l",xlab="Time (h)",ylab="Protein (μM)")
time=out[,1]
n=out[,2]
m=out[,3]
plot(time,m,col="red",type="l",lwd=3)
lines(time,n,col="blue",lwd=3)
legend("topright",c("Drug in gut","Drug in blood"),col=c("red","blue"),lty =1,lwd=3)
out[301,]
