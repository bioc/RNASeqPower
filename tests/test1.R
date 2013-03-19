library(RNASeqPower)

# Error message: "missing depth"
#rnapower(n=10, cv=.1, effect=2, alpha=.05, power=.9)

rnapower(10, cv=.5, effect=c(1.5, 1.75, 2), alpha=.05, power=c(.8, .9))
rnapower(10, n=20 , effect=c(1.5, 1.75, 2), alpha=.05, power=c(.8, .9))
rnapower(10, n= c(10,20), cv=.5, alpha=.05, power=c(.8, .9))
rnapower(10, n= c(10,20), cv=.5, effect=1.5, power=c(.8, .9))
rnapower(10, n= c(10,20), cv=.5, effect=c(1.5,2), alpha=.05)


# More formal tests: cv and n
t1 <- rnapower(15, n=20, effect=1.5, alpha=.05, power=.8)
t2 <- rnapower(15,  effect=1.5, alpha=.05, power=.8, cv=t1)
all.equal(t2, 20)


# effect and n
t3 <- rnapower(15, n=20, cv=.5, alpha=.05, power=.8)
t4 <- rnapower(15, cv=.5, alpha=.05, power=.8, effect=t3)
all.equal(t4, 20)

# power and n
t5 <- rnapower(15, n=20, effect=1.5, cv=.5, alpha=.05)
t6 <- rnapower(15, effect=1.5, cv=.5, alpha=.05, power=t5)
all.equal(t6, 20)

# alpha and n
t7 <- rnapower(15, n=20, effect=1.5, cv=.5, power=.8)
t8 <- rnapower(15, effect=1.5, cv=.5, power=.8, alpha=t7)
all.equal(t8, 20)

