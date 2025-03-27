#*****************************************************************
# PROJECT 3: HOPFIELD NETWORK - by Josh Goh 16 Feb 2007
#*****************************************************************

# Specify output file paths
# Simulation reports
filepath=sprintf("%s/HopfieldNetworkoutput_joshgoh.txt",getwd())
outputfile<-file(filepath,"a")

# Data for graphing
datafilepath=sprintf("%s/HopfieldNetworkdata_joshgoh.txt",getwd())
datafile<-file(datafilepath,"a")

# Energy for each simulation test (requested in the assignment 1.B.2)
efilepath=sprintf("%s/HopfieldNetworkenergy_joshgoh.txt",getwd())
efile<-file(efilepath,"a")


#-----------------------------------------------------------------
# READ FUNCTIONS INTO WORKSPACE
#-----------------------------------------------------------------

# Function 1: trainnet_hopfield
#-----------------------------------------------------------------
# Trains a weight matrix, W, based on a series of 
# training input vectors, a, using the Hopfield 
# learning rule W=(2a-1)*t(2a-1).
#
# Usage:
# W<-trainnet_hopfield(W,a)
#
# W - (input)initial weight matrix to be trained (must be square 
#     matrix of same dimension as input vector, a.
#   - if W is set to scalar, then an initial W of the scalar value
#     is created based on input vector, a, dimensions.
# a - matrix of training vectors with the units to be
#     train in rows, and the training pattern for each 
#     unit in columns.
#
# Note: So far this only does 2d weight matrices.
# Created by Josh Goh 13 Feb 2007

trainnet_hopfield<-function(W,a) {
	
	# Initialize weight matrix, W
	if (length(W)==1) {
	W<-matrix(W,dim(a)[1],dim(a)[1])
	}
	
	# Train W using training vectors
	for (i in 1:dim(a)[2]) {
		v<-matrix(a[,i],dim(a)[1],1)  # Current training vector
		dw<-(2*v-1)%*%t(2*v-1)        # Compute change in weights
		W<-W+dw			      # Update weight matrix
		}

	# Zero the diagonals (units are not connected to themselves).
	W<-W-((dim(a)[2])*diag(dim(a)[1]))
		
	# Return final weight matrix
	W

	}
# End function 1

# Function 2: testnet_hopfield
#--------------------------------------------------------------
# Tests a Hopfield trained weight matrix, W, with an
# input vector pattern, a, for t iterations of consecutive
# zero hamming distances between each activation state.
# Returns a matrix of unit activation changes, a, a vector
# of network energy level changes, e, and a vector of
# hamming distances changes, hd. Flags f=1 if network fails
# to settle, f=0 otherwise.
#
# Usage:
# testnet_hopfield(W,a,t)
#
# W - trained hopfield weight matrix to be tested
# a - test vector pattern
# t - number of iterations of consecutive zero hamming
#     distances desired
#
# Created by Josh Goh 13 Feb 2007

testnet_hopfield<-function(W,a,t) {

	# Initialize parameters
	u<-0    # activation threshold
	e<-c(0)  # energy zero start vector
	hd<-c(0) # hamming distance zero start vector
	s<-1    # consecutive no-hd-change counter
	r<-1    # run counter
	f<-0	# fail to settle flag, init 0
	
	# Start testing
	while (s<=t) {

		# Randomly select a unit (asynchronous)
		atemp<-matrix(a[,r],dim(a)[1],1)
		i<-round(runif(1,1,dim(a)[1]),0)
	
		# Update unit i using activation function
		if (W[i,]%*%a[,r]>u) {atemp[i,1]<-1}
		else {atemp[i,1]<-0}

		# Append activation
		a<-cbind(a,atemp)

		# Calculate energy
		e<-cbind(e,(-0.5)*(t(W%*%atemp)%*%a[,r]))

		# Calculate hamming distance
		hdtemp<-0
		for (n in 1:dim(a)[1]) {
			if (atemp[n]!=a[n,r]) {
				hdtemp<-hdtemp+1
				}
			}
		hd<-cbind(hd,hdtemp)

		# Reset change counter if hd changed
		if (hdtemp==0) {s<-s+1}
		else {
			if (r<1000) {s<-1}
			
			# Flag if fail to settle after 1000 runs
			else {f<-1}
			}
		
		# Update run counter every time
		r<-r+1
		
		}
		
	# Format and assign variables to global environment
	assign("a",a[,2:dim(a)[2]],envir=.GlobalEnv)
	
	e<-e[2:length(e)]
	dim(e)<-c(1,length(e))
	assign("e",e,envir=.GlobalEnv)
	
	hd<-hd[2:length(hd)]
	dim(hd)<-c(1,length(hd))
	assign("hd",hd,envir=.GlobalEnv)
	
	assign("f",f,envir=.GlobalEnv)

	}
# End function 2
# End read functions


#-----------------------------------------------------------------
# BEGIN WORKSPACE COMPUTATIONS FOR PROJECT SIMULATIONS
#-----------------------------------------------------------------

# Set simulation parameters
t<-30   # Number of iterations to use
runs<-20 # Number of runs when testing

# PART 1: WALSH FUNCTIONS
#-----------------------------------------------------------------

# 1A: Training patterns
# Create training patterns
tp<-c()
tp<-cbind(tp,matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16,1))
tp<-cbind(tp,matrix(c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0),16,1))
tp<-cbind(tp,matrix(c(1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0),16,1))
tp<-cbind(tp,matrix(c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0),16,1))

# Train network with training patterns
W<-trainnet_hopfield(0,tp)

# Print trained weights to data file
cat("PART 1A: Walsh-trained weight matrix\n",file=datafilepath,append=TRUE)
write.table(W,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# 1B: Test with distorted patterns
cat("PART 1B: Test Walsh-trained network with distorted training patterns\n",file=filepath,append=TRUE)

# Set test pattern parameters
	prob<-c(0,0.1,0.2,0.3,0.4,0.5) # List distortion probabilities

	# Initiate a zero test pattern matrix with activation placeholders for
	# each unit, training pattern, distortion probability, and run
	testp<-c(rep(0,(dim(tp)[1]*dim(tp)[2]*length(prob)*runs)))
	dim(testp)<-c(dim(tp)[1],dim(tp)[2],length(prob),runs)

# Create test patterns
for (p in 1:length(prob)) {									  # Distortion probability counter
	pdist<-c(rep(1,(100*prob[p])),rep(0,(100-(100*prob[p])))) # Create probability distribution
	for (pat in 1:dim(tp)[2]) {								  # Training pattern counter
		for (patv in 1:dim(tp)[1]) {						  # Unit counter
			for (nrun in 1:runs) {							  # Run counter
				pn<-sample(pdist,1) # Randomly select from distribution
				
				# Update test activation state
				if (pn==1) {
					if (tp[patv,pat]==1) {
						testp[patv,pat,p,nrun]<-0
					}
					else {
						testp[patv,pat,p,nrun]<-1
					}
				}
				else {
					testp[patv,pat,p,nrun]<-tp[patv,pat]
				} # End update state
				
			} # End run counter
		} # End unit counter
		
	#Print test pattern to data file
	cat(sprintf("1B: Training pattern %.0f with distortion probability: %.1f\n",pat,prob[p]),file=datafilepath,append=TRUE)
	write.table(testp[,pat,p,],sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
	cat("\n",file=datafilepath,append=TRUE)
		
	} # End training pattern counter
} # End distortion probability counter
	
# Initiate required report variables
ihdf<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ihdf)<-c(length(prob),dim(testp)[2],runs)	# Zero matrix for hamming distance

ni<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ni)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of iterations to settle

fi<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(fi)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of failures to settle


# Perform tests
for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Test each 'runs' times
		finala<-c()
		# Print header for energy state per iteration report
		cat(sprintf("\n1B.2: Energy for each iteration of each run for distortion probability %.1f of training pattern %.0f\n",prob[p],pat),file=efilepath,append=TRUE)
		
		for (r in 1:runs) {
			testnet_hopfield(W,matrix(testp[,pat,p,r],16,1),t)
			
			# Store final activations	
			finala<-cbind(finala,a[,dim(a)[2]])
			
			# Calculate hamming distance of final activation to original
			hdf<-0
			c<-dim(a)[2]
			for (n in 1:16) {
				if (a[n,c]!=tp[n,pat]) {
				hdf<-hdf+1
				}
			}
				
			# Store hamming distance of final activation to original
			ihdf[p,pat,r]<-hdf
			
			# Store number of iterations to settle
			ni[p,pat,r]<-dim(a)[2]
			
			# Store flag if failure to settle
			fi[p,pat,r]<-f
			
			# Print energy
			write.table(e,sep="\t",file=efilepath,append=TRUE,row.names=r,col.names=FALSE)
			#cat("\n",file=datafilepath,append=TRUE)
			
		}
		
		# Print final activation states
		cat(sprintf("1B.3: Final activation of units for distortion probability %.1f of training pattern %.0f, %.0f runs\n",prob[p],pat,runs),file=datafilepath,append=TRUE)
		write.table(finala,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
		cat("\n",file=datafilepath,append=TRUE)
	}
}

# Organize and calculate mean and standard deviations for hamming distance, iteration to settle, and failure to settle reports

mhd<-c(rep(0,length(prob)*dim(testp)[2])) # Mean hamming distance matrix
dim(mhd)<-c(length(prob),dim(testp)[2])

sdhd<-c(rep(0,length(prob)*dim(testp)[2])) # Standard deviation hamming distance matrix
dim(sdhd)<-c(length(prob),dim(testp)[2])

mi<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean iterations to settle matrix
dim(mi)<-c(length(prob),dim(testp)[2])

sdi<-c(rep(0,length(prob)*dim(testp)[2]))  # Standard deviation iterations to settle matrix
dim(sdi)<-c(length(prob),dim(testp)[2])

nf<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean failure to settle matrix
dim(nf)<-c(length(prob),dim(testp)[2])

for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Means
		mhd[p,pat]<-round(mean(ihdf[p,pat,]),1)
		mi[p,pat]<-round(mean(ni[p,pat,]),1)
		
		# Standard deviations
		sdhd[p,pat]<-round(sd(ihdf[p,pat,]),2)
		sdi[p,pat]<-round(sd(ni[p,pat,]),2)
		
		# Totals
		nf[p,pat]<-sum(fi[p,pat,])
	}
}

# Print reports
cat(sprintf("\nSimulation report for %.0f distortion probabilities (rows) of %.0f training patterns (col) across %.0f runs:\n",length(prob),dim(testp)[2],runs),file=filepath,append=TRUE)

cat("Mean hamming distance\n",file=filepath,append=TRUE)
write.table(mhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of hamming distance\n",file=filepath,append=TRUE)
write.table(sdhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Mean number of iterations to settle\n",file=filepath,append=TRUE)
write.table(mi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of number of iterations to settle\n",file=filepath,append=TRUE)
write.table(sdi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Number of failures to settle\n",file=filepath,append=TRUE)
write.table(nf,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)


#-----------------------------------------------------------------

# PART 2: REPETITION
#-----------------------------------------------------------------

# 2A: Training patterns
# Create training patterns
	tp<-c()
	tp<-cbind(tp,matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16,1)) # Pattern A
	tp<-cbind(tp,matrix(c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0),16,1)) # Pattern B
	tp<-cbind(tp,matrix(c(1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0),16,1)) # Pattern C

# Train network with training patterns
	
	W=0 # Start with zero weights
	
	# Train pattern A 20 times
	for (rp in 1:20) {
		W<-trainnet_hopfield(W,matrix(tp[,1],dim(tp)[1],1))
	}
	
	# Train pattern B 10 times
	for (rp in 1:10) {
		W<-trainnet_hopfield(W,matrix(tp[,2],dim(tp)[1],1))
	}
	
	# Train pattern C once
		W<-trainnet_hopfield(W,matrix(tp[,3],dim(tp)[1],1))

# Print trained weights to data file
cat("PART 2A: Repetition-trained weight matrix\n",file=datafilepath,append=TRUE)
write.table(W,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# 2B: Test with distorted patterns
cat("PART 2B: Test Repetition-trained network with distorted training patterns\n",file=filepath,append=TRUE)

# Set test pattern parameters
	prob<-c(0,0.1,0.2,0.3,0.4,0.5) # List distortion probabilities

	# Initiate a zero test pattern matrix with activation placeholders for
	# each unit, training pattern, distortion probability, and run
	testp<-c(rep(0,(dim(tp)[1]*dim(tp)[2]*length(prob)*runs)))
	dim(testp)<-c(dim(tp)[1],dim(tp)[2],length(prob),runs)

# Create test patterns
for (p in 1:length(prob)) {									  # Distortion probability counter
	pdist<-c(rep(1,(100*prob[p])),rep(0,(100-(100*prob[p])))) # Create probability distribution
	for (pat in 1:dim(tp)[2]) {								  # Training pattern counter
		for (patv in 1:dim(tp)[1]) {						  # Unit counter
			for (nrun in 1:runs) {							  # Run counter
				pn<-sample(pdist,1) # Randomly select from distribution
				
				# Update test activation state
				if (pn==1) {
					if (tp[patv,pat]==1) {
						testp[patv,pat,p,nrun]<-0
					}
					else {
						testp[patv,pat,p,nrun]<-1
					}
				}
				else {
					testp[patv,pat,p,nrun]<-tp[patv,pat]
				} # End update state
				
			} # End run counter
		} # End unit counter
		
	#Print test pattern to data file
	cat(sprintf("2B: Training pattern %.0f with distortion probability: %.1f\n",pat,prob[p]),file=datafilepath,append=TRUE)
	write.table(testp[,pat,p,],sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
	cat("\n",file=datafilepath,append=TRUE)
		
	} # End training pattern counter
} # End distortion probability counter
	
# Initiate required report variables
ihdf<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ihdf)<-c(length(prob),dim(testp)[2],runs)	# Zero matrix for hamming distance

ni<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ni)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of iterations to settle

fi<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(fi)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of failures to settle


# Perform tests
for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Test each 'runs' times
		finala<-c()
		# Print header for energy state per iteration report
		cat(sprintf("\n2B.2: Energy for each iteration of each run for distortion probability %.1f of training pattern %.0f\n",prob[p],pat),file=efilepath,append=TRUE)
		
		for (r in 1:runs) {
			testnet_hopfield(W,matrix(testp[,pat,p,r],16,1),t)
			
			# Store final activations	
			finala<-cbind(finala,a[,dim(a)[2]])
			
			# Calculate hamming distance of final activation to original
			hdf<-0
			c<-dim(a)[2]
			for (n in 1:16) {
				if (a[n,c]!=tp[n,pat]) {
				hdf<-hdf+1
				}
			}
				
			# Store hamming distance of final activation to original
			ihdf[p,pat,r]<-hdf
			
			# Store number of iterations to settle
			ni[p,pat,r]<-dim(a)[2]
			
			# Store flag if failure to settle
			fi[p,pat,r]<-f
			
			# Print energy
			write.table(e,sep="\t",file=efilepath,append=TRUE,row.names=r,col.names=FALSE)
			#cat("\n",file=datafilepath,append=TRUE)
			
		}
		
		# Print final activation states
		cat(sprintf("2B.3: Final activation of units for distortion probability %.1f of training pattern %.0f, %.0f runs\n",prob[p],pat,runs),file=datafilepath,append=TRUE)
		write.table(finala,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
		cat("\n",file=datafilepath,append=TRUE)
	}
}

# Organize and calculate mean and standard deviations for hamming distance, iteration to settle, and failure to settle reports

mhd<-c(rep(0,length(prob)*dim(testp)[2])) # Mean hamming distance matrix
dim(mhd)<-c(length(prob),dim(testp)[2])

sdhd<-c(rep(0,length(prob)*dim(testp)[2])) # Standard deviation hamming distance matrix
dim(sdhd)<-c(length(prob),dim(testp)[2])

mi<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean iterations to settle matrix
dim(mi)<-c(length(prob),dim(testp)[2])

sdi<-c(rep(0,length(prob)*dim(testp)[2]))  # Standard deviation iterations to settle matrix
dim(sdi)<-c(length(prob),dim(testp)[2])

nf<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean failure to settle matrix
dim(nf)<-c(length(prob),dim(testp)[2])

for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Means
		mhd[p,pat]<-round(mean(ihdf[p,pat,]),1)
		mi[p,pat]<-round(mean(ni[p,pat,]),1)
		
		# Standard deviations
		sdhd[p,pat]<-round(sd(ihdf[p,pat,]),2)
		sdi[p,pat]<-round(sd(ni[p,pat,]),2)
		
		# Totals
		nf[p,pat]<-sum(fi[p,pat,])
	}
}

# Print reports
cat(sprintf("\nSimulation report for %.0f distortion probabilities (rows) of %.0f training patterns (col) across %.0f runs:\n",length(prob),dim(testp)[2],runs),file=filepath,append=TRUE)

cat("Mean hamming distance\n",file=filepath,append=TRUE)
write.table(mhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of hamming distance\n",file=filepath,append=TRUE)
write.table(sdhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Mean number of iterations to settle\n",file=filepath,append=TRUE)
write.table(mi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of number of iterations to settle\n",file=filepath,append=TRUE)
write.table(sdi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Number of failures to settle\n",file=filepath,append=TRUE)
write.table(nf,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

#-----------------------------------------------------------------

# PART 3: RANDOM PATTERNS
#-----------------------------------------------------------------

# 3A: Training patterns
# Create training patterns
	
	# Initiate random parameters
	randt<-3 # Number of random training patterns to create
	tp<-matrix(0,16,randt)
	pdist<-c(rep(1,50),rep(0,50)) # 0.5 probability distribution of 1s and 0s
	
	# Randomly assign states to training patterns
	for (i in 1:randt) {
		for (j in 1:16) {
			# Randomly sample and assign unit value from probability distribution
			tp[j,i]<-sample(pdist,1)
		}
	}

# Print random training patterns to data file
cat("PART 3A: Random training patterns\n",file=datafilepath,append=TRUE)
write.table(tp,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# Train network with training patterns
W<-trainnet_hopfield(0,tp)

# Print trained weights to data file
cat("PART 3A: Random-trained weight matrix\n",file=datafilepath,append=TRUE)
write.table(W,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# 3B: Test with distorted patterns
cat("PART 3Ba: Test Random-trained network with distorted test patterns\n",file=filepath,append=TRUE)

# Set test pattern parameters
	prob<-c(0,0.1,0.2,0.3,0.4,0.5) # List distortion probabilities

	# Initiate a zero test pattern matrix with activation placeholders for
	# each unit, training pattern, distortion probability, and run
	testp<-c(rep(0,(dim(tp)[1]*dim(tp)[2]*length(prob)*runs)))
	dim(testp)<-c(dim(tp)[1],dim(tp)[2],length(prob),runs)

# Create test patterns
for (p in 1:length(prob)) {									  # Distortion probability counter
	pdist<-c(rep(1,(100*prob[p])),rep(0,(100-(100*prob[p])))) # Create probability distribution
	for (pat in 1:dim(tp)[2]) {								  # Training pattern counter
		for (patv in 1:dim(tp)[1]) {						  # Unit counter
			for (nrun in 1:runs) {							  # Run counter
				pn<-sample(pdist,1) # Randomly select from distribution
				
				# Update test activation state
				if (pn==1) {
					if (tp[patv,pat]==1) {
						testp[patv,pat,p,nrun]<-0
					}
					else {
						testp[patv,pat,p,nrun]<-1
					}
				}
				else {
					testp[patv,pat,p,nrun]<-tp[patv,pat]
				} # End update state
				
			} # End run counter
		} # End unit counter
		
	#Print test pattern to data file
	cat(sprintf("3Ba: Training pattern %.0f with distortion probability: %.1f\n",pat,prob[p]),file=datafilepath,append=TRUE)
	write.table(testp[,pat,p,],sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
	cat("\n",file=datafilepath,append=TRUE)
		
	} # End training pattern counter
} # End distortion probability counter
	
# Initiate required report variables
ihdf<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ihdf)<-c(length(prob),dim(testp)[2],runs)	# Zero matrix for hamming distance

ni<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(ni)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of iterations to settle

fi<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
dim(fi)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of failures to settle


# Perform tests
for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Test each 'runs' times
		finala<-c()
		# Print header for energy state per iteration report
		cat(sprintf("\n3Ba.2: Energy for each iteration of each run for distortion probability %.1f of training pattern %.0f\n",prob[p],pat),file=efilepath,append=TRUE)
		
		for (r in 1:runs) {
			testnet_hopfield(W,matrix(testp[,pat,p,r],16,1),t)
			
			# Store final activations	
			finala<-cbind(finala,a[,dim(a)[2]])
			
			# Calculate hamming distance of final activation to original
			hdf<-0
			c<-dim(a)[2]
			for (n in 1:16) {
				if (a[n,c]!=tp[n,pat]) {
				hdf<-hdf+1
				}
			}
				
			# Store hamming distance of final activation to original
			ihdf[p,pat,r]<-hdf
			
			# Store number of iterations to settle
			ni[p,pat,r]<-dim(a)[2]
			
			# Store flag if failure to settle
			fi[p,pat,r]<-f
			
			# Print energy
			write.table(e,sep="\t",file=efilepath,append=TRUE,row.names=r,col.names=FALSE)
			#cat("\n",file=datafilepath,append=TRUE)
			
		}
		
		# Print final activation states
		cat(sprintf("3Ba.3: Final activation of units for distortion probability %.1f of training pattern %.0f, %.0f runs\n",prob[p],pat,runs),file=datafilepath,append=TRUE)
		write.table(finala,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
		cat("\n",file=datafilepath,append=TRUE)
	}
}

# Organize and calculate mean and standard deviations for hamming distance, iteration to settle, and failure to settle reports

mhd<-c(rep(0,length(prob)*dim(testp)[2])) # Mean hamming distance matrix
dim(mhd)<-c(length(prob),dim(testp)[2])

sdhd<-c(rep(0,length(prob)*dim(testp)[2])) # Standard deviation hamming distance matrix
dim(sdhd)<-c(length(prob),dim(testp)[2])

mi<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean iterations to settle matrix
dim(mi)<-c(length(prob),dim(testp)[2])

sdi<-c(rep(0,length(prob)*dim(testp)[2]))  # Standard deviation iterations to settle matrix
dim(sdi)<-c(length(prob),dim(testp)[2])

nf<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean failure to settle matrix
dim(nf)<-c(length(prob),dim(testp)[2])

for (p in 1:length(prob)) {
	for (pat in 1: dim(testp)[2]) {
		# Means
		mhd[p,pat]<-round(mean(ihdf[p,pat,]),1)
		mi[p,pat]<-round(mean(ni[p,pat,]),1)
		
		# Standard deviations
		sdhd[p,pat]<-round(sd(ihdf[p,pat,]),2)
		sdi[p,pat]<-round(sd(ni[p,pat,]),2)
		
		# Totals
		nf[p,pat]<-sum(fi[p,pat,])
	}
}

# Print reports
cat(sprintf("\nSimulation report for %.0f distortion probabilities (rows) of %.0f training patterns (col) across %.0f runs:\n",length(prob),dim(testp)[2],runs),file=filepath,append=TRUE)

cat("Mean hamming distance\n",file=filepath,append=TRUE)
write.table(mhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of hamming distance\n",file=filepath,append=TRUE)
write.table(sdhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Mean number of iterations to settle\n",file=filepath,append=TRUE)
write.table(mi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of number of iterations to settle\n",file=filepath,append=TRUE)
write.table(sdi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Number of failures to settle\n",file=filepath,append=TRUE)
write.table(nf,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)


# 3B: Additional training patterns
# Create 4 more random patterns and train one at a time
for (k in 1:4) { 

	# Create training pattern
	# Initiate random parameters
	randt<-1 # Number of random training patterns to create
	tp<-cbind(tp,matrix(0,16,randt))
	pdist<-c(rep(1,50),rep(0,50)) # 0.5 probability distribution of 1s and 0s

	for (j in 1:16) {
		# Randomly sample and assign unit value from probability distribution
		tp[j,dim(tp)[2]]<-sample(pdist,1)
	}

	cat(sprintf("PART 3Bb: Additional random training pattern %.0f\n",k),file=datafilepath,append=TRUE)
	write.table(tp,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
	cat("\n",file=datafilepath,append=TRUE)

	# Train network with additional pattern
	W<-trainnet_hopfield(W,matrix(tp[,dim(tp)[2]],dim(tp)[1],1))
	
	cat(sprintf("PART 3Bb: Random-trained weight matrix for additional pattern %.0f\n",k),file=datafilepath,append=TRUE)
	write.table(W,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
	cat("\n",file=datafilepath,append=TRUE)

	# 3B: Test with distorted patterns
	cat(sprintf("PART 3Bb: Test network additionally trained on %.0f random pattern(s) with distorted training patterns\n",k),file=filepath,append=TRUE)

	# Set test pattern parameters
	prob<-c(0.2) # List distortion probabilities

	# Initiate a zero test pattern matrix with activation placeholders for
	# each unit, training pattern, distortion probability, and run
	testp<-c(rep(0,(dim(tp)[1]*dim(tp)[2]*length(prob)*runs)))
	dim(testp)<-c(dim(tp)[1],dim(tp)[2],length(prob),runs)

	# Create test patterns
	for (p in 1:length(prob)) {									  # Distortion probability counter
		pdist<-c(rep(1,(100*prob[p])),rep(0,(100-(100*prob[p])))) # Create probability distribution
		for (pat in 1:dim(tp)[2]) {								  # Training pattern counter
			for (patv in 1:dim(tp)[1]) {						  # Unit counter
				for (nrun in 1:runs) {							  # Run counter
					pn<-sample(pdist,1) # Randomly select from distribution
					
					# Update test activation state
					if (pn==1) {
						if (tp[patv,pat]==1) {
							testp[patv,pat,p,nrun]<-0
						}
						else {
							testp[patv,pat,p,nrun]<-1
						}
					}
					else {
						testp[patv,pat,p,nrun]<-tp[patv,pat]
					} # End update state
					
				} # End run counter
			} # End unit counter
			
		#Print test pattern to data file
		cat(sprintf("3Bb: Training pattern %.0f with distortion probability: %.1f\n",pat,prob[p]),file=datafilepath,append=TRUE)
		write.table(testp[,pat,p,],sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
		cat("\n",file=datafilepath,append=TRUE)
			
		} # End training pattern counter
	} # End distortion probability counter
		
	# Initiate required report variables
	ihdf<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
	dim(ihdf)<-c(length(prob),dim(testp)[2],runs)	# Zero matrix for hamming distance
	
	ni<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
	dim(ni)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of iterations to settle
	
	fi<-c(rep(0,length(prob)*dim(testp)[2]*runs)) 
	dim(fi)<-c(length(prob),dim(testp)[2],runs)		# Zero matrix for number of failures to settle
	
	
	# Perform tests
	for (p in 1:length(prob)) {
		for (pat in 1: dim(testp)[2]) {
			# Test each 'runs' times
			finala<-c()
			# Print header for energy state per iteration report
			cat(sprintf("\n3Bb.2: Energy for each iteration of each run for distortion probability %.1f of training pattern %.0f\n",prob[p],pat),file=efilepath,append=TRUE)
			
			for (r in 1:runs) {
				testnet_hopfield(W,matrix(testp[,pat,p,r],16,1),t)
				
				# Store final activations	
				finala<-cbind(finala,a[,dim(a)[2]])
				
				# Calculate hamming distance of final activation to original
				hdf<-0
				c<-dim(a)[2]
				for (n in 1:16) {
					if (a[n,c]!=tp[n,pat]) {
					hdf<-hdf+1
					}
				}
					
				# Store hamming distance of final activation to original
				ihdf[p,pat,r]<-hdf
				
				# Store number of iterations to settle
				ni[p,pat,r]<-dim(a)[2]
				
				# Store flag if failure to settle
				fi[p,pat,r]<-f
				
				# Print energy
				write.table(e,sep="\t",file=efilepath,append=TRUE,row.names=r,col.names=FALSE)
				#cat("\n",file=datafilepath,append=TRUE)
				
			}
			
			# Print final activation states
			cat(sprintf("3Bb.3: Final activation of units for distortion probability %.1f of training pattern %.0f, %.0f runs\n",prob[p],pat,runs),file=datafilepath,append=TRUE)
			write.table(finala,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
			cat("\n",file=datafilepath,append=TRUE)
		}
	}

	# Organize and calculate mean and standard deviations for hamming distance, iteration to settle, and failure to settle reports

	mhd<-c(rep(0,length(prob)*dim(testp)[2])) # Mean hamming distance matrix
	dim(mhd)<-c(length(prob),dim(testp)[2])

	sdhd<-c(rep(0,length(prob)*dim(testp)[2])) # Standard deviation hamming distance matrix
	dim(sdhd)<-c(length(prob),dim(testp)[2])

	mi<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean iterations to settle matrix
	dim(mi)<-c(length(prob),dim(testp)[2])

	sdi<-c(rep(0,length(prob)*dim(testp)[2]))  # Standard deviation iterations to settle matrix
	dim(sdi)<-c(length(prob),dim(testp)[2])

	nf<-c(rep(0,length(prob)*dim(testp)[2]))  # Mean failure to settle matrix
	dim(nf)<-c(length(prob),dim(testp)[2])

	for (p in 1:length(prob)) {
		for (pat in 1: dim(testp)[2]) {
			# Means
			mhd[p,pat]<-round(mean(ihdf[p,pat,]),1)
			mi[p,pat]<-round(mean(ni[p,pat,]),1)
			
			# Standard deviations
			sdhd[p,pat]<-round(sd(ihdf[p,pat,]),2)
			sdi[p,pat]<-round(sd(ni[p,pat,]),2)
			
			# Totals
			nf[p,pat]<-sum(fi[p,pat,])
		}
	}

	# Print reports
	cat(sprintf("\nSimulation report for %.0f distortion probabilities (rows) of %.0f training patterns (col) across %.0f runs:\n",length(prob),dim(testp)[2],runs),file=filepath,append=TRUE)

	cat("Mean hamming distance\n",file=filepath,append=TRUE)
	write.table(mhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
	cat("\n",file=filepath,append=TRUE)
	cat("Standard deviations of hamming distance\n",file=filepath,append=TRUE)
	write.table(sdhd,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
	cat("\n",file=filepath,append=TRUE)

	cat("Mean number of iterations to settle\n",file=filepath,append=TRUE)
	write.table(mi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
	cat("\n",file=filepath,append=TRUE)
	cat("Standard deviations of number of iterations to settle\n",file=filepath,append=TRUE)
	write.table(sdi,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
	cat("\n",file=filepath,append=TRUE)

	cat("Number of failures to settle\n",file=filepath,append=TRUE)
	write.table(nf,sep="\t",file=filepath,append=TRUE,row.names=prob,col.names=FALSE)
	cat("\n",file=filepath,append=TRUE)

}
#-----------------------------------------------------------------

# PART 4: BASE PATTERN
#-----------------------------------------------------------------
# 4A: Base training pattern
bp<-matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16,1)

# Create training patterns with .123 distortion probability
tp<-matrix(0,16,6)
pdist<-c(rep(1,125),rep(0,1000-125)) # 0.125 probability distribution of 1s

for (i in 1:6) {
	for (j in 1:16) {
		# Randomly sample and assign unit value from probability distribution
		pn<-sample(pdist,1)
		if (pn==1) {
				if (bp[j]==1) {tp[j,i]<-0}
				else {tp[j,i]<-1}
				}
			else {tp[j,i]<-bp[j]}
		}
	}

cat("PART 4A: Training patterns created from base pattern (0.125 switch probability)\n",file=datafilepath,append=TRUE)
write.table(tp,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# Train network with training patterns
W<-trainnet_hopfield(0,tp)
cat("PART 4A: Weight matrix from distorted base patterns\n",file=datafilepath,append=TRUE)
write.table(W,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# Test network with base pattern
cat("PART 4B: Test network with base pattern\n",file=filepath,append=TRUE)

testp<-bp

# Initiate required report variables
ihdf<-c(rep(0,dim(testp)[2]*runs)) # Zero matrix for hamming distance
dim(ihdf)<-c(dim(testp)[2],runs)

ni<-c(rep(0,dim(testp)[2]*runs)) # Zero matrix for number of iterations to settle
dim(ni)<-c(dim(testp)[2],runs)

fi<-c(rep(0,dim(testp)[2]*runs)) # Zero matrix for number of failures to settle
dim(fi)<-c(dim(testp)[2],runs)

# Test each 'runs' times
# Print header for energy state per iteration report
cat("\n4Bb.2: Energy for each iteration of each run\n",file=efilepath,append=TRUE)
finala<-c()
			
for (r in 1:runs) {
	testnet_hopfield(W,matrix(testp,16,1),t)
	
	# Store final activation state
	finala<-cbind(finala,a[,dim(a)[2]])
	
	# Calculate hamming distance of final activation to base pattern
	hdf<-0
	c<-dim(a)[2]
	for (n in 1:16) {
		if (a[n,c]!=bp[n]) {
			hdf<-hdf+1
		}
	}
		
	# Store hamming distance of final activation to base pattern
	ihdf[dim(testp)[2],r]<-hdf
	
	# Store number of iterations to settle
	ni[dim(testp)[2],r]<-dim(a)[2]
	
	# Store flag if failure to settle
	fi[dim(testp)[2],r]<-f
	
	# Print energy
	write.table(e,sep="\t",file=efilepath,append=TRUE,row.names=r,col.names=FALSE)
	#cat("\n",file=datafilepath,append=TRUE)
	
}

# Print activations	
cat(sprintf("4B.3: Final activation of units for run %.0f\n",r),file=datafilepath,append=TRUE)
write.table(finala,sep="\t",file=datafilepath,append=TRUE,col.names=FALSE)
cat("\n",file=datafilepath,append=TRUE)

# Calculate mean and standard deviations of hamming distance, iteration to settle, and failure to settle report
mhd<-round(mean(ihdf),1)
mi<-round(mean(ni),1)
sdhd<-round(sd(t(ihdf)),2)
sdi<-round(sd(t(ni)),2)
nf<-sum(fi)


# Print reports
cat(sprintf("\nSimulation report for base pattern across %.0f runs:\n",runs),file=filepath,append=TRUE)

cat("Mean hamming distance\n",file=filepath,append=TRUE)
write.table(mhd,sep="\t",file=filepath,append=TRUE,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of hamming distance\n",file=filepath,append=TRUE)
write.table(sdhd,sep="\t",file=filepath,append=TRUE,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Mean number of iterations to settle\n",file=filepath,append=TRUE)
write.table(mi,sep="\t",file=filepath,append=TRUE,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)
cat("Standard deviations of number of iterations to settle\n",file=filepath,append=TRUE)
write.table(sdi,sep="\t",file=filepath,append=TRUE,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

cat("Number of failures to settle\n",file=filepath,append=TRUE)
write.table(nf,sep="\t",file=filepath,append=TRUE,col.names=FALSE)
cat("\n",file=filepath,append=TRUE)

#-----------------------------------------------------------------


# Close files
close(outputfile)
close(datafile)
close(efile)

# Clean environment workspace
rm(list=ls())


