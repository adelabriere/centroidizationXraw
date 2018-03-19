if(!require(Rcpp)) stop("Rcpp required.")
if(!require(inline)) stop("inline required.")

src_include_centroid <- '#include <vector>
#include<math.h>
#include <iostream>
#include<Rcpp.h>
#define SIZE_BUFF 1000


int insertCentroid(std::vector<double>& BufferMz, std::vector<double>& BufferIntensity,std::vector<bool> DiffVec,
int& i, std::vector<double>& mz,std::vector<double>& intensity,
int b1, int b2, int b3, double itr)
{
	if(!((b1<=b2)&(b2<=b3)&(b3>b1)))
	{
	return 0;
	}

	//Checking that diff of masses are ok.
	for(int j = b1;j < b3;j++){
	if(DiffVec[j]) return 0;
	}
	//std::cout << b1 << " " << b2 << " " << b3 << std::endl;
	int range_c = (b2-b1)>(b3-b2) ? b3-b2 : b2-b1;
	//We chose to consider intensity as the max of the detected peaks.
	double max_int = 0;
	double accu_int = 0;
	//std::cout <<"accuint_"<<std::endl;
	for(int j = b1; j<=b3; j++)
	{
	if(intensity[j]>max_int) max_int = intensity[j];
	accu_int += intensity[j];
	}
	if(accu_int<itr)
	{
	return 1;
	}
	double accu_mz = 0;
	
	double accu_int_p = 0;
	for(int j = b2-range_c; j<=b2+range_c; j++)
	{
	accu_mz += mz[j]*intensity[j];
	accu_int_p += intensity[j];
	}
	accu_mz /= accu_int_p;
	BufferMz[i] = accu_mz;
	BufferIntensity[i] = accu_int;
	i++;
	return 1;
}

void checkBufferSize(std::vector<double>& BufferMz,std::vector<double>& BufferInt,int& i,std::vector<double>& cmz,std::vector<double>& cint)
{
	if(i<SIZE_BUFF)
	{
	return;
	}
	cmz.insert( cmz.end() , BufferMz.begin() , BufferMz.end());
	cint.insert( cint.end() , BufferInt.begin() , BufferInt.end());
	i=0;
}

double thresholdf(double x, double a, double b)
{
	return a*x+b;
}


SEXP centroidScan(std::vector<double>& mz,std::vector<double>& intensity,double intthreshold, std::vector<double>& cmz, std::vector<double>& cint, double a,double b)
{
	
	
	std::vector<double> BufferMz(SIZE_BUFF,0);
	std::vector<double> BufferInt(SIZE_BUFF,0);
	std::vector<bool> DiffVec(mz.size()-1,bool(0));

	//First pass on the difference.
	int i;
	for(i=1;i<(int)mz.size();i++){
		DiffVec[i-1]=(mz[i]-mz[i-1])>thresholdf(sqrt((mz[i]+mz[i-1])/2),a,b);
	}
	int pBuffer = 0;
	
	int peakBeginning = 0;
	int peakEnd = -1;
	int peakMax = 0;
	//std :: cout << std::endl;
	
	//A point is a centroid, if it is a local maxima, his diff with the previous point
	//are inferior 
	for(i=1; i<(int)mz.size()-1; i++){

		//We restart the peak everytime a peak cross the mz limit.
		if(DiffVec[i-1])
		{
			peakEnd = i-1;
			insertCentroid(BufferMz, BufferInt, DiffVec,pBuffer,mz,intensity, peakBeginning, peakMax,peakEnd, intthreshold);
			checkBufferSize(BufferMz,BufferInt,pBuffer,cmz,cint);
			peakBeginning = i;
			continue;
		}
		else if((intensity[i]<=intensity[i-1])&(intensity[i]<=intensity[i+1])) //Case where it a local minima
		{
			peakEnd = i;
			int boolp = insertCentroid(BufferMz, BufferInt, DiffVec, pBuffer,mz,intensity, peakBeginning, peakMax, peakEnd, intthreshold);
			checkBufferSize(BufferMz,BufferInt,pBuffer,cmz,cint);
			if(boolp)
			{
				peakBeginning = i;
			}
			continue;
		}
		
		if((intensity[i]>intensity[i-1])&(intensity[i]>intensity[i+1])) //Maximum of the peak.
		{
			peakMax = i;
		}
	}

	//We add the content of the last buffer to the scan.
	if(DiffVec[mz.size()-2]){
		peakEnd = mz.size()-2;
	}else{
		peakEnd = mz.size()-1;
	}
	insertCentroid(BufferMz, BufferInt, DiffVec, pBuffer,mz,intensity, peakBeginning, peakMax,peakEnd, intthreshold);
	cmz.insert( cmz.end() , BufferMz.begin() , BufferMz.begin()+pBuffer);
	cint.insert( cint.end() , BufferInt.begin() , BufferInt.begin()+pBuffer);

	//Return a data.frame.
	return Rcpp::DataFrame::create(Rcpp::Named("mz")=Rcpp::wrap(cmz),
	Rcpp::Named("intensity")=Rcpp::wrap(cint));
}'

	src_centroid <- 'std::vector<double> mz = Rcpp::as<std::vector<double> >(rmz);
	std::vector<double> intensity = Rcpp::as<std::vector<double> >(rintensity);
	std::vector<double> cmz;
	std::vector<double> cint;
	double intthreshold = Rcpp::as<double>(rthresholdint);
	double a = Rcpp::as<double>(ra);
	double b = Rcpp::as<double>(rb);
	return centroidScan(mz,intensity,intthreshold,cmz,cint, a, b);'
	
	cpp_centroid <- cxxfunction(signature(rmz="numeric",rintensity="numeric",rthresholdint="numeric",ra="numeric",rb="numeric"),plugin="Rcpp",
								includes = src_include_centroid,
								body=src_centroid)
	
	
	centroidXraw <- function(xraw, intthreshold = 200,graphical=TRUE,frac=2.2,levels=c("MS","MSn")) {
		cxraw <- NULL
		if("MS" %in% levels){
			##Mz threshold function is approximated on most intense scan.
			pmax <- which.max(xraw@tic)
			x <- xraw@env$mz[(xraw@scanindex[pmax] + 1):xraw@scanindex[pmax + 1]]
			dfmz <- diff(x)
			if (graphical)
				plot(
					sqrt((x[1:(length(x) - 1)] + x[2:length(x)]) / 2),
					dfmz,
					ylim = c(0, max(dfmz)),
					col = topo.colors(20)[1:10][ceiling(log10(xraw@env$intensity[(xraw@scanindex[pmax] + 1):xraw@scanindex[pmax + 1]]))],
					xlab = "sqrt(m/z)",
					ylab = "diff(m/z)"
				)
			
			vdiffg <- runmed(dfmz, 11)
			
			xabs <- sqrt((x[1:(length(x) - 1)] + x[2:length(x)]) / 2)
			
			#Liner regression to obtain the value of the parameter.
			reg <- lm(vdiffg ~ xabs, weights = 1 / vdiffg)
			sreg <- summary(reg)$coefficients[, 4]
			if (any(sreg > 0.1)) {
				print(summary(reg))
				warning(
					"Regression failed, centroidization may fail. Check if data are not already in centroid mode."
				)
			}
			
			cc <- coef(reg)
			if (graphical) {
				plot(
					xabs,
					vdiffg,
					ylim = c(0, max(dfmz)),
					col = topo.colors(20)[1:10][ceiling(log10(xraw@env$intensity[(xraw@scanindex[pmax] + 1):xraw@scanindex[pmax + 1]]))],
					xlab = "sqrt(m/z)",
					ylab = "diff(m/z)"
				)
				lines(xabs, xabs * cc[2] + cc[1], col = "red", lwd = 2)
			}
			a <- frac*cc[2]
			b <- cc[1]
			sci <- c(xraw@scanindex, length(xraw@env$mz))
			lcentroid <- list()
			for (i in seq_along(xraw@scanindex)) {
				lcentroid[[i]] <-
					cpp_centroid(xraw@env$mz[(sci[i] + 1):sci[i + 1]], xraw@env$intensity[(sci[i] +
																						   	1):sci[i + 1]],
								 intthreshold, a, b)
			}
			
			cxraw <- xraw
			cxraw@env <- as.environment(as.list(xraw@env, all.names = TRUE))
			putative_scanindex <- as.integer(cumsum(c(0, sapply(lcentroid, nrow))))
			cxraw@scanindex <-
				as.integer(cumsum(c(0, sapply(lcentroid, nrow)))[-(length(lcentroid)+1)])
			cxraw@env$mz <- unlist(sapply(lcentroid, function(x) {
				x[, 1]
			}))
			cxraw@env$intensity <- unlist(sapply(lcentroid, function(x) {
				x[, 2]
			}))
		}
		
		if("MSn" %in% levels & length(xraw@msnPrecursorIntensity)!=0){
			pmax <- which.max(xraw@msnPrecursorIntensity)
			sidx <- c(xraw@msnScanindex,length(xraw@env$msnMz))
			x <- xraw@env$msnMz[(sidx[pmax] + 1):sidx[pmax + 1]]
			dfmz <- diff(x)
			if (graphical)
				plot(
					sqrt((x[1:(length(x) - 1)] + x[2:length(x)]) / 2),
					dfmz,
					ylim = c(0, max(dfmz)),
					col = topo.colors(20)[1:10][ceiling(log10(xraw@env$msnIntensity[(sidx[pmax] + 1):sidx[pmax + 1]]))],
					xlab = "sqrt(m/z)",
					ylab = "diff(m/z)"
				)
			
			vdiffg <- runmed(dfmz, 11)
			
			xabs <- sqrt((x[1:(length(x) - 1)] + x[2:length(x)]) / 2)
			
			#Liner regression to obtain the value of the parameter.
			reg <- lm(vdiffg ~ xabs, weights = 1 / vdiffg)
			sreg <- summary(reg)$coefficients[, 4]
			if (any(sreg > 0.1)) {
				print(summary(reg))
				warning(
					"Regression failed, centroidization may fail. Check if data are not laready in centroid mode."
				)
			}
			
			cc <- coef(reg)
			if (graphical) {
				plot(
					xabs,
					vdiffg,
					ylim = c(0, max(dfmz)),
					col = topo.colors(20)[1:10][ceiling(log10(xraw@env$msnIntensity[(sidx[pmax] + 1):sidx[pmax + 1]]))],
					xlab = "sqrt(m/z)",
					ylab = "diff(m/z)"
				)
				lines(xabs, xabs * cc[2] + cc[1], col = "red", lwd = 2)
			}
			a <- frac*cc[2]
			b <- cc[1]
			sci <- sidx
			lcentroid <- list()
			for (i in seq_along(xraw@msnScanindex)) {
				
				lcentroid[[i]] <-
					cpp_centroid(xraw@env$msnMz[(sci[i] + 1):sci[i + 1]], xraw@env$msnIntensity[(sci[i] +
																								 	1):sci[i + 1]],
								 intthreshold, a, b)
			}
			if(!("MS" %in% levels)){
				cxraw <- xraw
				cxraw@env <- as.environment(as.list(xraw@env, all.names = TRUE))
			}
			putative_scanindex <- as.integer(cumsum(c(0, sapply(lcentroid, nrow))))
			cxraw@msnScanindex <-
				as.integer(cumsum(c(0, sapply(lcentroid, nrow)))[-(length(lcentroid)+1)])
			cxraw@env$msnMz <- unlist(sapply(lcentroid, function(x) {
				x[, 1]
			}))
			cxraw@env$msnIntensity <- unlist(sapply(lcentroid, function(x) {
				x[, 2]
			}))
		}
		
		message(paste("\nCentroidization for ",basename(xraw@filepath)," finished data passed from ",length(xraw@env$mz)," to ",length(cxraw@env$mz)," MS points.",
					  "\nand from ",length(xraw@env$msnMz)," to ",length(cxraw@env$msnMz)," MSn data points",sep=""))
		cxraw
	}
	