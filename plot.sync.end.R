#
# Author: ahfyth
# Date  : Dec 12, 2018
# R script for plotting the coverage and end signal around tissue-specific open chromatin regions
#

argv = commandArgs(T);
if( length(argv) < 1 )
{
	print( 'usage: R --slave --args <sid> < plot.R' );
	q();
}

enhancerVersion='AllOC';
sid=argv[1];

outfileName=paste( sid, "vs", enhancerVersion, "fragmentation.analysis.pdf", sep=".");
pdf( outfileName, width=12 );

## figure configurations
yrange=c(4, 6);
cs=c( 'black', 'red', 'blue' );
loessspan=0.05;

tissueList=c('Tcell', 'Liver', 'Placenta', 'Lung', 'Ovary', 'Breast', 'Intestine');

for( t in tissueList )
{
	cov=read.table( paste(sid, "vs", enhancerVersion, t, "sync.value", sep=".") );
	end=read.table( paste(sid, "vs", enhancerVersion, t, "sync.end",   sep=".") );

	## RAW SIGNALS FOLLOWED BY LOWESS FITTED
	plot( cov[,3]~cov[,1], type='l', lwd=4, col=cs[1], ylim=yrange,
			xlab='Relative position to center of open chromation regions (bp)',
			ylab='Normalized signal/end density',
			main=paste(t, "-spc enhancers, CaseID=", sid, sep="" ));
	lines( end[,3]~end[,1], col=cs[2] );
	lines( end[,5]~end[,1], col=cs[3] );
	legend( 'bottomleft', c('Coverage', 'U end', 'E end'), col=cs, lty=rep(1,3), bty='n', cex=1.1 );

	loess3=loess(end[,3]~end[,1], span=loessspan);
	pend3=predict(loess3, end[,1]);
	loess5=loess(end[,5]~end[,1], span=loessspan);
	pend5=predict(loess5, end[,1]);
	ymin=min(cov[,3], pend3, pend5);
	ymax=max(cov[,3], pend3, pend5);
	if( ymin>4 && ymax<6 )
	{
		ymin = 4;
		ymax = 6;
	} else {
		dist=max( 5-ymin, ymax-5 );
		ymin = 5-dist;
		ymax = 5+dist;
	}
	plot( cov[,3]~cov[,1], type='l', lwd=4, col=cs[1], ylim=c(ymin, ymax),
			xlab='Relative position to center of open chromation regions (bp)',
			ylab='Normalized signal/end density (LOWESS smoothed)',
			main=paste(t, "-spc enhancers, CaseID=", sid, sep="" ));
	lines( end[,1], pend3, col=cs[2], lwd=2 );
	lines( end[,1], pend5, col=cs[3], lwd=2 );
	legend( 'bottomleft', c('Coverage', 'U end', 'E end'), col=cs, lty=rep(1,3), bty='n', cex=1.1 );
}

