package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.ExpandingArrayList;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.VariantDataManager;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.VariantDatum;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.VariantRecalibratorEngine;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.mode.GaussianMixtureModel;

public class RScriptUtils {
	public static void createVisualizationScript(
			final String fileName,
			final VariantRecalibratorEngine engine,
			final VariantDataManager dataManager,
	        final List<VariantDatum> randomData,
	        final GaussianMixtureModel goodModel,
	        final GaussianMixtureModel badModel,
	        final double lodCutoff,
	        final String[] annotationKeys ) {
	    final PrintStream stream;
	    try {
	        stream = new PrintStream(fileName);
	    } catch (final FileNotFoundException e ) {
	        throw new UserException.CouldNotCreateOutputFile(fileName, e);
	    }

	    // We make extensive use of the ggplot2 R library: http://had.co.nz/ggplot2/
	    stream.println("library(ggplot2)");
	    // For compactPDF in R 2.13+
	    stream.println("library(tools)");
	    // For graphical functions R 2.14.2+
	    stream.println("library(grid)");

	    createArrangeFunction( stream );

	    stream.println("outputPDF <- \"" + fileName + ".pdf\"");
	    stream.println("pdf(outputPDF)"); // Unfortunately this is a huge pdf file, BUGBUG: need to work on reducing the file size

	    for(int iii = 0; iii < annotationKeys.length; iii++) {
	        for( int jjj = iii + 1; jjj < annotationKeys.length; jjj++) {

	            final List<VariantDatum> fakeData = new ExpandingArrayList<>();
	            double minAnn1 = 100.0, maxAnn1 = -100.0, minAnn2 = 100.0, maxAnn2 = -100.0;
	            for( final VariantDatum datum : randomData ) {
	                minAnn1 = Math.min(minAnn1, datum.annotations[iii]);
	                maxAnn1 = Math.max(maxAnn1, datum.annotations[iii]);
	                minAnn2 = Math.min(minAnn2, datum.annotations[jjj]);
	                maxAnn2 = Math.max(maxAnn2, datum.annotations[jjj]);
	            }
	            // Create a fake set of data which spans the full extent of these two annotation dimensions in order
	            // to calculate the model PDF projected to 2D
	            final double NUM_STEPS = 60.0;
	            for(double ann1 = minAnn1; ann1 <= maxAnn1; ann1+= (maxAnn1 - minAnn1) / NUM_STEPS) {
	                for(double ann2 = minAnn2; ann2 <= maxAnn2; ann2+= (maxAnn2 - minAnn2) / NUM_STEPS) {
	                    final VariantDatum datum = new VariantDatum();
	                    datum.prior = 0.0;
	                    datum.annotations = new double[randomData.get(0).annotations.length];
	                    datum.isNull = new boolean[randomData.get(0).annotations.length];
	                    for(int ann=0; ann< datum.annotations.length; ann++) {
	                        datum.annotations[ann] = 0.0;
	                        datum.isNull[ann] = true;
	                    }
	                    datum.annotations[iii] = ann1;
	                    datum.annotations[jjj] = ann2;
	                    datum.isNull[iii] = false;
	                    datum.isNull[jjj] = false;
	                    fakeData.add(datum);
	                }
	            }

	            engine.evaluateData( fakeData, goodModel, false );
	            engine.evaluateData( fakeData, badModel, true );

	            stream.print("surface <- c(");
	            for( final VariantDatum datum : fakeData ) {
	                stream.print(String.format("%.4f, %.4f, %.4f, ",
	                        dataManager.denormalizeDatum(datum.annotations[iii], iii),
	                        dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
	                        Math.min(4.0, Math.max(-4.0, datum.lod))));
	            }
	            stream.println("NA,NA,NA)");
	            stream.println("s <- matrix(surface,ncol=3,byrow=T)");

	            stream.print("data <- c(");
	            for( final VariantDatum datum : randomData ) {
	                stream.print(String.format("%.4f, %.4f, %.4f, %d, %d,",
	                        dataManager.denormalizeDatum(datum.annotations[iii], iii),
	                        dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
	                        (datum.lod < lodCutoff ? -1.0 : 1.0),
	                        (datum.atAntiTrainingSite ? -1 : (datum.atTrainingSite ? 1 : 0)), (datum.isKnown ? 1 : -1)));
	            }
	            stream.println("NA,NA,NA,NA,1)");
	            stream.println("d <- matrix(data,ncol=5,byrow=T)");

	            final String surfaceFrame = "sf." + annotationKeys[iii] + "." + annotationKeys[jjj];
	            final String dataFrame = "df." + annotationKeys[iii] + "." + annotationKeys[jjj];

	            stream.println(surfaceFrame + " <- data.frame(x=s[,1], y=s[,2], lod=s[,3])");
	            stream.println(dataFrame + " <- data.frame(x=d[,1], y=d[,2], retained=d[,3], training=d[,4], novelty=d[,5])");
	            stream.println("dummyData <- " + dataFrame + "[1,]");
	            stream.println("dummyData$x <- NaN");
	            stream.println("dummyData$y <- NaN");
	            stream.println("p <- ggplot(data=" + surfaceFrame + ", aes(x=x, y=y)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
	            stream.println("p1 = p +ggtitle(\"model PDF\") + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + geom_tile(aes(fill = lod)) + scale_fill_gradient(high=\"green\", low=\"red\", space=\"rgb\")");
	            stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=retained, alpha=I(1/7),legend=FALSE) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
	            stream.println("q <- geom_point(aes(x=x,y=y,color=retained),data=dummyData, alpha=1.0, na.rm=TRUE)");
	            stream.println("p2 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"outcome\", high=\"black\", low=\"red\",breaks=c(-1,1),guide=\"legend\",labels=c(\"filtered\",\"retained\"))");
	            stream.println("p <- qplot(x,y,data="+ dataFrame + "["+dataFrame+"$training != 0,], color=training, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
	            stream.println("q <- geom_point(aes(x=x,y=y,color=training),data=dummyData, alpha=1.0, na.rm=TRUE)");
	            stream.println("p3 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(high=\"green\", low=\"purple\",breaks=c(-1,1),guide=\"legend\", labels=c(\"neg\", \"pos\"))");
	            stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=novelty, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
	            stream.println("q <- geom_point(aes(x=x,y=y,color=novelty),data=dummyData, alpha=1.0, na.rm=TRUE)");
	            stream.println("p4 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"novelty\", high=\"blue\", low=\"red\",breaks=c(-1,1),guide=\"legend\", labels=c(\"novel\",\"known\"))");
	            stream.println("arrange(p1, p2, p3, p4, ncol=2)");
	        }
	    }
	    stream.println("dev.off()");

	    stream.println("if (exists(\"compactPDF\")) {");
	    stream.println("compactPDF(outputPDF)");
	    stream.println("}");

	    stream.close();
	}
	
	private static void createArrangeFunction( final PrintStream stream ) {
	    stream.println("vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)");
	    stream.println("arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {");
	    stream.println("dots <- list(...)");
	    stream.println("n <- length(dots)");
	    stream.println("if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}");
	    stream.println("if(is.null(nrow)) { nrow = ceiling(n/ncol)}");
	    stream.println("if(is.null(ncol)) { ncol = ceiling(n/nrow)}");
	    stream.println("grid.newpage()");
	    stream.println("pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )");
	    stream.println("ii.p <- 1");
	    stream.println("for(ii.row in seq(1, nrow)){");
	    stream.println("ii.table.row <- ii.row ");
	    stream.println("if(as.table) {ii.table.row <- nrow - ii.table.row + 1}");
	    stream.println("for(ii.col in seq(1, ncol)){");
	    stream.println("ii.table <- ii.p");
	    stream.println("if(ii.p > n) break");
	    stream.println("print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))");
	    stream.println("ii.p <- ii.p + 1");
	    stream.println("}");
	    stream.println("}");
	    stream.println("}");
	}
}
