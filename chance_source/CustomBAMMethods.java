package songlab;

import net.sf.samtools.*;

import java.io.File;
import java.util.*;
import javax.swing.JProgressBar;

public class CustomBAMMethods {
    /**
     * Read a SAM or BAM file; output phred_hist, d (# chromosome reads), and nucleotide frequencies
     */
    public static Object[] makeDensityFromFile(final File inputFile, int bins, final String[] chromosomeNames, final int[] chromosomeLengths, JProgressBar progressBar) {
		
		long fileSize = inputFile.length();
		long progress = 0, percentProgress;
		final SAMFileReader inputSam = new SAMFileReader(inputFile);
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
		final SAMRecordIterator recordIterator = inputSam.iterator();
		SAMRecord samRecord = recordIterator.next();
		final int sequenceLength = samRecord.getReadLength();
		int[][] phredHistogram = new int[127][sequenceLength]; // a new array's elements are initialized to 0 in Java
		List<SAMSequenceRecord> sequenceRecordList = samRecord.getHeader().getSequenceDictionary().getSequences();
		List<String> sequenceNames = new ArrayList<String>();
		for (SAMSequenceRecord sequenceRecord : sequenceRecordList) {
			sequenceNames.add(sequenceRecord.getSequenceIndex(), sequenceRecord.getSequenceName());
		}
		String[] sequenceNameArray = new String[sequenceNames.size()];
    	sequenceNameArray = sequenceNames.toArray(sequenceNameArray);
    	List<String> chromosomeNamesList = Arrays.<String>asList(chromosomeNames);
		int[][] alignmentDensities = new int[sequenceNameArray.length][]; //first index: chromosome; second index: bin
		for (int i=0; i<sequenceNameArray.length; i++) {
			alignmentDensities[i] = new int[(chromosomeLengths[chromosomeNamesList.indexOf(sequenceNameArray[i])] - 1) / bins + 1]; //number of bins per chromosome is chromosome length / bins + 1
		}
		int[][] baseCounts = new int[5][sequenceLength]; // Rows: A, T, C, G, and N
		int readCount = 0;
		byte[] baseReads, baseQualities;
		
        while (recordIterator.hasNext()) {
        	samRecord = recordIterator.next();
        	if (!samRecord.getReadUnmappedFlag()) { //ensure read is mapped
        		baseReads = samRecord.getReadBases();
        		baseQualities = samRecord.getBaseQualities();
        		if (!samRecord.getReadNegativeStrandFlag()) {
        			for (int i=0; i<baseReads.length; i++) {
        				if (baseReads[i] == 65) {
        					baseCounts[0][i]++; // A
        				}
        				else if (baseReads[i] == 84) {
        					baseCounts[1][i]++; // T
        				}
        				else if (baseReads[i] == 67) {
        					baseCounts[2][i]++; // C
        				}
        				else if (baseReads[i] == 71) {
        					baseCounts[3][i]++; // G
        				}
        				else if (baseReads[i] == 78) {
        					baseCounts[4][i]++; // N
        				}
        			}
        			for (int i=0; i<baseQualities.length; i++) {
        				phredHistogram[baseQualities[i]+33][i]++;
        			}
        		}
        		else { //reverse complement
        			for (int i=0; i<baseReads.length; i++) {
        				if (baseReads[i] == 65) { 
        					baseCounts[1][baseReads.length-i-1]++; // T
        				}
        				else if (baseReads[i] == 84) {
        					baseCounts[0][baseReads.length-i-1]++; // A
        				}
        				else if (baseReads[i] == 67) {
        					baseCounts[3][baseReads.length-i-1]++; // G
        				}
        				else if (baseReads[i] == 71) {
        					baseCounts[2][baseReads.length-i-1]++; // C
        				}
        				else if (baseReads[i] == 78) {
        					baseCounts[4][baseReads.length-i-1]++; // N
        				}
        			}
        			for (int i=0; i<baseQualities.length; i++) {
        				phredHistogram[baseQualities[i]+33][baseReads.length-i-1]++;
        			}
        		}
        		alignmentDensities[sequenceNames.indexOf(samRecord.getReferenceName())][(samRecord.getAlignmentStart()-1) / bins]++;
        		readCount++;
        	}
        	progress += samRecord.getAttributesBinarySize()*3;
        	percentProgress = progress*100/fileSize;
        	if (percentProgress > 100) {
        		percentProgress = 100;
        	}
        	progressBar.setValue((int) percentProgress);
        }

		return new Object[]{readCount, sequenceNameArray, alignmentDensities, baseCounts, phredHistogram};
    }
}