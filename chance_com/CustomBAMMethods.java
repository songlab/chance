import net.sf.samtools.*;
import java.io.File;
import java.util.*;

public class CustomBAMMethods {
    /**
     * Read a SAM or BAM file; output phred_hist, d (# chromosome reads), and nucleotide frequencies
     */
    public static Object[] makeDensityFromFile(final File inputFile, int bins, final String[] chromosomeNames, final int[] chromosomeLengths) {
		
		long fileSize = inputFile.length();
		final SAMFileReader inputSam = new SAMFileReader(inputFile);
		final SAMRecordIterator recordIterator = inputSam.iterator();
		SAMRecord samRecord = recordIterator.next();
		final int sequenceLength = samRecord.getReadLength();
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
		int readCount = 0;
		while (recordIterator.hasNext()) {
		    samRecord = recordIterator.next();
		    if (!samRecord.getReadUnmappedFlag()) { //ensure read is mapped
        		alignmentDensities[sequenceNames.indexOf(samRecord.getReferenceName())][(samRecord.getAlignmentStart()-1) / bins]++;
        		readCount++;
		    }
		}
		return new Object[]{readCount, sequenceNameArray, alignmentDensities};
    }
}