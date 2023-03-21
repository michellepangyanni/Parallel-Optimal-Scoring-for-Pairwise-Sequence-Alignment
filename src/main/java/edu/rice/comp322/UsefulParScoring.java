package edu.rice.comp322;

import edu.rice.hj.api.SuspendableException;
import edu.rice.hj.api.*;
import static edu.rice.hj.Module0.*;
import static edu.rice.hj.Module1.*;

/**
 * A scorer that works in parallel.
 */
public class UsefulParScoring extends AbstractDnaScoring {

    /**
     * <p>main.</p> Takes the names of two files, and in parallel calculates the sequence aligment scores of the two DNA
     * strands that they represent.
     *
     * @param args The names of two files.
     */

    /*
     * Fields.
     */
    /**
     * The length of the first sequence.
     */
    private final int xLength;
    /**
     * The length of the second sequence.
     */
    private final int yLength;
    /**
     * The Smith-Waterman matrix.
     */
    private final int[][] s;

    public static void main(final String[] args) throws Exception {
        final ScoringRunner scoringRunner = new ScoringRunner(UsefulParScoring::new);
        scoringRunner.start("UsefulParScoring", args);
    }

    /**
     * Creates a new UsefulParScoring.
     *
     * @param xLength length of the first sequence
     * @param yLength length of the second sequence
     */
    public UsefulParScoring(final int xLength, final int yLength) {
        if (xLength <= 0 || yLength <= 0) {
            throw new IllegalArgumentException("Lengths (" + xLength + ", " + yLength + ") must be positive!");
        }
        this.xLength = xLength;
        this.yLength = yLength;
        // Preallocate the matrix for alignment, dimension+1 to initialize the matrix
        this.s = new int[xLength + 1][yLength + 1];

        // Initialize the row of the matrix
        for (int ii = 1; ii < xLength + 1; ii++) {
            this.s[ii][0] = getScore(1, 0) * ii;
        }

        // Initialize the first column
        for (int jj = 1; jj < yLength + 1; jj++) {
            this.s[0][jj] = getScore(0, 1) * jj;
        }
        // Initialize top left corner
        this.s[0][0] = 0;
//        this.rows = Math.min(2 * numWorkerThreads() + 1, xLength);
//        this.cols = Math.min(2 * numWorkerThreads() + 1, yLength);
//        final HjDataDrivenFuture<Integer>[][] ddfMatrix;
//        for (int r = 1; r <= rows; )
    }

    /**
     * This is a sequential implementation of the Smith-Waterman alignment algorithm used for sub-matrices in the parallel
     * implementation of scoreSequences.
     * @param x The first sequence to be compared.
     * @param y The second sequence to be compared.
     * @param rowStart The first index in a row.
     * @param colStart The first index in a column.
     * @param rowEnd The last index in a row.
     * @param colEnd The last index in a column.
     * @return final score of the baby matrix
     * @throws SuspendableException
     */
    public int scoreBabySequences(final String x, final String y, int rowStart, int colStart, int rowEnd, int colEnd){
        for (int row = rowStart; row <= rowEnd; row++) {
            for (int col = colStart; col <= colEnd; col++) {
                // the two characters to be compared
                final char XChar = x.charAt(row - 1);
                final char YChar = y.charAt(col - 1);

                // find the largest point to jump from, and use it
                final int diagScore = this.s[row - 1][col - 1] + getScore(charMap(XChar), charMap(YChar));
                final int topColScore = this.s[row - 1][col] + getScore(charMap(XChar), 0);
                final int leftRowScore = this.s[row][col - 1] + getScore(0, charMap(YChar));
                int score = Math.max(diagScore, Math.max(leftRowScore, topColScore));
                doWork(1);
                this.s[row][col] = score;
            }
        }
        // return score of baby matrix
        return this.s[rowEnd - 1][colEnd - 1];
    }

    /**
     * This is a parallel implementation of the Smith-Waterman algorithm that
     * achieve the smallest execution time using 16 cores on a dedicated NOTS compute node with DDF and chunking.
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {
        if (xLength <= 0 || yLength <= 0) {
            throw new IllegalArgumentException("Lengths (" + xLength + ", " + yLength + ") must be positive!");
        }
        // 1. Initialization.
        final HjDataDrivenFuture<Integer>[][] ddfMatrix;
        final int xPartition;
        final int yPartition;

        // 2. Partition.
        xPartition = (int) Math.ceil((double) this.xLength / 80);
        yPartition = (int) Math.ceil((double) this.yLength / 80);
        int numXPartition = this.xLength / xPartition;
        int numYPartition = this.yLength / yPartition;

        // 3. Create DDF matrix.
        ddfMatrix = new HjDataDrivenFuture[numXPartition + 1][numYPartition + 1];
        // Initialize DDF matrix.
        for (int i = 0; i <= numXPartition; i++) {
            for (int j = 0; j <= numYPartition; j++) {
                ddfMatrix[i][j] = newDataDrivenFuture();
            }
        }
        // Initialize top corner.
        ddfMatrix[0][0].put(0);

        // Initialize first row.
        for (int row = 1; row <= numXPartition; row ++) {
            ddfMatrix[row][0].put(0);
        }
        // Initialize first column.
        for (int col = 1; col <= numYPartition; col ++) {
            ddfMatrix[0][col].put(0);
        }

        // 3. Scoring.
        finish(()->{
            // Iterate over DDF matrix.
//            forall(1, numXPartition, (i)->{
//
//            });
            for (int i = 1; i <= numXPartition; i++) {
//            forall(1, numXPartition, 1, numYPartition, (i, j)->{

                for (int j = 1; j <= numYPartition; j++) {
                    final int ii = i;
                    final int jj = j;
                    // Get this DDF's dependencies.
                    // For each column.
                    HjFuture<Integer> top = ddfMatrix[ii - 1][jj];
                    // For each diagonal.
                    HjFuture<Integer> diag = ddfMatrix[ii - 1][jj - 1];
                    // For each row.
                    HjFuture<Integer> left = ddfMatrix[ii][jj - 1];

                    // Wait for the dependencies
                    asyncAwait(top, left, diag, () -> {
                        // Convert DDF indices to scoring matrix indices.
                        int rowStart = ((ii - 1) * xPartition) + 1;
                        int colStart = ((jj - 1) * yPartition) + 1;
                        int rowEnd = Math.min(ii * xPartition, xLength);
                        int colEnd = Math.min(jj * yPartition, yLength);

                        // Score this sub-matrix.
                        int chunkScore = scoreBabySequences(x, y, rowStart, colStart, rowEnd, colEnd);
                        // Testing
                        ddfMatrix[ii][jj].put(chunkScore);
                    });

//            });
                }
            }
        });
        // Return the score in the bottom right corner of the scoring matrix.
        return this.s[this.xLength][this.yLength];
    }
}

