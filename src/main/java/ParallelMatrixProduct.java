public class ParallelMatrixProduct {
    private final int countThread;
    private final Matrix a;
    private final Matrix b;
    private final Matrix result;
    private final boolean bIsTranspose;

    ParallelMatrixProduct(Matrix a, Matrix b, int countThread,
                          boolean bIsTranspose) throws Exception {
        this.a = a;
        this.b = b;
        this.bIsTranspose = bIsTranspose;
        this.countThread = countThread;
        this.result = new Matrix(a.getRows(), b.getColumns());
    }

    public Matrix multiplyParallel() {
        assert countThread > 0;

        final int rowCount = a.getRows();
        final int colCount = b.getColumns();
        final int[][] result = new int[rowCount][colCount];

        final int cellsForThread = (rowCount * colCount) / countThread;
        int firstIdxCell = 0;
        final MyThread[] multiplierThreads = new MyThread[countThread];

        for (int threadIdx = countThread - 1; threadIdx >= 0; --threadIdx) {
            // Индекс последней вычисляемой ячейки.
            int lastIdxCell = firstIdxCell + cellsForThread;
            if (threadIdx == 0) {
                /* Один из потоков должен будет вычислить не только свой блок ячеек,
                   но и остаток, если число ячеек не делится нацело на число потоков. */
                lastIdxCell = rowCount * colCount;
            }
            multiplierThreads[threadIdx] = new MyThread(firstIdxCell, lastIdxCell);
            multiplierThreads[threadIdx].start();
            firstIdxCell = lastIdxCell;
        }

        // Ожидание завершения потоков.
        try {
            for (final MyThread multiplierThread : multiplierThreads)
                multiplierThread.join();
        }
        catch (InterruptedException e) {
            e.printStackTrace();
        }
        return new Matrix(result);
    }

    private class MyThread extends Thread {
        private final int startIdx;
        private final int endIdx;
        private final int sumLength;

        public MyThread(final int startIdx, final int endIdx) {
            this.startIdx = startIdx;
            this.endIdx = endIdx;

            sumLength = b.getRows();
        }

        private void calcValue(final int row, final int col) {
            int sum = 0;
            for (int i = 0; i < sumLength; ++i) {
                    sum += (bIsTranspose)
                            ? a.getMatrix()[row][i] * b.getMatrix()[col][i]
                            : a.getMatrix()[row][i] * b.getMatrix()[i][col];
            }
            result.getMatrix()[row][col] = sum;
        }

        @Override
        public void run() {
            final int colCount = b.getColumns();
            for (int index = startIdx; index < endIdx; ++index)
                calcValue(index / colCount, index % colCount);
        }
    }
}
