public class Main {
    public static void main(String[] args) throws Exception {
        Matrix aMtrx = new Matrix(1000,1000);
        aMtrx.fill(1);
        Matrix bMtrx = new Matrix(1000,1000);
        bMtrx.fill(1);

        long startTime = System.currentTimeMillis();
        aMtrx.multiplyByDefinition(bMtrx);
        long estimatedTime = System.currentTimeMillis() - startTime;
        System.out.println("Умножается за время: " + estimatedTime + "ms");

        ParallelMatrixProduct test = new ParallelMatrixProduct(aMtrx, bMtrx.transpose(),
                Runtime.getRuntime().availableProcessors(), false);
        long startTime1 = System.currentTimeMillis();
        test.multiplyParallel();
        long estimatedTime1 = System.currentTimeMillis() - startTime1;
        System.out.println("Умножается за время: " + estimatedTime1 + "ms");
    }

}