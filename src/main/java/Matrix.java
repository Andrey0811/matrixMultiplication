import java.util.Arrays;

public class Matrix {

    private final int[][] matrix;

    public Matrix(int rows, int columns) throws Exception {
        if (rows <= 0 || columns <= 0)
            throw new Exception("Заданы неверные размеры матрицы");

        this.matrix = new int[rows][columns];
    }

    public Matrix(int size) throws Exception {
        if (size <= 0)
            throw new Exception("Заданы неверные размеры матрицы");

        this.matrix = new int[size][size];
    }

    public Matrix(int[] arr, int size) throws Exception {
        if (size <= 0 && arr.length == size * size)
            throw new Exception("Заданы неверные размеры матрицы");

        this.matrix = new int[size][size];
        for (int i = 0; i < size; i++)
            System.arraycopy(arr, i * size, this.matrix[i], 0, size);

    }

    public Matrix(int[][] matrix) {
        this.matrix = matrix.clone();
    }

    public int getRows() {
        return matrix.length;
    }

    public int getColumns() {
        return matrix[0].length;
    }

    public int[][] getMatrix() {
        return matrix;
    }

    public Matrix sum(Matrix matrix2) throws Exception {
        if (this.getRows() != matrix2.getRows()
                || this.getColumns() != matrix2.getColumns())
            throw new Exception("Размеры матриц не совпадают");
        int[][] matrix = new int[this.getRows()][this.getColumns()];
        for (int i = 0; i < this.getRows(); i++)
            for (int j = 0; j < this.getColumns(); j++)
                matrix[i][j] = this.matrix[i][j] + matrix2.getMatrix()[i][j];
        return new Matrix(matrix);
    }

    public Matrix multiplyByDefinition(Matrix other) throws Exception {
        if (checkSizes(other))
            throw new Exception("Размеры матриц не совпадают");
        int[][] result = new int[getRows()][other.getColumns()];
        for (int i = 0; i < getRows(); i++)
            for (int j = 0; j < other.getColumns(); j++)
                for (int k = 0; k < getColumns(); k++)
                    result[i][j] +=
                            matrix[i][k] * other.getMatrix()[k][j];
        return new Matrix(result);
    }

    public Matrix multiplyByDefinitionWithOptimization(Matrix other) throws Exception {
        if (checkSizes(other))
            throw new Exception("Размеры матриц не совпадают");

        int columns = other.getColumns();
        Matrix temp = other.transpose();
        int[][] result = new int[getRows()][other.getColumns()];

        for (int i = 0; i < getRows(); i++)
            for (int j = 0; j < columns; j++)
                for (int k = 0; k < getColumns(); k++) {
                    result[i][j] +=
                            matrix[i][k] * temp.getMatrix()[j][k];
                }
        return new Matrix(result);
    }

    public Matrix multiplyBlock(Matrix other) throws Exception {
        if (!other.isSquare() && isSquare()
                && other.getRows() == getRows())
            throw new Exception();
        int n = getRows() * getRows();
        int[] mA = new int[n];
        int[] mB = new int[n];
        int[] mC = new int[n];
        for (int i = 0; i < getRows(); i++)
            for (int j = 0; j < getRows(); j++){
                int idx = i * getRows() + j;
                mA[idx] = matrix[i][j];
                mB[idx] = other.getMatrix()[i][j];
                mC[idx] = 0;
            }
        multiplyBlockIn(mC, 0, mA, 0, mB, 0, getRows(), getRows());
            return new Matrix(mC, getRows());
    }

    private void multiplyBlockIn(int[] result, int shiftResult, int[] a, int shiftA,
                                 int[] b, int shiftB, int tempSize, int rowSize) {
        int d11 = 0;
        // умеем умнажать матрицы 2 х 2
        if (tempSize == 2) {
            int d12 = 1;
            int d21 = rowSize;
            int d22 = rowSize + 1;

            result[d11] += a[d11] * b[d11] + a[d12] * b[d21];
            result[d12] += a[d11] * b[d12] + a[d12] * b[d22];
            result[d21] += a[d21] * b[d11] + a[d22] * b[d21];
            result[d22] += a[d21] * b[d12] + a[d22] * b[d22];
        }
        else {
            final int d12 = tempSize / 2;
            final int d21 = (tempSize / 2) * rowSize;
            final int d22 = (tempSize / 2) * (rowSize + 1);

            final int result11 = shiftResult + d11;
            final int a11 = shiftA + d11;
            final int b11 = shiftB + d11;

            final int result12 = shiftResult + d12;
            final int a12 = shiftA + d12;
            final int b12 = shiftB + d12;

            final int result21 = shiftResult + d21;
            final int a21 = shiftA + d21;
            final int b21 = shiftB + d21;

            final int result22 = shiftResult + d22;
            final int a22 = shiftA + d22;
            final int b22 = shiftB + d22;

            // C11 += A11 * B11
            multiplyBlockIn(result, result11, a, a11, b, b11, tempSize / 2, rowSize);
            // C11 += A12 * B21
            multiplyBlockIn(result, result11, a, a12, b, b21, tempSize / 2, rowSize);
            // C12 += A11 * B12
            multiplyBlockIn(result, result12, a, a11, b, b12, tempSize / 2, rowSize);
            // C12 += A12 * B22
            multiplyBlockIn(result, result12, a, a12, b, b22, tempSize / 2, rowSize);
            // C21 += A21 * B11
            multiplyBlockIn(result, result21, a, a21, b, b11, tempSize / 2, rowSize);
            // C21 += A22 * B21
            multiplyBlockIn(result, result21, a, a22, b, b21, tempSize / 2, rowSize);
            // C22 += A21 * B12
            multiplyBlockIn(result, result22, a, a21, b, b12, tempSize / 2, rowSize);
            // C22 += A22 * B22
            multiplyBlockIn(result, result22, a, a22, b, b22, tempSize / 2, rowSize);
        }
    }

    public Matrix multiplyWinograd(Matrix other) throws Exception {
        if (getColumns() < 2)
            throw new Exception("Размеры матрицы некорректны");
        
        int coeff = getColumns() / 2;
        int[] rowFactor = new int[getRows()];
        int[] columnFactor = new int[other.getColumns()];

        // вычисление rowFactors для this
        for (int i = 0; i < getRows(); i++) {
            rowFactor[i] = matrix[i][0] * matrix[i][1];
            for (int j = 2; j < coeff; j++)
                rowFactor[i] = rowFactor[i] + matrix[i][2 * j - 1] * matrix[i][2*j];
        }

        // вычисление	columnFactors для other
        for (int i = 0; i < other.getColumns(); i++) {
            columnFactor[i] = other.getMatrix()[0][i] * other.getMatrix()[1][i];
            for (int j = 2; j < coeff; j++)
                columnFactor[i] = columnFactor[i]
                        + other.getMatrix()[2 * j - 1][i] * other.getMatrix()[2 * j][i];
        }

        int[][] result = new int[getRows()][other.getColumns()];

        // вычисление матрицы result
        for (int i = 0; i < getRows(); i++)
            for (int j = 0; j < other.getColumns(); j++) {
                result[i][j] = - rowFactor[i] - columnFactor[j];
                for (int k = 0; k < coeff; k++)
                    result[i][j] = result[i][j]
                            + (matrix[i][2 * k] + other.getMatrix()[2 * k + 1][j])
                            * (matrix[i][2 * k + 1] + other.getMatrix()[2 * k][j]);
            }

        // прибавление членов в случае нечетной общей размерности
        if (2 * (getColumns()/2) != getColumns())
            for (int i = 0; i < getRows(); i++)
                for (int j = 0; j < other.getColumns(); j++)
                    result[i][j] = result[i][j]
                            + matrix[i][getColumns() - 1]
                            * other.getMatrix()[getColumns() - 1][j];
        return new Matrix(result);
    }

    public Matrix multiplyStrassen(Matrix other) throws Exception {
        if (isSquare()
                && other.isSquare()
                && isPowerOfTwo(getRows())
                && isPowerOfTwo(other.getRows()))
            return new Matrix(multiplyStrassenIn(matrix, other.getMatrix()));
        throw new Exception("Неподходящие матрицы");
    }

    private int[][] multiplyStrassenIn(int[][] first, int[][] second){
        int size = first.length;
        int[][] res = new int[size][size];

        // если размер матрицы 1x1
        if (size == 1) {
            res[0][0] = first[0][0] * second[0][0];
        } else {

            // 1 матрица
            int[][] a = new int[size / 2][size / 2];
            int[][] b = new int[size / 2][size / 2];
            int[][] c = new int[size / 2][size / 2];
            int[][] d = new int[size / 2][size / 2];

            // 2 матрица
            int[][] e = new int[size / 2][size / 2];
            int[][] f = new int[size / 2][size / 2];
            int[][] g = new int[size / 2][size / 2];
            int[][] h = new int[size / 2][size / 2];

            // деление матрицы a на 4 части
            divideArray(first, a, 0, 0);
            divideArray(first, b, 0, size / 2);
            divideArray(first, c, size / 2, 0);
            divideArray(first, d, size / 2, size / 2);

            // деление матрицы b на 4 части
            divideArray(second, e, 0, 0);
            divideArray(second, f, 0, size / 2);
            divideArray(second, g, size / 2, 0);
            divideArray(second, h, size / 2, size / 2);

            /*p1 = (a + d)(e + h)
             p2 = (c + d)e
             p3 = a(f - h)
             p4 = d(g - e)
             p5 = (a + b)h
             p6 = (c - a) (e + f)
             p7 = (b - d) (g + h)*/
            int[][] p1 = multiplyStrassenIn(addMatrices(a, d), addMatrices(e, h));
            int[][] p2 = multiplyStrassenIn(addMatrices(c, d), e);
            int[][] p3 = multiplyStrassenIn(a, subMatrices(f, h));
            int[][] p4 = multiplyStrassenIn(d, subMatrices(g, e));
            int[][] p5 = multiplyStrassenIn(addMatrices(a, b), h);
            int[][] p6 = multiplyStrassenIn(subMatrices(c, a), addMatrices(e, f));
            int[][] p7 = multiplyStrassenIn(subMatrices(b, d), addMatrices(g, h));

            /*C11 = p1 + p4 - p5 + p7
             C12 = p3 + p5
             C21 = p2 + p4
             C22 = p1 - p2 + p3 + p6 */
            int[][] C11 = addMatrices(subMatrices(addMatrices(p1, p4), p5), p7);
            int[][] C12 = addMatrices(p3, p5);
            int[][] C21 = addMatrices(p2, p4);
            int[][] C22 = addMatrices(subMatrices(addMatrices(p1, p3), p2), p6);

            // собираем результат
            copySubArray(C11, res, 0, 0);
            copySubArray(C12, res, 0, size / 2);
            copySubArray(C21, res, size / 2, 0);
            copySubArray(C22, res, size / 2, size / 2);
        }
        return res;
    }

    public int[][] addMatrices(int[][] first, int[][] second) {
        int n = first.length;
        int[][] res = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res[i][j] = first[i][j] + second[i][j];
            }
        }
        return res;
    }

    public int[][] subMatrices(int[][] first, int[][] second) {
        int n = first.length;
        int[][] res = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res[i][j] = first[i][j] - second[i][j];
            }
        }
        return res;
    }

    public void divideArray(int[][] first, int[][] second,
                            int startRowIdx, int startColumnIdx) {
        for (int i1 = 0, i2 = startRowIdx; i1 < second.length; i1++, i2++)
            for (int j1 = 0, j2 = startColumnIdx; j1 < second.length; j1++, j2++)
                second[i1][j1] = first[i2][j2];
    }

    public void copySubArray(int[][] first, int[][] second,
                             int startRowIdx, int startColumnIdx) {
        for (int i1 = 0, i2 = startRowIdx; i1 < first.length; i1++, i2++)
            for (int j1 = 0, j2 = startColumnIdx; j1 < first.length; j1++, j2++)
                second[i2][j2] = first[i1][j1];
    }

    private boolean checkSizes(Matrix matrix2) {
        return !((this.getRows() == matrix2.getRows())
                && (this.getColumns() == matrix2.getColumns())
                || (this.getColumns() == matrix2.getRows()));
    }

    public Matrix transpose() throws Exception {
        Matrix trans = new Matrix(getColumns(), getRows());
        for (int i = 0; i < getRows(); i++)
            for (int j = 0; j < getColumns(); j++)
                trans.getMatrix()[j][i] = this.matrix[i][j];
        return trans;
    }

    private boolean isPowerOfTwo(int n){
        double item = Math.log(n) / Math.log(2);
        return (int)(Math.ceil(item))
                == (int)(Math.floor(item));
    }

    @Override
    public String toString() {
        StringBuilder build = new StringBuilder();
        for (int[] ints : matrix) {
            for (int j = 0; j < matrix[0].length; ++j)
                build.append("[").append(ints[j]).append("]");
            build.append("\n");
        }
        return build.toString();
    }

    public void fill(int number) {
        for (int i = 0; i < getRows(); i++)
            Arrays.fill(matrix[i], number);
    }

    public boolean isSquare() {
        return getRows() == getColumns();
    }

    @Override
    public boolean equals(Object obj){
        if (this == obj)
            return true;
        if (obj instanceof Matrix){
            Matrix other = (Matrix) obj;
            if (other.getColumns() != getColumns()
                    || other.getRows() != getRows())
                return false;

            for (int i = 0; i < getRows(); i++)
                for (int j = 0; j < getColumns(); j++)
                    if (matrix[i][j] != other.getMatrix()[i][j])
                        return false;
                    return true;
        }
        return false;
    }
}

