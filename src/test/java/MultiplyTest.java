import org.junit.Before;
import org.junit.Test;

public class MultiplyTest {
    Matrix m1;
    Matrix m2;
    Matrix m3;
    Matrix m4;
    Matrix m5;
    Matrix m6;
    Matrix m7;
    Matrix m8;
    Matrix m9;
    Matrix m10;
    Matrix m11;
    Matrix m12;

    @Before
    public void setUp() {
        m1 = new Matrix(new int[][]{
                {4, 2},
                {9, 0}
        });
        m2 = new Matrix(new int[][]{
                {3, 1},
                {-3, 4}
        });
        m3 = new Matrix(new int[][]{
                {1, 0},
                {0, 1}
        });
        m4 = new Matrix(new int[][]{
                {2, 1},
                {-3, 0},
                {4, -1}
        });
        m5 = new Matrix(new int[][]{
                {5, -1, 6},
                {-3, 0, 7}
        });
        m6 = new Matrix(new int[][]{
                {7, -2, 19},
                {-15, 3, -18},
                {23, -4, 17}
        });
        m7 = new Matrix(new int[][]{
                {6, 12},
                {27, 9}
        });
        m8 = new Matrix(new int[][]{
                {4, -1, 20},
                {3, 56, 8},
                {0, -9, 2},
                {-3, 2, 23},
                {7, 1, 9},
                {15, 6, 12},
                {0, 99, 65}
        });
        m9 = new Matrix(new int[][]{
                {91, -16},
                {-130, -5},
                {35, -2},
                {80, -26},
                {47, -2},
                {60, 3},
                {-37, -65}
        });
        m10 = new Matrix(new int[][]{
                {1, 2, 3, 4, 5, 6, 7, 8},
                {12, 13, 14, 15, 16, 17, 18, 19},
                {23, 24, 25, 26, 27, 28, 29, 0},
                {-71, 35, 36, 37, 38, 39, -1, -33},
                {45, 46, 47, 48, 49, -42, 43, -41},
                {56, 57, 58, 59, 95, 85, 75, 65},
                {67, 68, 69, 96, 86, 76, -2, 0},
                {0, 43, 99, 96, 44, 10, 59, 8},
        });
        m11 = new Matrix(new int[][]{
                {1, 4, 65, 57, 83, 8, 9, 0},
                {-1, 2, 39, -1, -84, 739, 0, 32},
                {23, 17, 45, 24, 74, 54, 67, 90},
                {80, -23, 41, 16, 10, 91, 70, 2},
                {1, -6, 72, 60, 69, 21, 3, 6},
                {5, 8, 9, 7, 1, 1, 1, 1},
                {32, 42, 56, 87, 90, 23, 79, 33},
                {1, 2, 3, 0, 4, 5, 6, 7}
        });
        m12 = new Matrix(new int[][]{
                {655,	 295,	 1272,	 1142,	 1190,	 2324,	 1112,	  665},
                {2217,	 801,	 4902,	 3892,	 3907,	12686,	 3697,	 2546},
                {3749,	1247,	 8442,	 6642,	 6504,	22898,	 6102,	 4217},
                {3850,	-477,	 2819,	 -160,	-3360,	31257,	 4239,	 4437},
                {6094,	1061,	14237,	10802,	10874,	43031,	10170,	 7182},
                {9038,	3357,	22892,	18291,	18392,	55202,	15205,	10747},
                {9668,	-623,	20812,	12461,	11745,	65086,	12122,	 9104},
                {11904,	1871,	16654,	11712,	13062,	48190,	18204,	12755},

        });
    }

    @Test
    public void testTranspose() throws Exception {
        assert m1.transpose().equals(new Matrix(new int[][]{
                {4, 9},
                {2, 0}
        }));
        assert m4.transpose().equals(new Matrix(new int[][]{
                {2, -3, 4},
                {1, 0, -1}
        }));
    }

    @Test
    public void testSum() throws Exception {
        assert m1.sum(m2).equals(new Matrix(new int[][]{
                {7, 3},
                {6, 4}
        }));
    }

    @Test(expected = Exception.class)
    public void testException() throws Exception {
        m7.multiplyByDefinition(m8);
        m7.multiplyByDefinitionWithOptimization(m8);
        m7.multiplyStrassen(m8);
        m7.multiplyWinograd(m8);
        m1.sum(m7);
    }

    @Test
    public void test_1_2() throws Exception {
        Matrix temp = m1.multiplyByDefinition(m2);
        assert temp.equals(m7);
        temp = m1.multiplyByDefinitionWithOptimization(m2);
        assert temp.equals(m7);
        temp = m1.multiplyStrassen(m2);
        assert temp.equals(m7);
        temp = m1.multiplyWinograd(m2);
        assert temp.equals(m7);
        assert !m1.multiplyByDefinition(m2).equals(m2.multiplyByDefinition(m1));
        assert m1.multiplyByDefinition(m2).equals(m1.multiplyByDefinition(m2));
        temp = m1.multiplyBlock(m2);
        assert temp.equals(m7);
    }

    @Test
    public void test_1_3() throws Exception {
        Matrix temp = m1.multiplyByDefinition(m3);
        assert temp.equals(m1);
        temp = m1.multiplyByDefinitionWithOptimization(m3);
        assert temp.equals(m1);
        temp = m1.multiplyStrassen(m3);
        assert temp.equals(m1);
        temp = m1.multiplyWinograd(m3);
        assert temp.equals(m1);
    }

    @Test
    public void test_4_5() throws Exception {
        Matrix temp = m4.multiplyByDefinition(m5);
        assert temp.equals(m6);
        temp = m4.multiplyByDefinitionWithOptimization(m5);
        assert temp.equals(m6);
        temp = m4.multiplyWinograd(m5);
        assert temp.equals(m6);
    }

    @Test
    public void test_8_4() throws Exception {
        Matrix temp = m8.multiplyByDefinition(m4);
        assert temp.equals(m9);
        temp = m8.multiplyByDefinitionWithOptimization(m4);
        assert temp.equals(m9);
        temp = m8.multiplyWinograd(m4);
        assert temp.equals(m9);
    }

    @Test
    public void test_10_11() throws Exception {
        Matrix temp = m10.multiplyByDefinition(m11);
        assert temp.equals(m12);
        temp = m10.multiplyByDefinitionWithOptimization(m11);
        assert temp.equals(m12);
        temp = m10.multiplyStrassen(m11);
        assert temp.equals(m12);
//        temp = m10.multiplyWinograd(m11);
//        assert temp.equals(m12);
    }
}
