/**
 * Класс для выполнения матричных вычислений, специфичных для PCA.
 * @author	Kushal Ranjan
 * @version	051413
 */
class Matrix {
	
	static int numMults = 0; //Отслеживает количество выполненных умножений
	
	/**
	 * Test code for SVD. Uses example from MIT video: http://www.youtube.com/watch?v=cOUTpqlX-Xs
	 */
	public static void main(String[] args) {
		System.out.println("Original matrix:");
		double[][] test = {{5,-1},{5,7}}; //C
		Matrix.print(test);
		double[][][] SVD = Matrix.singularValueDecomposition(test);
		double[][] U = SVD[0];
		double[][] S = SVD[1];
		double[][] V = SVD[2];
		System.out.println("U-matrix:");
		Matrix.print(U);
		System.out.println("Sigma-matrix:");
		Matrix.print(S);
		System.out.println("V-matrix:");
		Matrix.print(V);
		System.out.println("Продукт разложения (C = US(V^T)):");
		Matrix.print(Matrix.multiply(U, Matrix.multiply(S, Matrix.transpose(V)))); //Should be C
	}
	
	
	/**
	 * Вычисляет сингулярное разложение (SVD)входной матрицы.
	 * @param input		the input matrix
	 * @return			the SVD of input, {U,S,V}, such that input = US(V^T). U and S are
	 * 					orthogonal matrix, and the non-zero entries of the diagonal matrix S are
	 * 					the 
	 */
	static double[][][] singularValueDecomposition(double[][] input) {
		double[][] C = Matrix.copy(input);
		double[][] CTC = multiply(transpose(C), C); //(C^T)C = V(S^T)S(V^T)
		EigenSet eigenC = eigenDecomposition(CTC);
		double[][] S = new double[C.length][C.length]; //Diagonal matrix
		for(int i = 0; i < S.length; i++) {
			S[i][i] = Math.sqrt(eigenC.values[i]); //Squareroots of eigenvalues are entries of S
		}
		double[][] V = eigenC.vectors;
		double[][] CV = multiply(C, V); //CV = US
		double[][] invS = copy(S); //Inverse of S
		for(int i = 0; i < invS.length; i++) {
			invS[i][i] = 1.0/S[i][i];
		}
		double[][] U = multiply(CV, invS); //U = CV(S^-1)
		return new double[][][] {U, S, V};
	}
	
	/**
	 * Определяет собственные значения и собственные векторы матрицы с помощью QR-алгоритма.
	 * Повторяется до тех пор, пока собственное значение не изменится более чем на 1/100000.
	 * @param	input	input matrix; must be square
	 * @return			an EigenSet containing the eigenvalues and corresponding eigenvectors of
	 * 					input
	 */
	static EigenSet eigenDecomposition(double[][] input) {
		if(input.length != input[0].length) {
			throw new MatrixException("Собственная композиция, не определенная на неквадратных матрицах.");
		}
		double[][] copy = copy(input);
		double[][] Q = new double[copy.length][copy.length];
		for(int i = 0; i < Q.length; i++) {
			Q[i][i] = 1; //Q starts as an identity matrix
		}
		boolean done = false;
		while(!done) {
			double[][][] fact = Matrix.QRFactorize(copy);
			double[][] newMat = Matrix.multiply(fact[1], fact[0]); //[A_k+1] := [R_k][Q_k]
			Q = Matrix.multiply(fact[0], Q);
			//Stop the loop if no eigenvalue changes by more than 1/100000
			for(int i = 0; i < copy.length; i++) {
				if(Math.abs(newMat[i][i] - copy[i][i]) > 0.00001) {
					copy = newMat;
					break;
				} else if(i == copy.length - 1) { //End of copy table
					done = true;
				}
			}
		}
		EigenSet ret = new EigenSet();
		ret.values = Matrix.extractDiagonalEntries(copy); //Eigenvalues lie on diagonal
		ret.vectors = Q; //Columns of Q converge to the eigenvectors
		return ret;
	}
	
	/**
	 * Создает массив диагональных записей во входной матрице.
	 * @param input	input matrix
	 * @return		the entries on the diagonal of input
	 */
	static double[] extractDiagonalEntries(double[][] input) {
		double[] out = new double[input.length];
		for(int i = 0; i<input.length; i++) {
			out[i] = input[i][i];
		}
		return out;
	}
	
	/**
	 * Выполняет QR факторизацию на входной матрице.
	 * @param input	input matrix
	 * @return		{Q, R}, the QR factorization of input.
	 */
	static double[][][] QRFactorize(double[][] input) {
		double[][][] out = new double[2][][];
		double[][] orthonorm = gramSchmidt(input);
		out[0] = orthonorm; //Q is the matrix of the orthonormal vectors formed by GS on input
		double[][] R = new double[orthonorm.length][orthonorm.length];
		for(int i = 0; i < R.length; i++) {
			for(int j = 0; j <= i; j++) {
				R[i][j] = dot(input[i], orthonorm[j]);
			}
		}
		out[1] = R;
		return out;
	}
	
	/**
	 * Преобразует входной список векторов в ортонормированный список с тем же диапазоном.
	 * @param input	list of vectors
	 * @return		orthonormal list with the same span as input
	 */
	static double[][] gramSchmidt(double[][] input) {
		double[][] out = new double[input.length][input[0].length];
		for(int outPos = 0; outPos < out.length; outPos++) {
			double[] v = input[outPos];
			for(int j = outPos - 1; j >= 0; j--) {
				double[] sub = proj(v, out[j]);
				v = subtract(v, sub); //Subtract off non-orthogonal components
			}
			out[outPos] = normalize(v); //return an orthonormal list
		}
		return out;
	}
	
	/**
	 * Возвращает заданную матрицу вращения с параметрами (i, j, th).
	 * @param size	total number of rows/columns in the matrix
	 * @param i		the first axis of the plane of rotation; i > j
	 * @param j		the second axis of the plane of rotation; i > j
	 * @param th	the angle of the rotation
	 * @return		the Givens rotation matrix G(i,j,th)
	 */
	static double[][] GivensRotation(int size, int i, int j, double th) {
		double[][] out = new double[size][size];
		double sine = Math.sin(th);
		double cosine = Math.cos(th);
		for(int x = 0; x < size; x++) {
			if(x != i && x != j) {
				out[x][x] = cosine;
			} else {
				out[x][x] = 1;
			}
		}
		out[i][j] = -sine;//ith column, jth row
		out[j][i] = sine;
		return out;
	}
	
	/**
	 * Возвращает транспонирование входной матрицы
	 * @param matrix	double[][] matrix of values
	 * @return			the matrix transpose of matrix
	 */
	static double[][] transpose(double[][] matrix) {
		double[][] out = new double[matrix[0].length][matrix.length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out[0].length; j++) {
				out[i][j] = matrix[j][i];
			}
		}
		return out;
	}
	
	/**
	 * Возвращает сумму a и b.
	 * @param a	double[][] matrix of values
	 * @param b	double[][] matrix of values
	 * @return	the matrix sum a + b
	 */
	static double[][] add(double[][] a, double[][] b) {
		if(a.length != b.length || a[0].length != b[0].length) {
			throw new MatrixException("Matrices not same size.");
		}
		double[][] out = new double[a.length][a[0].length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out[0].length; j++) {
				out[i][j] = a[i][j] + b[i][j];
			}
		}
		return out;
	}
	
	/**
	 * Возвращает разницу между a и b.
	 * @param a	double[][] matrix of values
	 * @param b	double[][] matrix of values
	 * @return	the matrix difference a - b
	 */
	static double[][] subtract(double[][] a, double[][] b) {
		if(a.length != b.length || a[0].length != b[0].length) {
			throw new MatrixException("Matrices not same size.");
		}
		double[][] out = new double[a.length][a[0].length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out[0].length; j++) {
				out[i][j] = a[i][j] - b[i][j];
			}
		}
		return out;
	}
	
	/**
	 * Возвращает сумму a и b.
	 * @param a	double[] vector of values
	 * @param b	double[] vector of values
	 * @return	the vector sum a + b
	 */
	static double[] add(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new MatrixException("Векторы не одинаковой длины.");
		}
		double[] out = new double[a.length];
		for(int i = 0; i < out.length; i++) {
			out[i] = a[i] + b[i];
		}
		return out;
	}
	
	/**
	 * Возвращает разницу между a и b.
	 * @param a	double[] vector of values
	 * @param b	double[] vector of values
	 * @return	the vector difference a - b
	 */
	static double[] subtract(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new MatrixException("Векторы не одинаковой длины.");
		}
		double[] out = new double[a.length];
		for(int i = 0; i < out.length; i++) {
			out[i] = a[i] - b[i];
		}
		return out;
	}
	
	/**
	 * Возвращает матричное произведение a и b; если горизонтальная длина a не равна вертикальной длине b, создает исключение.
	 * @param a	double[][] matrix of values
	 * @param b	double[][] matrix of values
	 * @return	the matrix product ab
	 */
	static double[][] multiply(double[][] a, double[][] b) {
		if(a.length != b[0].length) {
			throw new MatrixException("Матрицы, не совместимые для умножения.");
		}
		double[][] out = new double[b.length][a[0].length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out[0].length; j++) {
				double[] row = getRow(a, j);
				double[] column = getColumn(b, i);
				out[i][j] = dot(row, column);
			}
		}
		return out;
	}
	
	/**
	 * Возвращает версию,mat масштабируемую константой
	 * @param mat	input matrix
	 * @param coeff	constant by which to scale
	 * @return		mat scaled by coeff
	 */
	static double[][] scale(double[][] mat, double coeff) {
		double[][] out = new double[mat.length][mat[0].length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out[0].length; j++) {
				out[i][j] = mat[i][j] * coeff;
			}
		}
		return out;
	}
	
	/**
	 * Принимает точечное произведение двух векторов, {a[0]b[0], ..., a[n]b[n]}.
	 * @param a	double[] of values
	 * @param b	double[] of values
	 * @return	the dot product of a with b
	 */
	static double dot(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new MatrixException("Длины векторов не равны: " + a.length + "=/=" + b.length);
		}
		double sum = 0;
		for(int i = 0; i < a.length; i++) {
			numMults++;
			sum += a[i] * b[i];
		}
		return sum;
	}
	
	/**
	 * Возвращает копию входной матрицы.
	 * @param input	double[][] to be copied
	 */
	static double[][] copy(double[][] input) {
		double[][] copy = new double[input.length][input[0].length];
		for(int i = 0; i < copy.length; i++) {
			for(int j = 0; j < copy[i].length; j++) {
				copy[i][j] = input[i][j];
			}
		}
		return copy;
	}
	
	/**
	 * Возвращает I-й столбец входной матрицы.
	 */
	static double[] getColumn(double[][] matrix, int i) {
		return matrix[i];
	}
	
	/**
	 * Returns the ith row of the input matrix.
	 */
	static double[] getRow(double[][] matrix, int i) {
		double[] vals = new double[matrix.length];
		for(int j = 0; j < vals.length; j++) {
			vals[j] = matrix[j][i];
		}
		return vals;
	}
	
	/**
	 * Возвращает проекцию век на подпространство spanned by proj
	 * @param vec	вектор для проецирования
	 * @param proj	охватывающий вектор целевого подпространства
	 * @return		proj_proj(vec)
	 */
	static double[] proj(double[] vec, double[] proj) {
		double constant = dot(proj, vec)/dot(proj, proj);
		double[] projection = new double[vec.length];
		for(int i = 0; i < proj.length; i++) {
			projection[i] = proj[i]*constant;
		}
		return projection;
	}
	
	/**
	 * Возвращает нормализованную версию входного вектора, т. е. вектор масштабируется таким образом, что ||vec|| = 1.
	 * @return	vec/||vec||
	 */
	static double[] normalize(double[] vec) {
		double[] newVec = new double[vec.length];
		double norm = norm(vec);
		for(int i = 0; i < vec.length; i++) {
			newVec[i] = vec[i]/norm;
		}
		return newVec;
	}
	
	/**
	 * Вычисляет норму входного вектора
	 * @return ||vec||
	 */
	static double norm(double[] vec) {
		return Math.sqrt(dot(vec,vec));
	}
	
	/**
	 * Выводит входную матрицу с округлением каждого значения до 4 значащих цифр
	 */
	static void print(double[][] matrix) {
		for(int j = 0; j < matrix[0].length; j++) {
			for(int i = 0; i < matrix.length; i++) {
				double formattedValue = matrix[i][j];
				if(Math.abs(formattedValue) < 0.00001) { //Hide negligible values
					formattedValue = 0;
				}
				System.out.print(formattedValue + "\t");
			}
			System.out.print("\n");
		}
		System.out.println("");
	}
}

/**
 * Класс исключения, возникающий при попытке неверных вычислений матрицы
 */
class MatrixException extends RuntimeException {
	MatrixException(String string) {
		super(string);
	}
}

/**
 * Класс носителей данных, содержащий набор собственных значений и соответствующих им собственных векторов.
 */
class EigenSet {
	double[] values;
	double[][] vectors;
}