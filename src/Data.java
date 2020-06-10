/**
 * Содержит информацию о наборе данных. Каждая строка содержит одну точку данных.
 * Первичные расчеты РСА осуществляются данные объекта.
 * @author	Kushal Ranjan
 * @version	051313
 */
class Data {
	double[][] matrix; //матрица[i] - это i-я строка; матрица[i][j] - это i-я строка, j-й столбец
	
	/**
	 * Построение новой матрицы данных.
	 * @param vals	данные для нового объекта данных; измерения в виде столбцов, точки данных в виде строк.
	 */
	Data(double[][] vals) {
		matrix = Matrix.copy(vals);
	}
	
	/**
	 * Тестовый код. Строит произвольную таблицу данных из 5 точек данных с 3 переменными,
     * нормализует ее и вычисляет ковариационную матрицу и ее собственные значения
     * и ортонормированные собственные векторы.
     * Затем определяются два главных компонента.
	 */
	public static void main(String[] args) {
		double[][] data = {{4.0, 4.2, 3.9, 4.3, 4.1}, {2.0, 2.1, 2, 2.1, 2.2},
				{0.6, 0.59, 0.58, 0.62, 0.63}};
		System.out.println("Raw data:");
		Matrix.print(data);
		Data dat = new Data(data);
		dat.center();
		double[][] cov = dat.covarianceMatrix();
		System.out.println("Ковариационная матрица:");
		Matrix.print(cov);
		EigenSet eigen = dat.getCovarianceEigenSet();
		double[][] vals = {eigen.values};
		System.out.println("Собственное значение:");
		Matrix.print(vals);
		System.out.println("Соответствующие собственные векторы:");
		Matrix.print(eigen.vectors);
		System.out.println("Два главных компонента:");
		Matrix.print(dat.buildPrincipalComponents(2, eigen));
		System.out.println("Преобразование основного компонента:");
		Matrix.print(Data.principalComponentAnalysis(data, 2));
	}
	
	/**
	 * PCA реализован с использованием алгоритма NIPALS.
     * Возвращаемое значение-это double [] [], где каждый double []
     * j-это массив оценок j-й точки данных, соответствующий желаемому числу главных компонент.
	 * @param input			input raw data array
	 * @param numComponents	desired number of PCs
	 * @return				the scores of the data array against the PCS
	 */
	static double[][] PCANIPALS(double[][] input, int numComponents) {
		Data data = new Data(input);
		data.center();
		double[][][] PCA = data.NIPALSAlg(numComponents);
		double[][] scores = new double[numComponents][input[0].length];
		for(int point = 0; point < scores[0].length; point++) {
			for(int comp = 0; comp < PCA.length; comp++) {
				scores[comp][point] = PCA[comp][0][point];
			}
		}
		return scores;
	}
	
	/**
	 * Реализация нелинейного итерационного алгоритма частичных наименьших квадратов на матрице данных
     * для данного объекта данных. Количество возвращаемых компьютеров определяется пользователем.
	 * @param numComponents	желаемое количество основных компонентов
	 * @return				a double[][][] где i-й дубль [] [] содержит
     * ti и pi, баллы и нагрузки соответственно I-го главного компонента.
	 */
	double[][][] NIPALSAlg(int numComponents) {
		final double THRESHOLD = 0.00001;
		double[][][] out = new double[numComponents][][];
		double[][] E = Matrix.copy(matrix);
		for(int i = 0; i < out.length; i++) {
			double eigenOld = 0;
			double eigenNew = 0;
			double[] p = new double[matrix[0].length];
			double[] t = new double[matrix[0].length];
			double[][] tMatrix = {t};
			double[][] pMatrix = {p};
			for(int j = 0; j < t.length; j++) {
				t[j] = matrix[i][j];
			}
			do {
				eigenOld = eigenNew;
				double tMult = 1/Matrix.dot(t, t);
				tMatrix[0] = t;
				p = Matrix.scale(Matrix.multiply(Matrix.transpose(E), tMatrix), tMult)[0];
				p = Matrix.normalize(p);
				double pMult = 1/Matrix.dot(p, p);
				pMatrix[0] = p;
				t = Matrix.scale(Matrix.multiply(E, pMatrix), pMult)[0];
				eigenNew = Matrix.dot(t, t);
			} while(Math.abs(eigenOld - eigenNew) > THRESHOLD);
			tMatrix[0] = t;
			pMatrix[0] = p;
			double[][] PC = {t, p}; //{scores, loadings}
			E = Matrix.subtract(E, Matrix.multiply(tMatrix, Matrix.transpose(pMatrix)));
			out[i] = PC;
		}
		return out;
	}
	
	/**
	 * Предыдущие алгоритмы выполнения PCA
	 */
	
	/**
	 * Выполняет анализ главных компонентов с заданным числом главных компонентов.
	 * @param input			input data; каждый двойной[] вход представляет собой массив значений
     *                         одной переменной для каждой точки данных
	 * @param numComponents	number of components desired
	 * @return				the transformed data set
	 */
	static double[][] principalComponentAnalysis(double[][] input, int numComponents) {
		Data data = new Data(input);
		data.center();
		EigenSet eigen = data.getCovarianceEigenSet();
		double[][] featureVector = data.buildPrincipalComponents(numComponents, eigen);
		double[][] PC = Matrix.transpose(featureVector);
		double[][] inputTranspose = Matrix.transpose(input);
		return Matrix.transpose(Matrix.multiply(PC, inputTranspose));
	}
	
	/**
	 * Возвращает список, содержащий основные компоненты этого набора данных с указанным числом загрузок.
	 * @param numComponents	the number of principal components desired
	 * @param eigen			Собственный набор, содержащий собственные значения и собственные векторы
	 * @return				the numComponents most significant eigenvectors
	 */
	double[][] buildPrincipalComponents(int numComponents, EigenSet eigen) {
		double[] vals = eigen.values;
		if(numComponents > vals.length) {
			throw new RuntimeException("Нельзя производить больше основных компонентов, чем те, которые предусмотрены.");
		}
		boolean[] chosen = new boolean[vals.length];
		double[][] vecs = eigen.vectors;
		double[][] PC = new double[numComponents][];
		for(int i = 0; i < PC.length; i++) {
			int max = 0;
			while(chosen[max]) {
				max++;
			}
			for(int j = 0; j < vals.length; j++) {
				if(Math.abs(vals[j]) > Math.abs(vals[max]) && !chosen[j]) {
					max = j;
				}
			}
			chosen[max] = true;
			PC[i] = vecs[max];
		}
		return PC;
	}
	
	/**
	 * Использует QR-алгоритм для определения собственных значений и собственных векторов ковариационной матрицы
     * для этого набора данных.
     * Итерация продолжается до тех пор, пока собственное значение не изменится более чем на 1/10000.
	 * @return	an EigenSet containing the eigenvalues and eigenvectors of the covariance matrix
	 */
	EigenSet getCovarianceEigenSet() {
		double[][] data = covarianceMatrix();
		return Matrix.eigenDecomposition(data);
	}
	
	/**
	 * Строит ковариационную матрицу для этого набора данных.
	 * @return	the covariance matrix of this data set
	 */
	double[][] covarianceMatrix() {
		double[][] out = new double[matrix.length][matrix.length];
		for(int i = 0; i < out.length; i++) {
			for(int j = 0; j < out.length; j++) {
				double[] dataA = matrix[i];
				double[] dataB = matrix[j];
				out[i][j] = covariance(dataA, dataB);
			}
		}
		return out;
	}
	
	/**
	 * Возвращает ковариацию двух векторов данных.
	 * @param a	double[] of data
	 * @param b	double[] of data
	 * @return	the covariance of a and b, cov(a,b)
	 */
	static double covariance(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new MatrixException("Cannot take covariance of different dimension vectors.");
		}
		double divisor = a.length - 1;
		double sum = 0;
		double aMean = mean(a);
		double bMean = mean(b);
		for(int i = 0; i < a.length; i++) {
			sum += (a[i] - aMean) * (b[i] - bMean);
		}
		return sum/divisor;
	}
	
	/**
	 * Центрирует каждый столбец матрицы данных в его среднем значении.
	 */
	void center() {
		matrix = normalize(matrix);
	}
	
	
	/**
	 * Нормализует входную матрицу так, чтобы каждый столбец был центрирован на 0.
	 */
	double[][] normalize(double[][] input) {
		double[][] out = new double[input.length][input[0].length];
		for(int i = 0; i < input.length; i++) {
			double mean = mean(input[i]);
			for(int j = 0; j < input[i].length; j++) {
				out[i][j] = input[i][j] - mean;
			}
		}
		return out;
	}
	
	/**
	 * Вычисляет среднее значение массива двойников.
	 * @param entries	input array of doubles
	 */
	static double mean(double[] entries) {
		double out = 0;
		for(double d: entries) {
			out += d/entries.length;
		}
		return out;
	}
}
