import java.util.Random;

import Plotter.Plotter;

public class Sorting {

	final static int SELECTION_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int MERGE_VS_QUICK_SORTED_LENGTH = 12;
	final static int SELECT_VS_MERGE_LENGTH = 16;
	final static double T = 600.0;

	/**
	 * Sorts a given array using the quick sort algorithm. At each stage the
	 * pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity
	 * of O(n^2)
	 * 
	 * @param arr
	 *            - the array to be sorted
	 */
	public static void quickSort(double[] arr) {
		quickSort(arr, 0, arr.length - 1);
	}

	private static void quickSort(double[] arr, int p, int r) {
		if (p < r) {
			int q = partition(arr, p, r);
			quickSort(arr, p, q - 1);
			quickSort(arr, q + 1, r);
		}
	}

	private static int partition(double[] arr, int p, int r) {
		double x = arr[r];
		int i = p;
		for (int j = p; j < r; j++) {
			if (arr[j] <= x) {
				swap(arr, i, j);
				i++;
			}
		}
		swap(arr, i, r);
		return (i);
	}

	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr
	 *            - the array to be sorted
	 */
	public static void mergeSort(double[] arr) {
		mergeSort(arr, 0, arr.length - 1);
	}

	private static void mergeSort(double[] arr, int left, int right) {
		if (left < right) {
			int q = (left + right) / 2;
			mergeSort(arr, left, q);
			mergeSort(arr, q + 1, right);
			merge(arr, left, q, right);
		}
	}

	private static void merge(double[] arr, int left, int q, int right) {
		double[] Left = new double[q - left + 2];
		double[] Right = new double[right - q + 1];
		for (int i = 0; i < Left.length - 1; i++) {
			Left[i] = arr[left + i];
		}
		for (int j = 0; j < Right.length - 1; j++) {
			Right[j] = arr[q + 1 + j];
		}
		Right[Right.length - 1] = Double.MAX_VALUE;
		Left[Left.length - 1] = Double.MAX_VALUE;
		int rIndex = 0;
		int lIndex = 0;
		for (int k = left; k <= right; k++)
			if (Right[rIndex] <= Left[lIndex]) {
				arr[k] = Right[rIndex];
				rIndex++;
			} else {
				arr[k] = Left[lIndex];
				lIndex++;
			}
	}

	/**
	 * finds the i'th order statistic of a given array.
	 * 
	 * Should run in complexity O(n) in the worst case.
	 * 
	 * @param arr
	 *            - the array.
	 * @param i
	 *            - a number between 0 and arr.length - 1.
	 * @return the number which would be at index i, if the array was to be
	 *         sorted
	 */

	public static double select(double[] arr, int i) {
		return select(arr, 0, arr.length - 1, i);
	}

	private static double select(double[] arr, int left, int right, int i){
		if (left == right)
			return arr[left];
		int rest = (right - left + 1) % 5;
		double[] B = new double[(right - left)/5 + 1];
		if (rest == 0){
			for (int k =0; k < B.length; k ++){
				int j = k*5 + left;
				B[k]= median(arr[j], arr[j+1], arr[j+2], arr[j+3], arr[j+4]);
			}
		}
		else{
			double[] temp = new double[rest];
			for (int l = 0; l < rest; l++){
				temp[l] = arr[right - l];
			}
			quickSort(temp);
			B[B.length - 1] = temp[rest / 2];
			for (int m = 0; m < B.length - 1; m ++){
				int n = m*5 + left;
				B[m]= median(arr[n], arr[n+1], arr[n+2], arr[n+3], arr[n+4]);
			}
		}
		double pivot = select(B, 0, B.length - 1, B.length/2);
		int q = partition(arr,left, right, pivot);
		if(q == i)
			return arr[q];
		else if (q < i){ 
			return select(arr, q + 1, right, i);
		}
		else {
			return select(arr, left, q - 1, i);
		}
	}

	private static int partition(double[] arr, int left, int right, double pivot) {
		int location = 0;
		for (int i = left; i <= right; i++) {
			if (arr[i] == pivot) {
				location = i;
				break;
			}
		}
		swap(arr, location, right);
		return partition(arr, left, right);
	}

	private static double median(double a, double b, double c, double d, double e) {
		double[] arr = { a, b, c, d, e };
		quickSort(arr);
		return arr[2];

	}
	/**
	 * Sorts a given array using the selection sort algorithm.
	 * 
	 * Should run in complexity O(n^2) in the worst case.
	 * 
	 * @param arr
	 *            - the array to be sorted
	 */
	public static void selectionSort(double[] arr) {
		for (int i = 0; i < arr.length - 1; i++) {
			int min = i;
			for (int j = i + 1; j < arr.length; j++) {
				if (arr[j] < arr[min]) {
					min = j;
				}
			}
			if (min != i) {
				swap(arr, i, min);
			}
		}
	}

	private static void swap(double[] arr, int i, int j) {
		double temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
	
	}

	public static void main(String[] args) {
		selectionVsQuick();
		mergeVsQuick();
	    mergeVsQuickOnSortedArray();
	    selectVsMerge(); 
	}

	/**
	 * Compares the selection sort algorithm against quick sort on random arrays
	 */
	public static void selectionVsQuick() {
		double[] quickTimes = new double[SELECTION_VS_QUICK_LENGTH];
		double[] selectionTimes = new double[SELECTION_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < SELECTION_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumSelection = 0;
			for (int k = 0; k < T; k++) {
				int size = (int) Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				selectionSort(b);
				endTime = System.currentTimeMillis();
				sumSelection += endTime - startTime;
			}
			quickTimes[i] = sumQuick / T;
			selectionTimes[i] = sumSelection / T;
		}
		Plotter.plot("quick sort", quickTimes, "selection sort", selectionTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick() {
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int) Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick / T;
			mergeTimes[i] = sumMerge / T;
		}
		Plotter.plot("quick sort", quickTimes, "merge sort", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray() {
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int) Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick / T;
			mergeTimes[i] = sumMerge / T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}

	/**
	 * Compares the select algorithm against sorting an array.
	 */
	public static void selectVsMerge() {
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] selectTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		double x;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumMerge = 0;
			long sumSelect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int) Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				int index = (int) (Math.random() * size);
				startTime = System.currentTimeMillis();
				mergeSort(a);
				x = a[index];
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				x = select(b, index);
				endTime = System.currentTimeMillis();
				sumSelect += endTime - startTime;
			}
			mergeTimes[i] = sumMerge / T;
			selectTimes[i] = sumSelect / T;
		}
		Plotter.plot("merge sort and select", mergeTimes, "select", selectTimes);
	}
}
