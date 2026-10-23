package com.github.abundmax.seemnv.core.fieldArray;

import java.util.Arrays;
import java.util.StringJoiner;
import java.util.function.IntFunction;

/**
 * A FieldArray represents all values of all samples of ONE VCF entry (line) for ONE genotypeField (e.g. GT, DP, etc..)
 *
 *
 * @param <T>
 */
public class FieldArray<T extends Comparable<? super T>>{

    public final FieldConfig fieldConfig;
    public final T[][] values;
    public final int length;

    public final IntFunction<T[]> rowFactory;
    public final IntFunction<T[][]> matrixFactory;


    public FieldArray(FieldConfig fieldConfig, IntFunction<T[]> rowFactory, IntFunction<T[][]> matrixFactory, T[][] values){
        this.fieldConfig = fieldConfig;
        this.values = values;
        this.length = values.length;
        this.rowFactory = rowFactory;
        this.matrixFactory = matrixFactory;
    }

    public FieldArray(FieldConfig fieldConfig, IntFunction<T[]> rowFactory, IntFunction<T[][]> matrixFactory, String[] initialStringArray){
        this.fieldConfig = fieldConfig;
        this.rowFactory = rowFactory;
        this.matrixFactory = matrixFactory;
        this.values = parse(initialStringArray);
        this.length = values.length;
    }

    /**
     * Parses the input string array into a two-dimensional array of the generic type T.
     * Each string in the input array is split based on commas. The resulting 2D array
     * is constructed using the provided factories and the specific field type for parsing.
     *
     * @param initialArray an array of strings to be parsed and converted into a 2D array.
     *                     Each string is split into individual elements using a comma as the delimiter.
     * @return a two-dimensional array of type T[][] where each element is parsed and cast
     *         according to the field type configuration.
     */
    private T[][] parse(String[] initialArray) {
        T[][] parsedArray = matrixFactory.apply(initialArray.length);

        for (int i = 0; i < initialArray.length; i++) {

            // If there is more than 1 value given at the Field, it has to be split
            String[] stringFieldArray = initialArray[i].split(",");

            parsedArray[i] = rowFactory.apply(stringFieldArray.length);

            for (int j = 0; j < stringFieldArray.length; j++) {
                parsedArray[i][j] = fieldConfig.fieldType().parseAndCast(stringFieldArray[j]);
            }
        }

        return parsedArray;
    }

    /**
     * Converts this Array to an Array of Strings
     *
     * <p>Example:</p>
     * <pre>
     *      List of Integer FieldArray:
     *      [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
     *
     *      returns:
     *      ["1,2,3", "4,5,6", "7,8,9"]
     * </pre>
     *
     * @return A string array containing the values of this FieldArray as strings
     */
    public String[] toStringArray() {
        String[] result = new String[this.length];
        StringBuilder valueBuilder = new StringBuilder();

        for (int i = 0; i < result.length; i++) {
            valueBuilder.setLength(0); // clear StringBuilder
            T[] valueArray = this.values[i];

            for (T value : valueArray) {
                valueBuilder.append(value).append(",");
            }

            valueBuilder.deleteCharAt(valueArray.length);
            result[i] = valueBuilder.toString();
        }

        return result;
    }

    /**
     * Converts the array of values at the specified index into a single comma-separated string.
     *
     * @param sampleIdx the index of the desired array in the {@code values} field.
     * @return a comma-separated string representation of the values in the array at the specified index.
     */
    public String valueAtAsString(int sampleIdx) {
        T[] valueArray = this.values[sampleIdx];
        StringJoiner joiner = new StringJoiner(",");

        for (T value : valueArray) {
            joiner.add(String.valueOf(value));
        }

        return joiner.toString();
    }

    /**
     * Combines several FieldArrays of the same type using the specified combiner
     * @param fieldArrays
     * @return
     */
    @SafeVarargs
    private static <T extends Comparable<? super T>> FieldArray<T> combine(
            FieldArray<T>... fieldArrays
    ) {
        if (fieldArrays == null || fieldArrays.length == 0) {
            throw new IllegalArgumentException("At least one FieldArray is required");
        }

        FieldArray<T> first = fieldArrays[0];
        FieldConfig fieldConfig = first.fieldConfig;
        FieldAggregationStrategy strategy = fieldConfig.aggregationStrategy();
        IntFunction<T[][]> matrixFactory = first.matrixFactory;

        validateCompatible(fieldArrays);

        int length = first.length;
        T[][] combinedArray = Arrays.copyOf(first.values, first.length);

        for (int sample = 0; sample < length; sample++) {
            T[][] valuesOverArraysAtSample = matrixFactory.apply(fieldArrays.length);

            for (int fieldArrayIndex = 0; fieldArrayIndex < fieldArrays.length; fieldArrayIndex++) {
                valuesOverArraysAtSample[fieldArrayIndex] =
                        fieldArrays[fieldArrayIndex].values[sample];
            }

            combinedArray[sample] = strategy.aggregate(valuesOverArraysAtSample);
        }

        return new FieldArray<>(
                fieldConfig,
                first.rowFactory,
                first.matrixFactory,
                combinedArray
        );
    }

    public static FieldArray<?> combineWildcard(FieldArray<?>... fieldArrays) {
        validateCompatibleWildcards(fieldArrays);
        return combineUnchecked(fieldArrays);
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static FieldArray<?> combineUnchecked(FieldArray<?>[] fieldArrays) {
        return combine((FieldArray[]) fieldArrays);
    }

    private static void validateCompatibleWildcards(FieldArray<?>[] fieldArrays) {
        if (fieldArrays == null || fieldArrays.length == 0) {
            throw new IllegalArgumentException("At least one FieldArray is required");
        }

        int length = fieldArrays[0].length;
        FieldConfig config = fieldArrays[0].fieldConfig;
        Class<?> type = config.fieldType().getType();

        for (FieldArray<?> fieldArray : fieldArrays) {
            if (fieldArray.length != length) {
                throw new IllegalStateException("FieldArrays have to be the same size");
            }

            if (fieldArray.fieldConfig != config) {
                throw new IllegalArgumentException("FieldArrays must have the same FieldConfig");
            }

            if (fieldArray.fieldConfig.fieldType().getType() != type) {
                throw new IllegalArgumentException("FieldArrays must have the same value type");
            }
        }
    }

    private static <T extends Comparable<? super T>> void validateCompatible(
            FieldArray<T>[] fieldArrays
    ) {
        int length = fieldArrays[0].length;
        FieldConfig config = fieldArrays[0].fieldConfig;

        for (FieldArray<T> fieldArray : fieldArrays) {
            if (fieldArray.length != length) {
                throw new IllegalStateException("FieldArrays have to be the same size");
            }

            if (fieldArray.fieldConfig != config) {
                throw new IllegalArgumentException("FieldArrays must have the same FieldConfig");
            }
        }
    }
}
